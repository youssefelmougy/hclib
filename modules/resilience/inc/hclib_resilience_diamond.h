#ifndef HCLIB_RESILIENCE_DIAMOND_H
#define HCLIB_RESILIENCE_DIAMOND_H

#include <vector>
#include <mutex>

namespace ref_count = hclib::ref_count;

namespace hclib {
namespace resilience {
namespace diamond {

/*
vector that allows concurrent insertion
Read operations are not concurrent with insertion
*/
template<typename T>
class safe_vector {

    std::vector<T> vec;
    std::mutex mtx;

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
        std::lock_guard<std::mutex> lkg(mtx);
        vec.push_back(std::forward<T>(value));
    }

    auto begin() noexcept -> decltype(vec.begin()) {
        return vec.begin();
    }

    auto end() noexcept -> decltype(vec.end()) {
        vec.end();
    }

    auto size() const noexcept -> decltype(vec.size()){
        return vec.size();
    }

    auto data() noexcept -> decltype(vec.data()) {
        return vec.data();
    }
};

/*
Reslient Promise for pointer type data
*/

template<typename T>
class future_t: public ref_count::future_t<T> {

    void ref_count_decr();

  public:

    void release() {
      ref_count_decr();
    }
};

//TODO: To enable non pointer and reference types, 
//create template specialisation for pointer and reference
//and get the required data pointer as given in default impl.
//and store that data to tmp_data.
//the specialisations only override put().
//put_acutal() remains same
template<typename T, int N=3>
class promise_t: public ref_count::promise_t<T> {
    //replica i stores datum in tmp_data[i] when put is performed
    //and later at the end of diamond task it performs actual put
    T tmp_data[N];

    //since there are N replicas we need to reduce the actual
    //reference count when all the N replicas perform release
    //This variable keeps track of the number of replica release
    int *diamond_ref_count = nullptr;

  public:

    static const int TYPE = ref_count::promise_t<T>::TYPE+1;

    promise_t(int n, ref_count::DelType del_type = ref_count::DelType::NORMAL)
    : ref_count::promise_t<T>(n),
      diamond_ref_count(new int(0)) {
        hclib_promise_t::type = TYPE;
    }


    void ref_count_decr() {
        if(hclib_promise_t::type < TYPE) return;
        assert(diamond_ref_count != nullptr);

        //decrement the actual reference count when all the
        //N replicas preform the release
        int d_cnt = __sync_add_and_fetch(diamond_ref_count, 1);
        if(d_cnt % N == 0)
	    ref_count::promise_t<T>::ref_count_decr();
    }

    bool equal(int i, int j) {
        obj* obj1 = static_cast<obj*>(tmp_data[i]);
        obj* obj2 = static_cast<obj*>(tmp_data[j]);

        //Each replica should have allocated its own data
        //otherwise all replicas are modifying the same data
        assert(obj1 != obj2);

        return obj1->equals(obj2);
    }

    void put(T datum);

    void put_actual(int i) {
        hclib::promise_t<T>::put(tmp_data[i]);
    }
};

template<typename T>
void future_t<T>::ref_count_decr() {
    promise_t<T> * p = static_cast<promise_t<T>*>(hclib_future_t::owner);
    p->ref_count_decr();
}

/*
safe vector of promises
*/
template<typename T>
class promise_vector : public safe_vector<promise_t<T>*> {
};

/*
metadata that is passed between tasks created within each replica
*/
template<typename T>
struct diamond_task_params_t {
    //replica id
    int index;
    //vector which captures all the puts
    promise_vector<T> *put_vec;
};

/*
put operation only adds data to a temporary storage indexed 
by its replica id. And the first replica adds the promise 
to the vector later used to perform the actual put
*/
template<typename T, int N>
void promise_t<T,N>::put(T datum) {
    HASSERT_STATIC(std::is_pointer<T>::value,
        "resilient promise_t currently only supports pointers\n");

    void * task_local = *hclib_get_curr_task_local();

    //if called from non-diamond task, just perform the put
    if(task_local == nullptr)
        hclib::promise_t<T>::put(datum);

    //if called from diamond task, save to a temp and delay the put
    //also add it to the vector that captures all put opertations
    else {
        auto diamond_params = (diamond_task_params_t<T>*)(task_local);
        tmp_data[diamond_params->index] = datum;

        //TODO:makes assumption that all replicas call put() on same promises.
        //i.e. no replica misses to do the put and the error might only be inside the data.
        //To enable catching of missed puts by some replica, one option is to use
        //use atomics and make the first replica that calls the put perform the push_back
        if(diamond_params->index == 0)
            diamond_params->put_vec->push_back(this);
    }
}

/*
Resilient async_await
*/

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    //fetch the diamond_task_parameters from task local and pass is to the async
    auto dtp = (diamond_task_params_t<void*>*)(*hclib_get_curr_task_local());

    hclib::async_await( [=, lambda_mv = std::move(lambda)] () {
        *(hclib_get_curr_task_local()) = dtp;
	lambda_mv();
        if(future1 != nullptr)
            static_cast<future_t<void*>*>(future1)->release();
        if(future2 != nullptr)
            static_cast<future_t<void*>*>(future2)->release();
        if(future3 != nullptr)
            static_cast<future_t<void*>*>(future3)->release();
        if(future4 != nullptr)
            static_cast<future_t<void*>*>(future4)->release();
    }, future1, future2, future3, future4);
}

bool check_result_helper(promise_vector<void*>* put_vec,
        int replica1, int replica2) ;

//TODO: currently only works for triple redundancy
int check_result(promise_vector<void*>* put_vec, int N);

/*
async_await  triplicates the task and checks for equivalence of all puts
*/
template <typename T, int N=3>
void async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    HASSERT_STATIC(N==3, "Currently only supports triple redudancy\n");

    hclib::async_await([=, lambda_mv1 = std::move(lambda)]() {
      //dtp_arr is used to pass around required data of each replica
      auto dtp_arr = new diamond_task_params_t<void*>[N];
      //put_vec is used to collect all the put() operations
      //peformed within the finish scope
      auto put_vec = new promise_vector<void*>();

      hclib::finish([=, lambda_mv2 = std::move(lambda)]() {
        //create N tasks and give each of it dtp_arr[i]
        for(int i=0; i<N; i++) {
            dtp_arr[i].index = i;
            dtp_arr[i].put_vec = put_vec;
            *(hclib_get_curr_task_local()) = &dtp_arr[i];
            async_await(lambda_mv2, f1, f2, f3, f4);
        }
      });

      int result = check_result(put_vec, N);
      if(result < 0 )
          prom_check->put(0);
      else {
          prom_check->put(1);
          //perform the actual put
          for(auto && elem: *put_vec)
              elem->put_actual(result);
      }

      delete put_vec;
      delete[] dtp_arr;
    }, nullptr);
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif

