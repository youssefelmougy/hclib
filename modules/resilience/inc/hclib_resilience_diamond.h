#ifndef HCLIB_RESILIENCE_DIAMOND_H
#define HCLIB_RESILIENCE_DIAMOND_H

#define N_CNT 3

namespace ref_count = hclib::ref_count;

namespace hclib {
namespace resilience {
namespace diamond {

/*
Reslient Future and Promise for pointer type data
*/

template<typename T, int N=N_CNT>
class future_t: public ref_count::future_t<T> {

    void ref_count_decr();

  public:

    void release() {
      ref_count_decr();
    }
};

//template specialization for N_CNT
//when the replication count is 2 (even when replication count
//is 2, it needs to use promises with triple redundancy)
//add futures to a vector so that they can later be released.
template<typename T>
class future_t<T, N_CNT>: public ref_count::future_t<T> {

    void ref_count_decr();
    void add_future_vector();

  public:

    void release() {
      ref_count_decr();
      add_future_vector();
    }
};

//Currently only pointer type promise is allowed
template<typename T, int N=N_CNT>
class promise_t: public ref_count::promise_t<T> {
    promise_t() = delete;
    promise_t(const promise_t &) = delete;
    promise_t(promise_t &&) = delete;
};


//TODO: To enable non pointer and reference types, 
//create template specialisation for pointer and reference
//and get the required data pointer as given in default impl.
//and store that data to tmp_data.
//the specialisations only override put().
//put_acutal() remains same
template<typename T, int N>
class promise_t<T*, N>: public ref_count::promise_t<T*> {
    //replica i stores datum in tmp_data[i] when put is performed
    //and later at the end of diamond task it performs actual put
    T* tmp_data[N] = {nullptr};

    //since there are N replicas we need to reduce the actual
    //reference count when all the N replicas perform release
    //This variable keeps track of the number of replica release
    int *diamond_ref_count = nullptr;

  public:

    static const int TYPE = ref_count::promise_t<T*>::TYPE+1;

    promise_t(int ref_cnt, ref_count::DelType del_type = ref_count::DelType::NORMAL)
    : ref_count::promise_t<T*>(ref_cnt, del_type),
      diamond_ref_count(new int(0)) {
          HASSERT_STATIC(N>=3, "N should be greater than 3.\
              For 2-replica do not provide value of N.\
              Default value will be automically used.\n");

          hclib_promise_t::type = TYPE;
    }

    void ref_count_decr() {
        if(hclib_promise_t::type < TYPE) return;
        assert(diamond_ref_count != nullptr);

        auto task_local = * hclib_get_curr_task_local();
        //invoked from outside of diamond tasks
        if(task_local == nullptr) 
            ref_count::promise_t<T*>::ref_count_decr();
        //invoked from inside diamond tasks
        else{
            //decrement the actual reference count when all the
            //N replicas preform the release
            int d_cnt = __sync_add_and_fetch(diamond_ref_count, 1);
            if(d_cnt % N == 0)
	        ref_count::promise_t<T*>::ref_count_decr();
        }
    }

    bool equal(int i, int j) {
        obj* obj1 = static_cast<obj*>(tmp_data[i]);
        obj* obj2 = static_cast<obj*>(tmp_data[j]);

        //Each replica should have allocated its own data
        //otherwise all replicas are modifying the same data
        assert(obj1 != obj2);
        return obj1->equals(obj2);
    }

    void put(T* datum);

    future_t<T*,N> *get_future() {
        return static_cast<future_t<T*,N>*>( hclib::promise_t<T*>::get_future());
    }

    void put_actual(int replica_index) {
        for(int i=0; i<N; i++)
          if( i!=replica_index )
             (* ref_count::promise_t<T*>::deleter)(tmp_data[i]);
        hclib::promise_t<T*>::put(tmp_data[replica_index]);
    }
};

/*
safe vector of promises and futures
*/
template<typename T, int N>
using promise_vector = hclib::resilience::util::safe_vector<promise_t<T, N>*>;

template<typename T, int N>
using future_vector = hclib::resilience::util::safe_vector<future_t<T, N>*>;

/*
metadata that is passed between tasks created within each replica
*/

struct diamond_task_params_base_t {
    int index; //replica id
    int count; //replica count TODO: can we delete this field?
};

//all put operations should be added to a list, so that they
//can be later compared.
template<typename T, int N>
struct diamond_task_params_t : public diamond_task_params_base_t {
    //vector which captures all the puts
    promise_vector<T,N> *put_vec;
};

//when replica count is 2, the futures that are used as task dependencies
//should be added to a list, so that they can be released later
template<typename T>
struct diamond_task_params_t<T,2> : public diamond_task_params_base_t {
    //vector which captures all the puts
    promise_vector<T,N_CNT> *put_vec;
    //vector of futures that performed release
    future_vector<T,N_CNT> *rel_vec;
};

/*
definition of future_t and promise_t functions
*/
template<typename T, int N>
void future_t<T,N>::ref_count_decr() {
    promise_t<T,N> * p = static_cast<promise_t<T,N>*>(hclib_future_t::owner);
    p->ref_count_decr();
}

template<typename T>
void future_t<T,N_CNT>::ref_count_decr() {
    promise_t<T,N_CNT> * p = static_cast<promise_t<T,N_CNT>*>(hclib_future_t::owner);
    p->ref_count_decr();
}

template<typename T>
void future_t<T,N_CNT>::add_future_vector() {
    auto task_local = (diamond_task_params_base_t*)(*hclib_get_curr_task_local());
    if(task_local == nullptr) return;

    //when number of replica is two, the first replica will add the promise
    //to a vector so that it can be later released for the third time 
    if(task_local->count == 2 && task_local->index==0) {
        auto diamond_params = (diamond_task_params_t<T,2>*)(task_local);
        diamond_params->rel_vec->push_back(this);
    }
}

/*
put operation only adds data to a temporary storage indexed 
by its replica id. And the first replica adds the promise 
to the vector later used to perform the actual put
*/
template<typename T, int N>
void promise_t<T*,N>::put(T* datum) {
    //HASSERT_STATIC(std::is_pointer<T*>::value,
    //    "resilient promise_t currently only supports pointers\n");

    void * task_local = *hclib_get_curr_task_local();

    //if called from non-diamond task, just perform the put
    if(task_local == nullptr)
        hclib::promise_t<T*>::put(datum);

    //if called from diamond task, save to a temp and delay the put
    //also add it to the vector that captures all put opertations
    else {
        auto diamond_params = (diamond_task_params_t<T*,N>*)(task_local);
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

template<int N, typename T>
bool check_result_helper(promise_vector<T,N>* put_vec,
        int replica1, int replica2) {

    for(auto && elem: *put_vec) {
        if(!elem->equal(replica1, replica2))
            return false;
    }
    return true;
}

//TODO: currently only works for triple redundancy
template<int N, typename T>
int check_result(promise_vector<T,N>* put_vec) {
    if(check_result_helper<N>(put_vec, 0, 1))
        return 0;
    else if(check_result_helper<N>(put_vec, 1, 2))
        return 1;
    else if(check_result_helper<N>(put_vec, 0, 2))
        return 0;
    else
        return -1;
}

template <int N=N_CNT-1, typename T>
std::enable_if_t< N>=N_CNT-1, void>
async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    //fetch the diamond_task_parameters from task local and pass is to the async
    auto dtp = *hclib_get_curr_task_local();

    constexpr int N_NEW = (N == (N_CNT-1))? N+1: N;
    hclib::async_await( [=, lambda_mv = std::move(lambda)] () {
        *(hclib_get_curr_task_local()) = dtp;
	lambda_mv();
        if(future1 != nullptr)
            static_cast<future_t<void*, N_NEW>*>(future1)->release();
        if(future2 != nullptr)
            static_cast<future_t<void*, N_NEW>*>(future2)->release();
        if(future3 != nullptr)
            static_cast<future_t<void*, N_NEW>*>(future3)->release();
        if(future4 != nullptr)
            static_cast<future_t<void*, N_NEW>*>(future4)->release();
    }, future1, future2, future3, future4);
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif

