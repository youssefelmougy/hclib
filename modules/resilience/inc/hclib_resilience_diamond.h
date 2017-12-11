#ifndef HCLIB_RESILIENCE_DIAMOND_H
#define HCLIB_RESILIENCE_DIAMOND_H

namespace diamond {

/*
vector that allows concurrent insertion
*/
template<typename T>
class safe_vector {

    std::vector<T> vec;
    std::mutex mtx;

  public:
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
class future_t: public hclib::ref_count::future_t<T> {

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
class promise_t: public hclib::ref_count::promise_t<T> {
    T tmp_data[N];
    int *diamond_ref_count = nullptr;

  public:

    promise_t(int n) : hclib::ref_count::promise_t<T>(n) {
        diamond_ref_count = new int(0);
    }

    void ref_count_decr() {
        if(diamond_ref_count != nullptr) {
            int d_cnt = __sync_add_and_fetch(diamond_ref_count, 1);
            if(d_cnt % N == 0)
	        hclib::ref_count::promise_t<T>::ref_count_decr();
	}
	else
	    hclib::ref_count::promise_t<T>::ref_count_decr();
    }

    bool equal(int i, int j) {
        return !memcmp(tmp_data[i], tmp_data[j], sizeof(int));
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

template<typename T>
struct diamond_task_params_t {
    int index = 0;
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
    //auto put_vec = (safe_vector<promise_t*>*)(*hclib_get_curr_task_local());
    auto diamond_params = (diamond_task_params_t<T>*)(*hclib_get_curr_task_local());
    tmp_data[diamond_params->index] = datum;
    //TODO:makes assumption that all replicas call put() on same promises.
    //i.e. no replica misses to do the put and the error might only be inside the data.
    //To enable catching of missed puts by some replica, one option is to use
    //use atomics and make the first replica that calls the put perform the push_back
    if(diamond_params->index == 0)
        diamond_params->put_vec->push_back(this);
}

/*
Resilient async_await
*/
template <typename T>
class AsyncResilience  {
    T _lambda;
    hclib_future_t *future1, *future2, *future3, *future4;
    diamond_task_params_t<void*> *dtp;

    public:
        AsyncResilience(T lambda, hclib_future_t *f1, hclib_future_t *f2,
            hclib_future_t *f3, hclib_future_t *f4,
            diamond_task_params_t<void*>* tmp_dtp)
	: _lambda(lambda), future1(f1), future2(f2)
        , future3(f3), future4(f4), dtp(tmp_dtp) {}

    //set the task local to the diamond_task_parameters and invoke the lambda
    void operator() () {
        *(hclib_get_curr_task_local()) = dtp;
	_lambda();
        future_t<void*> *f = static_cast<future_t<void*>*>(future1);
        if(f != NULL)
            f->release();
        f = static_cast<future_t<void*>*>(future2);
        if(f != NULL)
            f->release();
        f = static_cast<future_t<void*>*>(future3);
        if(f != NULL)
            f->release();
        f = static_cast<future_t<void*>*>(future4);
        if(f != NULL)
            f->release();
    }
};

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    //fetch the diamond_task_parameters from task local and pass is to the async
    auto dtp = (diamond_task_params_t<void*>*)(*hclib_get_curr_task_local());
    AsyncResilience<T> async_res(lambda, future1, future2, future3, future4, dtp);
    hclib::async_await(std::move(async_res), future1, future2, future3, future4);
}

bool check_result_helper(promise_vector<void*>* put_vec, int index1, int index2) {
    for(auto && elem: *put_vec) {
        //printf("eq %p %d %d %d\n", elem, *(int*)elem->tmp_data[index1], 
        //*(int*)elem->tmp_data[index2], elem->equal(index1, index2));
        if(!elem->equal(index1, index2))
	    return false;
    }
    return true;
}

//TODO: currently only works for triple redundancy
int check_result(promise_vector<void*>* put_vec) {
    //check task 1 with 2
    if(check_result_helper(put_vec, 0, 1))
        return 0;
    else if(check_result_helper(put_vec, 1, 2))
        return 1;
    else if(check_result_helper(put_vec, 0, 2))
        return 0;
    else 
        return -1;
}

/*
async_await  triplicates the task and checks for equivalence of all puts
*/
template <typename T, int N=3>
void async_await_check(T&& lambda, 
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    HASSERT_STATIC(N==3, "Currently only supports triple redudancy\n");

    //dtp_arr is used to pass around required data of each replica
    auto dtp_arr = new diamond_task_params_t<void*>[N];
    //put_vec is used to collect all the put() operations
    //peformed within the finish scope
    auto put_vec = new promise_vector<void*>();

    hclib::finish([=]() {
      //create N tasks and give each of it dtp_arr[i]
      for(int i=0; i<N; i++) {
          dtp_arr[i].index = i;
          dtp_arr[i].put_vec = put_vec;
          *(hclib_get_curr_task_local()) = &dtp_arr[i];
          async_await(lambda, f1, f2, f3, f4);
      }
    });

    
    int result = check_result(put_vec);
    if(result < 0 ) {
        printf("Error discoverd in diamond task. unable to recover\n");
	exit(0);
    };
    int put_size = put_vec->size();
    auto put_vec_data = put_vec->data();
    
    for(int i=0; i<put_size; i++) {
        put_vec_data[i]->put_actual(result);
    }

    delete put_vec;
    delete[] dtp_arr;
}

} // namespace diamond

#endif

