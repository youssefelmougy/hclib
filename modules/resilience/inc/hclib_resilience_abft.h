
#ifndef HCLIB_RESILIENCE_ABFT_H
#define HCLIB_RESILIENCE_ABFT_H

namespace hclib {
namespace resilience {
namespace abft {

//TODO: make the check stronger
inline bool is_abft_task(void *ptr)
{
    return nullptr != ptr;
}

#ifndef USE_RESILIENT_PROMISE

/*
Reslient Future and Promise for pointer type data
*/

template<typename T>
class future_t: public ref_count::future_t<T> {
  public:
    void add_future_vector();
};

//Currently only pointer type promise is allowed
template<typename T>
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
template<typename T>
class promise_t<T*>: public ref_count::promise_t<T*> {

    static const int TYPE = ref_count::promise_t<T*>::TYPE+3;

    //datum is stored in tmp_data when put is performed
    //and later at the end of actual/abft task it performs actual put
    T* tmp_data = nullptr;

  public:

    promise_t(int ref_cnt, ref_count::DelType del_type = ref_count::DelType::NORMAL)
    : ref_count::promise_t<T*>(ref_cnt, del_type) {
          hclib_promise_t::type = TYPE;
    }

    void put(T* datum);

    future_t<T*>* get_future() {
        return static_cast<future_t<T*>*>( hclib::promise_t<T*>::get_future());
    }

    void put_actual(int index) {
        hclib::promise_t<T*>::put(tmp_data);
    }

    void tmp_delete() {
        (* ref_count::promise_t<T*>::deleter)(tmp_data);
    }
};

/*
safe vector of promises and futures
*/
template<typename T>
using promise_vector = hclib::resilience::util::safe_promise_vector<promise_t<T>*>;

template<typename T>
using future_vector = hclib::resilience::util::safe_future_vector<future_t<T>*>;

/*
metadata that is passed between tasks created within each run
*/
//all put operations should be added to a list, so that they
//can be later compared.
//all futures that are used as task dependencies should be
//added to a list, so that they can be later released.
template<typename T>
struct abft_task_params_t {
    int index; //normal or correction

    //vector which captures all the puts
    promise_vector<T> *put_vec;

    //vector which capture task dependencies
    future_vector<T> *rel_vec;
};

/*
definition of future_t and promise_t functions
*/

template<typename T>
void future_t<T>::add_future_vector() {
    auto task_local = static_cast<abft_task_params_t<T>*>(*hclib_get_curr_task_local());
    assert(is_abft_task(task_local));
    task_local->rel_vec->push_back(this);
}

/*
put operation only adds data to a temporary storage.
Also adds the promise to the vector which is
later used to perform the actual put
*/
template<typename T>
void promise_t<T*>::put(T* datum) {
    auto task_local = static_cast<abft_task_params_t<T*>*>(*hclib_get_curr_task_local());

    //if called from abft task, save to a temp and delay the put
    //also add it to the vector that captures all put opertations
    if(is_abft_task(task_local)){
        if(task_local->index == 0)
            task_local->put_vec->push_back(this);
        else if(task_local->index == 1)
            tmp_delete();
        else assert(false);

        tmp_data = datum;
    }
    //if called from non-abft task, just perform the put
    else {
        hclib::promise_t<T*>::put(datum);
    }
}

#endif // USE_RESILIENT_PROMISE

inline int get_index() {
  auto task_local = static_cast<abft_task_params_t<void*>*>(*hclib_get_curr_task_local());
  assert(is_abft_task(task_local));
  return task_local->index;
}

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    //fetch the abft_task_parameters from task local and pass is to the async
    auto atp = *hclib_get_curr_task_local();

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);
    hclib::async_await( [=] () {
        *(hclib_get_curr_task_local()) = atp;
        (*lambda_ptr)();
        delete lambda_ptr;
        auto task_local = static_cast<abft_task_params_t<void*>*>(*hclib_get_curr_task_local());
        //TODO: assuming the same future will be used in all replays
        if(task_local->index == 0) {
            if(future1 != nullptr)
                //static_cast<future_t<void*>*>(future1)->release();
                static_cast<future_t<void*>*>(future1)->add_future_vector();
            if(future2 != nullptr)
                //static_cast<future_t<void*>*>(future2)->release();
                static_cast<future_t<void*>*>(future2)->add_future_vector();
            if(future3 != nullptr)
                //static_cast<future_t<void*>*>(future3)->release();
                static_cast<future_t<void*>*>(future3)->add_future_vector();
            if(future4 != nullptr)
                //static_cast<future_t<void*>*>(future4)->release();
                static_cast<future_t<void*>*>(future4)->add_future_vector();
        }
    }, future1, future2, future3, future4);
}

template <typename T>
inline void async_await(T&& lambda, std::vector<hclib_future_t *> *futures) {

    //fetch the abft_task_parameters from task local and pass is to the async
    auto atp = *hclib_get_curr_task_local();

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);
    hclib::async_await( [=] () {
        *(hclib_get_curr_task_local()) = atp;
        (*lambda_ptr)();
        delete lambda_ptr;
        auto task_local = static_cast<abft_task_params_t<void*>*>(*hclib_get_curr_task_local());
        //TODO: assuming the same future will be used in all replays
        if(task_local->index == 0) {
            for(auto && elem: *(futures))
                static_cast<future_t<void*>*>(elem)->add_future_vector();
        }
    }, futures);
}

template <typename T1, typename T2>
void async_await_check(T1&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        T2&& abft_lambda,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    typedef typename std::remove_reference<T1>::type U1;
    U1* lambda_ptr = new U1(lambda);
    typedef typename std::remove_reference<T2>::type U2;
    U2* abft_lambda_ptr = new U2(abft_lambda);

    hclib::async_await([=]() {
        auto atp = new abft_task_params_t<void*>();;
        atp->put_vec = new promise_vector<void*>();
        atp->rel_vec = new future_vector<void*>();
        bool result = false;
        int index = 0;

        assert(*(hclib_get_curr_task_local()) ==  nullptr);
        hclib::finish([=]() {
            atp->index = 0;
            *(hclib_get_curr_task_local()) = atp;
            async_await(*lambda_ptr, f1, f2, f3, f4);
        });
        result = error_check_fn(params);
        //ran with errors
        if(result == 0){
            atp->index = 1;
            index = 1;
            (*abft_lambda_ptr)();
            result = error_check_fn(params);
        }
        delete lambda_ptr;
        delete abft_lambda_ptr;

        *(hclib_get_curr_task_local()) = nullptr;

        if(result) {
            atp->put_vec->do_puts(index);
            atp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete atp->put_vec;
        delete atp->rel_vec;
        delete atp;;
    }, f1, f2, f3, f4);
}

template <typename T1, typename T2>
void async_await_check(T1&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        T2&& abft_lambda,
        std::vector<hclib_future_t *> *futures) {

    typedef typename std::remove_reference<T1>::type U1;
    U1* lambda_ptr = new U1(lambda);
    typedef typename std::remove_reference<T2>::type U2;
    U2* abft_lambda_ptr = new U2(abft_lambda);

    hclib::async_await([=]() {
        auto atp = new abft_task_params_t<void*>();;
        atp->put_vec = new promise_vector<void*>();
        atp->rel_vec = new future_vector<void*>();
        bool result = false;
        int index = 0;

        assert(*(hclib_get_curr_task_local()) ==  nullptr);
        hclib::finish([=]() {
            atp->index = 0;
            *(hclib_get_curr_task_local()) = atp;
            async_await(*lambda_ptr, futures);
        });
        result = error_check_fn(params);
        //ran with errors
        if(result == 0){
            atp->index = 1;
            index = 1;
            (*abft_lambda_ptr)();
            result = error_check_fn(params);
        }

        delete lambda_ptr;
        delete abft_lambda_ptr;

        *(hclib_get_curr_task_local()) = nullptr;

        if(result) {
            atp->put_vec->do_puts(index);
            atp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete atp->put_vec;
        delete atp->rel_vec;
        delete atp;;
    }, futures);
}

} // namespace abft
} // namespace resilience
} // namespace hclib

#endif
