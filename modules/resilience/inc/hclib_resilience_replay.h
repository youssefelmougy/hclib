
#ifndef HCLIB_RESILIENCE_REPLAY_H
#define HCLIB_RESILIENCE_REPLAY_H

namespace hclib {
namespace resilience {
namespace replay {


bool is_replay_task(void *ptr)
{
    return nullptr != ptr;
}

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

    static const int TYPE = ref_count::promise_t<T*>::TYPE+2;

    //datum is stored in tmp_data when put is performed
    //and later at the end of replay task it performs actual put
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

    void put_actual() {
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
using promise_vector = hclib::resilience::util::safe_vector<promise_t<T>*>;

template<typename T>
using future_vector = hclib::resilience::util::safe_vector<future_t<T>*>;

/*
metadata that is passed between tasks created within each replay
*/
//all put operations should be added to a list, so that they
//can be later compared.
//all futures that are used as task dependencies should be
//added to a list, so that they can be later released.
template<typename T>
struct replay_task_params_t {
    int index; //replay id

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
    auto task_local = static_cast<replay_task_params_t<T>*>(*hclib_get_curr_task_local());
    assert(is_replay_task(task_local));
    task_local->rel_vec->push_back(this);
}

/*
put operation only adds data to a temporary storage.
Also adds the promise to the vector which is
later used to perform the actual put
*/
template<typename T>
void promise_t<T*>::put(T* datum) {
    auto task_local = static_cast<replay_task_params_t<T*>*>(*hclib_get_curr_task_local());

    //if called from replay task, save to a temp and delay the put
    //also add it to the vector that captures all put opertations
    if(is_replay_task(task_local)){
        tmp_data = datum;
        task_local->put_vec->push_back(this);
    }
    //if called from non-replay task, just perform the put
    else {
        hclib::promise_t<T*>::put(datum);
    }
}

template <typename T>
void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    //fetch the replay_task_parameters from task local and pass is to the async
    auto rtp = *hclib_get_curr_task_local();

    hclib::async_await( [=, lambda_mv = std::move(lambda)] () {
        *(hclib_get_curr_task_local()) = rtp;
        lambda_mv();
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
    }, future1, future2, future3, future4);
}

int get_replay_index() {
  auto task_local = static_cast<replay_task_params_t<void*>*>(*hclib_get_curr_task_local());
  assert(is_replay_task(task_local));
  return task_local->index;
}

template <int N=2, typename T>
void async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    hclib::async_await([=, lambda_mv1 = std::move(lambda)]() {
        auto rtp = new replay_task_params_t<void*>();;
	rtp->put_vec = new promise_vector<void*>();
	rtp->rel_vec = new future_vector<void*>();
        bool result = false;

        assert(*(hclib_get_curr_task_local()) ==  nullptr);
	for(int i=0; i<N; i++) {
	    hclib::finish([=]() {
              rtp->index = i;
	      *(hclib_get_curr_task_local()) = rtp;
	      async_await(lambda_mv1, f1, f2, f3, f4);
	    });
            //ran successfully without errors
	    if(error_check_fn(params) == 1) {
                result = true;
                break;
            }
            for(auto && elem: *(rtp->put_vec))
               elem->tmp_delete();
            rtp->put_vec->clear();
            rtp->rel_vec->clear();
	}

        *(hclib_get_curr_task_local()) = nullptr;

        if(result) {
            //perform the actual put
            for(auto && elem: *(rtp->put_vec))
                elem->put_actual();
            //perform the release
            for(auto && elem: *(rtp->rel_vec))
                elem->release();
            prom_check->put(1);
        }
	else {
	    prom_check->put(0);
	}
	delete rtp->put_vec;
	delete rtp->rel_vec;
        delete rtp;;
    }, f1, f2, f3, f4);
}

} // namespace replay
} // namespace resilience
} // namespace hclib

#endif
