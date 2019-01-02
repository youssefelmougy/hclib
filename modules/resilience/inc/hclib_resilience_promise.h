#ifndef HCLIB_RESILIENCE_PROMISE_H
#define HCLIB_RESILIENCE_PROMISE_H

#define N_CNT 3
#define FINISH_WORKAROUND

namespace hclib {
namespace resilience {

#ifdef USE_RESILIENT_PROMISE

inline bool is_resilient_task(void *ptr)
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

    static const int TYPE = ref_count::promise_t<T*>::TYPE+1;

    //replica i stores datum in tmp_data[i] when put is performed
    //and later at the end of resilient task it performs actual put
    T* tmp_data[N_CNT] = {nullptr};

  public:

    promise_t(int ref_cnt, ref_count::DelType del_type = ref_count::DelType::NORMAL)
    : ref_count::promise_t<T*>(ref_cnt, del_type) {
          hclib_promise_t::type = TYPE;
    }

    void put(T* datum);

    future_t<T*> *get_future() {
        return static_cast<future_t<T*>*>( hclib::promise_t<T*>::get_future());
    }

    void put_actual(const int index) {
    for(int i=index+1; i<index+N_CNT; i++)
        (* ref_count::promise_t<T*>::deleter)(tmp_data[i%N_CNT]);
        hclib::promise_t<T*>::put(tmp_data[index]);
    }

    void tmp_delete() {
    }

    bool equal(int i, int j) {
        obj* obj1 = static_cast<obj*>(tmp_data[i]);
        obj* obj2 = static_cast<obj*>(tmp_data[j]);

        //Each replica should have allocated its own data
        //otherwise all replicas are modifying the same data
        assert(obj1 != obj2);
        return obj1->equals(obj2);
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
metadata that is passed between tasks created within each replica

All put operations should be added to a list, so that they
can be later compared.
All futures that are used as task dependencies should be
added to a list, so that they can be later released.
*/
template<typename T>
struct resilient_task_params_t {
    int index; //replica id

    //vector which captures all the puts
    promise_vector<T> *put_vec;

    //vector which capture task dependencies
    future_vector<T> *rel_vec;

#ifdef FINISH_WORKAROUND
    volatile int *count;
    hclib::promise_t<int> *finish_prom;
#endif

};

/*
definition of future_t and promise_t functions
*/
template<typename T>
inline void future_t<T>::add_future_vector() {
    auto task_local = static_cast<resilient_task_params_t<T>*>(*hclib_get_curr_task_local());
    assert(is_resilient_task(task_local));
    task_local->rel_vec->push_back(this);
}

/*
put operation only adds data to a temporary storage indexed
by its replica id. And the first replica adds the promise
to the vector later used to perform the actual put
*/
template<typename T>
inline void promise_t<T*>::put(T* datum) {
    //HASSERT_STATIC(std::is_pointer<T*>::value,
    //    "resilient promise_t currently only supports pointers\n");
    auto task_local = static_cast<resilient_task_params_t<T*>*>(*hclib_get_curr_task_local());

    //if called from resilient task, save to a temp and delay the put
    //also add it to the vector that captures all put opertations
    if(is_resilient_task(task_local)){
        tmp_data[task_local->index] = datum;

        //TODO: assuming the same promise will be called for put in all replays
        if(task_local->index==0)
            task_local->put_vec->push_back(this);
    }
    //if called from non-resilient task, just perform the put
    else {
        hclib::promise_t<T*>::put(datum);
    }
}

inline int get_index() {
    auto task_local = static_cast<resilient_task_params_t<void*>*>(*hclib_get_curr_task_local());
    assert(is_resilient_task(task_local));
    return task_local->index;
}

template<typename T>
using diamond_task_params_t = resilient_task_params_t<T>;

template<typename T>
using replay_task_params_t = resilient_task_params_t<T>;

template<typename T>
using abft_task_params_t = resilient_task_params_t<T>;

namespace diamond {
template<typename T>
using promise_t = ::hclib::resilience::promise_t<T>;

template<typename T>
using future_t = ::hclib::resilience::future_t<T>;
}

namespace replay {
template<typename T>
using promise_t = ::hclib::resilience::promise_t<T>;

template<typename T>
using future_t = ::hclib::resilience::future_t<T>;
}

namespace abft {
template<typename T>
using promise_t = ::hclib::resilience::promise_t<T>;

template<typename T>
using future_t = ::hclib::resilience::future_t<T>;
}
#endif // USE_RESILIENT_PROMISE

} // namespace resilience
} // namespace hclib

#endif

