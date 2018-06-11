
#ifndef HCLIB_RESILIENCE_CHECKPOINT_H
#define HCLIB_RESILIENCE_CHECKPOINT_H

namespace hclib {
namespace resilience {
namespace checkpoint {

/*
The object that will be saved to the archive.
*/
struct archive_obj {

    //size of object blob
    int size = -1;
    //blob of object to be archived
    void *data = nullptr;

    archive_obj() {}

    archive_obj(archive_obj *ptr) {
        assert(ptr->size > 0);
        size = ptr->size;
        data = malloc(size);
        memcpy(data, ptr->data, size); 
    }

    ~archive_obj() {
        assert(size > 0);
        free(data);
    }
};

/*
The user object should include a serialize method that
will convert the current objct to an archive object
and a constructor that can create the current object
baack from the archive object
*/
class obj: public hclib::resilience::obj {
  public:
    //obj(archive_obj* ptr) { assert(false); }
    virtual archive_obj* serialize() {}
};

class archive_store {
  public:
    virtual ~archive_store() {}

    virtual void save(void* key, archive_obj* data) {}

    virtual archive_obj* retrieve(void* key) {}
};

static archive_store* store_var = nullptr;

inline void set_archive_store(archive_store *ptr) {
    assert(store_var == nullptr && ptr != nullptr);
    store_var = ptr;
}

inline archive_store* get_archive_store() {
    return store_var;
}

//TODO: make the check stronger
inline bool is_checkpoint_task(void *ptr)
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
    T get();
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

    static const int TYPE = ref_count::promise_t<T*>::TYPE+4;

    //datum is stored in tmp_data when put is performed
    //and later at the end of actual/replay task it performs actual put
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
        archive_obj* ar_ptr =  ((obj*) tmp_data)->serialize();
        get_archive_store()->save(this, ar_ptr);
        delete ar_ptr;
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
struct checkpoint_task_params_t {
    int index; //replay index

    int checkpoint_run; //use checkpoint data

    //vector which captures all the puts
    promise_vector<T> *put_vec;

    //vector which capture task dependencies
    future_vector<T> *rel_vec;
};

/*
definition of future_t and promise_t functions
*/
template<typename T>
T future_t<T>::get() {
    auto task_local = static_cast<checkpoint_task_params_t<T>*>(*hclib_get_curr_task_local());
    if(is_checkpoint_task(task_local) && task_local->checkpoint_run == 1) {
        archive_obj *ptr = get_archive_store()->retrieve(hclib_future_t::owner);
        if(ptr == nullptr)
            return static_cast<T>(hclib_future_get(this));
        else {
            typedef typename std::remove_pointer<T>::type U;
            //create compile time error if U* != T, i.e. T has to be a pointer type
            T data = new U(ptr);
            delete ptr;
            return data;
        }
    }
    else {
        return static_cast<T>(hclib_future_get(this));
    }
}

template<typename T>
void future_t<T>::add_future_vector() {
    auto task_local = static_cast<checkpoint_task_params_t<T>*>(*hclib_get_curr_task_local());
    assert(is_checkpoint_task(task_local));
    task_local->rel_vec->push_back(this);
}

/*
put operation only adds data to a temporary storage.
Also adds the promise to the vector which is
later used to perform the actual put
*/
template<typename T>
void promise_t<T*>::put(T* datum) {
    auto task_local = static_cast<checkpoint_task_params_t<T*>*>(*hclib_get_curr_task_local());

    //if called from replay  task, save to a temp and delay the put
    //also add it to the vector that captures all put opertations
    if(is_checkpoint_task(task_local)){
        tmp_data = datum;
        //TODO: assuming the same promise will be used in all replays
        if(task_local->index == 0)
            task_local->put_vec->push_back(this);
    }
    //if called from non-checkpoint task, just perform the put
    else {
        hclib::promise_t<T*>::put(datum);
    }
}

#endif // USE_RESILIENT_PROMISE

inline int get_index() {
  auto task_local = static_cast<checkpoint_task_params_t<void*>*>(*hclib_get_curr_task_local());
  assert(is_checkpoint_task(task_local));
  return task_local->index;
}

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    //fetch the checkpoint_task_parameters from task local and pass is to the async
    auto ctp = *hclib_get_curr_task_local();

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);
    hclib::async_await( [=] () {
        *(hclib_get_curr_task_local()) = ctp;
        (*lambda_ptr)();
        delete lambda_ptr;
        auto task_local = static_cast<checkpoint_task_params_t<void*>*>(*hclib_get_curr_task_local());
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

    //fetch the checkpoint_task_parameters from task local and pass is to the async
    auto ctp = *hclib_get_curr_task_local();

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);
    hclib::async_await( [=] () {
        *(hclib_get_curr_task_local()) = ctp;
        (*lambda_ptr)();
        delete lambda_ptr;
        auto task_local = static_cast<checkpoint_task_params_t<void*>*>(*hclib_get_curr_task_local());
        //TODO: assuming the same future will be used in all replays
        if(task_local->index == 0) {
            for(auto && elem: *(futures))
                static_cast<future_t<void*>*>(elem)->add_future_vector();
        }
    }, futures);
}


template <typename T>
inline void async_await_sync(T* lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    (*lambda)();
    auto task_local = static_cast<checkpoint_task_params_t<void*>*>(*hclib_get_curr_task_local());
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
}

template <typename T>
inline void async_await_sync(T* lambda, std::vector<hclib_future_t *> *futures) {

    (*lambda)();
    auto task_local = static_cast<checkpoint_task_params_t<void*>*>(*hclib_get_curr_task_local());

    //TODO: assuming the same future will be used in all replays
    if(task_local->index == 0) {
        for(auto && elem: *(futures))
            static_cast<future_t<void*>*>(elem)->add_future_vector();
    }
}

template <int LEAF=1, int N=1, typename T>
std::enable_if_t< LEAF==1, void>
async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await([=]() {
        auto ctp = new checkpoint_task_params_t<void*>();;
        ctp->put_vec = new promise_vector<void*>();
        ctp->rel_vec = new future_vector<void*>();
        ctp->checkpoint_run = 0;
        bool result = false;
        int index = 0;

        assert(*(hclib_get_curr_task_local()) ==  nullptr);
        for(int i=0; i<N; i++) {
            ctp->index = i;
            *(hclib_get_curr_task_local()) = ctp;
            async_await_sync(lambda_ptr, f1, f2, f3, f4);
            assert(ctp->checkpoint_run == 0);

            //ran successfully without errors
            if(error_check_fn(params) == 1)
            {
                result = true;
                index = ctp->index;
                break;
            }

            auto data = ctp->put_vec->data();
            auto size = ctp->put_vec->size();
            for(int i=0; i<size; i++)
                data[i]->tmp_delete();
        }

        //ran with errors and therefore access data from checkpoint if available
        if(result == false){
            assert(ctp->checkpoint_run == 0);
            ctp->index = N;
            ctp->checkpoint_run = 1;
            index = N;
            async_await_sync(lambda_ptr, f1, f2, f3, f4);
            result = error_check_fn(params);
        }

        delete lambda_ptr;

        *(hclib_get_curr_task_local()) = nullptr;

        if(result) {
            ctp->put_vec->do_puts(index);
            ctp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete ctp->put_vec;
        delete ctp->rel_vec;
        delete ctp;;
    }, f1, f2, f3, f4);
}

} // namespace checkpoint
} // namespace resilience
} // namespace hclib

#endif
