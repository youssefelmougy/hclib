#define NAPI_EXPERIMENTAL
#include "hclib_node.js_openshmem.h"
#include "addon_hclib_node.js_openshmem.h"
#include "util.h"
#include <utility>
#include <unordered_map>
#include <deque>
#include <shmem.h>
#include <uv.h>

//#define USE_MIX_SYNC_ASYNC_STYLE
//#define USE_UV_LOOP

namespace hclib {

template<typename T>
using promise_data = std::pair<Napi::Promise::Deferred *, T>;

template<typename T>
using callback_data = std::pair<Napi::FunctionReference *, T>;

std::unordered_map<long*, long> lock_count;
std::unordered_map<long*, std::deque<Napi::Promise::Deferred *>*> lock_queue_promise;
std::unordered_map<long*, std::deque<Napi::FunctionReference *>*> lock_queue_callback;
std::mutex lock_queue_mtx;
#ifdef USE_MIX_SYNC_ASYNC_STYLE
std::unordered_map<long*, bool> lock_is_async;
#endif

#ifndef USE_UV_LOOP
napi_threadsafe_function long_promise_ts_fn, long_callback_ts_fn;
napi_threadsafe_function double_promise_ts_fn, double_callback_ts_fn;
napi_threadsafe_function lock_promise_ts_fn, lock_callback_ts_fn;

template<typename T>
void promise_resolve_fn(napi_env env, napi_value js_callback, void* ctx, void* data) {
    auto elem = (promise_data<T>*)data;
    elem->first->Resolve( Napi::Number::New(env,elem->second) );
    delete elem->first;
    delete elem;
}

template<typename T>
void callback_fn(napi_env env, napi_value js_callback, void* ctx, void* data) {
    auto elem = (callback_data<T>*)data;
    elem->first->Call({Napi::Number::New(env,elem->second)});
    elem->first->Reset();
    delete elem;
}

void promise_resolve_lock_fn(napi_env env, napi_value js_callback, void* ctx, void* data) {
    auto elem = (Napi::Promise::Deferred *)data;
    napi_value undefined;
    napi_get_undefined(env, &undefined);
    elem->Resolve(undefined);
    delete elem;
}

void callback_lock_fn(napi_env env, napi_value js_callback, void* ctx, void* data) {
    auto elem = (Napi::FunctionReference *)data;
    napi_value undefined;
    napi_get_undefined(env, &undefined);
    elem->Call({undefined});
    delete elem;
}

void init_async_helper(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::Function cb_blank = info[0].As<Napi::Function>();
    napi_value async_name_blank;

    napi_status status;
    status = napi_create_string_utf8(env, "Blank Callback", NAPI_AUTO_LENGTH, &async_name_blank);
    status = napi_create_threadsafe_function(env, cb_blank, nullptr, async_name_blank,
            0, 1, nullptr, nullptr, nullptr, promise_resolve_fn<long>, &long_promise_ts_fn);

    status = napi_create_threadsafe_function(env, cb_blank, nullptr, async_name_blank,
            0, 1, nullptr, nullptr, nullptr, promise_resolve_fn<double>, &double_promise_ts_fn);

    status = napi_create_threadsafe_function(env, cb_blank, nullptr, async_name_blank,
            0, 1, nullptr, nullptr, nullptr, callback_fn<long>, &long_callback_ts_fn);

    status = napi_create_threadsafe_function(env, cb_blank, nullptr, async_name_blank,
            0, 1, nullptr, nullptr, nullptr, callback_fn<double>, &double_callback_ts_fn);

    status = napi_create_threadsafe_function(env, cb_blank, nullptr, async_name_blank,
            0, 1, nullptr, nullptr, nullptr, promise_resolve_lock_fn, &lock_promise_ts_fn);

    status = napi_create_threadsafe_function(env, cb_blank, nullptr, async_name_blank,
            0, 1, nullptr, nullptr, nullptr, callback_lock_fn, &lock_callback_ts_fn);
}

void finalize_async_helper(const Napi::CallbackInfo& info) {
    napi_release_threadsafe_function(long_promise_ts_fn, napi_tsfn_release);
    napi_release_threadsafe_function(double_promise_ts_fn, napi_tsfn_release);
    napi_release_threadsafe_function(long_callback_ts_fn, napi_tsfn_release);
    napi_release_threadsafe_function(double_callback_ts_fn, napi_tsfn_release);
    napi_release_threadsafe_function(lock_promise_ts_fn, napi_tsfn_release);
    napi_release_threadsafe_function(lock_callback_ts_fn, napi_tsfn_release);
}

#else // USE_UV_LOOP

safe_vector<promise_data<long>> prom_vec_long;
safe_vector<promise_data<double>> prom_vec_double;
safe_vector<callback_data<long>> callback_vec_long;
safe_vector<callback_data<double>> callback_vec_double;
safe_vector<Napi::Promise::Deferred *> prom_vec_lock;
safe_vector<Napi::FunctionReference *> callback_vec_lock;

template<typename T>
inline safe_vector<promise_data<T>>* get_prom_vec() {assert(false); return &prom_vec_long; }

template<>
inline safe_vector<promise_data<long>>* get_prom_vec() { return &prom_vec_long; }

template<>
inline safe_vector<promise_data<double>>* get_prom_vec() { return &prom_vec_double; }

template<typename T>
inline safe_vector<callback_data<T>>* get_callback_vec() {assert(false); return &callback_vec_long; }

template<>
inline safe_vector<callback_data<long>>* get_callback_vec() { return &callback_vec_long; }

template<>
inline safe_vector<callback_data<double>>* get_callback_vec() { return &callback_vec_double; }

uv_async_t async_var_prom_long, async_var_prom_double;
uv_async_t async_var_callback_long, async_var_callback_double;
uv_async_t async_var_prom_lock, async_var_callback_lock;
Napi::ObjectReference env_ref;

//From https://github.com/nodejs/node-addon-api/issues/223
struct CallbackScope {
    napi_env env;
    napi_handle_scope handle_scope;
    napi_callback_scope callback_scope;

    CallbackScope(napi_env env) : env(env) {
        napi_status status;
        napi_async_context context;
        napi_value resource_name, resource_object;

        status = napi_open_handle_scope(env, &handle_scope);
        status = napi_create_string_utf8(env, "threaded_callback_resolve", NAPI_AUTO_LENGTH, &resource_name);
        status = napi_create_object(env, &resource_object);
        status = napi_async_init(env, resource_object, resource_name, &context);
        status = napi_open_callback_scope(env, resource_object, context, &callback_scope);
    }

    ~CallbackScope() {
        napi_status status;
        status = napi_close_callback_scope(env, callback_scope);
        status = napi_close_handle_scope(env, handle_scope);
    }
};

template<typename T>
void promise_resolve_fn_helper(Napi::Env env) {
    auto prom_vec = get_prom_vec<T>();
    std::lock_guard<std::mutex> lkg(prom_vec->mtx);
    for(auto && elem:*prom_vec) {
        elem.first->Resolve(Napi::Number::New(env,elem.second));
        delete elem.first;
    }
    prom_vec->clear();
}

template<typename T>
void promise_resolve_fn(uv_async_t* handle) {
    Napi::Env env = env_ref.Env();
    //Napi::HandleScope scope(env);
    CallbackScope scope(env);
    promise_resolve_fn_helper<T>(env);
}

template<typename T>
void callback_fn_helper(Napi::Env env) {
    auto callback_vec = get_callback_vec<T>();
    std::lock_guard<std::mutex> lkg(callback_vec->mtx);
    for(auto && elem:*callback_vec) {
        elem.first->Call({Napi::Number::New(env,elem.second)});
        delete elem.first;
    }
    callback_vec->clear();
}

template<typename T>
void callback_fn(uv_async_t* handle) {
    Napi::Env env = env_ref.Env();
    Napi::HandleScope scope(env);
    //CallbackScope scope(env);
    callback_fn_helper<T>(env);
}

void promise_resolve_lock_fn_helper(Napi::Env env) {
    //Napi::Value undefined = env.Undefined();
    std::lock_guard<std::mutex> lkg(prom_vec_lock.mtx);
    for(auto && elem: prom_vec_lock) {
        elem->Resolve(Napi::Number::New(env,1));
        delete elem;
    }
    prom_vec_lock.clear();
}

void promise_resolve_lock_fn(uv_async_t* handle) {
    Napi::Env env = env_ref.Env();
    //Napi::HandleScope scope(env);
    CallbackScope scope(env);
    promise_resolve_lock_fn_helper(env);
}

void callback_lock_fn_helper(Napi::Env env) {
    Napi::Value undefined = env.Undefined();
    std::lock_guard<std::mutex> lkg(callback_vec_lock.mtx);
    for(auto && elem: callback_vec_lock) {
        elem->Call({undefined});
        delete elem;
    }
    callback_vec_lock.clear();
}

void callback_lock_fn(uv_async_t* handle) {
    Napi::Env env = env_ref.Env();
    Napi::HandleScope scope(env);
    //CallbackScope scope(env);
    callback_lock_fn_helper(env);
}

void init_async_helper(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    env_ref = Napi::Persistent(Napi::Object::New(env));
    env_ref.SuppressDestruct();
    uv_async_init(uv_default_loop(), &async_var_prom_long, promise_resolve_fn<long>);
    uv_async_init(uv_default_loop(), &async_var_prom_double, promise_resolve_fn<double>);
    uv_async_init(uv_default_loop(), &async_var_callback_long, callback_fn<long>);
    uv_async_init(uv_default_loop(), &async_var_callback_double, callback_fn<double>);
    uv_async_init(uv_default_loop(), &async_var_prom_lock, promise_resolve_lock_fn);
    uv_async_init(uv_default_loop(), &async_var_callback_lock, callback_lock_fn);
}

void finalize_async_helper(const Napi::CallbackInfo& info) {
    uv_close((uv_handle_t*) &async_var_prom_long, nullptr);
    uv_close((uv_handle_t*) &async_var_prom_double, nullptr);
    uv_close((uv_handle_t*) &async_var_callback_long, nullptr);
    uv_close((uv_handle_t*) &async_var_callback_double, nullptr);
    uv_close((uv_handle_t*) &async_var_prom_lock, nullptr);
    uv_close((uv_handle_t*) &async_var_callback_lock, nullptr);
}

Napi::Value clear_put_promises_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    promise_resolve_fn_helper<long>(env);
    promise_resolve_fn_helper<double>(env);
    callback_fn_helper<long>(env);
    callback_fn_helper<double>(env);
    promise_resolve_lock_fn_helper(env);
    callback_lock_fn_helper(env);
    return Napi::Value();
}

#endif // USE_UV_LOOP

Napi::Promise shmem_long_g_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *src = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        long val = shmem_long_g(src, pe);
#ifndef USE_UV_LOOP
        auto ret = new promise_data<long>(prom_ptr, val);
        const napi_status status = napi_call_threadsafe_function(long_promise_ts_fn, ret, napi_tsfn_nonblocking);
        assert(status == napi_ok);
#else
        prom_vec_long.push_back(std::make_pair(prom_ptr, val));
        uv_async_send(&async_var_prom_long);
#endif
    }, nic);

    return prom_ptr->Promise();
}

Napi::Value shmem_long_g_async_callback_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *src = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    auto cb_ptr = new Napi::FunctionReference(Napi::Persistent(info[2].As<Napi::Function>()));

    hclib::async_nb_at([=] () {
        long val = shmem_long_g(src, pe);
#ifndef USE_UV_LOOP
        auto ret = new callback_data<long>(cb_ptr, val);
        const napi_status status = napi_call_threadsafe_function(long_callback_ts_fn, ret, napi_tsfn_nonblocking);
        assert(status == napi_ok);
#else
        callback_vec_long.push_back(std::make_pair(cb_ptr, val));
        uv_async_send(&async_var_callback_long);
#endif
    }, nic);

    return Napi::Value();
}

Napi::Promise shmem_double_g_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *src = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        double val = shmem_double_g(src, pe);
#ifndef USE_UV_LOOP
        auto ret = new promise_data<double>(prom_ptr, val);
        const napi_status status = napi_call_threadsafe_function(double_promise_ts_fn, ret, napi_tsfn_nonblocking);
        assert(status == napi_ok);
#else
        prom_vec_double.push_back(std::make_pair(prom_ptr, val));
        uv_async_send(&async_var_prom_double);
#endif
    }, nic);

    return prom_ptr->Promise();
}

Napi::Value shmem_double_g_async_callback_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *src = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    auto cb_ptr = new Napi::FunctionReference(Napi::Persistent(info[2].As<Napi::Function>()));

    hclib::async_nb_at([=] () {
        double val = shmem_double_g(src, pe);
#ifndef USE_UV_LOOP
        auto ret = new callback_data<double>(cb_ptr, val);
        const napi_status status = napi_call_threadsafe_function(double_callback_ts_fn, ret, napi_tsfn_nonblocking);
        assert(status == napi_ok);
#else
        callback_vec_double.push_back(std::make_pair(cb_ptr, val));
        uv_async_send(&async_var_callback_double);
#endif
    }, nic);

    return Napi::Value();
}

Napi::Promise shmem_set_lock_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    lock_queue_mtx.lock();
    const long count = lock_count[ptr];
    lock_count[ptr] = count+1;
    assert(count >= 0);

    if(count == 0) {
        lock_queue_mtx.unlock();
        hclib::async_nb_at([=] () {
            shmem_set_lock(ptr);
#ifndef USE_UV_LOOP
            const napi_status status = napi_call_threadsafe_function(lock_promise_ts_fn, prom_ptr, napi_tsfn_nonblocking);
            assert(status == napi_ok);
#else
            prom_vec_lock.push_back((Napi::Promise::Deferred*)prom_ptr);
            uv_async_send(&async_var_prom_lock);
#endif
        }, nic);
#ifdef USE_MIX_SYNC_ASYNC_STYLE
        lock_is_async[ptr] = true;
#endif
    }
    else {
        std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue_promise[ptr];
        if(q_ptr == nullptr) {
            q_ptr = new std::deque<Napi::Promise::Deferred*>();
            lock_queue_promise[ptr] = q_ptr;
        }
        q_ptr->push_back(prom_ptr);
        lock_queue_mtx.unlock();
    }
    return prom_ptr->Promise();
}

Napi::Value shmem_set_lock_async_callback_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    auto cbr = new Napi::FunctionReference(Napi::Persistent(info[1].As<Napi::Function>()));

    lock_queue_mtx.lock();
    const long count = lock_count[ptr];
    lock_count[ptr] = count+1;
    assert(count >= 0);

    if(count == 0) {
        lock_queue_mtx.unlock();
        hclib::async_nb_at([=] () {
            shmem_set_lock(ptr);
#ifndef USE_UV_LOOP
            const napi_status status = napi_call_threadsafe_function(lock_callback_ts_fn, cbr, napi_tsfn_nonblocking);
            assert(status == napi_ok);
#else
            callback_vec_lock.push_back((Napi::FunctionReference*)cbr);
            uv_async_send(&async_var_callback_lock);
#endif
        }, nic);
    }
    else {
        auto q_ptr = lock_queue_callback[ptr];
        if(q_ptr == nullptr) {
            q_ptr = new std::deque<Napi::FunctionReference*>();
            lock_queue_callback[ptr] = q_ptr;
        }
        q_ptr->push_back(cbr);
        lock_queue_mtx.unlock();
    }
    return Napi::Value();
}

Napi::Promise shmem_clear_lock_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        shmem_clear_lock(ptr);
#ifndef USE_UV_LOOP
        const napi_status status = napi_call_threadsafe_function(lock_promise_ts_fn, prom_ptr, napi_tsfn_nonblocking);
        assert(status == napi_ok);
#else
        prom_vec_lock.push_back((Napi::Promise::Deferred*)prom_ptr);
        uv_async_send(&async_var_prom_lock);
#endif

        //if there is additional set_lock_async calls, then they will be queued.
        //invoke set_lock if there are additional calls made from JS.
        lock_queue_mtx.lock();
        const long count = lock_count[ptr];
        lock_count[ptr] = count-1;

        if(count > 1) {
            std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue_promise[ptr];
            assert(q_ptr != nullptr && !(q_ptr->empty()));
            Napi::Promise::Deferred* prom_ptr_next = q_ptr->front();
            q_ptr->pop_front();
            lock_queue_mtx.unlock();

            shmem_set_lock(ptr);
#ifndef USE_UV_LOOP
            const napi_status status = napi_call_threadsafe_function(lock_promise_ts_fn, prom_ptr_next, napi_tsfn_nonblocking);
            assert(status == napi_ok);
#else
            prom_vec_lock.push_back((Napi::Promise::Deferred*)prom_ptr_next);
            uv_async_send(&async_var_prom_lock);
#endif
        }
        else {
            lock_queue_mtx.unlock();
            assert(count == 1);
        }
    }, nic);
    return prom_ptr->Promise();
}

Napi::Value shmem_clear_lock_async_noreturn_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::Number ret = info[0].As<Napi::Number>();
    long *ptr = (long*) ret.Int64Value();

    hclib::async_nb_at([=] () {
        shmem_clear_lock(ptr);

        //if there is additional set_lock_async calls, then they will be queued.
        //invoke set_lock if there are additional calls made from JS.
        lock_queue_mtx.lock();
        const long count = lock_count[ptr];
        lock_count[ptr] = count-1;

        if(count > 1) {
            std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue_promise[ptr];
            assert(q_ptr != nullptr && !(q_ptr->empty()));
            Napi::Promise::Deferred* prom_ptr_next = q_ptr->front();
            q_ptr->pop_front();
            lock_queue_mtx.unlock();

            shmem_set_lock(ptr);
#ifndef USE_UV_LOOP
            const napi_status status = napi_call_threadsafe_function(lock_promise_ts_fn, prom_ptr_next, napi_tsfn_nonblocking);
            assert(status == napi_ok);
#else
            prom_vec_lock.push_back((Napi::Promise::Deferred*)prom_ptr_next);
            uv_async_send(&async_var_prom_lock);
#endif
        }
        else {
            lock_queue_mtx.unlock();
            assert(count == 1);
        }
    }, nic);
    return Napi::Value();
}

Napi::Value shmem_clear_lock_async_finish_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::Number ret = info[0].As<Napi::Number>();
    long *ptr = (long*) ret.Int64Value();

    const long count = lock_count[ptr];
    lock_count[ptr] = count-1;

#ifndef USE_MIX_SYNC_ASYNC_STYLE
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_clear_lock(ptr);
            if(count > 1) {
                shmem_set_lock(ptr);
            }
        }, nic);
    });
#else
    bool is_async = lock_is_async[ptr];
    if(is_async) {
        hclib::finish([=] () {
            hclib::async_nb_at([=] () {
                shmem_clear_lock(ptr);
            }, nic);
        });
        lock_is_async[ptr] = false;
    }
    else {
        shmem_clear_lock(ptr);
    }
    if(count > 1)
        shmem_set_lock(ptr);
#endif

    if(count > 1) {
        std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue_promise[ptr];
        assert(q_ptr != nullptr && !(q_ptr->empty()));
        Napi::Promise::Deferred* prom_ptr = q_ptr->front();
        q_ptr->pop_front();
        prom_ptr->Resolve(env.Undefined());
    }
    else {
        assert(count == 1);
    }
    return Napi::Value();
}

Napi::Value shmem_clear_lock_async_callback_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::Number ret = info[0].As<Napi::Number>();
    long *ptr = (long*) ret.Int64Value();

    hclib::async_nb_at([=] () {
        shmem_clear_lock(ptr);

        //if there is additional set_lock_async calls, then they will be queued.
        //invoke set_lock if there are additional calls made from JS.
        lock_queue_mtx.lock();
        const long count = lock_count[ptr];
        lock_count[ptr] = count-1;
        auto q_ptr = lock_queue_callback[ptr];

        if(count > 1) {
            assert(q_ptr != nullptr && !(q_ptr->empty()));
            auto callback_next = q_ptr->front();
            q_ptr->pop_front();
            lock_queue_mtx.unlock();

            shmem_set_lock(ptr);
#ifndef USE_UV_LOOP
            const napi_status status = napi_call_threadsafe_function(lock_callback_ts_fn, callback_next, napi_tsfn_nonblocking);
            assert(status == napi_ok);
#else
            callback_vec_lock.push_back((Napi::FunctionReference*)callback_next);
            uv_async_send(&async_var_callback_lock);
#endif
        }
        else {
            lock_queue_mtx.unlock();
            assert(count == 1);
        }
    }, nic);
    return Napi::Value();
}

Napi::Value Init_async(Napi::Env env, Napi::Object exports) {

    exports.Set(Napi::String::New(env, "long_g_async"),
            Napi::Function::New(env, shmem_long_g_async_fn));
    exports.Set(Napi::String::New(env, "long_g_async_callback"),
            Napi::Function::New(env, shmem_long_g_async_callback_fn));
    exports.Set(Napi::String::New(env, "double_g_async"),
            Napi::Function::New(env, shmem_double_g_async_fn));
    exports.Set(Napi::String::New(env, "double_g_async_callback"),
            Napi::Function::New(env, shmem_double_g_async_callback_fn));
    exports.Set(Napi::String::New(env, "set_lock_async"),
            Napi::Function::New(env, shmem_set_lock_async_fn));
    exports.Set(Napi::String::New(env, "set_lock_async_callback"),
            Napi::Function::New(env, shmem_set_lock_async_callback_fn));
    exports.Set(Napi::String::New(env, "clear_lock_async"),
            Napi::Function::New(env, shmem_clear_lock_async_fn));
    exports.Set(Napi::String::New(env, "clear_lock_async_noreturn"),
            Napi::Function::New(env, shmem_clear_lock_async_noreturn_fn));
    exports.Set(Napi::String::New(env, "clear_lock_async_finish"),
            Napi::Function::New(env, shmem_clear_lock_async_finish_fn));
    exports.Set(Napi::String::New(env, "clear_lock_async_callback"),
            Napi::Function::New(env, shmem_clear_lock_async_callback_fn));

#ifdef USE_UV_LOOP
    exports.Set(Napi::String::New(env, "clear_put_promises"),
            Napi::Function::New(env, clear_put_promises_fn));
#endif
    return Napi::Number::New(env, 1);

}

}
