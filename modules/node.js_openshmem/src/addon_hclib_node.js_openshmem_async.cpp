
#include "hclib_node.js_openshmem.h"
#include "addon_hclib_node.js_openshmem.h"
#include "util.h"
#include <utility>
#include <unordered_map>
#include <deque>
#include <shmem.h>
#include <uv.h>

#define USE_MIX_SYNC_ASYNC_STYLE

namespace hclib {

using promise_long = std::pair<Napi::Promise::Deferred *, long>;
using promise_double = std::pair<Napi::Promise::Deferred *, double>;

safe_vector<promise_long> prom_vec_long;
safe_vector<promise_double> prom_vec_double;
safe_vector<Napi::Promise::Deferred *> prom_vec_lock_promise;
safe_vector<Napi::FunctionReference *> prom_vec_lock_callback;

std::unordered_map<long*, long> lock_count;
std::unordered_map<long*, std::deque<Napi::Promise::Deferred *>*> lock_queue;
std::unordered_map<long*, std::deque<Napi::FunctionReference *>*> lock_queue_callback;
std::mutex lock_queue_mtx;
#ifdef USE_MIX_SYNC_ASYNC_STYLE
std::unordered_map<long*, bool> lock_is_async;
#endif

uv_async_t async_var;
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
        status = napi_create_string_utf8(env, "threaded_promise_resolve", NAPI_AUTO_LENGTH, &resource_name);
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

void clear_callback_fn_helper(Napi::Env env) {
    Napi::Number napi_one = Napi::Number::New(env,1);
    std::lock_guard<std::mutex> lkg(prom_vec_lock_callback.mtx);
    for(auto && elem: prom_vec_lock_callback) {
        elem->Call({napi_one});
        delete elem;
    }
    prom_vec_lock_callback.clear();
}

void clear_put_promises_fn_helper(Napi::Env env) {
    {
      std::lock_guard<std::mutex> lkg(prom_vec_long.mtx);
      for(auto && elem: prom_vec_long) {
        elem.first->Resolve( Napi::Number::New(env,elem.second) );
        delete elem.first;
      }
      prom_vec_long.clear();
    }
    {
      std::lock_guard<std::mutex> lkg(prom_vec_double.mtx);
      for(auto && elem: prom_vec_double) {
        elem.first->Resolve( Napi::Number::New(env,elem.second) );
        delete elem.first;
      }
      prom_vec_double.clear();
    }
    {
      Napi::Number napi_one = Napi::Number::New(env,1);
      std::lock_guard<std::mutex> lkg(prom_vec_lock_promise.mtx);
      for(auto && elem: prom_vec_lock_promise) {
        elem->Resolve(napi_one);
        delete elem;
      }
      prom_vec_lock_promise.clear();
    }
}

void async_promise_resolve(uv_async_t* handle) {
    Napi::Env env = env_ref.Env();
    //Napi::HandleScope scope(env);
    CallbackScope scope(env);
    clear_callback_fn_helper(env);
    clear_put_promises_fn_helper(env);
}

void init_async_helper(Napi::Env env) {
    env_ref = Napi::Persistent(Napi::Object::New(env));
    env_ref.SuppressDestruct();
    uv_async_init(uv_default_loop(), &async_var, async_promise_resolve);
}

void finalize_async_helper(Napi::Env env) {
    uv_close((uv_handle_t*) &async_var, nullptr);
}

Napi::Number clear_put_promises_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    clear_callback_fn_helper(env);
    clear_put_promises_fn_helper(env);
    return Napi::Number::New(env, 1);
}

Napi::Promise shmem_long_g_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *src = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        long val = shmem_long_g(src, pe);
        prom_vec_long.push_back(std::make_pair(prom_ptr, val));
        uv_async_send(&async_var);
    }, nic);

    return prom_ptr->Promise();
}

Napi::Promise shmem_double_g_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *src = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        double val = shmem_double_g(src, pe);
        prom_vec_double.push_back(std::make_pair(prom_ptr, val));
        uv_async_send(&async_var);
    }, nic);

    return prom_ptr->Promise();
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
            prom_vec_lock_promise.push_back((Napi::Promise::Deferred*)prom_ptr);
            uv_async_send(&async_var);
        }, nic);
#ifdef USE_MIX_SYNC_ASYNC_STYLE
        lock_is_async[ptr] = true;
#endif
    }
    else {
        std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue[ptr];
        if(q_ptr == nullptr) {
            q_ptr = new std::deque<Napi::Promise::Deferred*>();
            lock_queue[ptr] = q_ptr;
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
            prom_vec_lock_callback.push_back((Napi::FunctionReference*)cbr);
            uv_async_send(&async_var);
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
        prom_vec_lock_promise.push_back((Napi::Promise::Deferred*)prom_ptr);
        uv_async_send(&async_var);

        //if there is additional set_lock_async calls, then they will be queued.
        //invoke set_lock if there are additional calls made from JS.
        lock_queue_mtx.lock();
        const long count = lock_count[ptr];
        lock_count[ptr] = count-1;
        std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue[ptr];

        if(count > 1) {
            assert(q_ptr != nullptr && !(q_ptr->empty()));
            Napi::Promise::Deferred* prom_ptr_next = q_ptr->front();
            q_ptr->pop_front();
            lock_queue_mtx.unlock();

            shmem_set_lock(ptr);
            prom_vec_lock_promise.push_back((Napi::Promise::Deferred*)prom_ptr_next);
            uv_async_send(&async_var);
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
        std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue[ptr];

        if(count > 1) {
            assert(q_ptr != nullptr && !(q_ptr->empty()));
            Napi::Promise::Deferred* prom_ptr_next = q_ptr->front();
            q_ptr->pop_front();
            lock_queue_mtx.unlock();

            shmem_set_lock(ptr);
            prom_vec_lock_promise.push_back((Napi::Promise::Deferred*)prom_ptr_next);
            uv_async_send(&async_var);
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
        std::deque<Napi::Promise::Deferred*>* q_ptr = lock_queue[ptr];
        assert(q_ptr != nullptr && !(q_ptr->empty()));
        Napi::Promise::Deferred* prom_ptr = q_ptr->front();
        q_ptr->pop_front();
        prom_ptr->Resolve(ret);
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
            prom_vec_lock_callback.push_back((Napi::FunctionReference*)callback_next);
            uv_async_send(&async_var);
        }
        else {
            lock_queue_mtx.unlock();
            assert(count == 1);
        }
    }, nic);
    return Napi::Value();
}

Napi::Value Init_async(Napi::Env env, Napi::Object exports) {

    exports.Set(Napi::String::New(env, "clear_put_promises"),
            Napi::Function::New(env, clear_put_promises_fn));
    exports.Set(Napi::String::New(env, "long_g_async"),
            Napi::Function::New(env, shmem_long_g_async_fn));
    exports.Set(Napi::String::New(env, "double_g_async"),
            Napi::Function::New(env, shmem_double_g_async_fn));
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
    return Napi::Number::New(env, 1);

}

}
