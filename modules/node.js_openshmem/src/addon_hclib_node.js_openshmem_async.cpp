
#include "hclib_node.js_openshmem.h"
#include "addon_hclib_node.js_openshmem.h"
#include "util.h"
#include <utility>
#include <shmem.h>

namespace hclib {

using promise_long = std::pair<Napi::Promise::Deferred *, long>;
using promise_double = std::pair<Napi::Promise::Deferred *, double>;

safe_vector<promise_long> prom_vec_long;
safe_vector<promise_double> prom_vec_double;

Napi::Number clear_put_promises_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
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
    return Napi::Number::New(env, 1);
}

Napi::Promise shmem_long_g_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *src = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    long *val = (long*)malloc(sizeof(long));
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        *val = shmem_long_g(src, pe);
        prom_vec_long.push_back(std::make_pair(prom_ptr, *val));
        free(val);
    }, nic);

    return prom_ptr->Promise();
}

Napi::Promise shmem_double_g_async_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *src = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    double *val = (double*)malloc(sizeof(double));
    Napi::Promise::Deferred * prom_ptr = new Napi::Promise::Deferred(env);

    hclib::async_nb_at([=] () {
        *val = shmem_double_g(src, pe);
        prom_vec_double.push_back(std::make_pair(prom_ptr, *val));
        free(val);
    }, nic);

    return prom_ptr->Promise();
}

Napi::Value Init_async(Napi::Env env, Napi::Object exports) {

    exports.Set(Napi::String::New(env, "clear_put_promises"),
            Napi::Function::New(env, clear_put_promises_fn));
    exports.Set(Napi::String::New(env, "long_g_async"),
            Napi::Function::New(env, shmem_long_g_async_fn));
    exports.Set(Napi::String::New(env, "double_g_async"),
            Napi::Function::New(env, shmem_double_g_async_fn));
    return Napi::Number::New(env, 1);

}

}
