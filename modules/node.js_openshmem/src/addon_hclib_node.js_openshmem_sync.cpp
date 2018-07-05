#include "hclib_node.js_openshmem.h"
#include "addon_hclib_node.js_openshmem.h"
#include <shmem.h>

namespace hclib {

Napi::Number shmem_malloc_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    int64_t size = info[0].As<Napi::Number>().Int64Value();
    void **out_alloc = (void **)malloc(sizeof(void *));
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            *out_alloc = shmem_malloc(size);
        }, nic);
    });
    void *ptr = *out_alloc;
    free(out_alloc);
    return Napi::Number::New(env, (uintptr_t)ptr);
}

Napi::Number shmem_free_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *ptr = (void*) info[0].As<Napi::Number>().Int64Value();
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_free(ptr);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_long_g_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *src = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    long *val = (long*)malloc(sizeof(long));
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            *val = shmem_long_g(src, pe);
        }, nic);
    });
    long val_tmp = *val;
    free(val);
    return Napi::Number::New(env, val_tmp);
}

Napi::Number shmem_long_p_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *dest = (long*) info[0].As<Napi::Number>().Int64Value();
    long val = info[1].As<Napi::Number>().Int64Value();
    int64_t pe = info[2].As<Napi::Number>().Int64Value();
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_long_p(dest, val, pe);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_double_g_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *src = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    double *val = (double*)malloc(sizeof(long));
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            *val = shmem_double_g(src, pe);
        }, nic);
    });
    double val_tmp = *val;
    free(val);
    return Napi::Number::New(env, val_tmp);
}

Napi::Number shmem_double_p_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *dest = (double*) info[0].As<Napi::Number>().Int64Value();
    double val = info[1].As<Napi::Number>().DoubleValue();
    int64_t pe = info[2].As<Napi::Number>().Int64Value();
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_double_p(dest, val, pe);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_barrier_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    hclib::finish([] () {
        hclib::async_nb_at([] () {
            shmem_barrier_all();
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_set_lock_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

    //TODO: only one task per rank can lock as of now. otherwise it will cause deadlock.
    //Use promise/future list as done in openshmem module.

    void *ptr = (void*) info[0].As<Napi::Number>().Int64Value();
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_set_lock((long*)ptr);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_clear_lock_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *ptr = (void*) info[0].As<Napi::Number>().Int64Value();
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_clear_lock((long*)ptr);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_broadcast64_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *dest = (void*) info[0].As<Napi::Number>().Int64Value();
    void *src = (void*) info[1].As<Napi::Number>().Int64Value();
    size_t nelems = info[2].As<Napi::Number>().Int64Value();
    int PE_root = info[3].As<Napi::Number>().Int64Value();
    int PE_start = info[4].As<Napi::Number>().Int64Value();
    int logPE_stride = info[5].As<Napi::Number>().Int64Value();
    int PE_size = info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
          shmem_broadcast64(dest, src, nelems, PE_root, PE_start, logPE_stride, PE_size, pSync);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_long_sum_to_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *dest = (long*) info[0].As<Napi::Number>().Int64Value();
    long *src  = (long*) info[1].As<Napi::Number>().Int64Value();
    int nreduce = info[2].As<Napi::Number>().Int64Value();
    int PE_start = info[3].As<Napi::Number>().Int64Value();
    int logPE_stride = info[4].As<Napi::Number>().Int64Value();
    int PE_size = info[5].As<Napi::Number>().Int64Value();
    long *pWrk = (long*) info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_long_sum_to_all(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_double_sum_to_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *dest = (double*) info[0].As<Napi::Number>().Int64Value();
    double *src  = (double*) info[1].As<Napi::Number>().Int64Value();
    int nreduce = info[2].As<Napi::Number>().Int64Value();
    int PE_start = info[3].As<Napi::Number>().Int64Value();
    int logPE_stride = info[4].As<Napi::Number>().Int64Value();
    int PE_size = info[5].As<Napi::Number>().Int64Value();
    double *pWrk = (double*) info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_double_sum_to_all(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Value Init_sync(Napi::Env env, Napi::Object exports) {

  exports.Set(Napi::String::New(env, "malloc_sync"),
              Napi::Function::New(env, shmem_malloc_sync_fn));
  exports.Set(Napi::String::New(env, "free_sync"),
              Napi::Function::New(env, shmem_free_sync_fn));
  exports.Set(Napi::String::New(env, "long_g_sync"),
              Napi::Function::New(env, shmem_long_g_sync_fn));
  exports.Set(Napi::String::New(env, "long_p_sync"),
              Napi::Function::New(env, shmem_long_p_sync_fn));
  exports.Set(Napi::String::New(env, "double_g_sync"),
              Napi::Function::New(env, shmem_double_g_sync_fn));
  exports.Set(Napi::String::New(env, "double_p_sync"),
              Napi::Function::New(env, shmem_double_p_sync_fn));
  exports.Set(Napi::String::New(env, "barrier_all_sync"),
              Napi::Function::New(env, shmem_barrier_all_sync_fn));
  exports.Set(Napi::String::New(env, "set_lock_sync"),
              Napi::Function::New(env, shmem_set_lock_sync_fn));
  exports.Set(Napi::String::New(env, "clear_lock_sync"),
              Napi::Function::New(env, shmem_clear_lock_sync_fn));
  exports.Set(Napi::String::New(env, "broadcast64_sync"),
              Napi::Function::New(env, shmem_broadcast64_sync_fn));
  exports.Set(Napi::String::New(env, "long_sum_to_all_sync"),
              Napi::Function::New(env, shmem_long_sum_to_all_sync_fn));
  exports.Set(Napi::String::New(env, "double_sum_to_all_sync"),
              Napi::Function::New(env, shmem_double_sum_to_all_sync_fn));
  return Napi::Number::New(env, 1);
}

}
