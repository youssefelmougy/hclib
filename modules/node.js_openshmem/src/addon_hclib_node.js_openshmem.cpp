
#include "hclib_node.js_openshmem.h"
#include <napi.h>
#include <shmem.h>

namespace hclib {

Napi::Number get_SHMEM_BCAST_SYNC_SIZE_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, SHMEM_BCAST_SYNC_SIZE);
}

Napi::Number get_SHMEM_BARRIER_SYNC_SIZE_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, SHMEM_BARRIER_SYNC_SIZE);
}

Napi::Number get_SHMEM_REDUCE_SYNC_SIZE_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, SHMEM_REDUCE_SYNC_SIZE);
}

Napi::Number get_SHMEM_REDUCE_MIN_WRKDATA_SIZE_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, SHMEM_REDUCE_MIN_WRKDATA_SIZE);
}

Napi::Number get_SHMEM_SYNC_VALUE_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, SHMEM_SYNC_VALUE);
}

Napi::Number init_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    //shmem_init();
    return Napi::Number::New(env, 1);
}

Napi::Number finalize_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    //shmem_finalize();
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_my_rank_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, shmem_my_pe());
}

Napi::Number shmem_num_ranks_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, shmem_n_pes());
}

Napi::Number get_value_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    int index = info[1].As<Napi::Number>().Int64Value();
    return Napi::Number::New(env, ptr[index]);
}

Napi::Number put_value_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    long val = info[1].As<Napi::Number>().Int64Value();
    int index = info[2].As<Napi::Number>().Int64Value();
    ptr[index] = val;
    return Napi::Number::New(env, 1);
}

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
    void *src = (void*) info[0].As<Napi::Number>().Int64Value();
    int pe = info[1].As<Napi::Number>().Int64Value();
    long *val = (long *)malloc(sizeof(long));
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            *val = shmem_long_g((long*)src, pe);
        }, nic);
    });
    long val_tmp = *val;
    free(val);
    return Napi::Number::New(env, val_tmp);
}

Napi::Number shmem_long_p_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *dest = (void*) info[0].As<Napi::Number>().Int64Value();
    long val = info[1].As<Napi::Number>().Int64Value();
    int pe = info[2].As<Napi::Number>().Int64Value();
    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_long_p((long*)dest, val, pe);
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

Napi::Number shmem_int_sum_to_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    int *dest = (int*) info[0].As<Napi::Number>().Int64Value();
    int *src = (int*) info[1].As<Napi::Number>().Int64Value();
    int nreduce = info[2].As<Napi::Number>().Int64Value();
    int PE_start = info[3].As<Napi::Number>().Int64Value();
    int logPE_stride = info[4].As<Napi::Number>().Int64Value();
    int PE_size = info[5].As<Napi::Number>().Int64Value();
    int *pWrk = (int*) info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_int_sum_to_all(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_longlong_sum_to_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long long *dest = (long long*) info[0].As<Napi::Number>().Int64Value();
    long long *src = (long long*) info[1].As<Napi::Number>().Int64Value();
    int nreduce = info[2].As<Napi::Number>().Int64Value();
    int PE_start = info[3].As<Napi::Number>().Int64Value();
    int logPE_stride = info[4].As<Napi::Number>().Int64Value();
    int PE_size = info[5].As<Napi::Number>().Int64Value();
    long long *pWrk = (long long*) info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    hclib::finish([=] () {
        hclib::async_nb_at([=] () {
            shmem_longlong_sum_to_all(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
        }, nic);
    });
    return Napi::Number::New(env, 1);
}

Napi::Object Init(Napi::Env env, Napi::Object exports) {
  exports.Set(Napi::String::New(env, "get_SHMEM_BCAST_SYNC_SIZE"),
              Napi::Function::New(env, get_SHMEM_BCAST_SYNC_SIZE_fn));
  exports.Set(Napi::String::New(env, "get_SHMEM_BARRIER_SYNC_SIZE"),
              Napi::Function::New(env, get_SHMEM_BARRIER_SYNC_SIZE_fn));
  exports.Set(Napi::String::New(env, "get_SHMEM_REDUCE_SYNC_SIZE"),
              Napi::Function::New(env, get_SHMEM_REDUCE_SYNC_SIZE_fn));
  exports.Set(Napi::String::New(env, "get_SHMEM_REDUCE_MIN_WRKDATA_SIZE"),
              Napi::Function::New(env, get_SHMEM_REDUCE_MIN_WRKDATA_SIZE_fn));
  exports.Set(Napi::String::New(env, "get_SHMEM_SYNC_VALUE"),
              Napi::Function::New(env, get_SHMEM_SYNC_VALUE_fn));
  exports.Set(Napi::String::New(env, "init"),
              Napi::Function::New(env, init_fn));
  exports.Set(Napi::String::New(env, "finalize"),
              Napi::Function::New(env, finalize_fn));
  exports.Set(Napi::String::New(env, "my_rank"),
              Napi::Function::New(env, shmem_my_rank_fn));
  exports.Set(Napi::String::New(env, "num_ranks"),
              Napi::Function::New(env, shmem_num_ranks_fn));
  exports.Set(Napi::String::New(env, "get_value"),
              Napi::Function::New(env, get_value_fn));
  exports.Set(Napi::String::New(env, "put_value"),
              Napi::Function::New(env, put_value_fn));
  exports.Set(Napi::String::New(env, "malloc_sync"),
              Napi::Function::New(env, shmem_malloc_sync_fn));
  exports.Set(Napi::String::New(env, "free_sync"),
              Napi::Function::New(env, shmem_free_sync_fn));
  exports.Set(Napi::String::New(env, "long_g_sync"),
              Napi::Function::New(env, shmem_long_g_sync_fn));
  exports.Set(Napi::String::New(env, "long_p_sync"),
              Napi::Function::New(env, shmem_long_p_sync_fn));
  exports.Set(Napi::String::New(env, "barrier_all_sync"),
              Napi::Function::New(env, shmem_barrier_all_sync_fn));
  exports.Set(Napi::String::New(env, "set_lock_sync"),
              Napi::Function::New(env, shmem_set_lock_sync_fn));
  exports.Set(Napi::String::New(env, "clear_lock_sync"),
              Napi::Function::New(env, shmem_clear_lock_sync_fn));
  exports.Set(Napi::String::New(env, "broadcast64_sync"),
              Napi::Function::New(env, shmem_broadcast64_sync_fn));
  exports.Set(Napi::String::New(env, "int_sum_to_all_sync"),
              Napi::Function::New(env, shmem_int_sum_to_all_sync_fn));
  exports.Set(Napi::String::New(env, "longlong_sum_to_all_sync"),
              Napi::Function::New(env, shmem_longlong_sum_to_all_sync_fn));
  return exports;
}

NODE_API_MODULE(NODE_GYP_MODULE_NAME, Init)

}

