
#include "addon_hclib_node.js_openshmem.h"
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
    init_async_helper(env);
    return Napi::Number::New(env, 1);
}

Napi::Number finalize_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    //shmem_finalize();
    finalize_async_helper(env);
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

Napi::Number get_long_value_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t index = info[1].As<Napi::Number>().Int64Value();
    return Napi::Number::New(env, ptr[index]);
}

Napi::Number put_long_value_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *ptr = (long*) info[0].As<Napi::Number>().Int64Value();
    long val = info[1].As<Napi::Number>().Int64Value();
    int64_t index = info[2].As<Napi::Number>().Int64Value();
    ptr[index] = val;
    return Napi::Number::New(env, 1);
}

Napi::Number get_double_value_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *ptr = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t index = info[1].As<Napi::Number>().Int64Value();
    return Napi::Number::New(env, ptr[index]);
}

Napi::Number put_double_value_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *ptr = (double*) info[0].As<Napi::Number>().Int64Value();
    double val = info[1].As<Napi::Number>().DoubleValue();
    int64_t index = info[2].As<Napi::Number>().Int64Value();
    ptr[index] = val;
    return Napi::Number::New(env, 1);
}

Napi::Value Init_util(Napi::Env env, Napi::Object exports) {
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
  exports.Set(Napi::String::New(env, "get_long_value"),
              Napi::Function::New(env, get_long_value_fn));
  exports.Set(Napi::String::New(env, "put_long_value"),
              Napi::Function::New(env, put_long_value_fn));
  exports.Set(Napi::String::New(env, "get_double_value"),
              Napi::Function::New(env, get_double_value_fn));
  exports.Set(Napi::String::New(env, "put_double_value"),
              Napi::Function::New(env, put_double_value_fn));
  return Napi::Number::New(env, 1);
}

}
