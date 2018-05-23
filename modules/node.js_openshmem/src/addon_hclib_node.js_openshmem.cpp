
#include "hclib_node.js_openshmem.h"
#include <napi.h>
#include <shmem.h>

namespace hclib {

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

Napi::Number my_rank_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, shmem_my_pe());
}

Napi::Number num_ranks_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    return Napi::Number::New(env, shmem_n_pes());
}

Napi::Object Init(Napi::Env env, Napi::Object exports) {
  exports.Set(Napi::String::New(env, "init"),
              Napi::Function::New(env, init_fn));
  exports.Set(Napi::String::New(env, "finalize"),
              Napi::Function::New(env, finalize_fn));
  exports.Set(Napi::String::New(env, "my_rank"),
              Napi::Function::New(env, my_rank_fn));
  exports.Set(Napi::String::New(env, "num_ranks"),
              Napi::Function::New(env, num_ranks_fn));
  return exports;
}

NODE_API_MODULE(NODE_GYP_MODULE_NAME, Init)

}
