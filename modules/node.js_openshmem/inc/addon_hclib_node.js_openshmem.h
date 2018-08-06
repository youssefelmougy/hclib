#ifndef ADDON_HCLIB_NODEJS_OPENSHMEM_H
#define ADDON_HCLIB_NODEJS_OPENSHMEM_H
#include <napi.h>

namespace hclib {

Napi::Value Init_util(Napi::Env env, Napi::Object exports);

Napi::Value Init_sync(Napi::Env env, Napi::Object exports);

Napi::Value Init_async(Napi::Env env, Napi::Object exports);

void init_async_helper(Napi::Env env);
void finalize_async_helper(Napi::Env env);

}

#endif //ADDON_HCLIB_NODEJS_OPENSHMEM_H

