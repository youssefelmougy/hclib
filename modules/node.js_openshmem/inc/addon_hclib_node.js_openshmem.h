#ifndef ADDON_HCLIB_NODEJS_OPENSHMEM_H
#define ADDON_HCLIB_NODEJS_OPENSHMEM_H
#include <napi.h>

namespace hclib {

Napi::Value Init_util(Napi::Env env, Napi::Object exports);

Napi::Value Init_sync(Napi::Env env, Napi::Object exports);

Napi::Value Init_async(Napi::Env env, Napi::Object exports);
}

#endif //ADDON_HCLIB_NODEJS_OPENSHMEM_H

