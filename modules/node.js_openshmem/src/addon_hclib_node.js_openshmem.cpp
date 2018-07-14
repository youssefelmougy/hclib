
#include "addon_hclib_node.js_openshmem.h"

namespace hclib {

Napi::Object Init(Napi::Env env, Napi::Object exports) {

  Init_util(env, exports);
  Init_sync(env, exports);
  Init_async(env, exports);
  return exports;
}

NODE_API_MODULE(NODE_GYP_MODULE_NAME, Init)

}

