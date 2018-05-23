#ifndef HCLIB_NODE_JS_OPENSHMEM_H
#define HCLIB_NODE_JS_OPENSHMEM_H

#include "hclib-module.h"
#include "hclib_cpp.h"

extern int nic_locale_id;
extern hclib::locale_t *nic;

namespace hclib {

HCLIB_MODULE_INITIALIZATION_FUNC(nodejs_openshmem_pre_initialize);
HCLIB_MODULE_INITIALIZATION_FUNC(nodejs_openshmem_post_initialize);
HCLIB_MODULE_INITIALIZATION_FUNC(nodejs_openshmem_finalize);

int shmem_my_pe();
int shmem_n_pes();

}

#endif
