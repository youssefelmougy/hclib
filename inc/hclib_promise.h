#ifndef HCLIB_PROMISE_H
#define HCLIB_PROMISE_H

#include "hclib-promise.h"
#include "hclib_future.h"

namespace hclib {

class promise_t {
    public:
        hclib_promise_t internal;

        promise_t() {
            hclib_promise_init(&internal);
        }
        ~promise_t() { }

        void put(void *datum) {
            hclib_promise_put(&internal, datum);
        }

        future_t *get_future() {
            return new future_t(&internal.future);
        }
};

}

#endif
