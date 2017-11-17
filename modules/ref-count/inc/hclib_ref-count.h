#ifndef HCLIB_REF_COUNT_H
#define HCLIB_REF_COUNT_H

#include "hclib_cpp.h"

namespace hclib {
namespace ref_count {

/*
Future
*/

template<typename T>
class future_t: public hclib::future_t<T> {
};

template<typename T>
class future_t<T*>: public hclib::future_t<T*> {

    void ref_count_decr();

  public:

    void release() {
      ref_count_decr();
    }
};

/*
Promise
*/

template<typename T>
class promise_t: public hclib::promise_t<T> {
};

// Specialized for pointers
template<typename T>
class promise_t<T*>: public hclib::promise_t<T*> {

    int *ref_count = NULL;

  public:

    promise_t(int n) {
      ref_count = new int(n);
    }

    void ref_count_decr() {
        int cnt = __sync_sub_and_fetch(ref_count, 1);
	HASSERT(cnt >= 0);
	if(cnt == 0) {
	    free(hclib_promise_t::datum);
            hclib_promise_t::datum = NULL;
	    delete ref_count;
	}
    }

    future_t<T*> *get_future() {
        return static_cast<future_t<T*>*>( hclib::promise_t<T*>::get_future());
    }
};

template<typename T>
void future_t<T*>::ref_count_decr() {
    promise_t<T*> * p = static_cast<promise_t<T*>*>(hclib_future_t::owner);
    p->ref_count_decr();
}

} // namespace ref_count
} // namespace hclib

#endif
