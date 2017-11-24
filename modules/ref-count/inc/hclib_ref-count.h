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
        if(ref_count == NULL) return;
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

/*
Async
*/

template <typename T>
class AsyncRefCount {
    T _lambda;
    hclib_future_t *future1, *future2, *future3, *future4;

  public:
    AsyncRefCount(T lambda, hclib_future_t *f1, hclib_future_t *f2,
        hclib_future_t *f3, hclib_future_t *f4)
        : _lambda(lambda), future1(f1), future2(f2)
        , future3(f3), future4(f4) {}

    void operator() () {
        _lambda();
	future_t<void*> *f = static_cast<future_t<void*>*>(future1);
	if(f != NULL)
	    f->release();
        f = static_cast<future_t<void*>*>(future2);
	if(f != NULL)
	    f->release();
        f = static_cast<future_t<void*>*>(future3);
	if(f != NULL)
	    f->release();
        f = static_cast<future_t<void*>*>(future4);
	if(f != NULL)
	    f->release();
    }
};

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future) {
    AsyncRefCount<T> async_ref_count(lambda, future, NULL, NULL, NULL);
    hclib::async_await(async_ref_count, future);
}

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2) {
    AsyncRefCount<T> async_ref_count(lambda, future1, future2, NULL, NULL);
    hclib::async_await(async_ref_count, future1, future2);
}

template <typename T>
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2, hclib_future_t *future3, hclib_future_t *future4) {
    AsyncRefCount<T> async_ref_count(lambda, future1, future2, future3, future4);
    hclib::async_await(async_ref_count, future1, future2, future3, future4);
}

} // namespace ref_count
} // namespace hclib

#endif
