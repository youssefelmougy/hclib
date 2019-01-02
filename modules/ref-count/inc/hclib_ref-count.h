#ifndef HCLIB_REF_COUNT_H
#define HCLIB_REF_COUNT_H

#include "hclib_cpp.h"

namespace hclib {
namespace ref_count {

/*
Future
*/

enum class DelType {NORMAL, ARR, FREE};
enum PromiseType {REF_COUNT=2};

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
Type erased deleter
*/

class deleter_t {
  public:
   virtual ~deleter_t() {}
   virtual void operator()( void* ) = 0;
};

template <typename T>
class deleter_norm : public deleter_t {
  public:
    void operator()( void* p ) {
      delete static_cast<T*>(p);
    }
};

template <typename T>
class deleter_arr : public deleter_t {
  public:
    void operator()( void* p ) {
      delete[] static_cast<T*>(p);
    }
};

template <typename T>
class deleter_free : public deleter_t {
  public:
    void operator()( void* p ) {
      free(p);
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

    int *ref_count = nullptr;
  protected:
    deleter_t* deleter = nullptr;

  public:
    static const int TYPE = PromiseType::REF_COUNT;

    promise_t(int n, DelType d_type=DelType::NORMAL) 
    : ref_count(new int(n)) {

      hclib_promise_t::type = TYPE;

      switch(d_type) {
        case DelType::NORMAL:
	  deleter = new deleter_norm<T>();
	  break;
	case DelType::ARR:
	  deleter = new deleter_arr<T>();
	  break;
	case DelType::FREE:
	  deleter = new deleter_free<T>();
	  break;
      }
    }

    ~promise_t() {
      delete ref_count;
      delete deleter;
    }

    void ref_count_decr() {
        if(hclib_promise_t::type < TYPE) return;

        HASSERT(ref_count != nullptr);
        int cnt = __sync_sub_and_fetch(ref_count, 1);
	HASSERT(cnt >= 0);
	if(cnt == 0) {
            (*deleter)(hclib_promise_t::datum);
            hclib_promise_t::datum = nullptr;
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
inline void async_await(T&& lambda, hclib_future_t *future1,
        hclib_future_t *future2=nullptr, hclib_future_t *future3=nullptr,
        hclib_future_t *future4=nullptr) {

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await( [=]() {
        (*lambda_ptr)();
        delete lambda_ptr;
	if(future1 != nullptr)
	    static_cast<future_t<void*>*>(future1)->release();
	if(future2 != nullptr)
	    static_cast<future_t<void*>*>(future2)->release();
	if(future3 != nullptr)
	    static_cast<future_t<void*>*>(future3)->release();
	if(future4 != nullptr)
	    static_cast<future_t<void*>*>(future4)->release();
    }, future1, future2, future3, future4);
}

template <typename T>
inline void async_await(T&& lambda, std::vector<hclib_future_t *> *futures){

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await( [=]() {
        (*lambda_ptr)();
        delete lambda_ptr;

        for(auto && elem: *(futures))
            static_cast<future_t<void*>*>(elem)->release();
    }, futures);
}

} // namespace ref_count
} // namespace hclib

#endif
