#ifndef HCLIB_RESILIENCE_H
#define HCLIB_RESILIENCE_H

#include "hclib_cpp.h"
#include "hclib_ref-count.h"

#if defined(USE_STD_VEC) || defined(USE_C_ARRAY_WITH_ATOMIC) || defined(USE_C_ARRAY_WITH_LOCK)
#else
#define USE_C_ARRAY_WITH_ATOMIC
#endif

//#define USE_STD_VEC
//#define USE_C_ARRAY_WITH_ATOMIC
//#define USE_C_ARRAY_WITH_LOCK

//Hack to remove atomics/lock when leaf task are used
#define ONLY_LEAF

#ifdef USE_STD_VEC
#include <vector>
#include <mutex>
#elif defined USE_C_ARRAY_WITH_ATOMIC
#include "hclib-atomics.h"
#elif defined USE_C_ARRAY_WITH_LOCK
#include <mutex>
#else
#error No data container for resilient promises specified
#endif

//TODO: Assumes lambdas provided to tasks are non mutable.
//Need to find out what happens it not?

namespace hclib {
namespace resilience {

/*
Resilient object.
Users should extend this class to provide their own
equals method to be used for comparison at the end
of resilient task
*/
class obj {
  public:
    virtual ~obj() {}

    //return true if objs are equal else false
    virtual bool equals(obj* obj2) { assert(false); return false; }
};

namespace util {
/*
vector that allows concurrent insertion
Read operations are not concurrent with insertion
*/
template<typename T>
class safe_vector {
#ifdef USE_STD_VEC
    std::vector<T> vec;
    std::mutex mtx;

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
        std::lock_guard<std::mutex> lkg(mtx);
        vec.push_back(std::forward<T>(value));
    }

    size_t size() const noexcept {
        return vec.size();
    }

    T* data() noexcept {
        return vec.data();
    }

    void clear() noexcept {
        vec.clear();
    }

    //auto begin() noexcept -> decltype(vec.begin()) {
    //    return vec.begin();
    //}

    //auto end() noexcept -> decltype(vec.end()) {
    //    return vec.end();
    //}
#elif defined (USE_C_ARRAY_WITH_ATOMIC) || defined (USE_C_ARRAY_WITH_LOCK)
#ifndef SAFE_VECTOR_CAPACITY
#define SAFE_VECTOR_CAPACITY 16
#endif
    T vec[SAFE_VECTOR_CAPACITY];
#ifdef ONLY_LEAF
    int pos = -1;
#else
    volatile int pos = -1;
#endif
#ifdef USE_C_ARRAY_WITH_LOCK
    std::mutex mtx;
#endif

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
#ifdef ONLY_LEAF
        assert(pos < SAFE_VECTOR_CAPACITY-1);
        vec[++pos] = value;
#elif defined (USE_C_ARRAY_WITH_LOCK)
        std::lock_guard<std::mutex> lkg(mtx);
        pos++;
        assert(pos < SAFE_VECTOR_CAPACITY);
        vec[pos] = value;
#elif defined (USE_C_ARRAY_WITH_ATOMIC)
        const int pos1 = hc_atomic_inc(&pos);
        assert(pos1 < SAFE_VECTOR_CAPACITY);
        vec[pos1] = std::forward<T>(value);
#endif
    }

    size_t size() const noexcept {
        return pos + 1;
    }

    T* data() noexcept {
        return &vec[0];
    }

    void clear() noexcept {
    // TODO: mutual exclusion?
        pos = -1;
    }
#endif
};

template<typename T>
class safe_promise_vector : public safe_vector<T> {
  public:
    void do_puts(int index) {
        //perform the actual put
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++)
            data_var[i]->put_actual(index);
    }
};

template<typename T>
class safe_future_vector : public safe_vector<T> {
  public:
    void do_releases() {
        //perform the release
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++)
            data_var[i]->release();
    }
};

} // namespace util

} // namespace resilience
} // namespace hclib

#include "hclib_resilience_promise.h"
#include "hclib_resilience_diamond.h"
#include "hclib_resilience_replay.h"
#include "hclib_resilience_abft.h"
#include "hclib_resilience_checkpoint.h"

#endif

