#ifndef HCLIB_RESILIENCE_H
#define HCLIB_RESILIENCE_H

#include "hclib_cpp.h"
#include "hclib_ref-count.h"
#include <vector>
#include <mutex>

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

    std::vector<T> vec;
    std::mutex mtx;

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
        std::lock_guard<std::mutex> lkg(mtx);
        vec.push_back(std::forward<T>(value));
    }

    auto begin() noexcept -> decltype(vec.begin()) {
        return vec.begin();
    }

    auto end() noexcept -> decltype(vec.end()) {
        return vec.end();
    }

    auto size() const noexcept -> decltype(vec.size()){
        return vec.size();
    }

    auto data() noexcept -> decltype(vec.data()) {
        return vec.data();
    }

    auto clear() noexcept -> decltype(vec.clear()) {
        vec.clear();
    }
};

} // namespace util

} // namespace resilience
} // namespace hclib

#include "hclib_resilience_diamond2.h"
#include "hclib_resilience_diamond3.h"
#include "hclib_resilience_replay.h"

#endif

