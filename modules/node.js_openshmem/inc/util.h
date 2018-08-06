#ifndef ADDON_HCLIB_NODEJS_UTIL_H
#define ADDON_HCLIB_NODEJS_UTIL_H

#include<vector>
#include<mutex>

template<typename T>
class safe_vector {

    std::vector<T> vec;

  public:
    std::mutex mtx;

    //TODO: can this take both lvalue and rvalue?
    void push_back(T &&value) {
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

#endif //ADDON_HCLIB_NODEJS_UTIL_H
