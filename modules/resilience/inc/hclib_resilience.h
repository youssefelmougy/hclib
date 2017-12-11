#ifndef HCLIB_RESILIENCE_H
#define HCLIB_RESILIENCE_H

#include "hclib_cpp.h"
#include "hclib_ref-count.h"
#include <vector>
#include <mutex>

namespace hclib {
namespace resilience {

//template<typename T>
//class obj {
//    T *data = nullptr;
//    int size = 0;
//
//  public:
//    obj(T *ptr, int sz) : data(ptr), size(sz) {}
//
//    bool equals(obj & obj2) {
//        return size==obj2.size && memcmp(data, obj2.data, size)==0;
//    }
//};

#include "hclib_resilience_diamond.h"

}
}

#endif

