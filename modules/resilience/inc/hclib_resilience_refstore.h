#ifndef HCLIB_RESILIENCE_REFSTORE_H
#define HCLIB_RESILIENCE_REFSTORE_H

#include <unordered_map>
#include "hclib_resilience_checkpoint.h"

namespace hclib {
namespace resilience {
namespace checkpoint {

//performs a copy store of the objects in memory
class ref_store: public archive_store {

  std::unordered_map<void*, archive_obj*> store;

  public:
    ~ref_store() {
        for (auto iter : store) 
            delete iter.second;
    }

    void save(void* key, archive_obj* data) {
        assert(store.count(key) == 0);
        auto ptr = new archive_obj(data);
        store[key] = ptr;
    }

    archive_obj* retrieve(void* key) {
        if(store.count(key) == 1)
            return new archive_obj(store[key]);
        else
            return nullptr;
    }
};

} // namespace checkpoint
} // namespace resilience
} // namespace hclib

#endif 

