
#include "hclib_resilience.h"

namespace hclib {
namespace resilience{

namespace diamond {

bool check_result_helper(promise_vector<void*>* put_vec,
        int replica1, int replica2) {

    for(auto && elem: *put_vec) {
        if(!elem->equal(replica1, replica2))
            return false;
    }
    return true;
}

//TODO: currently only works for triple redundancy
int check_result(promise_vector<void*>* put_vec, int N) {
    if(check_result_helper(put_vec, 0, 1))
        return 0;
    else if(check_result_helper(put_vec, 1, 2))
        return 1;
    else if(check_result_helper(put_vec, 0, 2))
        return 0;
    else
        return -1;
}

} // namespace diamond

} // namespace resilience
} // namespace hclib
