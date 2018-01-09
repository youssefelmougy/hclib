#ifndef HCLIB_RESILIENCE_DIAMOND3_H
#define HCLIB_RESILIENCE_DIAMOND3_H

#include "hclib_resilience_diamond.h"

namespace hclib {
namespace resilience {
namespace diamond {

/*
triplicates the task and checks for equivalence of all puts
*/
template <int N, typename T>
std::enable_if_t< N>=N_CNT, void>
//typename std::enable_if< N>=3, void>::type
async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    HASSERT_STATIC(N<=N_CNT, "Currently only supports double and triple redudancy\n");

    hclib::async_await([=, lambda_mv1 = std::move(lambda)]() {
      //dtp_arr is used to pass around required data of each replica
      auto dtp_arr = new diamond_task_params_t<void*,N>[N];
      //put_vec is used to collect all the put() operations
      //peformed within the finish scope
      auto put_vec = new promise_vector<void*,N>();

      hclib::finish([=, lambda_mv2 = std::move(lambda_mv1)]() {
        //create N tasks and give each of it dtp_arr[i]
        for(int i=0; i<N; i++) {
            dtp_arr[i].index = i;
            dtp_arr[i].count = N;
            dtp_arr[i].put_vec = put_vec;
            *(hclib_get_curr_task_local()) = &dtp_arr[i];
            async_await<N>(lambda_mv2, f1, f2, f3, f4);
        }
      });

      int result = check_result<N>(put_vec);
      if(result < 0 )
          prom_check->put(0);
      else {
          //perform the actual put
          for(auto && elem: *put_vec)
              elem->put_actual(result);
          prom_check->put(1);
      }

      delete put_vec;
      delete[] dtp_arr;
    }, nullptr);
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif

