#ifndef HCLIB_RESILIENCE_DIAMOND2_H
#define HCLIB_RESILIENCE_DIAMOND2_H

namespace hclib {
namespace resilience {
namespace diamond {

/*
duplicates the task and checks for equivalence of all puts
*/
template <int N=N_CNT-1, typename T>
std::enable_if_t< N==N_CNT-1, void>
async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await([=]() {
      //dtp_arr is used to pass around required data of each replica
      auto dtp_arr = new diamond_task_params_t<void*>[N+1];
      //put_vec is used to collect all the put() operations
      //peformed within the finish scope
      auto put_vec = new promise_vector<void*>();
      //rel_vec is used to collect all futures used as dependences
      //so that they can be released later if third task is not created
      auto rel_vec = new future_vector<void*>();

      hclib::finish([=]() {
        //create N tasks and give each of it dtp_arr[i]
        for(int i=0; i<N; i++) {
            dtp_arr[i].index = i;
            dtp_arr[i].put_vec = put_vec;
            dtp_arr[i].rel_vec = rel_vec;
            *(hclib_get_curr_task_local()) = &dtp_arr[i];
            async_await(*lambda_ptr, f1, f2, f3, f4);
        }
      });

      bool result = check_result_helper(put_vec, 0, 1);
      //it all puts in both replica are same then do the actual put
      if(result) {
          *(hclib_get_curr_task_local()) = nullptr;
          //perform the actual put
          for(auto && elem: *put_vec)
              elem->put_actual(0);
          //perform the release
          for(auto && elem: *rel_vec)
              elem->release();
          prom_check->put(1);
      }
      //if there is error in put, start a third task and then compare output
      else {
          hclib::finish([=]() {
              dtp_arr[N].index = N;
              dtp_arr[N].put_vec = put_vec;
              dtp_arr[N].rel_vec = rel_vec;
              *(hclib_get_curr_task_local()) = &dtp_arr[N];
              async_await(*lambda_ptr, f1, f2, f3, f4);
          });

          *(hclib_get_curr_task_local()) = nullptr;
          int index = -1;
	  result = check_result_helper(put_vec, 0, 2);
	  if(result)
	    index = 0;
	  else {
	    result = check_result_helper(put_vec, 1, 2);
	    if(result)
	      index = 1;
	  }

	  if(index < 0) 
              prom_check->put(0);
	  else {
	      //perform the actual put
	      for(auto && elem: *put_vec)
	          elem->put_actual(index);
              //perform the release
              for(auto && elem: *rel_vec)
                  elem->release();
	      prom_check->put(1);
	  }
      }

      delete lambda_ptr;
      delete put_vec;
      delete rel_vec;
      delete[] dtp_arr;
    }, f1, f2, f3, f4);
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif
