#ifndef HCLIB_RESILIENCE_DIAMOND2_LEAF_H
#define HCLIB_RESILIENCE_DIAMOND2_LEAF_H

namespace hclib {
namespace resilience {
namespace diamond {

#define release_futures(f1, f2, f3, f4) \
do {\
    if(f1 != nullptr) \
        static_cast<future_t<void*>*>(f1)->release();\
    if(f2 != nullptr) \
        static_cast<future_t<void*>*>(f2)->release();\
    if(f3 != nullptr) \
        static_cast<future_t<void*>*>(f3)->release();\
    if(f4 != nullptr) \
        static_cast<future_t<void*>*>(f4)->release();\
} while(0)

#define release_futures_vec(futures) \
do {\
    auto data_var = futures->data();\
    auto size_var = futures->size();\
    for(int i=0; i<size_var; i++)\
      static_cast<future_t<void*>*>(data_var[i])->release();\
} while(0)

template <typename T>
inline void async_await_at_leaf(T&& lambda, hclib_locale_t *locale) {
    auto dtp = *hclib_get_curr_task_local();
    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);
    hclib::async_nb_await_at( [=]() {
        *(hclib_get_curr_task_local()) = dtp;
        (*lambda_ptr)();
        delete lambda_ptr;
     }, future_nullptr, locale);
}

/*
duplicates the task and checks for equivalence of all puts
*/
template <int LEAF=1, int N=N_CNT-1, typename T>
std::enable_if_t< LEAF==1 && N==N_CNT-1, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        hclib_future_t *f1, hclib_future_t *f2,
        hclib_future_t *f3, hclib_future_t *f4,
        hclib_locale_t *locale) {

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        //dtp_arr is used to pass around required data of each replica
        auto dtp_arr = new diamond_task_params_t<void*>[N+1];
        //put_vec is used to collect all the put() operations
        //peformed within the finish scope
        auto put_vec = new promise_vector<void*>();
        //rel_vec is used to collect all futures used as dependences
        //so that they can be released later if third task is not created
        //auto rel_vec = new future_vector<void*>();

        hclib::finish([=]() {
            //create N tasks and give each of it dtp_arr[i]
            for(int i=0; i<N; i++) {
                dtp_arr[i].index = i;
                dtp_arr[i].put_vec = put_vec;
                *(hclib_get_curr_task_local()) = &dtp_arr[i];
                async_await_at_leaf(*lambda_ptr, locale);
            }
        });

        *(hclib_get_curr_task_local()) = nullptr;
        bool result = check_result_helper(put_vec, 0, 1);
        //it all puts in both replica are same then do the actual put
        if(result) {
            put_vec->do_puts(0);
            release_futures(f1, f2, f3, f4);
            //rel_vec->do_releases();
            prom_check->put(1);
        }
        //if there is error in put, start a third task and then compare output
        else {
                dtp_arr[N].index = N;
                dtp_arr[N].put_vec = put_vec;
                *(hclib_get_curr_task_local()) = &dtp_arr[N];
                (*lambda_ptr)();

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
                put_vec->do_puts(index);
                release_futures(f1, f2, f3, f4);
                //rel_vec->do_releases();
                prom_check->put(1);
            }
        }
        delete lambda_ptr;
        delete put_vec;
        //delete rel_vec;
        delete[] dtp_arr;
    }, f1, f2, f3, f4, locale);
}

template <int LEAF=1, int N=N_CNT-1, typename T>
std::enable_if_t< LEAF==1 && N==N_CNT-1, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::vector<hclib_future_t *> *futures, hclib_locale_t *locale) {

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        //dtp_arr is used to pass around required data of each replica
        auto dtp_arr = new diamond_task_params_t<void*>[N+1];
        //put_vec is used to collect all the put() operations
        //peformed within the finish scope
        auto put_vec = new promise_vector<void*>();
        //rel_vec is used to collect all futures used as dependences
        //so that they can be released later if third task is not created
        //auto rel_vec = new future_vector<void*>();

        hclib::finish([=]() {
            //create N tasks and give each of it dtp_arr[i]
            for(int i=0; i<N; i++) {
                dtp_arr[i].index = i;
                dtp_arr[i].put_vec = put_vec;
                *(hclib_get_curr_task_local()) = &dtp_arr[i];
		async_await_at_leaf(*lambda_ptr, locale);
            }
        });

        *(hclib_get_curr_task_local()) = nullptr;
        bool result = check_result_helper(put_vec, 0, 1);
        //it all puts in both replica are same then do the actual put
        if(result) {
            put_vec->do_puts(0);
	    release_futures_vec(futures);
            //rel_vec->do_releases();
            prom_check->put(1);
        }
        //if there is error in put, start a third task and then compare output
        else {
                dtp_arr[N].index = N;
                dtp_arr[N].put_vec = put_vec;
                *(hclib_get_curr_task_local()) = &dtp_arr[N];
		(*lambda_ptr)();

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
              put_vec->do_puts(index);
	      release_futures_vec(futures);
              //rel_vec->do_releases();
              prom_check->put(1);
            }
        }

        delete lambda_ptr;
        delete put_vec;
        //delete rel_vec;
        delete[] dtp_arr;
    }, futures, locale);
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif
