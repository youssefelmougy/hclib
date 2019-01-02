#ifndef HCLIB_RESILIENCE_DIAMOND3_H
#define HCLIB_RESILIENCE_DIAMOND3_H

namespace hclib {
namespace resilience {
namespace diamond {

/*
triplicates the task and checks for equivalence of all puts
*/
//TODO: For now only both leaf task and non leaf task behave the same.
//Create an optimized version for leaf tasks.
template <int LEAF, int N, typename T>
std::enable_if_t< N>=N_CNT, void>
//typename std::enable_if< N>=3, void>::type
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        hclib_future_t *f1, hclib_future_t *f2,
        hclib_future_t *f3, hclib_future_t *f4,
        hclib_locale_t *locale) {

    HASSERT_STATIC(N<=N_CNT, "Currently only supports double and triple redudancy\n");

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        //dtp_arr is used to pass around required data of each replica
        auto dtp_arr = new diamond_task_params_t<void*>[N];
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
                async_await_at(*lambda_ptr, f1, f2, f3, f4, locale);
            }
        });
        delete lambda_ptr;

        *(hclib_get_curr_task_local()) = nullptr;
        int index = check_result(put_vec);
        if(index < 0 )
            prom_check->put(0);
        else {
            put_vec->do_puts(index);
            rel_vec->do_releases();
            prom_check->put(1);
        }

        delete put_vec;
        delete rel_vec;
        delete[] dtp_arr;
    }, f1, f2, f3, f4, locale);
}

template <int LEAF, int N, typename T>
std::enable_if_t< N>=N_CNT, void>
//typename std::enable_if< N>=3, void>::type
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::vector<hclib_future_t *> *futures, hclib_locale_t *locale) {

    HASSERT_STATIC(N<=N_CNT, "Currently only supports double and triple redudancy\n");

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        //dtp_arr is used to pass around required data of each replica
        auto dtp_arr = new diamond_task_params_t<void*>[N];
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
                async_await_at(*lambda_ptr, futures, locale);
            }
        });
        delete lambda_ptr;

        *(hclib_get_curr_task_local()) = nullptr;
        int index = check_result(put_vec);
        if(index < 0 )
            prom_check->put(0);
        else {
            put_vec->do_puts(index);
            rel_vec->do_releases();
            prom_check->put(1);
        }

        delete put_vec;
        delete rel_vec;
        delete[] dtp_arr;
    }, futures, locale);
}

template <int LEAF=1, int N=N_CNT-1, typename T>
inline void async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {
    async_await_check_at<LEAF, N>(lambda, prom_check, f1, f2, f3, f4, nullptr);
}

template <int LEAF=1, int N=N_CNT-1, typename T>
inline void async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        std::vector<hclib_future_t *> *futures) {
    async_await_check_at<LEAF, N>(lambda, prom_check, futures, nullptr);
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif

