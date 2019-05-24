
#ifndef HCLIB_RESILIENCE_REPLAYN_H
#define HCLIB_RESILIENCE_REPLAYN_H

namespace hclib {
namespace resilience {
namespace replay {

template <int LEAF, int N=2, typename T>
std::enable_if_t< LEAF==0, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        hclib_future_t *f1, hclib_future_t *f2,
        hclib_future_t *f3, hclib_future_t *f4,
        hclib_locale_t *locale) {


#ifdef USE_RESILIENT_PROMISE
    HASSERT_STATIC(N>=2 && N<=N_CNT, "Currently only supports double and triple replay\n");
#endif

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        auto rtp = new replay_task_params_t<void*>();;
        rtp->put_vec = new promise_vector<void*>();
        rtp->rel_vec = new future_vector<void*>();
        bool result = false;
        int index = -1;

        assert(*(hclib_get_curr_task_local()) ==  nullptr);
        for(int i=0; i<N; i++) {
            hclib::finish([=]() {
                rtp->index = i;
                *(hclib_get_curr_task_local()) = rtp;
                async_await_at(*lambda_ptr, f1, f2, f3, f4, locale);
            });
            *(hclib_get_curr_task_local()) = nullptr;

            //ran successfully without errors
            if(error_check_fn(params) == 1)
            {
                result = true;
                index = rtp->index;
                break;
            }
            auto data = rtp->put_vec->data();
            auto size = rtp->put_vec->size();
            for(int i=0; i<size; i++)
                data[i]->tmp_delete();
            //rtp->put_vec->clear();
            //rtp->rel_vec->clear();
        }
        delete lambda_ptr;

        if(result) {
            rtp->put_vec->do_puts(index);
            rtp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete rtp->put_vec;
        delete rtp->rel_vec;
        delete rtp;;
    }, f1, f2, f3, f4, locale);
}

template <int LEAF=1, int N=2, typename T>
std::enable_if_t< LEAF==1, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        hclib_future_t *f1, hclib_future_t *f2,
        hclib_future_t *f3, hclib_future_t *f4,
        hclib_locale_t *locale) {

#ifdef USE_RESILIENT_PROMISE
    HASSERT_STATIC(N>=2 && N<=N_CNT, "Currently only supports double and triple replay\n");
#endif

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        auto rtp = new replay_task_params_t<void*>();;
        rtp->put_vec = new promise_vector<void*>();
        rtp->rel_vec = new future_vector<void*>();
        bool result = false;
        int index = -1;

        //assert(*(hclib_get_curr_task_local()) ==  nullptr);
        for(int i=0; i<N; i++) {
            rtp->index = i;
            *(hclib_get_curr_task_local()) = rtp;
            async_await_sync(lambda_ptr, f1, f2, f3, f4);

            *(hclib_get_curr_task_local()) = nullptr;

            //ran successfully without errors
            if(error_check_fn(params) == 1)
            {
                result = true;
                index = rtp->index;
                break;
            }

            auto data = rtp->put_vec->data();
            auto size = rtp->put_vec->size();
            for(int i=0; i<size; i++)
                data[i]->tmp_delete();
            //rtp->put_vec->clear();
            //rtp->rel_vec->clear();
        }
        delete lambda_ptr;

        if(result) {
            rtp->put_vec->do_puts(index);
            rtp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete rtp->put_vec;
        delete rtp->rel_vec;
        delete rtp;;
    }, f1, f2, f3, f4, locale);
}

template <int LEAF=1, int N=2, typename T>
inline void async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        hclib_future_t *f1, hclib_future_t *f2=nullptr,
        hclib_future_t *f3=nullptr, hclib_future_t *f4=nullptr) {
    async_await_check_at<LEAF, N>(lambda, prom_check,
        error_check_fn, params, f1, f2, f3, f4, nullptr);
}

template <int LEAF, int N=2, typename T>
std::enable_if_t< LEAF==0, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        std::vector<hclib_future_t *> *futures, hclib_locale_t *locale) {

#ifdef USE_RESILIENT_PROMISE
    HASSERT_STATIC(N>=2 && N<=N_CNT, "Currently only supports double and triple replay\n");
#endif

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        auto rtp = new replay_task_params_t<void*>();;
        rtp->put_vec = new promise_vector<void*>();
        rtp->rel_vec = new future_vector<void*>();
        bool result = false;
        int index = -1;

        assert(*(hclib_get_curr_task_local()) ==  nullptr);
        for(int i=0; i<N; i++) {
            hclib::finish([=]() {
                rtp->index = i;
                *(hclib_get_curr_task_local()) = rtp;
                async_await_at(*lambda_ptr, futures, locale);
            });

            *(hclib_get_curr_task_local()) = nullptr;

            ///ran successfully without errors
            if(error_check_fn(params) == 1)
            {
                 result = true;
                 index = rtp->index;
                 break;
            }
            auto data = rtp->put_vec->data();
            auto size = rtp->put_vec->size();
            for(int i=0; i<size; i++)
                data[i]->tmp_delete();
            //rtp->put_vec->clear();
            //rtp->rel_vec->clear();
        }
        delete lambda_ptr;

        if(result) {
            rtp->put_vec->do_puts(index);
            rtp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete rtp->put_vec;
        delete rtp->rel_vec;
        delete rtp;;
    }, futures, locale);
}

template <int LEAF=1, int N=2, typename T>
std::enable_if_t< LEAF==1, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        std::vector<hclib_future_t *> *futures, hclib_locale_t *locale) {

#ifdef USE_RESILIENT_PROMISE
    HASSERT_STATIC(N>=2 && N<=N_CNT, "Currently only supports double and triple replay\n");
#endif

    typedef typename std::remove_reference<T>::type U;
    U* lambda_ptr = new U(lambda);

    hclib::async_await_at([=]() {
        auto rtp = new replay_task_params_t<void*>();;
        rtp->put_vec = new promise_vector<void*>();
        rtp->rel_vec = new future_vector<void*>();
        bool result = false;
        int index = -1;

        //assert(*(hclib_get_curr_task_local()) ==  nullptr);
        for(int i=0; i<N; i++) {
            rtp->index = i;
            *(hclib_get_curr_task_local()) = rtp;
            async_await_sync(lambda_ptr, futures);

            *(hclib_get_curr_task_local()) = nullptr;

            //ran successfully without errors
            if(error_check_fn(params) == 1)
            {
                result = true;
                index = rtp->index;
                break;
            }

            auto data = rtp->put_vec->data();
            auto size = rtp->put_vec->size();
            for(int i=0; i<size; i++)
                data[i]->tmp_delete();
            //rtp->put_vec->clear();
            //rtp->rel_vec->clear();
        }
        delete lambda_ptr;

        if(result) {
            rtp->put_vec->do_puts(index);
            rtp->rel_vec->do_releases();
            prom_check->put(1);
        }
        else {
            prom_check->put(0);
        }
        delete rtp->put_vec;
        delete rtp->rel_vec;
        delete rtp;;
    }, futures, locale);
}

template <int LEAF=1, int N=2, typename T>
inline void async_await_check(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(void*)> error_check_fn, void * params,
        std::vector<hclib_future_t *> *futures) {
    async_await_check_at<LEAF, N>(lambda, prom_check,
        error_check_fn, params, futures, nullptr);
}

} // namespace replay
} // namespace resilience
} // namespace hclib

#endif
