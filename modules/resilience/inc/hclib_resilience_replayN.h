
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

template <int LEAF=1, int N=2, typename T, typename V>
std::enable_if_t< LEAF==1, void>
async_await_check_at(T&& lambda, hclib::promise_t<int> *prom_check,
        std::function<int(V)> error_check_fn, V params,
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

template <int LEAF=1, int N=2, typename T>
inline void forasync1D_await_check(loop_domain_1d *loop, T&& lambda, hclib::promise_t<int> *prom_check,
                                   std::function<int(int)> error_check_fn, hclib_future_t *future) {
    std::vector<hclib_future_t *> result_futures;
    const hclib_loop_domain_t* loop0 = loop->get_internal();
    for (int i = loop0->low; i < loop0->high; i++) {
        hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
        replay::async_await_check_at<LEAF, N>([=]() { lambda(i); }, prom_res, error_check_fn, i, future, nullptr, nullptr, nullptr, nullptr);
        result_futures.push_back(prom_res->get_future());
    }
    hclib::async_await([=]() {
            bool result = true;
            for (int i = loop0->low; i < loop0->high; i++) {
                int res = ((hclib::future_t<int>*)result_futures[i])->get();
                if (!res) result = false;
            }
            if (result) {
                prom_check->put(1);
            } else {
                prom_check->put(0);
            }
        }, result_futures);
}

template <int LEAF=1, int N=2, typename T>
inline void forasync1D_await_bulkcheck(loop_domain_1d *loop, T&& lambda, hclib::promise_t<int> *prom_check,
                                       std::function<int()> error_check_fn, hclib_future_t *future) {
    hclib::async_await([=]() {
            auto rtps = new std::vector<replay_task_params_t<void*>*>();
            bool result = false;
            int index = -1;
            assert(*(hclib_get_curr_task_local()) ==  nullptr);
            for (int i = 0; i<N; i++) {
                hclib::finish([=] {
                        //                        hclib::forasync1D(loop, lambda, false, FORASYNC_MODE_FLAT, future);
                        const hclib_loop_domain_t* loop0 = loop->get_internal();
                        const int low = loop0->low, high = loop0->high, stride = loop0->stride, tile = loop0->tile;
                        const loop_dist_func func = hclib_lookup_dist_func(HCLIB_DEFAULT_LOOP_DIST);

                        int nb_chunks = (int) (high/tile);
                        int size = tile * nb_chunks;
                        int low0;
                        for (low0 = low; low0 < size; low0 += tile) {
                            auto rtp = new replay_task_params_t<void*>();;
                            rtp->put_vec = new promise_vector<void*>();
                            rtp->rel_vec = new future_vector<void*>();
                            rtp->index = i;
                            rtps->push_back(rtp);
                            hclib_loop_domain_t ld = {low0, low0 + tile, stride, tile};
                            auto lambda_wrapper = [=]() {
                                *(hclib_get_curr_task_local()) = rtp;
                                forasync1D_runner<T>(&ld, lambda);
                            };

                            hclib_locale_t *locale = func(1, &ld, loop0, FORASYNC_MODE_FLAT);
                            hclib::async_await_at(lambda_wrapper, future, locale);
                        }
                        // handling leftover
                        if (size < high) {
                            auto rtp = new replay_task_params_t<void*>();;
                            rtp->put_vec = new promise_vector<void*>();
                            rtp->rel_vec = new future_vector<void*>();
                            rtp->index = i;
                            rtps->push_back(rtp);
                            hclib_loop_domain_t ld = {low0, high, stride, tile};
                            auto lambda_wrapper = [=]() {
                                *(hclib_get_curr_task_local()) = rtp;
                                forasync1D_runner<T>(&ld, lambda);
                            };
                            hclib_locale_t *locale = func(1, &ld, loop0, FORASYNC_MODE_FLAT);
                            hclib::async_await_at(lambda_wrapper, future, locale);
                        }
                    });
                *(hclib_get_curr_task_local()) = nullptr;

                //ran successfully without errors
                if(error_check_fn() == true)
                    {
                        result = true;
                        index = (*rtps)[0]->index;
                        break;
                    }
                for (int i = 0; i < (*rtps).size(); i++) {
                    auto rtp = (*rtps)[i];
                    auto data = rtp->put_vec->data();
                    auto size = rtp->put_vec->size();
                    for(int i=0; i<size; i++)
                        data[i]->tmp_delete();
                }
            }
            //delete lambda_ptr;

            if(result) {
                for (int i = 0; i < (*rtps).size(); i++) {
                    auto rtp = (*rtps)[i];
                    rtp->put_vec->do_puts(index);
                    rtp->rel_vec->do_releases();
                }
                prom_check->put(1);
            }
            else {
                prom_check->put(0);
            }
            for (int i = 0; i < (*rtps).size(); i++) {
                auto rtp = (*rtps)[i];
                delete rtp->put_vec;
                delete rtp->rel_vec;
                delete rtp;
            }
        }, future);
}

} // namespace replay
} // namespace resilience
} // namespace hclib

#endif
