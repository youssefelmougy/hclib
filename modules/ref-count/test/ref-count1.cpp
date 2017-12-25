
#include "hclib_cpp.h"
#include "hclib_ref-count.h"
#include <unistd.h>

int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system"};
    hclib::launch(deps, 1, [=]() {
        hclib::ref_count::promise_t<int*> *prom = new hclib::ref_count::promise_t<int*>(2,
                                                      hclib::ref_count::DelType::ARR);
        hclib::finish([=]() {

            hclib::ref_count::async_await([=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value task1 %d\n", *signal);
                }, prom->get_future());

            hclib::async_await([=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value task2 %d\n", *signal);
                    prom->get_future()->release();
                }, prom->get_future());

            hclib::async_await([=]() {
                    sleep(1);
                    int *n= new int[1];
                    n[0] = SIGNAL_VALUE;
                    prom->put(n);
                }, nullptr);
        });
        assert(prom->get_future()->get() == NULL);
    });
    printf("Exiting...\n");
    return 0;
}

