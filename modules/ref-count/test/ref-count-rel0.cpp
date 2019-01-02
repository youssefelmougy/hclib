
#include "hclib_cpp.h"
#include "hclib_ref-count.h"
#include <unistd.h>

int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system"};
    hclib::launch(deps, 1, [=]() {
        hclib::finish([=]() {
            hclib::ref_count::promise_t<int*> *prom = new hclib::ref_count::promise_t<int*>(1);

            hclib::async_await([=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value %d\n", *signal);
		    prom->get_future()->release();
		    assert(prom->get_future()->get() == NULL);
                }, prom->get_future());
            hclib::async_await([=]() {
                    sleep(1);
		    int *n = new int(SIGNAL_VALUE);
                    prom->put(n);
                }, nullptr);
        });
    });
    printf("Exiting...\n");
    return 0;
}

