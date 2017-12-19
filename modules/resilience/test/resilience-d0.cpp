
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <unistd.h>

namespace diamond = hclib::resilience::diamond;

class int_obj : public hclib::resilience::obj {
  public:
    int n;

    bool equals(obj* obj2) {
      return n == ((int_obj*)obj2)->n;
    }
};


int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system"};
    hclib::launch(deps, 1, [=]() {
        hclib::finish([=]() {
            int a_ref=10;
            hclib::promise_t<int*> *prom = new hclib::promise_t<int*>();
            diamond::promise_t<int_obj*> *prom1 = new diamond::promise_t<int_obj*>(1);

            hclib::async( [=]() {
                   sleep(1);
                   int *n= (int*) malloc(sizeof(int));
                   *n = SIGNAL_VALUE;
                   prom->put(n);
            });

	    hclib::async_await( [=]() {
                    int_obj *n2_tmp = prom1->get_future()->get();
                    printf("Value2 %d\n", n2_tmp->n);
            }, prom1->get_future());
           
            diamond::async_await_check( [=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value1 %d\n", *signal);

		    diamond::async_await( [=]() {
                        int_obj *n2 = new int_obj();
                        n2->n = 22;
		        prom1->put(n2);
		    }, nullptr);
            }, prom->get_future());
        });
    });
    printf("Exiting...\n");
    return 0;
}

