
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <unistd.h>

namespace ref_count = hclib::ref_count;
namespace diamond = hclib::resilience::diamond;

enum TASK_STATE {NON_LEAF, LEAF};

class int_obj : public hclib::resilience::obj {
  public:
    int n;

    int_obj() {
      printf("Creating int_obj\n");
    }

    ~int_obj() {
        printf("Deleting int_obj\n");
    }

    bool equals(obj* obj2) {
      return n == ((int_obj*)obj2)->n;
    }
};

int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system"};
    hclib::launch(deps, 1, [=]() {
        hclib::finish([=]() {

            hclib::promise_t<int*> *prom = new hclib::promise_t<int*>();
            diamond::promise_t<int_obj*> *prom1 = new diamond::promise_t<int_obj*>(1);
            diamond::promise_t<int*> *prom2 = new diamond::promise_t<int*>(2, ref_count::DelType::ARR);
            hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();

            hclib::async_await( [=]() {
                   sleep(1);
                   int *n= new int[1];
                   n[0] = SIGNAL_VALUE;
                   prom->put(n);
		   prom2->put(n);
            }, future_nullptr);

	    hclib::ref_count::async_await( [=]() {
                    int_obj *n2_tmp = prom1->get_future()->get();
                    printf("Value2 %d\n", n2_tmp->n);
            }, prom1->get_future(), prom2->get_future());
           
            diamond::async_await_check<NON_LEAF, 3>( [=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value1 %d replica %d\n", *signal, diamond::get_index());

		    diamond::async_await( [=]() {
                        int_obj *n2 = new int_obj();
			n2->n = 22;
		        prom1->put(n2);
		    }, prom2->get_future());
            }, prom_res, prom->get_future());

            hclib::async_await( [=]() {
                int res = prom_res->get_future()->get();
                printf("result %d\n", res);
                if(res == 0) exit(0);
            }, prom_res->get_future());

        });
    });
    printf("Exiting...\n");
    return 0;
}

