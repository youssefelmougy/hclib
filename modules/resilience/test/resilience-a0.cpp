
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <unistd.h>

namespace abft = hclib::resilience::abft;

//int_obj is not required for abft promises, base types can be used.
//Here it is just used to print inside constructor/destructor
class int_obj {
  public:
    int n;
    int_obj() { printf("creating int_obj\n"); }
    ~int_obj() { printf("deleting int_obj\n"); }
};

int check(void *args) {
    int *ptr = (int*)(args);
    if(*ptr == 22)
      return 1;
    else return 0;
}

int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system"};
    hclib::launch(deps, 1, [=]() {
        hclib::finish([=]() {
            hclib::promise_t<int*> *prom = new hclib::promise_t<int*>();
            abft::promise_t<int_obj*> *prom1 = new abft::promise_t<int_obj*>(1);
            hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();

            hclib::async_await( [=]() {
                   sleep(1);
                   int *n= (int*) malloc(sizeof(int));
                   *n = SIGNAL_VALUE;
                   prom->put(n);
            }, future_nullptr);

	    hclib::async_await( [=]() {
                    int_obj *n2_tmp = prom1->get_future()->get();
                    printf("Value2 %d\n", n2_tmp->n);
                    prom1->get_future()->release();
            }, prom1->get_future());

            //This could be an array, vector, hash table or anything
            int *args = (int*)malloc(sizeof(int)*1);
            abft::async_await_check( [=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value1 %d index %d\n", *signal, abft::get_index());

		    abft::async_await( [=]() {
                        int_obj *n2 = new int_obj();
                        n2->n = 21;
		        prom1->put(n2);
                        *args = n2->n;
		    }, future_nullptr);
            }, prom_res, check, args,
            [=]() {
                (*args)++;
                int_obj *n3 = new int_obj();
                n3->n = *args;
                prom1->put(n3);
            },
            prom->get_future());

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

