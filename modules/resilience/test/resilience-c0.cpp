
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include "hclib_resilience_refstore.h"
#include <unistd.h>

namespace checkpoint = hclib::resilience::checkpoint;

enum TASK_STATE {NON_LEAF, LEAF};

class int_obj : public checkpoint::obj {
  public:

    int n;

    int_obj() { printf("creating int_obj\n"); }

    int_obj(checkpoint::archive_obj* ar_ptr) {
        n = *(int*)(ar_ptr->data);
    }

    ~int_obj() { printf("deleting int_obj\n"); }

    checkpoint::archive_obj* serialize() {
        auto ar_ptr = new checkpoint::archive_obj();
        ar_ptr->size = sizeof(int);
        ar_ptr->data = malloc(ar_ptr->size);
        memcpy(ar_ptr->data, &n, ar_ptr->size);
        return ar_ptr;
    }
};

//set count to 0 or less to checkpoint
int count = 1;

int check(void *args) {
    if(count <= 0 ) {
      count++;
      return 0;
    }
    int *ptr = (int*)(args);
    if(*ptr == 22)
      return 1;
    else return 0;
}

int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system"};
    hclib::launch(deps, 1, [=]() {

        checkpoint::set_archive_store(new checkpoint::ref_store());

        hclib::finish([=]() {
            hclib::promise_t<int*> *prom = new hclib::promise_t<int*>();
            checkpoint::promise_t<int_obj*> *prom1 = new checkpoint::promise_t<int_obj*>(2);
            hclib::promise_t<int*> *prom2 = new hclib::promise_t<int*>();
            hclib::promise_t<int>* prom_res1 = new hclib::promise_t<int>();
            hclib::promise_t<int>* prom_res2 = new hclib::promise_t<int>();

            hclib::async_await( [=]() {
                   sleep(1);
                   int *n= (int*) malloc(sizeof(int));
                   *n = SIGNAL_VALUE;
                   prom->put(n);
            }, future_nullptr);

            hclib::async_await( [=]() {
                int_obj *n2_tmp = prom1->get_future()->get();

                //introduce error
                n2_tmp->n++;

                prom2->put(&(n2_tmp->n));
                prom1->get_future()->release();
            }, prom1->get_future());

            //This could be an array, vector, hash table or anything
            int *args1 = (int*)malloc(sizeof(int)*1);
            checkpoint::async_await_check<LEAF,2>( [=]() {
                int* signal = prom->get_future()->get();
                assert(*signal == SIGNAL_VALUE);
                printf("Value1 %d replay %d\n", *signal, checkpoint::get_index());

                int_obj *n2 = new int_obj();
                n2->n = 22;
                prom1->put(n2);
                *args1 = n2->n;
            }, prom_res1, check, args1, prom->get_future());

            int *args2 = (int*)malloc(sizeof(int)*1);
            checkpoint::async_await_check( [=]() {
                int_obj *n2_tmp1 = prom1->get_future()->get();
                int *n2_tmp2 = prom2->get_future()->get();
                printf("Value2 tmp1 %d tmp2 %d\n", n2_tmp1->n, *n2_tmp2);
                *args2 = n2_tmp1->n;
            }, prom_res2, check, args2, prom1->get_future(), prom2->get_future());
 
            hclib::async_await( [=]() {
                int res = prom_res1->get_future()->get();
                printf("result1 %d\n", res);
                if(res == 0) exit(0);
            }, prom_res1->get_future());

            hclib::async_await( [=]() {
                int res = prom_res2->get_future()->get();
                printf("result2 %d\n", res);
                if(res == 0) exit(0);
            }, prom_res2->get_future());
 
        });

        delete(checkpoint::get_archive_store());

    });
    printf("Exiting...\n");
    return 0;
}

