#define MPI_COMMUNICATION

#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <unistd.h>

namespace replay = hclib::resilience::replay;

enum TASK_STATE {NON_LEAF, LEAF};

int check(void *args) {
    return true;
}

int main(int argc, char ** argv) {
    int SIGNAL_VALUE = 42;
    const char *deps[] = { "system", "resilience"};
    //MPI_Init(&argc, &argv);
    //int provided;
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    //assert(provided == MPI_THREAD_MULTIPLE);
    hclib::launch(deps, 2, [=]() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        hclib::finish([=]() {
             auto prom = new hclib::promise_t<int*>();
             auto prom_res = new hclib::promise_t<int>();
             auto prom_allred1 = new hclib::promise_t<long*>();
             auto prom_allred2 = new hclib::promise_t<long*>();

            prom->put(new int(2));

            hclib::async_await( [=]() {
                    long *red = prom_allred1->get_future()->get();
                    printf("Reduce from replay task %ld\n", *red);
            }, prom_allred1->get_future());

            hclib::async_await( [=]() {
                    long *red = prom_allred2->get_future()->get();
                    printf("Reduce from non resilient task %ld\n", *red);
            }, prom_allred2->get_future());
 
            //This could be an array, vector, hash table or anything
            int *args = (int*)malloc(sizeof(int)*1);
            replay::async_await_check<NON_LEAF>( [=]() {
                    long *send_data = new long (10);
                    //reduction from replay task
                    communication::Iallreduce_tmp(send_data, MPI_LONG, MPI_SUM, 0, prom_allred1);
            }, prom_res, check, args, prom->get_future());

            long *send_data = new long (15);
            //reduction from non resilient task
            communication::Iallreduce_tmp(send_data, MPI_LONG, MPI_SUM, 0, prom_allred2);

            hclib::async_await( [=]() {
                    int res = prom_res->get_future()->get();
                    printf("result %d\n", res);
                    if(res == 0) exit(0);
            }, prom_res->get_future());
        });
    });
    printf("Going to Exiting\n");
    //MPI_Finalize();
    printf("Exiting...\n");
    return 0;
}

