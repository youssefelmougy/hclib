
#define USE_RESILIENT_PROMISE

#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include "hclib_resilience_replay_mpi.h"
#include <unistd.h>

namespace replay = hclib::resilience::replay;

enum TASK_STATE {NON_LEAF, LEAF};

//int_obj is not required for replay promises, base types can be used.
//Here it is just used to print inside constructor/destructor
class int_obj {
  public:
    int n;
    int_obj() { printf("creating int_obj\n"); }
    ~int_obj() { printf("deleting int_obj\n"); }
};

//set count to 0 or less to replay
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
    //MPI_Init(&argc, &argv);
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    assert(provided == MPI_THREAD_MULTIPLE);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    hclib::launch(deps, 1, [=]() {
        hclib::finish([=]() {
            hclib::promise_t<int*> *prom = new hclib::promise_t<int*>();
            replay::promise_t<int_obj*> *prom1 = new replay::promise_t<int_obj*>(1);
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
            replay::async_await_check<NON_LEAF>( [=]() {
                    int* signal = prom->get_future()->get();
                    assert(*signal == SIGNAL_VALUE);
                    printf("Value1 %d replay %d\n", *signal, replay::get_index());

                    replay::async_await( [=]() {
                        int_obj *n2 = new int_obj();
                        n2->n = 22;
                        prom1->put(n2);
                        *args = n2->n;
                        if(rank == 0) {
                            replay::MPI_Send(&(n2->n), 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
                        }
                    }, future_nullptr);

                    if(rank == 1) {
                        int buf;
                        replay::MPI_Recv(&buf, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        printf("Received %d in rank %d\n", buf, rank);
                    }
            }, prom_res, check, args, prom->get_future());

            hclib::async_await( [=]() {
                    int res = prom_res->get_future()->get();
                    printf("result %d\n", res);
                    if(res == 0) exit(0);
            }, prom_res->get_future());
        });
    });
    printf("Going to Exiting %d\n", rank); fflush(stdout);
    MPI_Finalize();
    printf("Exiting...\n");
    return 0;
}

