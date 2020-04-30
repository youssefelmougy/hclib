#define MPI_COMMUNICATION
#define USE_FENIX

#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <unistd.h>
#include <signal.h>

const int kKillID = 1;

class int_obj : public communication::obj {
  public:
    int n;
    int_obj() { printf("creating int_obj\n"); }
    ~int_obj() { printf("deleting int_obj\n"); }

    void deserialize(communication::archive_obj ar_ptr) {
        n = *(int*)(ar_ptr.data);
    }

    communication::archive_obj serialize() {
        communication::archive_obj ar_ptr;
        ar_ptr.size = sizeof(int);
        ar_ptr.data = malloc(ar_ptr.size);
        memcpy(ar_ptr.data, &n, ar_ptr.size);
        return ar_ptr;
    }
};

int main(int argc, char **argv) {
    const char *deps[] = { "system", "resilience" };
    communication::launch(deps, 2, [] () {

        int recovered = !communication::is_initial_role();
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD_RES, &rank);

        printf("in rank %d recovered %d pid %d\n", rank, recovered, getpid());

        if (rank == kKillID &&  recovered == 0) {
          pid_t pid = getpid();
          kill(pid, SIGTERM);
        }

        if(rank == kKillID) {
            auto prom_send = new hclib::promise_t<int_obj*>();
            int_obj *a = new int_obj();
            a->n = 23;
            printf("trying to send in rank %d\n", rank);
            communication::Isend(a, kKillID-1, 1, 0, prom_send, MPI_COMM_WORLD_RES);
            printf("trying to send1 in rank %d\n", rank);

            hclib::async_await( [=]() {
                printf("send %d in rank %d\n", prom_send->get_future()->get()->n, rank);
            }, prom_send->get_future());
        }
        if(rank == kKillID-1) {
            auto prom_recv = new hclib::promise_t<int_obj*>();
            printf("trying to recv in rank %d\n", rank);
            communication::Irecv((int)sizeof(int), kKillID, 1, prom_recv, nullptr, MPI_COMM_WORLD_RES);
            printf("trying to recv1 in rank %d\n", rank);
            hclib::async_await( [=]() {
                printf("recv %d in rank %d\n", prom_recv->get_future()->get()->n, rank);
            }, prom_recv->get_future());
        }
    });
    return 0;
}
