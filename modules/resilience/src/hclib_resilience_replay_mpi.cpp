
#define USE_RESILIENT_PROMISE

#include "hclib_resilience.h"
#include "hclib_resilience_replay_mpi.h"

namespace hclib {
namespace resilience{
namespace replay{

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    auto task_local = static_cast<replay_task_params_t<void*>*>(*hclib_get_curr_task_local());
    assert(is_replay_task(task_local));
    auto temp = new mpi_data(MPI_Send_lbl, buf, count, datatype, dest, tag, comm);
    task_local->mpi_send_vec->push_back(temp);
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    ::MPI_Recv(buf, count, datatype, source, tag, comm, status);
}

} // namespace replay
} // namespace resilience
} // namespace hclib
