#ifndef HCLIB_RESILIENCE_REPLAY_MPI_H
#define HCLIB_RESILIENCE_REPLAY_MPI_H

#include<mpi.h>

namespace hclib {
namespace resilience {
namespace replay {

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

} // namespace replay
} // namespace resilience
} // namespace hclib

#endif

