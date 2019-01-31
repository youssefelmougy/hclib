#ifndef HCLIB_RESILIENCE_REPLAY_MPI_H
#define HCLIB_RESILIENCE_REPLAY_MPI_H

#include<mpi.h>

namespace hclib {
namespace resilience {
namespace replay {

/*
//int Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, int do_free, MPI_Comm comm=MPI_COMM_WORLD);
void Send(void *buf, int count, int dest, int tag, int do_free, MPI_Comm comm=MPI_COMM_WORLD);

//int Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Status *status, MPI_Comm comm=MPI_COMM_WORLD);
void Recv(void *buf, int count, int source, int tag, MPI_Status *status, MPI_Comm comm=MPI_COMM_WORLD);
*/

void Isend(void *buf, int count, int dest, int tag, int do_free, hclib::promise_t<void*> *prom, MPI_Comm comm=MPI_COMM_WORLD);

void Irecv(int count, int source, int tag, hclib::promise_t<void*> *prom, MPI_Comm comm=MPI_COMM_WORLD);

} // namespace replay
} // namespace resilience
} // namespace hclib

#endif

