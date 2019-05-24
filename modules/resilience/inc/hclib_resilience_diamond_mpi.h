#ifndef HCLIB_RESILIENCE_REPLAY_MPI_H
#define HCLIB_RESILIENCE_REPLAY_MPI_H

#include <mpi.h>
#include <hclib_resilience_checkpoint.h>
#include "hclib-module-common.h"

namespace checkpoint = hclib::resilience::checkpoint;

struct pending_mpi_op {
    MPI_Request req;
    hclib_promise_t *prom;
    communication_obj *data;
    checkpoint::archive_obj *serialize;
    pending_mpi_op *next;
};

pending_mpi_op *pending = nullptr;
hclib::locale_t *nic = nullptr;

bool test_mpi_completion(void *generic_op);

int Isend_helper(communication_obj *data, MPI_Datatype datatype, int dest, int tag, hclib_promise_t *prom, MPI_Comm comm);

namespace hclib {
namespace resilience {
namespace diamond {

/*
//int Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, int do_free, MPI_Comm comm=MPI_COMM_WORLD);
void Send(void *buf, int count, int dest, int tag, int do_free, MPI_Comm comm=MPI_COMM_WORLD);

//int Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Status *status, MPI_Comm comm=MPI_COMM_WORLD);
void Recv(void *buf, int count, int source, int tag, MPI_Status *status, MPI_Comm comm=MPI_COMM_WORLD);
*/

//Assumption: the serialized archive_obj can be deleted after data is send to remote node
template<class COMMUNICATION_OBJ>
void Isend(COMMUNICATION_OBJ *data, int dest, int tag, int do_free, hclib::promise_t<COMMUNICATION_OBJ*> *prom, MPI_Comm comm=MPI_COMM_WORLD) {
    auto task_local = static_cast<diamond_task_params_t<void*>*>(*hclib_get_curr_task_local());
    //assert(is_replay_task(task_local));
    if(is_resilient_task(task_local)) {
        auto temp = new mpi_data(MPI_Isend_lbl, data, MPI_BYTE, dest, tag, do_free, prom, comm);
        task_local->mpi_send_vec->push_back(temp);
    }
    else
        Isend_helper(data, MPI_BYTE, dest, tag, prom, comm);
}

template<class COMMUNICATION_OBJ>
void Irecv(int count, int source, int tag, hclib::promise_t<COMMUNICATION_OBJ*> *prom, MPI_Comm comm=MPI_COMM_WORLD) {
    if(resilience::get_index() == 0) {
        hclib::async_nb_await_at([=] {
            auto ar_ptr = new checkpoint::archive_obj();
            ar_ptr->size = count;
            ar_ptr->data = malloc(count);
            MPI_Request req;
            ::MPI_Irecv(ar_ptr->data, count, MPI_BYTE, source, tag, comm, &req);

            pending_mpi_op *op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
            assert(op);
            op->req = req;
            op->prom = prom;
            op->data = new COMMUNICATION_OBJ();
            op->serialize = ar_ptr;
            hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
        }, nullptr, nic);
    }
    else {
        //do nothing on replays since the mpi communication is already scheduled during first execution
    }
}

} // namespace diamond
} // namespace resilience
} // namespace hclib

#endif
