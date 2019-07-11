
#define USE_RESILIENT_PROMISE

#include "hclib_resilience.h"
#include "hclib_resilience_replay_mpi.h"
#include "hclib-module.h"

static int nic_locale_id;

HCLIB_MODULE_INITIALIZATION_FUNC(mpi_pre_initialize) {
    nic_locale_id = hclib_add_known_locale_type("Interconnect");
}

HCLIB_MODULE_INITIALIZATION_FUNC(mpi_post_initialize) {
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    assert(provided == MPI_THREAD_FUNNELED);

    int n_nics;
    hclib::locale_t **nics = hclib::get_all_locales_of_type(nic_locale_id,
            &n_nics);
    HASSERT(n_nics == 1);
    HASSERT(nics);
    HASSERT(nic == NULL);
    nic = nics[0];

    hclib_locale_mark_special(nic, "COMM");
}

HCLIB_MODULE_INITIALIZATION_FUNC(mpi_finalize) {
    MPI_Finalize();
}

HCLIB_REGISTER_MODULE("resilience_mpi", mpi_pre_initialize, mpi_post_initialize, mpi_finalize)

bool test_mpi_completion(void *generic_op) {
    pending_mpi_op *op = (pending_mpi_op *)generic_op;

    int complete;
    ::MPI_Test(&op->req, &complete, MPI_STATUS_IGNORE);

    if (complete) {
        return true;
    } else {
        return false;
    }
}

int Isend_helper(communication_obj *data, MPI_Datatype datatype, int dest, int tag, hclib_promise_t *prom, MPI_Comm comm) {
    hclib::async_nb_await_at([=] {
        MPI_Request req;
        pending_mpi_op *op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
        assert(op);
        //TODO: since serialized object it not saved, they cannot be deleted
        //Need to improve garbage collection
        op->serialize = nullptr;
        auto ar_ptr = data->serialize();
        ::MPI_Isend(ar_ptr->data, ar_ptr->size, MPI_BYTE, dest, tag, comm, &req);

        op->req = req;
        op->prom = prom;
        op->data = data;
        hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
    }, nullptr, nic);
    return 0;
}

int Iallreduce_helper(void *data, MPI_Datatype datatype, int mpi_op, hclib_promise_t *prom, MPI_Comm comm) {
     hclib::async_nb_await_at([=] {
        MPI_Request req;
        pending_mpi_op *op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
        assert(op);
        op->serialize = nullptr;
        void *recv_data = malloc(sizeof(double));
        ::MPI_Iallreduce(data, recv_data, 1, datatype, mpi_op, comm, &req);

        op->req = req;
        op->prom = prom;
        op->data = (communication_obj*)recv_data;
        hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
    }, nullptr, nic);
    return 0;
}

#if 0
namespace hclib {
namespace resilience{
namespace replay{

/*
//blocking APIs are used only for debugging non-error executions

void Send(void *buf, int count, int dest, int tag, int do_free, MPI_Comm comm) {
    auto task_local = static_cast<replay_task_params_t<void*>*>(*hclib_get_curr_task_local());
    assert(is_replay_task(task_local));
    auto temp = new mpi_data(MPI_Send_lbl, buf, count, MPI_BYTE, dest, tag, do_free, nullptr, comm);
    task_local->mpi_send_vec->push_back(temp);
}

void Recv(void *buf, int count, int source, int tag, MPI_Status *status, MPI_Comm comm) {
    if(resilience::get_index() == 0) {
        ::MPI_Recv(buf, count, MPI_BYTE, source, tag, comm, status);
    }
    else
    {
        assert(false);
    }
}
*/


} // namespace replay
} // namespace resilience
} // namespace hclib
#endif

