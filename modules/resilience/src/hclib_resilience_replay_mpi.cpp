
#define USE_RESILIENT_PROMISE

#include "hclib_resilience.h"
#include "hclib_resilience_replay_mpi.h"
#include "hclib-module.h"
#include "hclib-module-common.h"

struct pending_mpi_op {
    MPI_Request req;
    hclib::promise_t<void*> *prom;
    void *buf;
    pending_mpi_op *next;
};

pending_mpi_op *pending = NULL;

static int nic_locale_id;
static hclib::locale_t *nic = NULL;

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

static bool test_mpi_completion(void *generic_op) {
    pending_mpi_op *op = (pending_mpi_op *)generic_op;

    int complete;
    ::MPI_Test(&op->req, &complete, MPI_STATUS_IGNORE);

    if (complete) {
        return true;
    } else {
        return false;
    }
}

int Isend_helper(void *buf, int count, MPI_Datatype datatype, int dest, int tag, hclib::promise_t<void*> *prom, MPI_Comm comm) {
    hclib::async_nb_await_at([=] {
        MPI_Request req;
        ::MPI_Isend(buf, count, MPI_BYTE, dest, tag, comm, &req);

        pending_mpi_op *op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
        assert(op);
        op->req = req;
        op->prom = prom;
        op->buf = buf;
        hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
    }, nullptr, nic);
    return 0;
}

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

void Isend(void *buf, int count, int dest, int tag, int do_free, hclib::promise_t<void*> *prom, MPI_Comm comm) {
    auto task_local = static_cast<replay_task_params_t<void*>*>(*hclib_get_curr_task_local());
    assert(is_replay_task(task_local));
    auto temp = new mpi_data(MPI_Isend_lbl, buf, count, MPI_BYTE, dest, tag, do_free, prom, comm);
    task_local->mpi_send_vec->push_back(temp);
}

void Irecv(int count, int source, int tag, hclib::promise_t<void*> *prom, MPI_Comm comm) {
    if(resilience::get_index() == 0) {
        hclib::async_nb_await_at([=] {
            void *buf = (void*)malloc(count);
            MPI_Request req;
            ::MPI_Irecv(buf, count, MPI_BYTE, source, tag, comm, &req);

            pending_mpi_op *op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
            assert(op);
            op->req = req;
            op->prom = prom;
            op->buf = buf;
            hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
        }, nullptr, nic);
    }
    else {
        //do nothing on replays since the mpi communication is already scheduled during first execution
    }
}

} // namespace replay
} // namespace resilience
} // namespace hclib
