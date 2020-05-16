
#ifndef MPI_COMMUNICATION
#define MPI_COMMUNICATION
#endif

#include "hclib_resilience.h"
#include "hclib-module.h"

#ifdef USE_FENIX
#include<fenix.h>
#endif

namespace hclib {
namespace resilience{
namespace communication {

pending_mpi_op *pending = nullptr;
hclib::locale_t *nic = nullptr;

static int nic_locale_id;

#ifdef USE_FENIX
MPI_Comm new_comm;
static int fenix_status;
#endif

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

#ifdef COMM_PROFILE
int64_t test_count=0,send_count=0,recv_count=0,send_size=0,recv_size=0;
#endif

HCLIB_MODULE_INITIALIZATION_FUNC(mpi_finalize) {
#ifdef COMM_PROFILE
    printf("test %ld, send/recv count %ld / %ld, send/recv bytes %ld / %ld\n",test_count, send_count, recv_count, send_size, recv_size);
#endif
#ifdef USE_FENIX
    Fenix_Finalize();
#endif
    MPI_Finalize();
}

HCLIB_REGISTER_MODULE("resilience_mpi", mpi_pre_initialize, mpi_post_initialize, mpi_finalize)

bool test_mpi_completion(void *generic_op) {
    auto op = (pending_mpi_op *)generic_op;

    int complete;
    ::MPI_Test(&op->req, &complete, MPI_STATUS_IGNORE);
#ifdef COMM_PROFILE
    test_count++;
#endif

    if (complete) {
        return true;
    } else {
        return false;
    }
}

int Isend_helper(obj *data, MPI_Datatype datatype, int dest, int64_t tag, hclib_promise_t *prom, MPI_Comm comm) {
    hclib::async_nb_await_at([=] {
        MPI_Request req;
        auto op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
        assert(op);
        archive_obj ar_ptr = data->serialize();
        ::MPI_Isend(ar_ptr.data, ar_ptr.size, MPI_BYTE, dest, tag, comm, &req);
#ifdef COMM_PROFILE
        send_count++;
        send_size+=ar_ptr.size;
#endif

        op->serialized = ar_ptr;
        op->serialized.size = 0;
        op->req = req;
        op->prom = prom;
        op->data = data;
#ifdef USE_FENIX
        op->neighbor = dest;
        op->tag = tag;
        op->comm  = comm;
#endif
        hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
    }, nullptr, nic);
    return 0;
}

int Iallreduce_helper(void *data, MPI_Datatype datatype, int64_t mpi_op, hclib_promise_t *prom, MPI_Comm comm) {
     hclib::async_nb_await_at([=] {
        MPI_Request req;
        auto op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
        assert(op);
        op->serialized.size = 0;
        op->serialized.do_free = 0;
        void *recv_data = malloc(sizeof(double));
        ::MPI_Iallreduce(data, recv_data, 1, datatype, (MPI_Op)mpi_op, comm, &req);

        op->req = req;
        op->prom = prom;
        op->data = (obj*)recv_data;
        hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
    }, nullptr, nic);
    return 0;
}

#ifdef USE_FENIX

bool is_initial_role() {
    return fenix_status == FENIX_ROLE_INITIAL_RANK;
}

void hclib_launch(generic_frame_ptr fct_ptr, void *arg, const char **deps,
        int ndeps) {
    unsigned long long start_time = 0;
    unsigned long long end_time;

    const int instrument = (getenv("HCLIB_INSTRUMENT") != NULL);
    const int profile_launch_body = (getenv("HCLIB_PROFILE_LAUNCH_BODY") != NULL);

    hclib_init(deps, ndeps, instrument);

    if (profile_launch_body) {
        start_time = hclib_current_time_ns();
    }

    //Fenix setup {
        int old_world_size = -1, new_world_size = -1;
        int old_rank = -1, new_rank = -1;
        int spare_ranks=1;
        MPI_Comm world_comm;
        MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
        MPI_Comm_size(world_comm, &old_world_size);
        MPI_Comm_rank(world_comm, &old_rank);

        int error;
        Fenix_Init(&fenix_status, world_comm, &new_comm, nullptr, nullptr, spare_ranks, 0, MPI_INFO_NULL, &error);
        MPI_Comm_rank(new_comm, &new_rank);
#if COMM_PROFILE
        printf("init old rank: %d, new rank: %d\n", old_rank, new_rank);
#endif
    //} Fenix setup

    if(is_initial_role()) {
        //hclib_async(fct_ptr, arg, NULL, 0, hclib_get_closest_locale());

        //hclib_future_t *fut = hclib_async_future((future_fct_t)fct_ptr, arg, NULL, 0, hclib_get_closest_locale());
        //hclib::async_nb_await_at([=]() {
        //    printf("Barrier rank %d\n", new_rank);
        //    MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
        //}, fut, nic);

        hclib::async_at([=]() {
            hclib::finish([=]() {
              fct_ptr(arg);
            });
            //Perform barrier in communication worker
            hclib::async_at([=]() {
                MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
            }, nic);
        }, hclib_get_closest_locale());
    }

    //Recovering
    else {
        int *fails, num_fails;
        num_fails = Fenix_Process_fail_list(&fails);
        assert(num_fails == 1);
#if COMM_PROFILE
        printf("Recovered rank %d, failed rank %d\n", new_rank, fails[0]);
#endif
        //Failed rank restarts from beginning
        if(fails[0] == new_rank)
            hclib_async(fct_ptr, arg, NULL, 0, hclib_get_closest_locale());
        //Non failed tasks rexecutes the pending communication to all tasks
        //and all communication to failed task
        else {
            //Iterate though the pending list and enqueue back the communication operations
            //with the new status variable
            //TODO: save all completed pending ops and reuse them here to failed task
            pending_mpi_op *op = pending;
            while(op) {
                pending_mpi_op *next = op->next;

                //if(op->neighbor == fails[0]) {
                  archive_obj ar_ptr = op->serialized;

                  //TODO: map from old communicator to new communicator
                  //rather than setting it to default COMM_WORLD
                  //printf("rank %d %p %p\n", new_rank, op->comm, MPI_COMM_WORLD_DEFAULT);
                  op->comm = MPI_COMM_WORLD_DEFAULT;

                  //Isend
                  if(op->serialized.size == 0) {
                      ::MPI_Isend(ar_ptr.data, ar_ptr.size, MPI_BYTE, op->neighbor, op->tag, op->comm, &(op->req));
                  }
                  //Irecv
                  else if (op->serialized.size > 0) {
                      ::MPI_Irecv(ar_ptr.data, ar_ptr.size, MPI_BYTE, op->neighbor, op->tag, op->comm, &(op->req));
                  }
                  else {
                      assert(false);
                  }
                //}

                op = next;
            }

            //Restart the communication worker
            //Remove one item and put it back to the queue
            if(pending) {
                pending_mpi_op *op = pending;
                pending = pending->next;
                hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
            }

            //decrease the finish counter by one since the kill exited the task without cleanup
            check_out_finish_worker();
        }
    }

    hclib_finalize(instrument);
    if (profile_launch_body) {
        end_time = hclib_current_time_ns();
        printf("\nHCLIB TIME %llu ns\n", end_time - start_time);
    }
}

#endif // USE_FENIX

} // namespace communication
} // namespace resilience
} // namespace hclib

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

