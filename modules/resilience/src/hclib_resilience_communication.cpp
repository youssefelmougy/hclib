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
pending_mpi_op **completed = nullptr;
hclib::locale_t *nic = nullptr;

static int nic_locale_id;

#ifdef USE_FENIX
//TODO: see which of these variables can be made local
comm_op *comm_operations = nullptr;
bool just_recovered = false;
int old_rank = -1, new_rank = -1;
MPI_Comm new_comm;
static int fenix_status;
int callback_source = 0;
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

//TODO:Can we retain Fenix_Finalize in this module-finalize routine?
//#ifdef USE_FENIX
//    MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
//    Fenix_Finalize();
//#endif
    MPI_Finalize();
}

HCLIB_REGISTER_MODULE("resilience_mpi", mpi_pre_initialize, mpi_post_initialize, mpi_finalize)

bool test_mpi_completion(void *generic_op) {
    auto op = (pending_mpi_op *)generic_op;

    //callback_source = 0;
    int complete;
    ::MPI_Test(&op->req, &complete, MPI_STATUS_IGNORE);
#ifdef COMM_PROFILE
    test_count++;
#endif

#ifdef USE_FENIX
    //if just_recovered is true then re_enqueue has just happened,
    //which means the MPI_Test failed, in which case 'complete' might
    //be set to true instead of false due to failure. Since the
    //communication is not finished, we need to return false even it
    //'complete' is true. re_enqueue has re-executed the communication
    //operation related to op->req
    if(just_recovered) {
        just_recovered = false;
        return false;
    }
#endif // USE_FENIX

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
        //callback_source = 1;
        archive_obj ar_ptr = data->serialize();
        //TODO: should we use MPI_Isend or MPI_Issend?
        ::MPI_Isend(ar_ptr.data, ar_ptr.size, MPI_BYTE, dest, tag, comm, &req);
#ifdef COMM_PROFILE
        send_count++;
        send_size+=ar_ptr.size;
#endif

        op->serialized = ar_ptr;
        op->req = req;
        op->prom = prom;
        op->data = data;
#ifdef USE_FENIX
        op->msg_type = ISEND;
        op->neighbor = dest;
        op->tag = tag;
        op->comm  = comm;
#else
        op->serialized.size = 0;
#endif
        hclib::append_to_pending(op, &pending, completed, test_mpi_completion, nic);
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
        hclib::append_to_pending(op, &pending, completed, test_mpi_completion, nic);
    }, nullptr, nic);
    return 0;
}

#ifdef USE_FENIX

bool is_initial_role() {
    return fenix_status == FENIX_ROLE_INITIAL_RANK;
}

void Cart_create(MPI_Comm comm, int dim, int dims[], int periods[], int reorder, MPI_Comm *cart_comm) {
    comm_op *op = new comm_op(comm, dim, dims, periods, reorder, cart_comm);

    //TODO: For now we only allow one Cart_create before all MPI operations.
    //Allow multiple Cart_create in future.
    assert(comm_operations == nullptr);
    comm_operations = op;
    *cart_comm = comm;

    //MPI_Cart_create(comm, dim, dims, periods, reorder, cart_comm);
}

void Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest) {
    //TODO: For now we only allow one Cart_create before all MPI operations.
    //Allow multiple Cart_create in future.
    assert(comm == comm_operations->comm);

    int factor, thisdirection = 0, thisperiod = 0, ord;
    int srcord, destord, i, *d, *q;
    MPI_Comm_rank(comm, &ord);
    MPI_Comm_size(comm, &factor);

    d = comm_operations->dims;
    q = comm_operations->periods;

    for (i = 0; (i < comm_operations->ndims) && (i <= direction); ++i, ++d, ++q) {
        thisdirection = *d;
        thisperiod = *q;

        ord %= factor;
        factor /= thisdirection;
    }
    ord /= factor;

    *rank_source = *rank_dest = MPI_UNDEFINED;

    srcord = ord - disp;
    destord = ord + disp;

    if ( ((destord < 0) || (destord >= thisdirection)) && (!thisperiod) ) {
        *rank_dest = MPI_PROC_NULL;
    } else {
        destord %= thisdirection;
        if (destord < 0) destord += thisdirection;
        MPI_Comm_rank(comm, rank_dest);
        *rank_dest += ((destord - ord) * factor);
    }
    if ( ((srcord < 0) || (srcord >= thisdirection)) && (!thisperiod) ) {
        *rank_source = MPI_PROC_NULL;
    } else {
        srcord %= thisdirection;
        if (srcord < 0) srcord += thisdirection;
        MPI_Comm_rank(comm, rank_source);
        *rank_source += ((srcord - ord) * factor);
    }
}

void re_execute_comm(comm_op *comm_operations) {
    comm_op *op=comm_operations;
    while(op) {
        //printf("re_execute_comm in %d\n", new_rank);
        //TODO: Allow other communicators in future.
        op->comm = MPI_COMM_WORLD_DEFAULT;
        *(op->cart_comm) = MPI_COMM_WORLD_DEFAULT;

        //MPI_Cart_create(op->comm, op->dim, op->dims, op->periods, op->reorder, op->cart_comm);
        op = op->next;

        //TODO: Allow multiple communicator operations in future
        assert(op == nullptr);
    }
}

void re_enqueue(pending_mpi_op* list, bool isCompletedList, int fail_rank, pending_mpi_op** pending_ptr=nullptr) {
    pending_mpi_op *op = list;
    //printf("re_enqueue in %d with isCompletedList %d\n", new_rank, isCompletedList);
    while(op) {
        //printf("re_enqueue loop in %d with isCompletedList %d\n", new_rank, isCompletedList);
        pending_mpi_op *next = op->next;

        //if(op->neighbor == fails[0]) {
          //archive_obj ar_ptr = op->data->serialize();
          //op->serialized = ar_ptr;
          //TODO: use a seperate field to distinguish between send/recv rather
          //than using op->serialized.size() so that the previous op->serialized
          //data can be reused instead of calling op->data->serialize() again
          archive_obj ar_ptr = op->serialized;

          if(isCompletedList) {
            //if(op->neighbor == fail_rank) {
            //printf("re_enqueue out :isCompletedList in %d neighbor %d type %d tag %d\n", new_rank, op->neighbor, op->msg_type, op->tag);
            if(op->neighbor == fail_rank || op->msg_type == ISEND) {
                //printf("re_enqueue yes :isCompletedList in %d neighbor %d type %d tag %d\n", new_rank, op->neighbor, op->msg_type, op->tag);
                //put the op into pending list and re-enqueue the operation since it is to a failed rank
                //TODO:For now we enqueue sends to any rank. This should be removed.
                op->prom = nullptr;
                op->next = *pending_ptr;
                *pending_ptr = op;
            }
            else {
                //printf("re_enqueue no :isCompletedList in %d neighbor %d type %d tag %d\n", new_rank, op->neighbor, op->msg_type, op->tag);
                op = next;
                continue;;
            }
          }
          else {
            //if pending operation was completed to a non failed rank then no need to enqueue again
            int test_flag;
            //printf("try MPI_Test in\n");
            //MPI_Test( &(op->req), &test_flag, MPI_STATUS_IGNORE);
            //printf("try MPI_Test out\n");
            if(MPI_Test( &(op->req), &test_flag, MPI_STATUS_IGNORE) != FENIX_ERROR_CANCELLED
              && op->neighbor != fail_rank) {
                //TODO: for now we enqueue all sends. But once MPI_Issend completion can
                //be detected we can remove the check and ignore all completed operations
                //to non failed ranks.
                if(op->msg_type == IRECV){
                    //printf("re_enqueue no :MPI_Test in %d neighbor %d type %d tag %d\n", new_rank, op->neighbor, op->msg_type, op->tag);
                    op = next;
                    continue;
                }
            }
          }

          //TODO: map from old communicator to new communicator
          //rather than setting it to default COMM_WORLD
          //printf("rank %d %p %p\n", new_rank, op->comm, MPI_COMM_WORLD_DEFAULT);
          op->comm = MPI_COMM_WORLD_DEFAULT;

          //Enqueue incomplete operations using new communicator
          //Isend
          if(op->msg_type == ISEND) {
          //TODO: remove the additional serialization by adding a new field for send/recv to op
          //ar_ptr = op->data->serialize();
          //op->serialized = ar_ptr;
              //printf("re_enqueue yes :Isend in %d neighbor %d type %d tag %d\n", new_rank, op->neighbor, op->msg_type, op->tag);
              ::MPI_Isend(ar_ptr.data, ar_ptr.size, MPI_BYTE, op->neighbor, op->tag, op->comm, &(op->req));
          //op->serialized.size = 0;
          }
          //Irecv
          else if (op->msg_type == IRECV) {
              //printf("re_enqueue yes :Irecv in %d neighbor %d type %d tag %d\n", new_rank, op->neighbor, op->msg_type, op->tag);
              ::MPI_Irecv(ar_ptr.data, ar_ptr.size, MPI_BYTE, op->neighbor, op->tag, op->comm, &(op->req));
          }
          else {
              assert(false);
          }
        //}

        op = next;
    }
    //printf("re_enqueue out in %d\n", new_rank);
}

void my_recover_callback( MPI_Comm new_comm_world, int error, void *callback_data) {
    just_recovered = true;
    //printf("recovering in rank %d with callback_source %d\n", new_rank, *(int*)callback_data);

    //Recovering
    //else {
        int *fails, num_fails;
        num_fails = Fenix_Process_fail_list(&fails);
        //printf("num_fails %d\n", num_fails);fflush(stdout);
        assert(num_fails == 1);
#if COMM_PROFILE
        printf("Recovered rank %d, failed rank %d\n", new_rank, fails[0]);
#endif
        //Failed rank restarts from beginning
        //if(fails[0] == new_rank) {
        //    hclib::finish([=]() {
        //    completed = new pending_mpi_op*[new_world_size]{};
        //    hclib_async(fct_ptr, arg, NULL, 0, hclib_get_closest_locale());
        //    });
        //}

        //Non failed tasks rexecutes the pending communication to all tasks
        //and all communication to failed task
        //else {
            //check_out_finish_worker();
            //check_out_finish_worker();
            //hclib::finish([=]() {

            //re-execute communicator operations since
            //the old communicator is invalid
            re_execute_comm(comm_operations);

            //Iterate though the pending list and enqueue back the
            //incomplete communication operations with the
            //new communicator and status variable
            bool isCompletedList = false;
            re_enqueue(pending, isCompletedList, fails[0]);

            //add completed_ops to the failed rank in the pending ops list and re-enqueue them
            pending_mpi_op *completed_ops = completed[fails[0]];
            isCompletedList = true;
            re_enqueue(completed_ops, isCompletedList, fails[0], &pending);

            //Restart the communication worker
            //Remove one item and put it back to the queue
            //need to enqueue only if from the last MPI_Barrier
            //if(pending) {
            //    //pending_mpi_op *op = pending;
            //    //pending = pending->next;
            //    hclib::poll_on_pending(&pending, completed, test_mpi_completion, nic);
            //}

            //decrease the finish counter by one since the kill exited the task without cleanup
            //check_out_finish_worker();
            //});
        //}
        //hclib_end_finish();
    //}
    //printf("MPI_Barrier start new rank %d old rank %d\n", new_rank, old_rank);fflush(stdout);
    //MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
    //printf("MPI_Barrier end new rank %d old rank %d\n", new_rank, old_rank);fflush(stdout);
    //Fenix_Finalize();
    //printf("Fenix_Finalize end new rank %d old rank %d\n", new_rank, old_rank);fflush(stdout);
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
        //int old_rank = -1, new_rank = -1;
        int spare_ranks=1;
        MPI_Comm world_comm;
        MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
        MPI_Comm_size(world_comm, &old_world_size);
        MPI_Comm_rank(world_comm, &old_rank);

        MPI_Info info;
        MPI_Info_create(&info);
        MPI_Info_set(info, "FENIX_RESUME_MODE", "NO_JUMP");

        int error;
        Fenix_Init(&fenix_status, world_comm, &new_comm, nullptr, nullptr, spare_ranks, 0, info, &error);
        assert(fenix_status == FENIX_ROLE_INITIAL_RANK || fenix_status == FENIX_ROLE_RECOVERED_RANK);
        assert(fenix_status != FENIX_ROLE_SURVIVOR_RANK);
        MPI_Comm_rank(MPI_COMM_WORLD_DEFAULT, &new_rank);
        MPI_Comm_size(MPI_COMM_WORLD_DEFAULT, &new_world_size);

        Fenix_Callback_register(my_recover_callback, &callback_source);
#if COMM_PROFILE
        printf("init old rank: %d, new rank: %d\n", old_rank, new_rank);
#endif
    //} Fenix setup

    hclib::async_at([=](){
        printf("init old rank: %d, new rank: %d\n", old_rank, new_rank);
    //if(is_initial_role()) {
        completed = new pending_mpi_op*[new_world_size]{};
        //hclib_async(fct_ptr, arg, NULL, 0, hclib_get_closest_locale());

        //hclib::finish([=]() {
        //hclib_future_t *fut = hclib_async_future((future_fct_t)fct_ptr, arg, NULL, 0, hclib_get_closest_locale());
        //hclib::async_nb_await_at([=]() {
        //    printf("Barrier rank %d\n", new_rank);
        //    MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
        //}, fut, nic);
        //});

        //hclib::finish([=]() {
        //hclib::async_at([=]() {
            hclib::finish([=]() {
              fct_ptr(arg);
            });
            //Perform barrier in communication worker
            hclib::async_at([=]() {
                printf("invoking is_initial_role barrier in %d\n", new_rank);
                MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
                printf("done is_initial_role barrier in %d\n", new_rank);

                if(just_recovered == true) {
                    just_recovered = false;
                    if(pending) {
                        //pending_mpi_op *op = pending;
                        //pending = pending->next;
                        hclib::poll_on_pending(&pending, completed, test_mpi_completion, nic);
                    }
                    MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
                }
                Fenix_Finalize();
            }, nic);
        //}, hclib_get_closest_locale());
        //});
#if 0
    }
    else {
        completed = new pending_mpi_op*[new_world_size]{};
        hclib::finish([=]() {
            //hclib_async(fct_ptr, arg, NULL, 0, hclib_get_closest_locale());
            fct_ptr(arg);
        });
        hclib::async_at([=]() {
            printf("invoking next barrier in %d\n", new_rank);
            MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
            printf("done next barrier in %d\n", new_rank);
            //Fenix_Finalize();
        }, nic);
    }

    //Recovering
    else {
        int *fails, num_fails;
        num_fails = Fenix_Process_fail_list(&fails);
        printf("num_fails %d\n", num_fails);fflush(stdout);
        assert(num_fails == 1);
#if COMM_PROFILE
        printf("Recovered rank %d, failed rank %d\n", new_rank, fails[0]);
#endif
        //Failed rank restarts from beginning
        if(fails[0] == new_rank) {
            hclib::finish([=]() {
            completed = new pending_mpi_op*[new_world_size]{};
            hclib_async(fct_ptr, arg, NULL, 0, hclib_get_closest_locale());
            });
        }

        //Non failed tasks rexecutes the pending communication to all tasks
        //and all communication to failed task
        else {
            //check_out_finish_worker();
            check_out_finish_worker();
            hclib::finish([=]() {
            //Iterate though the pending list and enqueue back the
            //incomplete communication operations with the
            //new communicator and status variable
            bool isCompletedList = false;
            re_enqueue(pending, isCompletedList, fails[0]);

            //add completed_ops to the failed rank in the pending ops list and re-enqueue them
            pending_mpi_op *completed_ops = completed[fails[0]];
            isCompletedList = true;
            re_enqueue(completed_ops, isCompletedList, fails[0], &pending);

            //Restart the communication worker
            //Remove one item and put it back to the queue
            if(pending) {
                pending_mpi_op *op = pending;
                pending = pending->next;
                hclib::append_to_pending(op, &pending, completed, test_mpi_completion, nic);
            }

            //decrease the finish counter by one since the kill exited the task without cleanup
            //check_out_finish_worker();
            });
        }
        //hclib_end_finish();
    }
    printf("MPI_Barrier start new rank %d old rank %d\n", new_rank, old_rank);fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD_DEFAULT);
    printf("MPI_Barrier end new rank %d old rank %d\n", new_rank, old_rank);fflush(stdout);
    Fenix_Finalize();
    printf("Fenix_Finalize end new rank %d old rank %d\n", new_rank, old_rank);fflush(stdout);
#endif
    }, nic);

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

