
#ifndef HCLIB_RESILIENCE_COMMUNICATION_H
#define HCLIB_RESILIENCE_COMMUNICATION_H

#ifdef MPI_COMMUNICATION

#include "hclib-module-common.h"

namespace hclib {
namespace resilience {
namespace communication {

#ifdef USE_FENIX

extern MPI_Comm new_comm;
#define MPI_COMM_WORLD_RES communication::new_comm

bool is_initial_role();
void hclib_launch(generic_frame_ptr fct_ptr, void *arg, const char **deps, int ndeps);

template <typename T>
inline void launch(const char **deps, int ndeps, T &&lambda) {
    typedef typename std::remove_reference<T>::type U;
    hclib_task_t *user_task = _allocate_async(new U(lambda));
    hclib::resilience::communication::hclib_launch((generic_frame_ptr)spawn, user_task, deps, ndeps);
}

#endif //USE_FENIX

/*
The object that will be saved to the archive.
*/
struct archive_obj {

    //size of object blob
    int size = 0;
    //delete data after usage is over
    bool do_free = 0;
    //blob of object to be archived
    void *data = nullptr;
};

/*
The user object should include a serialize method that
will convert the current objct to an archive object
and a constructor that can create the current object
baack from the archive object
*/
class obj: public hclib::resilience::obj {
  public:
    //obj(archive_obj* ptr) { assert(false); }
    virtual void* allocate_buffer(size_t) { return nullptr; }
    virtual void deserialize(archive_obj) { assert(false); }
    virtual archive_obj serialize() { assert(false); return archive_obj(); }
};

struct pending_mpi_op {
    MPI_Request req;
    hclib_promise_t *prom;
    obj *data;
    archive_obj serialized;
#ifdef USE_FENIX
    int neighbor;
    int64_t tag;
    MPI_Comm comm;
#endif
    pending_mpi_op *next;
};

extern pending_mpi_op *pending;
extern hclib::locale_t *nic;

bool test_mpi_completion(void *generic_op);

int Isend_helper(obj *data, MPI_Datatype datatype, int dest, int64_t tag, hclib_promise_t *prom, MPI_Comm comm);
int Iallreduce_helper(obj *data, MPI_Datatype datatype, int64_t mpi_op, hclib_promise_t *prom, MPI_Comm comm);

//Assumption: the serialized archive_obj can be deleted after data is send to remote node
template<class COMMUNICATION_OBJ>
void Isend(COMMUNICATION_OBJ *data, int dest, int64_t tag, int do_free, hclib::promise_t<COMMUNICATION_OBJ*> *prom, MPI_Comm comm=MPI_COMM_WORLD) {
    auto task_local = static_cast<resilient_task_params_t<void*>*>(*hclib_get_curr_task_local());
    //assert(is_replay_task(task_local));
    if(is_resilient_task(task_local)) {
        auto temp = new mpi_data(MPI_Isend_lbl, data, MPI_BYTE, dest, tag, do_free, prom, comm);
        task_local->mpi_send_vec->push_back(temp);
    }
    else
        Isend_helper(data, MPI_BYTE, dest, tag, prom, comm);
}

#ifdef COMM_PROFILE
extern int64_t recv_count, recv_size;
#endif

template<class COMMUNICATION_OBJ>
void Irecv(int count, int source, int64_t tag, hclib::promise_t<COMMUNICATION_OBJ*> *prom, hclib_future_t *fut = nullptr, MPI_Comm comm=MPI_COMM_WORLD) {

    assert(count>0);
    if(resilience::get_index() == 0) {
        hclib::async_nb_await_at([=] {
            archive_obj ar_ptr;
            ar_ptr.size = count;

            pending_mpi_op *op = (pending_mpi_op *)malloc(sizeof(pending_mpi_op));
            assert(op);
            op->data = new COMMUNICATION_OBJ();
            ar_ptr.data = op->data->allocate_buffer(count);
            if(ar_ptr.data == nullptr)
                ar_ptr.data = malloc(count);

            MPI_Request req;
            ::MPI_Irecv(ar_ptr.data, count, MPI_BYTE, source, tag, comm, &req);
#ifdef COMM_PROFILE
            recv_count++;
            recv_size+=count;
#endif

            op->req = req;
            op->prom = prom;
            op->serialized = ar_ptr;
#ifdef USE_FENIX
            op->neighbor = source;
            op->tag = tag;
            op->comm = comm;
#endif
            hclib::append_to_pending(op, &pending, test_mpi_completion, nic);
        }, fut, nic);
    }
    else {
        //do nothing on replays since the mpi communication is already scheduled during first execution
    }
}

template<class COMMUNICATION_OBJ>
void Iallreduce_tmp(void *data, MPI_Datatype datatype, MPI_Op op, int do_free, hclib::promise_t<COMMUNICATION_OBJ*> *prom, MPI_Comm comm=MPI_COMM_WORLD) {
    assert(sizeof(MPI_Op) <= sizeof(int64_t));
    auto task_local = static_cast<replay_task_params_t<void*>*>(*hclib_get_curr_task_local());
    //assert(is_replay_task(task_local));
    if(is_resilient_task(task_local)) {
        auto temp = new mpi_data(MPI_Iallreduce_lbl, (obj*)data, datatype, -1, (int64_t)op, do_free, prom, comm);
        task_local->mpi_send_vec->push_back(temp);
    }
    else
        Iallreduce_helper(data, datatype, (int64_t)op, prom, comm);
}

} // namespace communication
} // namespace resilience
} // namespace hclib

#endif // MPI_COMMUNICATION

#endif
