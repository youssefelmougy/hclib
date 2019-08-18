#ifndef HCLIB_RESILIENCE_H
#define HCLIB_RESILIENCE_H

#include "hclib_cpp.h"
#include "hclib_ref-count.h"

#define MPI_COMMUNICATION
#define USE_STD_VEC
//#define USE_C_ARRAY_WITH_ATOMIC
//#define USE_C_ARRAY_WITH_LOCK

#ifdef USE_STD_VEC
#include <vector>
#include <mutex>
#elif defined USE_C_ARRAY_WITH_ATOMIC
#include "hclib-atomics.h"
#elif defined USE_C_ARRAY_WITH_LOCK
#include <mutex>
#else
#error No data container for resilient promises specified
#endif

#ifdef MPI_COMMUNICATION
#include<mpi.h>
enum MPI_FUNC_LABELS {
    MPI_Send_lbl,
    MPI_Isend_lbl,
    MPI_Iallreduce_lbl
};

//Forward declare communication object and communication helpers
namespace hclib { namespace resilience { namespace communication {
    class obj;

    int Isend_helper(communication::obj* buf, MPI_Datatype datatype, int dest, int tag, hclib_promise_t *prom, MPI_Comm comm);
    int Iallreduce_helper(void *buf, MPI_Datatype datatype, int mpi_op, hclib_promise_t *prom, MPI_Comm comm);
};};};

//short name for communication namespace
namespace communication = hclib::resilience::communication;

struct mpi_data {
    MPI_FUNC_LABELS type;
    communication::obj *buf;
    MPI_Datatype datatype;
    int src_dest;
    int tag;
    int do_free;
    hclib_promise_t *prom;
    MPI_Comm comm;

    mpi_data(MPI_FUNC_LABELS type, communication::obj *buf, MPI_Datatype datatype, int src_dest, int tag, int do_free, hclib_promise_t *prom, MPI_Comm comm)
        :type(type), buf(buf),  datatype(datatype), src_dest(src_dest), tag(tag), do_free(do_free), prom(prom), comm(comm) {}

    inline int send() {
        //assert(type == MPI_Isend_lbl);
        //blocking APIs are used only for debugging non-error executions
        //if(type == MPI_Send_lbl)
        //    return ::MPI_Send(buf, count, datatype, src_dest, tag, comm);
        //else
        switch(type) {
            case MPI_Isend_lbl:
                return communication::Isend_helper(buf, datatype, src_dest, tag, prom, comm);
            case MPI_Iallreduce_lbl:
                return communication::Iallreduce_helper(buf, datatype, tag, prom, comm);
            default:
                assert(false);
        }
        return -1;
    }

    inline void tmp_delete() {
        if(do_free == 1) free(buf);
    }
};
#endif

//TODO: Assumes lambdas provided to tasks are non mutable.
//Need to find out what happens it not?

namespace hclib {
namespace resilience {

/*
Resilient object.
Users should extend this class to provide their own
equals method to be used for comparison at the end
of resilient task
*/
class obj {
  public:
    virtual ~obj() {}

    //return true if objs are equal else false
    virtual bool equals(obj* obj2) { assert(false); return false; }
};

namespace util {
/*
vector that allows concurrent insertion
Read operations are not concurrent with insertion
*/
template<typename T>
class safe_vector {
#ifdef USE_STD_VEC
    std::vector<T> vec;
    std::mutex mtx;

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
        std::lock_guard<std::mutex> lkg(mtx);
        vec.push_back(std::forward<T>(value));
    }

    void push_back(T& value) {
        std::lock_guard<std::mutex> lkg(mtx);
        vec.push_back(value);
    }

    size_t size() const noexcept {
        return vec.size();
    }

    T* data() noexcept {
        return vec.data();
    }

    void clear() noexcept {
        vec.clear();
    }

    auto begin() noexcept -> decltype(vec.begin()) {
        return vec.begin();
    }

    auto end() noexcept -> decltype(vec.end()) {
        return vec.end();
    }
#elif defined (USE_C_ARRAY_WITH_ATOMIC) || defined (USE_C_ARRAY_WITH_LOCK)
#ifndef SAFE_VECTOR_CAPACITY
#define SAFE_VECTOR_CAPACITY 128
#endif
    T vec[SAFE_VECTOR_CAPACITY];
    volatile int pos = -1;
#ifdef USE_C_ARRAY_WITH_LOCK
    std::mutex mtx;
#endif

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
#if defined (USE_C_ARRAY_WITH_LOCK)
	std::lock_guard<std::mutex> lkg(mtx);
        pos++;
        assert(pos < SAFE_VECTOR_CAPACITY);
        vec[pos] = value;
#elif defined (USE_C_ARRAY_WITH_ATOMIC)
	hc_atomic_inc(&pos);
	assert(pos < SAFE_VECTOR_CAPACITY);
	vec[pos] = std::forward<T>(value);
#endif
    }

    void push_back(T& value) {
#if defined (USE_C_ARRAY_WITH_LOCK)
	std::lock_guard<std::mutex> lkg(mtx);
        pos++;
        assert(pos < SAFE_VECTOR_CAPACITY);
        vec[pos] = value;
#elif defined (USE_C_ARRAY_WITH_ATOMIC)
	hc_atomic_inc(&pos);
	assert(pos < SAFE_VECTOR_CAPACITY);
	vec[pos] = value;
#endif
    }

    size_t size() const noexcept {
        return pos + 1;
    }

    T* data() noexcept {
        return &vec[0];
    }

    void clear() noexcept {
	// TODO: mutual exclusion?
        pos = -1;
    }

    auto begin() noexcept -> decltype(&vec[0]) {
        return &vec[0];
    }

    auto end() noexcept -> decltype(&vec[pos]) {
        return &vec[pos];
    }
#endif
};

template<typename T>
class safe_promise_vector : public safe_vector<T> {
  public:
    void do_puts(int index) {
        //perform the actual put
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++)
            data_var[i]->put_actual(index);
    }
};

template<typename T>
class safe_future_vector : public safe_vector<T> {
  public:
    void do_releases() {
        //perform the release
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++)
            data_var[i]->release();
    }
};

#ifdef MPI_COMMUNICATION
template<typename T>
class safe_mpi_data_vector : public safe_vector<T> {
  public:
    void do_sends() {
        //perform the actual send
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++) {
            data_var[i]->send();
            delete data_var[i];
        }
    }
};
#endif

} // namespace util

} // namespace resilience
} // namespace hclib

#include "hclib_resilience_promise.h"
#include "hclib_resilience_diamond.h"
#include "hclib_resilience_replay.h"
#include "hclib_resilience_abft.h"
#include "hclib_resilience_checkpoint.h"
#include "hclib_resilience_communication.h"

#endif
