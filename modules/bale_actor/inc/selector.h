
#include "safe_buffer.h"
#include "hclib_bale_actor.h"
extern "C" {
#include "convey.h"
}

#define DONE_MARK -1
#define BUFFER_SIZE 10000001
#ifndef ELASTIC_BUFFER_SIZE
#define ELASTIC_BUFFER_SIZE 48
#endif

namespace hclib {

#ifdef USE_LAMBDA
class BaseLambdaPacket {
  public:
    virtual void invoke() = 0;
    virtual size_t get_bytes() = 0;
    virtual ~BaseLambdaPacket() {}
};

template<typename L>
class LambdaPacket : public BaseLambdaPacket {
    L lambda;

  public:
    LambdaPacket(L lambda) : lambda(lambda) {}

    void invoke() {
        lambda();
    }

    size_t get_bytes() {
        return sizeof(*this);
    }
};

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;
    BaseLambdaPacket* lambda_pkt;

    BufferPacket() : rank(0) {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
    BufferPacket(int64_t rank, BaseLambdaPacket* lambda_pkt) : rank(rank), lambda_pkt(lambda_pkt) {}
};

#else // USE_LAMBDA

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;

    BufferPacket() : rank(0) {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
};

#endif // USE_LAMBDA

template<typename T, int SIZE>
class Mailbox {

    hclib::conveyor::safe_buffer<BufferPacket<T>> *buff=nullptr;
    convey_t* conv=nullptr;
    BufferPacket<T> done_mark;
    hclib::promise_t<int> worker_loop_end;

  public:

    Mailbox() {
        buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
#ifdef SELECTOR_DEBUG
        printf("Creating Mailbox\n");
#endif
    }

    ~Mailbox() {
#ifdef SELECTOR_DEBUG
        printf("Deleting Mailbox\n");
#endif
        //delete buff;
        //convey_free(conv);
    }

    std::function<void (T, int)> process;

    hclib::future_t<int>* get_worker_loop_finish() {
        return worker_loop_end.get_future();
    }

    void start() {
#ifdef USE_LAMBDA
        conv = convey_new_elastic(1, ELASTIC_BUFFER_SIZE, SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
#else
        //buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
        //conv = convey_new(SIZE_MAX, 0, NULL, 0);
        conv = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
#endif
        assert( conv != nullptr );
        convey_begin(conv, sizeof(T));
        done_mark.rank = DONE_MARK;
    }

    void end() {
        delete buff;
        convey_free(conv);
    }

#ifdef USE_LAMBDA
    template<typename L>
    void send(int rank, L lambda) {
        assert(!buff->full());
        buff->push_back(BufferPacket<T>(rank, new LambdaPacket<L>(lambda)));
    }
#else
    void send(T pkt, int rank) {
        assert(!buff->full());
        buff->push_back(BufferPacket<T>(pkt, rank));
    }
#endif // USE_LAMBDA

    void done() {
        assert(!buff->full());
        buff->push_back(done_mark);
    }

    int start_worker_loop(int status=0) {

#ifndef YIELD_LOOP
        assert(status == 0);
        hclib::async_at([=]{
#endif
          while(true) {
              size_t buff_size = buff->size();
              if(buff_size > 0) break;
#ifndef YIELD_LOOP
              hclib::yield_at(nic);
#else
              if(status == 1)
                  return 1;
              else {
                  assert(status == 2);
                  break;
              }
#endif
          }

          BufferPacket<T> bp;
          if(buff->size() > 0)
            bp = buff->at(0);

          //Assumes once 'advance' is called with done=true, the conveyor
          //enters endgame and subsequent value of 'done' is ignored
          while(convey_advance(conv, bp.rank == DONE_MARK)) {
              int i;
              size_t buff_size = buff->size();
              for(i=0;i<buff_size; i++){
                  bp = buff->operator[](i);
                  if( bp.rank == DONE_MARK) break;
#ifdef USE_LAMBDA
                  //printf("size %d\n", bp.lambda_pkt->get_bytes());
                  if( !convey_epush(conv, bp.lambda_pkt->get_bytes(), bp.lambda_pkt, bp.rank)) break;
                  //delete bp.lambda_pkt;
#else
                  if( !convey_push(conv, &(bp.data), bp.rank)) break;
#endif //  USE_LAMBDA
              }

	          if(i>0)
              {
#ifdef USE_LOCK
                  std::lock_guard<std::mutex> lg(buff->get_mutex());
#endif
                  buff->erase_begin(i);
              }
              int64_t from;
#ifdef USE_LAMBDA
              convey_item_t item;
              while( convey_epull(conv, &item) == convey_OK) {
                  BaseLambdaPacket* data = (BaseLambdaPacket*)item.data;
                  data->invoke();
              }
#else

              T pop;
              while( convey_pull(conv, &pop, &from) == convey_OK) {
                  //hclib::async([=]() { process(pop, from); });
                  process(pop, from);
              }
#endif // USE_LAMBDA

#ifndef YIELD_LOOP
              hclib::yield_at(nic);
#else
              return 2;
#endif
          }
          worker_loop_end.put(1);
#ifndef YIELD_LOOP
        }, nic);
#endif
        return 0;
    }
};

template<int N, typename T=int64_t, int SIZE=BUFFER_SIZE>
class Selector {

  private:
#ifndef YIELD_LOOP
    void start_worker_loop() {
        for(int i=0;i<N;i++) {
            mb[i].start_worker_loop();
        }
    }
#else
    void start_worker_loop() {
        hclib::async_at([=]{
            int loop_stat[N];
            std::fill_n(loop_stat, N, 1);
            int finish_count = 0;

            while(finish_count < N) {
              for(int i=0;i<N;i++) {
                if(loop_stat[i] != 0) {
                  loop_stat[i] = mb[i].start_worker_loop(loop_stat[i]);
                  if(loop_stat[i] == 0)
                    finish_count++;
                }
              }
              hclib::yield_at(nic);
            }
        }, nic);
    }
#endif

    hclib::promise_t<int> end_prom;
    int num_work_loop_end = 0;

  protected:

  public:

    Mailbox<T, SIZE> mb[N];

    Selector(bool is_start = false) {
        if(is_start) {
            start();
        }
    }

    ~Selector() {
        for(int i=0; i<N; i++) {
            mb[i].end();
        }
    }

    void start() {
        for(int i=0; i<N; i++) {
            mb[i].start();
        }
        start_worker_loop();
    }

#ifdef USE_LAMBDA
    template<typename L>
    void send(int mb_id, int rank, L lambda) {
        mb[mb_id].send(rank, lambda);
    }

    template<typename L>
    void send(int rank, L lambda) {
        assert(N==1);
        send(0, rank, lambda);
    }
#else
    void send(int mb_id, T pkt, int rank) {
        mb[mb_id].send(pkt, rank);
    }

    void send(T pkt, int rank) {
        assert(N==1);
        send(0, pkt, rank);
    }
#endif // USE_LAMBDA

    void done(int mb_id) {
        mb[mb_id].done();
        hclib::async_await_at([=]() {
            num_work_loop_end++;
            if(num_work_loop_end < N) {
                done((mb_id+1)%N);
            }
            else {
                assert(num_work_loop_end == N);
                end_prom.put(1);
            }
        }, mb[mb_id].get_worker_loop_finish(), nic);
    }

    void done() {
        assert(N==1);
        done(0);
    }

    hclib::future_t<int>* get_future() {
        return end_prom.get_future();
    }
};

template<typename T=int64_t>
using Actor = Selector<1,T>;

}; // namespace hclib


