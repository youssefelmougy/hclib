
#include "safe_buffer.h"
#include "hclib_cpp.h"
extern "C" {
#include "convey.h"
}

#define DONE_MARK -1
#define BUFFER_SIZE 1000001

namespace hclib {

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;

    BufferPacket() : rank(0) {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
};

template<typename T, int SIZE>
class Mailbox {

    hclib::conveyor::safe_buffer<BufferPacket<T>> *buff=nullptr;
    convey_t* conv=nullptr;
    BufferPacket<T> done_mark;
    hclib::promise_t<int> worker_loop_end;

  public:

    Mailbox() {
        //buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
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
        buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
        //conv = convey_new(SIZE_MAX, 0, NULL, 0);
        conv = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
        assert( conv != nullptr );
        convey_begin(conv, sizeof(T));
        done_mark.rank = DONE_MARK;
    }

    void end() {
        delete buff;
        convey_free(conv);
    }

    void send(T pkt, int rank) {
        buff->push_back(BufferPacket<T>(pkt, rank));
    }

    void done() {
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
                  if( !convey_push(conv, &(bp.data), bp.rank)) break;
              }

	          if(i>0)
              {
#ifdef USE_LOCK
                  std::lock_guard<std::mutex> lg(buff->get_mutex());
#endif
                  buff->erase_begin(i);
              }
              T pop;
              int64_t from;
              while( convey_pull(conv, &pop, &from) == convey_OK) {
                  //hclib::async([=]() { process(pop, from); });
                  process(pop, from);
              }
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

template<int N, typename T, int SIZE=BUFFER_SIZE>
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

    void send(int mb_id, T pkt, int rank) {
        mb[mb_id].send(pkt, rank);
    }

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

    hclib::future_t<int>* get_future() {
        return end_prom.get_future();
    }
};

namespace selector {

template<typename S, typename T>
void finish(S slr, T lambda) {
    slr->start();
    lambda();
    //hclib::yield_at(nic);
}

}; // namespace selector

}; // namespace hclib


