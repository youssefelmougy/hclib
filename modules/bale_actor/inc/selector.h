
#include "safe_buffer.h"
#include "hclib_cpp.h"
extern "C" {
#include "convey.h"
}

#define DONE_MARK -1
#define BUFFER_SIZE 100000

namespace hclib {

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;

    BufferPacket() {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
};

template<typename T, int SIZE>
class Mailbox {

    hclib::conveyor::safe_buffer<BufferPacket<T>> *buff;
    convey_t* conv;
    BufferPacket<T> done_mark;

  public:

    Mailbox() {
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

    void start() {
        buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
        conv = convey_new(SIZE_MAX, 0, NULL, 0);
        assert( conv != nullptr );
        convey_begin(conv, sizeof(T));
        done_mark.rank = DONE_MARK;
    }

    void send(T pkt, int rank) {
        buff->push_back(BufferPacket<T>(pkt, rank));
    }

    void done() {
        buff->push_back(done_mark);
    }

    void start_worker_loop() {
        hclib::async_at([=]{
          while(true) {
              size_t buff_size = buff->size();
              if(buff_size > 0) break;
              //hclib::yield_at(nic);
          }

          BufferPacket<T> bp = buff->at(0);
          while(convey_advance(conv, bp.rank == DONE_MARK)) {
              int i;
              size_t buff_size = buff->size();
              for(i=1; i<=buff_size && bp.rank!=DONE_MARK; i++){
                  if( !convey_push(conv, &(bp.data), bp.rank)) break;
                  bp = buff->operator[](i);
              }
              {
                  std::lock_guard<std::mutex> lg(buff->get_mutex());
                  buff->erase_begin(i-1);
              }
              T pop;
              int64_t from;
              while( convey_pull(conv, &pop, &from) == convey_OK) {
                  hclib::async([=]() { process(pop, from); });
              }
              //hclib::yield_at(nic);
          }
        }, nic);
}
};



template<int N, typename T, int SIZE=BUFFER_SIZE>
class Selector {

  protected:

  public:

    Mailbox<T, SIZE> mb[N];

    void start() {
        for(int i=0; i<N; i++) {
            mb[i].start();
            mb[i].start_worker_loop();
        }
    }

    void send(int mb_id, T pkt, int rank) {
        mb[mb_id].send(pkt, rank);
    }

    void done(int mb_id) {
        mb[mb_id].done();
    }
};

namespace selector {

template<typename S, typename T>
void finish(S slr, T lambda) {
    slr->start();
    lambda();
    //TODO: add details about waiting for the end of all communication
}

}; // namespace selector

}; // namespace hclib


