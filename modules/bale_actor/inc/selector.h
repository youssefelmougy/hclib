
namespace selector {

template<typename S, typename T>
void finish(S slr, T lambda) {
    //TODO: add details
}

};

template<typename T>
class Mailbox {

  public:

    std::function<void (T, int)> process;

    void send(T pkt, int rank) {}

    void done() {}

    //TODO: add details
};

template<int N, typename T>
class Selector {

 protected:

    Mailbox<T> mb[N];

  public:

    void start() {}

    void send(int mb_id, T pkt, int rank) {
        mb[mb_id].send(pkt, rank);
    }

    void done(int mb_id) {
        mb[mb_id].done();
    }

    //TODO: add details
};

