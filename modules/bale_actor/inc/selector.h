
template<typename T>
class Mailbox {

  public:

    std::function<void (T, int)> process;

    void send(T pkt, int rank) {}

    void end() {}

    //TODO: add details
};

template<int N, typename T>
class Selector {

 protected:

    Mailbox<T> mb[N];

  public:

    void send(int mb_id, T pkt, int rank) {
        mb[mb_id].send(pkt, rank);
    }

    void end() {
        for(int i=0;i<N;i++) {
            mb[i].end();
        }
    }

    //TODO: add details
};

