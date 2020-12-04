
#ifndef SELECTOR_AGI_H
#define SELECTOR_AGI_H

#include<selector.h>

struct PutPkt {
    int64_t *loc;
    int64_t val;
};

struct GetPkt {
    int64_t *dest;
    union {
        int64_t *src;
        int64_t val;
    };
};

enum MailBoxType{REQUEST, RESPONSE};

class Put : public hclib::Actor<PutPkt> {

  void process(PutPkt pkt, int sender_rank) {
      *(pkt.loc) = pkt.val;
  }

  public:
    Put() {
        mb[0].process = [this](PutPkt pkt, int sender_rank) { this->process(pkt, sender_rank);};
        start();
    }

    void operator()(int64_t *loc, int64_t val, int pe) {
        send({loc, val}, pe);
    }
};

class Get : public hclib::Selector<2, GetPkt>{

  void req_process(GetPkt pkt, int sender_rank) {
      pkt.val = *(pkt.src);
      send(RESPONSE, pkt, sender_rank);
  }

  void resp_process(GetPkt pkt, int sender_rank) {
      *(pkt.dest) = pkt.val;
  }

  public:

    Get() {
        mb[REQUEST].process = [this](GetPkt pkt, int sender_rank) { this->req_process(pkt, sender_rank); };
        mb[RESPONSE].process = [this](GetPkt pkt, int sender_rank) { this->resp_process(pkt, sender_rank); };
        start();
    }

    void operator()(int64_t *dest, int64_t *src, int pe) {
       send(REQUEST, {dest, src}, pe); 
    }

    void done() {
        hclib::Selector<2, GetPkt>::done(REQUEST);
    }
};

#endif

