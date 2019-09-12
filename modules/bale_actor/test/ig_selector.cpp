#include <shmem.h>
#include "hclib_bale_actor.h"

struct IgPkt {
    int idx;
    int val;
};

enum MailBoxType{REQUEST, RESPONSE};

class IgSelector: public Selector<2, IgPkt> {
  
  // Allocate and populate the shared table array
  int ltab_siz, l_num_req;
  int64_t * table, *tgt;
  int num_processed = 0;

  void req_process(IgPkt pkt, int sender_rank) {
      pkt.val = table[pkt.val];
      send(RESPONSE, pkt, sender_rank);
      
      num_processed ++;
      if(num_processed == l_num_req)
        end();
  }

  void resp_process(IgPkt pkt, int sender_rank) {
      tgt[pkt.idx] = pkt.val;
  }

  public:

    IgSelector() {
        mb[REQUEST].process = [this](IgPkt pkt, int sender_rank) { this->req_process(pkt, sender_rank); };
        mb[RESPONSE].process = [this](IgPkt pkt, int sender_rank) { this->resp_process(pkt, sender_rank); };
    }

    void init_ig(int tab_siz, int num_req) {
        ltab_siz = tab_siz;
        l_num_req = num_req;
        table   =  (int64_t*)shmem_malloc(ltab_siz * sizeof(int64_t));
        tgt  =  (int64_t*)malloc(l_num_req* sizeof(int64_t));
    }

};

int main() {

  const char *deps[] = { "system", "bale_actor" };
  hclib::launch(deps, 2, [=] {

    int tab_siz = 10, num_req = 100;
    IgSelector igs;
    igs.init_ig(tab_siz, num_req);

    int64_t *index   =  (int64_t*)malloc(num_req* sizeof(int64_t));
    //TODO: populate index array with some indexes
    
    
    int num_ranks = shmem_n_pes();
    for(int i=0;i<num_req;i++) {
        IgPkt pkt;
        pkt.val = index[i] % num_ranks;
        pkt.idx = i;
        int dest_rank = index[i] / num_ranks;
        igs.send(REQUEST, pkt, dest_rank);
    }

    //Do we need to notify that send from here is over?

  });
  return 0;
}

