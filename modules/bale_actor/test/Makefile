include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

TARGETS=ig_conveyor ig_selector ig_agi test_conveyor test_selector histo_conveyor histo_selector histo_agi toposort_agi toposort_conveyor

all: $(TARGETS)

%: %.cpp
	$(CXX) -g -O3 -std=c++11 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm
	#$(CXX) -g -O3 -std=c++11 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lconvey
	#$(CXX) -g -O3 -std=c++11 -I $(HPCTOOLKIT_HOME)/include $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -L $(HPCTOOLKIT_HOME)/lib/hpctoolkit/ -lconvey -lhclib_bale_actor -lhpctoolkit
	#$(CXX) -g -O3 -std=c++11 -DUSE_LOCK $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@_lock $^ $(HCLIB_LDLIBS) -lconvey -lhclib_bale_actor_lock

clean:
	rm -f $(TARGETS)
