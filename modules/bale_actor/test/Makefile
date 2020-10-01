include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

TARGETS=test_conveyor test_selector test_selector3 \
        histo_agi histo_conveyor histo_selector histo_selector_lambda\
        ig_agi ig_conveyor ig_selector ig_selector_lambda \
        toposort_agi toposort_conveyor toposort_selector \
		triangle_agi triangle_conveyor triangle_selector \
		randperm_agi randperm_conveyor randperm_selector \
        permute_agi permute_conveyor permute_selector \
		transpose_agi transpose_conveyor transpose_selector

all: $(TARGETS)

%: %.cpp
	$(CXX) -g -O3 -std=c++11 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm
	#$(CXX) -g -O3 -std=c++11 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -llibgetput
	#$(CXX) -g -O3 -std=c++11 -I $(HPCTOOLKIT_HOME)/include $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -L $(HPCTOOLKIT_HOME)/lib/hpctoolkit/ -lspmat -lconvey -lexstack -llibgetput  -lhclib_bale_actor -lm -lhpctoolkit
	#$(CXX) -g -O3 -std=c++11 -DUSE_LOCK $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@_lock $^ $(HCLIB_LDLIBS) -lconvey -lhclib_bale_actor_lock

clean:
	rm -f $(TARGETS)
