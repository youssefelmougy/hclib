include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

TARGETS=cg_baseline cg_replay cg_replication

FLAGS=-O3 -std=c++14 -Wall -DMAX_NUM_WAITS=256 -DUSE_STD_VEC

all: $(TARGETS) targets.txt

.PHONY: targets.txt # always update
targets.txt:
	@echo "$(TARGETS)" > $@

%: %.cpp
	$(CXX) ${FLAGS} $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

clean:
	rm -f $(TARGETS)

