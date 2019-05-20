include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

TARGETS=heat_3d_simple_nonres heat_3d_simple_replay heat_3d_simple_replic heat_3d_simple_mix

FLAGS=-g -O3 -std=c++14 -Wall

all: $(TARGETS) targets.txt

.PHONY: targets.txt # always update
targets.txt:
	@echo "$(TARGETS)" > $@

%_nonres: %.cpp
	$(CXX) ${FLAGS} $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

%_replay: %.cpp
	$(CXX) ${FLAGS} -DUSE_REPLAY $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

%_replic: %.cpp
	$(CXX) ${FLAGS} -DUSE_REPLICATION $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

%: %.cpp
	$(CXX) ${FLAGS} $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

clean:
	rm -f $(TARGETS)

