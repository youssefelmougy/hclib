include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

TARGETS=resilience-d0 resilience-d1 resilience-r0 resilience-r1 resilience-a0 resilience-mix resilience-c0

all: $(TARGETS)

%: %.cpp
	g++ -g -O0 -std=c++14 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

test:
	for target in $(TARGETS); do ./$$target ; done

clean:
	rm -f $(TARGETS)
