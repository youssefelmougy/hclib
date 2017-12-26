include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

TARGETS=resilience-d0 resilience-d1

all: $(TARGETS)

%: %.cpp
	g++ -g -O0 -std=c++14 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

clean:
	rm -f $(TARGETS)
