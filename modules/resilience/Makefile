include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak

DEFINES=-DHC_ASSERTION_CHECK

OBJS=obj/hclib_resilience.o

TARGETS=lib/libhclib_resilience.so lib/libhclib_resilience.a

all: $(TARGETS)

lib/libhclib_resilience.a: $(OBJS)
	ar rvs $@ $(OBJS)

lib/libhclib_resilience.so: $(OBJS)
	g++ -shared $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -o $@ $(OBJS)
	cp $@ $(HCLIB_ROOT)/lib

obj/hclib_resilience.o: src/hclib_resilience.cpp
	g++ -g -O0 -std=c++14 -Iinc $^ -o $@ $(HCLIB_CFLAGS) --compile -fPIC $(DEFINES)

clean:
	rm -f obj/*.o lib/*
