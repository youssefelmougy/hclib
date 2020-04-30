include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak

#ifeq (1, 0)

DEFINES=-DHC_ASSERTION_CHECK
#DEFINES=-DHC_ASSERTION_CHECK -DUSE_FENIX
#DEFINES=-DHC_ASSERTION_CHECK -DCOMM_PROFILE

OBJS=obj/hclib_resilience.o obj/hclib_resilience_communication.o

TARGETS=lib/libhclib_resilience.so lib/libhclib_resilience.a

all: $(TARGETS)

lib/libhclib_resilience.a: $(OBJS)
	$(AR) rvs $@ $(OBJS)

lib/libhclib_resilience.so: $(OBJS)
	$(MPICXX) -shared $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -o $@ $(OBJS)
	cp $@ $(HCLIB_ROOT)/lib

obj/hclib_resilience.o: src/hclib_resilience.cpp
	$(CXX) -g -O3 -std=c++14 -Iinc $^ -o $@ $(HCLIB_CFLAGS) -c -fPIC $(DEFINES)

obj/hclib_resilience_communication.o: src/hclib_resilience_communication.cpp
	$(MPICXX) -g -O3 -std=c++14 -Iinc $^ -o $@ $(HCLIB_CFLAGS) -I $(FENIX)/include -c -fPIC $(DEFINES)

clean:
	rm -f obj/*.o lib/*

#endif
