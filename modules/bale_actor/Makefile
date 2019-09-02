include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

DEFINES=-DHC_ASSERTION_CHECK
# DEFINES=-DVERBOSE

OBJS=obj/hclib_bale_actor.o

TARGETS=lib/libhclib_bale_actor.so lib/libhclib_bale_actor.a

all: $(TARGETS)

lib/libhclib_bale_actor.a: $(OBJS)
	$(AR) rvs $@ $(OBJS)

lib/libhclib_bale_actor.so: $(OBJS)
	$(CXX) -shared -fPIC $(HCLIB_LDLIBS) $(HCLIB_LDFLAGS) -o $@ $(OBJS)

obj/hclib_bale_actor.o: src/hclib_bale_actor.cpp
	$(CXX) -g -c -fPIC -std=c++11 -Iinc $^ -o $@ $(HCLIB_CFLAGS) $(DEFINES)

clean:
	rm -f obj/*.o lib/*

install: $(TARGETS)
	cp lib/*.a lib/*.so $(HCLIB_ROOT)/lib
	cp inc/*.h $(HCLIB_ROOT)/include/
