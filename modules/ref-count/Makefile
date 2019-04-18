include $(HCLIB_ROOT)/include/hclib.mak

ifeq (1, 0)

DEFINES=-DHC_ASSERTION_CHECK

OBJS=obj/hclib_ref-count.o

TARGETS=lib/libhclib_ref-count.so lib/libhclib_ref-count.a

all: $(TARGETS)

lib/libhclib_ref-count.a: $(OBJS)
	$(AR) rvs $@ $(OBJS)

lib/libhclib_ref-count.so: $(OBJS)
	$(CXX) -shared $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -o $@ $(OBJS)
	cp $@ $(HCLIB_ROOT)/lib

obj/hclib_ref-count.o: src/hclib_ref-count.cpp
	$(CXX) -g -O3 -std=c++14 -Iinc $^ -o $@ $(HCLIB_CFLAGS) --compile -fPIC $(DEFINES)

clean:
	rm -f obj/*.o lib/*

endif
