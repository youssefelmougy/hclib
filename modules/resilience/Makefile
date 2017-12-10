include $(HCLIB_ROOT)/include/hclib.mak

DEFINES=-DHC_ASSERTION_CHECK

OBJS=obj/hclib_resilience.o

TARGETS=lib/libhclib_resilience.so lib/libhclib_resilience.a

all: $(TARGETS)

lib/libhclib_resilience.a: $(OBJS)
	ar rvs $@ $(OBJS)

lib/libhclib_resilience.so: $(OBJS)
	g++ -shared $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -o $@ $(OBJS)

obj/hclib_resilience.o: src/hclib_resilience.cpp
	g++ -g -O0 -std=c++11 -Iinc $^ -o $@ $(HCLIB_CFLAGS) --compile -fPIC $(DEFINES)

clean:
	rm -f obj/*.o lib/*
