include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

EXE=smith_waterman

FLAGS=-g -O3 -std=c++14 -Wall

all: clean $(EXE) clean-obj

smith_waterman: smith_waterman.cpp
	$(CXX) $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@.baseline.out $^ $(HCLIB_LDLIBS)
	$(CXX) $(FLAGS) -DUSE_REPLAY $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@.replay.out $^ $(HCLIB_LDLIBS)
	$(CXX) $(FLAGS) -DUSE_REPLICATION $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@.replication.out $^ $(HCLIB_LDLIBS)

clean-obj:
	rm -rf *.o *.dSYM

clean:
	rm -rf *.o $(EXE) *.dSYM smith_waterman.*.out
