include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

IXX=icpc

TARGETS=resilient_cholesky_nonres resilient_cholesky_replay resilient_cholesky_replic resilient_cholesky_abft

all: $(TARGETS)

%_nonres: %.cpp
	$(IXX) -g -O3 -std=c++14 $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

%_refct: %.cpp
	$(IXX) -g -O3 -std=c++14 -DUSE_REF_CT $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

%_replay: %.cpp
	$(IXX) -g -O3 -std=c++14 -DUSE_REPLAY $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

%_replic: %.cpp
	$(IXX) -g -O3 -std=c++14 -DUSE_REPLICATION $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

%_abft: %.cpp
	$(IXX) -g -O3 -std=c++14 -DUSE_ABFT $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

%_checkpoint: %.cpp
	$(IXX) -g -O3 -std=c++14 -DUSE_CHECKPOINT $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

%: %.cpp
	$(IXX) -g -O3 -std=c++14 $(FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -mkl

test:
	for target in $(TARGETS); do ./$$target ; done

clean:
	rm -f $(TARGETS)
