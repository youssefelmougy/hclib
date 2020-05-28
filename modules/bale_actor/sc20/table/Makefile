
OSHCC ?=oshcc
OSHRUN ?=oshrun

TARGETS=histo_shmem histo_conveyor\
		ig_shmem ig_nbi_shmem ig_conveyor

all: $(TARGETS)

%_shmem: %_shmem.cpp
	$(OSHCC) -g -O3 -o $@ $^

%_conveyor: %_conveyor.cpp
	$(OSHCC) -g -O3 -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput

test:
	for target in $(TARGETS); do $(OSHRUN) ./$$target ; done
clean:
	rm -f $(TARGETS)
