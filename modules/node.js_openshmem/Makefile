include $(HCLIB_ROOT)/include/hclib.mak

TARGETS=flags lib/libhclib_node.js_openshmem.so node_run node_addon

all: $(TARGETS)

flags:
ifeq ($(findstring tr7, tr${CRAY_SHMEM_VER}),tr7)
  SHMEM_FLAG += -DUSE_CRAY_SHMEM_7
endif

lib/libhclib_node.js_openshmem.so: obj/hclib_node.js_openshmem.o
	$(OSHCC) -O3 -shared $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -o $@ $^
	cp $@ $(HCLIB_ROOT)/lib

obj/hclib_node.js_openshmem.o: src/hclib_node.js_openshmem.cpp
	$(OSHCC) -O3 -c $(SHMEM_FLAG) -std=c++11 -fPIC -Iinc $^ -o $@ $(HCLIB_CFLAGS)

node_run:
	cd bin && $(MAKE) $@

node_addon: src/addon_hclib_node.js_openshmem.cpp
	CXX=$(OSHCC) npm install --prefix obj 
	cp obj/build/Release/hclib_openshmem.node lib/hclib_node.js_openshmem.node

clean:
	rm -rf lib/*.so lib/*.node obj/*.o obj/build obj/etc obj/node_modules obj/package-lock.json
	cd bin && $(MAKE) $@
