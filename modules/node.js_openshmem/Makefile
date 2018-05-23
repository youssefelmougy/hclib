include $(HCLIB_ROOT)/include/hclib.mak

TARGETS=lib/libhclib_node.js_openshmem.so node_run node_addon

all: $(TARGETS)

lib/libhclib_node.js_openshmem.so: obj/hclib_node.js_openshmem.o
	$(OPENSHMEM_INSTALL)/bin/oshcc -shared $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -o $@ $^
	cp $@ $(HCLIB_ROOT)/lib

obj/hclib_node.js_openshmem.o: src/hclib_node.js_openshmem.cpp
	$(OPENSHMEM_INSTALL)/bin/oshcc -c -std=c++11 -fPIC -Iinc $^ -o $@ $(HCLIB_CFLAGS)

node_run:
	cd bin && $(MAKE) $@

node_addon: src/addon_hclib_node.js_openshmem.cpp
	CXX=$(OPENSHMEM_INSTALL)/bin/oshcc npm install --prefix obj 
	cp obj/build/Release/hclib_openshmem.node lib/hclib_node.js_openshmem.node

clean:
	rm -rf lib/*.so lib/*.node obj/*.o obj/build obj/etc obj/node_modules obj/package-lock.json
	cd bin && $(MAKE) $@
