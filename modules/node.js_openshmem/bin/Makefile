include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/node.js_openshmem/inc/hclib_node.js_openshmem.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/node.js_openshmem/inc/hclib_node.js_openshmem.post.mak

NODEINC=$(NODE_HOME)/include/node
NODELIB=$(NODE_HOME)/lib/node
NODELIBS=-lnode

node_run: ../src/node_run.cpp
	g++ -std=c++11 $^ -I$(HCLIB_ROOT)/include \
		-I$(HCLIB_ROOT)/../modules/node.js_openshmem/inc \
		-I$(NODEINC) -L$(NODELIB) $(NODELIBS) \
		-lstdc++ -lhclib -lhclib_node.js_openshmem \
		-L$(HCLIB_ROOT)/lib -L../lib $(JSMN_HOME)/libjsmn.a -o $@

clean:
	rm -f node_run

