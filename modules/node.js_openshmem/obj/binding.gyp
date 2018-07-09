{
  "targets": [
    {
      "target_name": "hclib_openshmem",
      "cflags!": [ "-fno-exceptions" ],
      "cflags_cc!": [ "-fno-exceptions" ],
      "cflags": [ "-fPIC", "-DUSE_OFFLOAD" ],
      "sources": [
        "../src/addon_hclib_node.js_openshmem_util.cpp",
        "../src/addon_hclib_node.js_openshmem_sync.cpp",
        "../src/addon_hclib_node.js_openshmem.cpp"
      ],
      "include_dirs": [
        "<!@(node -p \"require('node-addon-api').include\")",
        "<!(echo $HCLIB_ROOT)/include",
        "<!(echo $HCLIB_HOME)/modules/system/inc",
        "<!(echo $HCLIB_HOME)/modules/node.js_openshmem/inc",
#        "<!(echo $OPENSHMEM_INSTALL)/include",
        "<!(echo $JSMN_HOME)",
      ],
      "link_settings": {
        "libraries": [
          "-lhclib -ljsmn -ldl -lhclib_system -lhclib_node.js_openshmem",
        ],
        "library_dirs": [
          "<!(echo $HCLIB_ROOT)/lib",
          "<!(echo $HCLIB_HOME)/modules/system/lib",
          "<!(echo $HCLIB_HOME)/modules/node.js_openshmem/lib",
#          "<!(echo $OPENSHMEM_INSTALL)/lib",
          "<!(echo $JSMN_HOME)",
        ],
      },
      'defines': [ 'NAPI_DISABLE_CPP_EXCEPTIONS' ],
    }
  ]
}

