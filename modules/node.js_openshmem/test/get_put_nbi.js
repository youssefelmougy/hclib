
const addon = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/hclib');

addon.init_hc(process.argv[2])
console.log("inside js pid "+process.pid);
const rank = addon.my_rank();
const NumProcs = addon.num_ranks();
console.log("rank " + rank + " num procs " + NumProcs);

var buffer1 = addon.malloc_ab_sync(addon.SIZEOF_64_t);
var arr1 = new Float64Array(buffer1);
var buffer2 = new ArrayBuffer(8);
var arr2 = new Float64Array(buffer2);

if(rank == 0) {
    arr1[0] = 29;
    console.log("local value get " + arr1[0]);
    addon.barrier_all_sync();
    addon.finalize_hc();
}
else if(rank == 1) {
    addon.barrier_all_sync();
    addon.double_g_nbi_sync(buffer2, 0, buffer1, 0, 0);
    console.log("value from 0 on 1 before quiet is " + arr2[0]);
    addon.quiet_sync();
    console.log("value from 0 on 1 after quiet is " + arr2[0]);
    addon.finalize_hc();
}
else {
    addon.barrier_all_sync();
    addon.double_g_nbi_sync(buffer2, 0, buffer1, 0, 0);
    console.log("value from 0 on 2 before quiet is " + arr2[0]);
    var prom = addon.quiet_async();
    prom.then(function (fulfilled) {
        console.log("value from 0 on 2 after quiet is " + arr2[0]);
        addon.finalize_hc();
    })
}

