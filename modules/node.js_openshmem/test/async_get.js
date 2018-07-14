
const addon = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/hclib');

addon.init_hc(process.argv[2])
console.log("inside js pid "+process.pid);
const rank = addon.my_rank();
const NumProcs = addon.num_ranks();
console.log("rank " + rank + " num procs " + NumProcs);

var ptr = addon.malloc_sync(addon.SIZEOF_64_t);

if(rank == 0) {
    addon.put_long_value(ptr, 29);
    var val1 = addon.get_long_value(ptr);
    console.log("local value get " + val1);
    addon.barrier_all_sync();
    addon.finalize_hc();
}
else if(rank == 1) {
    addon.barrier_all_sync();
    var prom = addon.long_g_async(ptr, 0);
    prom.then(function (fulfilled) {
        console.log("value from 0 is " + fulfilled);
        addon.finalize_hc();
    }).catch(function (error) {
        console.log(error.message);
    });
}
else {
    addon.barrier_all_sync();
    addon.finalize_hc();
}

