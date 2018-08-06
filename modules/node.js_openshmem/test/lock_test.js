
const addon = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/hclib');

addon.init_hc(process.argv[2])
console.log("inside js pid "+process.pid);
const rank = addon.my_rank();
const NumProcs = addon.num_ranks();
console.log("rank " + rank + " num procs " + NumProcs);

var ptr = addon.malloc_sync(addon.SIZEOF_64_t);
var lock_ptr = addon.malloc_sync(addon.SIZEOF_64_t);

let use_lock = (prom, i, last) => {
    prom.then(function (fulfilled) {
        console.log("value from lock " + i + " is " + fulfilled);
        addon.clear_lock_async(lock_ptr);
        //addon.clear_lock_async_finish(lock_ptr);
        if(last == true)
            addon.finalize_hc();
    }).catch(function (error) {
        console.log(error.message);
    });
}

if(rank == 0) {
    var prom_lock1 = addon.set_lock_async(lock_ptr);
    var prom_lock2 = addon.set_lock_async(lock_ptr);
    var prom_lock3 = addon.set_lock_async(lock_ptr);
    var prom_lock4 = addon.set_lock_async(lock_ptr);
    addon.put_long_value(ptr, 29);
    var val1 = addon.get_long_value(ptr);
    console.log("local value get " + val1);

    use_lock(prom_lock3, 3, false);
    use_lock(prom_lock4, 4, true);
    use_lock(prom_lock2, 2, false);
    use_lock(prom_lock1, 1, false);
}
else if(rank == 1) {
    addon.set_lock_sync(lock_ptr);
    var prom = addon.long_g_async(ptr, 0);
    addon.clear_lock_sync(lock_ptr);
    prom.then(function (fulfilled) {
        console.log("value from 0 is " + fulfilled);
        addon.finalize_hc();
    }).catch(function (error) {
        console.log(error.message);
    });
    addon.set_lock_sync(lock_ptr);
    addon.clear_lock_sync(lock_ptr);
    addon.finalize_hc();
}
else {
    addon.finalize_hc();
}

