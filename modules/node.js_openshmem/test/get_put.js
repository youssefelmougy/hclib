
const addon = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/hclib');

addon.init_hc(process.argv[2])
console.log("inside js pid "+process.pid);
const rank = addon.my_rank();
const NumProcs = addon.num_ranks();
console.log("rank " + rank + " num procs " + NumProcs);

console.log("SHMEM_BCAST_SYNC_SIZE "+ addon.SHMEM_BCAST_SYNC_SIZE);

const sizeof_long = 8;
var ptr = addon.malloc_sync(sizeof_long);
var lock_ptr = addon.malloc_sync(sizeof_long);
var dest_bcast = addon.malloc_sync(sizeof_long);
var src_bcast = addon.malloc_sync(sizeof_long);
var pSync_bcast = addon.malloc_sync(sizeof_long * addon.SHMEM_BCAST_SYNC_SIZE);
var pSync_reduce = addon.malloc_sync(sizeof_long * addon.SHMEM_REDUCE_SYNC_SIZE);
var ipWrk = addon.malloc_sync(sizeof_long * addon.SHMEM_REDUCE_MIN_WRKDATA_SIZE);

for (i = 0; i < addon.SHMEM_BCAST_SYNC_SIZE; i++){
    addon.put_value(pSync_bcast, addon.SHMEM_SYNC_VALUE, i);
}

for (i = 0; i < addon.SHMEM_REDUCE_SYNC_SIZE; i++){
    addon.put_value(pSync_reduce, addon.SHMEM_SYNC_VALUE, i);
}

if(rank == 0) {
    addon.long_p_sync(ptr, 33, 1);
    addon.long_p_sync(src_bcast, 37, 0);
    addon.barrier_all_sync();
    var val_recv = addon.long_g_sync(ptr, 1);
    console.log("value got on 0 from 1 is "+ val_recv);
    //addon.set_lock_sync(lock_ptr);
    addon.broadcast64_sync(dest_bcast, src_bcast, 1, 0, 0, 0, NumProcs, pSync_bcast);
    val_recv = addon.long_g_sync(dest_bcast, 0);
    console.log("bcast value " + val_recv + " at rank " + rank);
    addon.int_sum_to_all_sync(dest_bcast, src_bcast, 1, 0, 0, NumProcs, ipWrk, pSync_reduce);
    val_recv = addon.long_g_sync(dest_bcast, 0);
    console.log("sum_to_all value " + val_recv + " at rank " + rank);
}
else if(rank == 1) {
    addon.barrier_all_sync();
    var val_recv = addon.long_g_sync(ptr, 1);
    console.log("value got on 1 from 0 is "+ val_recv);
    //addon.set_lock_sync(lock_ptr);
    addon.broadcast64_sync(dest_bcast, src_bcast, 1, 0, 0, 0, NumProcs, pSync_bcast);
    val_recv = addon.long_g_sync(dest_bcast, 1);
    console.log("bcast value " + val_recv + " at rank " + rank);
    addon.int_sum_to_all_sync(dest_bcast, src_bcast, 1, 0, 0, NumProcs, ipWrk, pSync_reduce);
    val_recv = addon.long_g_sync(dest_bcast, 1);
    console.log("sum_to_all value " + val_recv + " at rank " + rank);
}
else {
    addon.barrier_all_sync();
    //addon.clear_lock_sync(lock_ptr);
    addon.broadcast64_sync(dest_bcast, src_bcast, 1, 0, 0, 0, NumProcs, pSync_bcast);
    var val_recv = addon.long_g_sync(dest_bcast, 2);
    console.log("bcast value " + val_recv + " at rank " + rank);
    addon.int_sum_to_all_sync(dest_bcast, src_bcast, 1, 0, 0, NumProcs, ipWrk, pSync_reduce);
    val_recv = addon.long_g_sync(dest_bcast, 2);
    console.log("sum_to_all value " + val_recv + " at rank " + rank);
}

addon.finalize_hc();

