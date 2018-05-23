
const addon = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/hclib');

addon.init_hc(process.argv[2])
console.log("inside js pid "+process.pid);
console.log("rank "+addon.my_rank()); 
addon.finalize_hc();
