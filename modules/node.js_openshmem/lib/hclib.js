
//const oshmem = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/node_hclib_openshmem.node');
const addon = require('./hclib_node.js_openshmem.node');
//const yaml = require('js-yaml');
const fs = require('fs');

var promFlush;
var ptr_array = new Object();

module.exports = {

init_hc : function (filename) {

  addon.init();
  //const type_arr = yaml.safeLoad(fs.readFileSync(filename, 'utf8'));
  //const indentedJson = JSON.stringify(type_arr, null, 4);
  //console.log(indentedJson);

  //for (const one_type of type_arr.arrays) {
    //console.log(one_type);
    //var size = one_type.size * one_type.count;
    //console.log("size " +size);
    //ptr_array[one_type.name] = addon.shmem_malloc(size);
  //}

  promFlush = setInterval(() => {
    //addon.clear_put_promises();
  }, 500);
},

finalize_hc : function () {
  addon.finalize();
  clearInterval(promFlush);
},

get_ptr_array : function(name) {
    return ptr_array[name];
},

my_rank: function(){
    return addon.my_rank();
},

num_ranks: function(){
    return addon.num_ranks();
}

};
