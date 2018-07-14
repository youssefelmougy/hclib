#include "hclib_cpp.h"
#include "hclib_node.js_openshmem.h"
#include "node.h"

#include <iostream>
#include <unistd.h>
//#include <shmem.h>

int main(int argc, char **argv) {
    const char *deps[] = { "system", "node.js_openshmem" };
    //shmem_init();
    hclib::launch(deps, 2, [=] {
     // int nic_locale_id = hclib_add_known_locale_type("sysmem");
     // int n_nics;
     // hclib::locale_t **nics = hclib::get_all_locales_of_type(nic_locale_id, &n_nics);
     // HASSERT(n_nics == 1);
     // HASSERT(nics);

     // hclib::async_at([=](){
        std::cout << "Hello world from rank " << hclib::shmem_my_pe() << " pid "<< getpid() << std::endl;
        node::Start(argc, argv);
      //}, nics[0]);
    });
    //shmem_finalize();
    return 0;
}
