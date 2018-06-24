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
      hclib::finish([=]() {
        std::cout << "Hello world from rank " << hclib::shmem_my_pe() << " pid "<< getpid() << std::endl;
        node::Start(argc, argv);
      });
    });
    //shmem_finalize();
    return 0;
}
