
//mpirun -n 4 ./hello-3test-mpif

#include <fenix.h>

const int kKillID = 0;
const int send_rank = 1;
const int recv_rank = 2;
MPI_Request r1, r2, r3;
int d1=-1, d2=-2, d3=-3;

int main(int argc, char **argv) {
        //int rank;
        //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Init(&argc, &argv);

        int old_world_size, new_world_size = - 1;
        int old_rank = 1, new_rank = - 1;
        int spare_ranks=1;

        MPI_Comm world_comm;
        MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
        MPI_Comm_size(world_comm, &old_world_size);
        MPI_Comm_rank(world_comm, &old_rank);

        int fenix_status;
        int recovered = 0;
        MPI_Comm new_comm;
        int error;
        Fenix_Init(&fenix_status, world_comm, &new_comm, nullptr, nullptr, spare_ranks, 0, MPI_INFO_NULL, &error);
          MPI_Comm_rank(new_comm, &new_rank);
        printf("inited rank %d %d\n", old_rank, new_rank);

        if (fenix_status != FENIX_ROLE_INITIAL_RANK) {
          MPI_Comm_size(new_comm, &new_world_size);
          MPI_Comm_rank(new_comm, &new_rank);
          recovered = 1;
            printf("inside FENIX_ROLE_INITIAL_RANK rank %d\n", new_rank);
        }

        if (old_rank == kKillID &&  recovered == 0) {
          MPI_Barrier(new_comm);
          pid_t pid = getpid();
          kill(pid, SIGTERM);
        }
        
        if(new_rank==send_rank) {
            //MPI_Request r1, r2, r3;
            if(recovered == 0) {
                d1=23; d2=37; d3=43;
                MPI_Isend(&d1, 1, MPI_INT, recv_rank, 1, new_comm, &r1);
                //MPI_Isend(&d2, 1, MPI_INT, recv_rank, 2, new_comm, &r2);
                MPI_Isend(&d3, 1, MPI_INT, recv_rank, 3, new_comm, &r3);

                int c1=0, c2=0, c3=0;
                while(!c1)
                    MPI_Test(&r1, &c1, MPI_STATUS_IGNORE);
                printf("send in rank %d data %d\n", new_rank, d1);
                MPI_Barrier(new_comm);
                MPI_Barrier(new_comm);
                while(!c2)
                    MPI_Test(&r2, &c2, MPI_STATUS_IGNORE);
                printf("send in rank %d data %d\n", new_rank, d2);
                while(!c3)
                    MPI_Test(&r3, &c3, MPI_STATUS_IGNORE);
                printf("send in rank %d data %d\n", new_rank, d3);
            }
            else {
                int c1=0, c2=0, c3=0;
                MPI_Test(&r2, &c2, MPI_STATUS_IGNORE);
                MPI_Test(&r3, &c3, MPI_STATUS_IGNORE);
                printf("recovered send in rank %d status %d %d\n", new_rank, c2, c3);
            }
        }
        else if(new_rank==recv_rank) {
            //MPI_Request r1, r2, r3;
            if(recovered == 0) {
                d1=-1; d2=-2; d3=-3;
                MPI_Irecv(&d1, 1, MPI_INT, send_rank, 1, new_comm, &r1);
                MPI_Irecv(&d2, 1, MPI_INT, send_rank, 2, new_comm, &r2);
                MPI_Irecv(&d3, 1, MPI_INT, send_rank, 3, new_comm, &r3);

                int c1=0, c2=0, c3=0;
                while(!c3)
                    MPI_Test(&r3, &c3, MPI_STATUS_IGNORE);
                printf("recieved in rank %d data %d\n", new_rank, d3);
                MPI_Barrier(new_comm);
                MPI_Barrier(new_comm);
                while(!c2)
                    MPI_Test(&r2, &c2, MPI_STATUS_IGNORE);
                printf("recieved in rank %d data %d\n", new_rank, d2);
                while(!c1)
                    MPI_Test(&r1, &c1, MPI_STATUS_IGNORE);
                printf("recieved in rank %d data %d\n", new_rank, d1);
            }
            else {
                int c1=0, c2=0, c3=0;
                MPI_Test(&r2, &c2, MPI_STATUS_IGNORE);
                MPI_Test(&r2, &c1, MPI_STATUS_IGNORE);
                printf("recovered recv in rank %d data %d(%d) %d(%d)\n", new_rank,d1, c1, d2, c2);
            }
        }

        //if(new_rank == kKillID) {
        //    int a=23;
        //    MPI_Send(&a, 1, MPI_INT, kKillID-1, 0, new_comm);
        //}
        //if(new_rank == kKillID-1) {
        //    int b =22;
        //    MPI_Status status;
        //    printf("trying to recv in rank %d\n", new_rank);
        //    MPI_Recv(&b, 1, MPI_INT, kKillID, 0, new_comm, &status);
        //    printf("recv %d in rank %d\n", b, new_rank);
        //}

        MPI_Barrier(new_comm);

        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

        printf("hello world: %s, old rank (MPI_COMM_WORLD): %d, new rank: %d, active ranks: %d, ranks before process failure: %d\n",
               processor_name, old_rank, new_rank, new_world_size, old_world_size);

        Fenix_Finalize();
        MPI_Finalize();

        //printf("Hello world from rank %d\n", rank);
    return 0;
}
