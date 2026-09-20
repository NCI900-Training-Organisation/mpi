/* Intentional unsafe ordering with 16 MiB messages.
 * Compile: mpicc -Wall -Wextra -o mpi_deadlock mpi_deadlock.c
 * Run: timeout --kill-after=2s 5s mpiexec -np 2 ./mpi_deadlock
 */
#include <stdio.h>
#include <mpi.h>

#define MESSAGE_BYTES (16 * 1024 * 1024)

/* Static buffers are zero-initialised and do not use the stack. */
static char outgoing[MESSAGE_BYTES];
static char incoming[MESSAGE_BYTES];

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        if (rank == 0) fprintf(stderr, "Run with exactly 2 ranks.\n");
        MPI_Finalize();
        return 1;
    }

    int peer = 1 - rank;
    printf("Rank %d: sending 16 MiB to rank %d\n", rank, peer);
    fflush(stdout);

    /* Both sends can wait for receives that neither rank has posted. */
    MPI_Send(outgoing, MESSAGE_BYTES, MPI_BYTE, peer, 0, MPI_COMM_WORLD);
    MPI_Recv(incoming, MESSAGE_BYTES, MPI_BYTE, peer, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("Rank %d: received the message\n", rank);
    MPI_Finalize();
    return 0;
}
