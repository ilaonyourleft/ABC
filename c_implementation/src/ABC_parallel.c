/*
 ============================================================================
 Name        : ABC_parallel.c
 Author      : Ilaria Malinconico
 Version     :
 Copyright   : Your copyright notice
 Description : ABC parallel implementation
 ============================================================================
 */

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int i;
    int id; // process rank
    int p; // number of processes

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /*CODE HERE*/

    printf("=== Process %d is done.\n", id);
    fflush(stdout);
    MPI_Finalize();
    return 0;
}