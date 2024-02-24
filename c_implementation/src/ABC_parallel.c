/*
 ============================================================================
 Name        : ABC_parallel.c
 Author      : Ilaria Malinconico
 Version     :
 Copyright   : Your copyright notice
 Description : ABC sequential implementation
 ============================================================================
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ABC_lib.h"

int main(int argc, char *argv[]) {
    // Open input and output files
    FILE *file = fopen("../data/dataset10000.csv", "r");
    FILE *output = fopen("../results/results_parallel_10000.txt", "w");

    // Variable declarations
	int i = 0, j = 0, g = 0, d = 0, e = 0, localCounter = 0, counter = 0, factor = BETA * N, otherFactor = 0;
	float x, y = 0.0;

    // Variable declarations for parallelization
    int rank, size;

    // Error handling for file opening
    if (file == NULL || output == NULL) {
        printf("Could not open file.\n");
        exit(OPEN_FILE_ERROR);
    }

    // MPI initialization
    MPI_Init(&argc, &argv);
    // Process rank in a communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Number of processes in a communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    float elapsed_time = - MPI_Wtime();

    MPI_Datatype mpi_double_float_type, mpi_triple_float_type;
    MPI_Type_contiguous(2, MPI_FLOAT, &mpi_double_float_type);
    MPI_Type_commit(&mpi_double_float_type);
    MPI_Type_contiguous(3, MPI_FLOAT, &mpi_triple_float_type);
    MPI_Type_commit(&mpi_triple_float_type);

    int chunkSize = N / size;
    int start = rank * chunkSize;
    int end = (rank == size - 1) ? N : start + chunkSize;
    int remainder = N % size;

    // --------------------------------------------------- MEMORY ALLOCATIONS
    struct double_float *ptrPoints;
    struct triple_float *ptrBorderPointsAll;
    struct point_label *borderPointsAndLabels;

    ptrPoints = calloc(N, sizeof(struct double_float));
    if (ptrPoints == NULL) {
        printErrorAllocation();
    }

    int *displs = calloc(size, sizeof(int));
    if (displs == NULL) {
        printErrorAllocation();
    }

    int *counts = calloc(size, sizeof(int));
    if (counts == NULL) {
        printErrorAllocation();
    }

    int local_N = 0;
    if (rank < remainder) {
        local_N = chunkSize + 1;
    } else {
        local_N = chunkSize;
    }

    int displacement = 0;
    for (int h = 0; h < size; h++) {
        counts[h] = (h < remainder) ? (chunkSize + 1) : chunkSize;
        displs[h] = displacement;
        displacement += counts[h];
    }

    if (rank == 0) {
        /*ptrPoints = calloc(N, sizeof(struct double_float));
        if (ptrPoints == NULL) {
            printErrorAllocation();
        }*/

        ptrBorderPointsAll = calloc(N, sizeof(struct triple_float));
        if (ptrBorderPointsAll == NULL) {
            printErrorAllocation();
        }

        // Create an array of points
        for (i = 0; i < N; i++) {
            fscanf(file, "%f,%f", &x, &y);
            ptrPoints[i].x = x;
            ptrPoints[i].y = y;
        }

        borderPointsAndLabels = calloc(factor, sizeof(struct point_label));
        if (borderPointsAndLabels == NULL) {
            printErrorAllocation();
        }
    }

    MPI_Bcast(ptrPoints, N, mpi_double_float_type, 0, MPI_COMM_WORLD);

    struct double_float *localPtrPoints = calloc(local_N, sizeof(struct double_float));
    if (localPtrPoints == NULL) {
        printErrorAllocation();
    }

    MPI_Scatterv(ptrPoints, counts, displs, mpi_double_float_type, localPtrPoints, local_N, mpi_double_float_type, 0, MPI_COMM_WORLD);

    struct double_float *localPtrKnnPoint = calloc(K, sizeof(struct double_float));
    if (localPtrKnnPoint == NULL) {
        printErrorAllocation();
    }

    float *localPtrMeanPoint = calloc(2, sizeof(float));
    if (localPtrMeanPoint == NULL) {
        printErrorAllocation();
    }

    float *localPtrDirectionalAnglesPoint = calloc(K, sizeof(float));
    if (localPtrDirectionalAnglesPoint == NULL) {
        printErrorAllocation();
    }

    struct triple_float *localPtrBorderPointsAll = calloc(local_N, sizeof(struct triple_float));
    if (localPtrBorderPointsAll == NULL) {
        printErrorAllocation();
    }

    // --------------------------------------------------- BORDER POINTS LOOP

    // Read data from the input file and create an array of points
    for (i = 0; i < counts[rank]; i++) {
        // Find k nearest neighbors for each point and the mean point
        getNeighbors(ptrPoints, localPtrPoints[i].x, localPtrPoints[i].y, localPtrKnnPoint, localPtrMeanPoint);

        // Find directional angles between the center, its k nearest neighbors, and the mean point
        for (j = 0; j < K; j++) {
            localPtrDirectionalAnglesPoint[j] = getDirectionalAngle(localPtrPoints[i], localPtrMeanPoint, localPtrKnnPoint[j]);
        }

        // Find the border points with enclosing angle for each point and border degree
        if (isBorderPoint(getEnclosingAngle(localPtrDirectionalAnglesPoint)) == 1) {
            localPtrBorderPointsAll[localCounter].x = localPtrPoints[i].x;
            localPtrBorderPointsAll[localCounter].y = localPtrPoints[i].y;
            localPtrBorderPointsAll[localCounter].z = getBorderDegree(localPtrDirectionalAnglesPoint);
            ++localCounter;
        }
    }

    MPI_Gatherv(localPtrBorderPointsAll, local_N, mpi_triple_float_type, ptrBorderPointsAll, counts, displs, mpi_triple_float_type, 0, MPI_COMM_WORLD);

    MPI_Reduce(&localCounter, &counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    free(localPtrPoints);
    free(displs);
    free(counts);
    free(localPtrKnnPoint);
    free(localPtrMeanPoint);
    free(localPtrDirectionalAnglesPoint);
    free(localPtrBorderPointsAll);

    elapsed_time += MPI_Wtime();
    if (rank == 0) {
        printf("Elapsed time (parallel) = %f\n", elapsed_time);
    }

    if (rank == 0) {
        // Get factor border points
        getBorderPoints(ptrBorderPointsAll, counter, borderPointsAndLabels);
        free(ptrBorderPointsAll);

        // Check if counter is less than factor
        if (counter < factor) {
            factor = counter;
        }

        // Get label for each border point
        getLabelsBorderPoints(factor, 1, 3, borderPointsAndLabels);

        // Write border points and labels to the output file
        for (d = 0; d < factor; d++) {
            if (borderPointsAndLabels[d].x != 0 && borderPointsAndLabels[d].y != 0) {
                fprintf(output, "%f,%f,%d\n", borderPointsAndLabels[d].x, borderPointsAndLabels[d].y, borderPointsAndLabels[d].label);
            }
        }

        otherFactor = N - factor;

        struct point_label *internalPointsAndLabels = calloc(otherFactor, sizeof(struct point_label));
        if (internalPointsAndLabels == NULL) {
            printErrorAllocation();
        }

        // Get non-border points
        getNonBorderPoints(ptrPoints, borderPointsAndLabels, factor, internalPointsAndLabels);
        free(ptrPoints);

        // Get labels for non-border points
        getLabelsNonBorderPoints(factor, borderPointsAndLabels, internalPointsAndLabels, otherFactor);

        // Write non border points and labels to the output file
        for (e = 0; e < otherFactor; e++) {
            if (internalPointsAndLabels[e].x != 0 && internalPointsAndLabels[e].y != 0) {
                fprintf (output, "%f,%f,%d\n", internalPointsAndLabels[e].x, internalPointsAndLabels[e].y, internalPointsAndLabels[e].label);
            }
        }
        
        // Free remaining allocated memory
        free(borderPointsAndLabels);
        free(internalPointsAndLabels);
    }

	// Close output file
	fclose(output);

    // MPI finalization
    MPI_Finalize();

    return 0;
}
