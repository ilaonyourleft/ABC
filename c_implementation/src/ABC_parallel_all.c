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
#include <stdlib.h>
#include <math.h>
#include "ABC_lib.h"

int main(int argc, char *argv[]) {
    // Open input and output files
    FILE *file = fopen("../data/dataset10000.csv", "r");
    FILE *output = fopen("../results/results_parallel_10000.txt", "w");

    // Variable declarations
	int i = 0, j = 0, g = 0, d = 0, e = 0, localCounter = 0, counter = 0, initialFactor = BETA * N, factor = BETA * N, otherFactor = 0;
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

    MPI_Datatype mpi_double_float_type, mpi_triple_float_type, mpi_point_label_type;
    MPI_Type_contiguous(2, MPI_FLOAT, &mpi_double_float_type);
    MPI_Type_commit(&mpi_double_float_type);
    MPI_Type_contiguous(3, MPI_FLOAT, &mpi_triple_float_type);
    MPI_Type_commit(&mpi_triple_float_type);
    int blocklengths[3] = {1, 1, 1};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(struct point_label, x);
    offsets[1] = offsetof(struct point_label, y);
    offsets[2] = offsetof(struct point_label, label);
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_INT};
    MPI_Type_create_struct(3, blocklengths, offsets, types, &mpi_point_label_type);
    MPI_Type_commit(&mpi_point_label_type);

    int chunkSize = N / size;
    int start = rank * chunkSize;
    int end = (rank == size - 1) ? N : start + chunkSize;
    int remainder = N % size;

    struct double_float *ptrPoints;
    struct triple_float *ptrBorderPointsAll;
    struct point_label *borderPointsAndLabels, *internalPointsAndLabels;

    // --------------------------------------------------- MEMORY ALLOCATIONS - PHASE 1
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
    }

    borderPointsAndLabels = calloc(initialFactor, sizeof(struct point_label));
    if (borderPointsAndLabels == NULL) {
        printErrorAllocation();
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

    // --------------------------------------------------- BORDER POINTS LOOP - PHASE 1

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

    // --------------------------------------------------- MEMORY ALLOCATION - PHASE 2

    if (rank == 0) {
        // Get factor border points
        getBorderPoints(ptrBorderPointsAll, counter, borderPointsAndLabels);
        free(ptrBorderPointsAll);

        // Check if counter is less than factor
        if (counter < factor) {
            factor = counter;
        }
    }

    MPI_Bcast(&factor, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Factor: %d, Rank: %d\n", factor, rank);

    MPI_Bcast(borderPointsAndLabels, initialFactor, mpi_point_label_type, 0, MPI_COMM_WORLD);

    int *borderDispls = calloc(size, sizeof(int));
    if (borderDispls == NULL) {
        printErrorAllocation();
    }

    int *borderCounts = calloc(size, sizeof(int));
    if (borderCounts == NULL) {
        printErrorAllocation();
    }

    int local_factor = 0;
    if (rank < remainder) {
        local_factor = chunkSize + 1;
    } else {
        local_factor = chunkSize;
    }

    int borderDisplacement = 0;
    for (int h = 0; h < size; h++) {
        borderCounts[h] = (h < remainder) ? (chunkSize + 1) : chunkSize;
        borderDispls[h] = borderDisplacement;
        borderDisplacement += borderCounts[h];
    }

    struct point_label *localBorderPointsAndLabels = calloc(local_factor, sizeof(struct point_label));
    if (localBorderPointsAndLabels == NULL) {
        printErrorAllocation();
    }

    MPI_Scatterv(borderPointsAndLabels, borderCounts, borderDispls, mpi_point_label_type, localBorderPointsAndLabels, local_factor, mpi_point_label_type, 0, MPI_COMM_WORLD);

    // --------------------------------------------------- BORDER POINTS CLUSTERING - PHASE 2

    // Get label for each border point
    getLabelsBorderPoints(factor, 1, 3, borderPointsAndLabels);

    for (d = 0; d < factor; d++) {
        localBorderPointsAndLabels[d].label = borderPointsAndLabels[d].label;
    }

    MPI_Gatherv(localBorderPointsAndLabels, local_factor, mpi_point_label_type, borderPointsAndLabels, borderCounts, borderDispls, mpi_point_label_type, 0, MPI_COMM_WORLD);

    free(localBorderPointsAndLabels);
    free(borderDispls);
    free(borderCounts);

    if (rank == 0) {
        // Write border points and labels to the output file
        for (d = 0; d < factor; d++) {
            if (borderPointsAndLabels[d].x != 0 && borderPointsAndLabels[d].y != 0) {
                fprintf(output, "%f,%f,%d\n", borderPointsAndLabels[d].x, borderPointsAndLabels[d].y, borderPointsAndLabels[d].label);
            }
        }

        otherFactor = N - factor;

        internalPointsAndLabels = calloc(otherFactor, sizeof(struct point_label));
        if (internalPointsAndLabels == NULL) {
            printErrorAllocation();
        }

        // Get non-border points
        getNonBorderPoints(ptrPoints, borderPointsAndLabels, factor, internalPointsAndLabels);
        free(ptrPoints);
    }

    MPI_Bcast(&otherFactor, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Other factor: %d, Rank: %d\n", otherFactor, rank);

    MPI_Bcast(internalPointsAndLabels, otherFactor, mpi_point_label_type, 0, MPI_COMM_WORLD);

    int *innerDispls = calloc(size, sizeof(int));
    if (innerDispls == NULL) {
        printErrorAllocation();
    }

    int *innerCounts = calloc(size, sizeof(int));
    if (innerCounts == NULL) {
        printErrorAllocation();
    }

    int local_other_factor = 0;
    if (rank < remainder) {
        local_other_factor = chunkSize + 1;
    } else {
        local_other_factor = chunkSize;
    }

    int innerDisplacement = 0;
    for (int h = 0; h < size; h++) {
        innerCounts[h] = (h < remainder) ? (chunkSize + 1) : chunkSize;
        innerDispls[h] = innerDisplacement;
        innerDisplacement += innerCounts[h];
    }

    struct point_label *localInternalPointsAndLabels = calloc(local_other_factor, sizeof(struct point_label));
    if (localInternalPointsAndLabels == NULL) {
        printErrorAllocation();
    }

    MPI_Scatterv(internalPointsAndLabels, innerCounts, innerDispls, mpi_point_label_type, localInternalPointsAndLabels, local_other_factor, mpi_point_label_type, 0, MPI_COMM_WORLD);

    // --------------------------------------------------- INNER POINTS CLUSTERING - PHASE 3

    // Get labels for non-border points
    getLabelsNonBorderPoints(factor, borderPointsAndLabels, internalPointsAndLabels, otherFactor);

    MPI_Gatherv(localInternalPointsAndLabels, local_other_factor, mpi_point_label_type, internalPointsAndLabels, innerCounts, innerDispls, mpi_point_label_type, 0, MPI_COMM_WORLD);

    free(localInternalPointsAndLabels);
    free(innerDispls);
    free(innerCounts);

    if (rank == 0) {
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
    
    elapsed_time += MPI_Wtime();
    if (rank == 0) {
        printf("Elapsed time (parallel) = %f\n", elapsed_time);
    }

	// Close output file
	fclose(output);

    // MPI finalization
    MPI_Finalize();

    return 0;
}
