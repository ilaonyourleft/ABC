/*
 ============================================================================
 Name        : ABC_sequential.c
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
#include <string.h>
#include "ABC_sequential_lib.h"

int main(int argc, char *argv[]) {
    // Open input and output files
    FILE *file = fopen("../data/dataset1.csv", "r");
    FILE *output = fopen("../results/results_parallel.txt", "w");

    // Variable declarations
	int i = 0, j = 0, g = 0, d = 0, e = 0, counter = 0, factor = BETA * N, otherFactor = 0;
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

    MPI_Datatype mpi_double_float_type, mpi_triple_float_type;
    MPI_Type_contiguous(2, MPI_FLOAT, &mpi_double_float_type);
    MPI_Type_commit(&mpi_double_float_type);
    MPI_Type_contiguous(3, MPI_FLOAT, &mpi_triple_float_type);
    MPI_Type_commit(&mpi_triple_float_type);

    // --------------------------------------------------- MEMORY ALLOCATIONS
    struct double_float *ptrPoints, *ptrKnnPoint, *ptrPointsAll;
    struct triple_float *ptrBorderPointsAll;
    struct point_label *borderPointsAndLabels;
    float *ptrMeanPoint, *ptrDirectionalAnglesPoint;

    if (rank == 0) {
        ptrPoints = calloc(N, sizeof(struct double_float));
        if (ptrPoints == NULL) {
            printErrorAllocation();
        }

        ptrPointsAll = calloc(N, sizeof(struct double_float));
        if (ptrPointsAll == NULL) {
            printErrorAllocation();
        }

        ptrKnnPoint = calloc(K, sizeof(struct double_float));
        if (ptrKnnPoint == NULL) {
            printErrorAllocation();
        }

        ptrMeanPoint = calloc(2, sizeof(float));
        if (ptrMeanPoint == NULL) {
            printErrorAllocation();
        }

        ptrDirectionalAnglesPoint = calloc(K, sizeof(float));
        if (ptrDirectionalAnglesPoint == NULL) {
            printErrorAllocation();
        }

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

        memcpy(ptrPointsAll, ptrPoints, N * sizeof(struct double_float));
    }

    // Broadcast the ptrPointsAll array to all processes
    MPI_Bcast(ptrPointsAll, N, mpi_double_float_type, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        for (i = 0; i < N; i++) {
            printf("rank: %d, x: %f, y: %f\n", rank, ptrPointsAll[i].x, ptrPointsAll[i].y);
        }
    }
    
    int local_n = N / size;
    int remainder = N % size;
    int *send_counts = calloc(size, sizeof(int));
    int *displs = calloc(size, sizeof(int));
    int displacement = 0;
    for (int i = 0; i < size; i++) {
        send_counts[i] = (i < remainder) ? local_n + 1 : local_n;
        displs[i] = displacement;
        displacement += send_counts[i];
    }

    printf("Rank: %d, send_counts: %d, displs: %d\n", rank, send_counts[rank], displs[rank]);

    struct double_float *local_ptrPoints = calloc(send_counts[rank], sizeof(struct double_float));
    if (local_ptrPoints == NULL) {
        printErrorAllocation();
    }

    MPI_Scatterv(ptrPoints, send_counts, displs, mpi_double_float_type, local_ptrPoints, send_counts[rank], mpi_double_float_type, 0, MPI_COMM_WORLD);

    for (int i = 0; i < send_counts[rank]; i++) {
        printf("Rank %d: %f, %f\n", rank, local_ptrPoints[i].x, local_ptrPoints[i].y);
    }

    struct triple_float *local_ptrBorderPointsAll = calloc(send_counts[rank], sizeof(struct triple_float));
    if (local_ptrBorderPointsAll == NULL) {
        printErrorAllocation();
    }

    // --------------------------------------------------- BORDER POINTS LOOP

    // Read data from the input file and create an array of points
    for (i = 0; i < send_counts[rank]; i++) {
        printf("Rank: %d, x: %f, y: %f\n", rank, local_ptrPoints[i].x, local_ptrPoints[i].y);
        // Find k nearest neighbors for each point and the mean point
        getNeighbors(ptrPointsAll, local_ptrPoints[i].x, local_ptrPoints[i].y, ptrKnnPoint, ptrMeanPoint);

        // Find directional angles between the center, its k nearest neighbors, and the mean point
        for (j = 0; j < K; j++) {
            ptrDirectionalAnglesPoint[j] = getDirectionalAngle(local_ptrPoints[i], ptrMeanPoint, ptrKnnPoint[j]);
        }

        // Find the border points with enclosing angle for each point and border degree
        if (isBorderPoint(getEnclosingAngle(ptrDirectionalAnglesPoint)) == 1) {
            local_ptrBorderPointsAll[counter].x = local_ptrPoints[i].x;
            local_ptrBorderPointsAll[counter].y = local_ptrPoints[i].y;
            local_ptrBorderPointsAll[counter].z = getBorderDegree(ptrDirectionalAnglesPoint);
            ++counter;
        }
    }

    MPI_Gatherv(local_ptrBorderPointsAll, send_counts[rank], mpi_triple_float_type, ptrBorderPointsAll, send_counts, displs, mpi_triple_float_type, 0, MPI_COMM_WORLD);

    free(local_ptrPoints);
    free(local_ptrBorderPointsAll);
        
    if (rank == 0) {
        // Free allocated memory
        free(ptrKnnPoint);
        free(ptrMeanPoint);
        free(ptrDirectionalAnglesPoint);

        // Get factor border points
        getBorderPoints(ptrBorderPointsAll, counter, borderPointsAndLabels);
        free(ptrBorderPointsAll);

        // Check if counter is less than factor
        if (counter < factor) {
            factor = counter;
        }

        // Get label for each border point
        getLabelsBorderPoints(factor, 19000, 3, borderPointsAndLabels);

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
