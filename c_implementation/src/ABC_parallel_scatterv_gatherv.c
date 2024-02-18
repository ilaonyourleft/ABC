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

    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_INT};
    MPI_Aint offsets[3];

    offsets[0] = offsetof(struct point_label, x);
    offsets[1] = offsetof(struct point_label, y);
    offsets[2] = offsetof(struct point_label, label);

    MPI_Datatype mpi_point_label_type;
    MPI_Type_create_struct(3, blocklengths, offsets, types, &mpi_point_label_type);
    MPI_Type_commit(&mpi_point_label_type);

    // --------------------------------------------------- MEMORY ALLOCATIONS
    struct double_float *ptrPoints, *ptrKnnPoint;
    struct triple_float *ptrBorderPointsAll;
    struct point_label *borderPointsAndLabels;
    float *ptrMeanPoint, *ptrDirectionalAnglesPoint;
    int *sendcounts, *displs;

    if (rank == 0) {
        ptrPoints = calloc(N, sizeof(struct double_float));
        if (ptrPoints == NULL) {
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

        // --------------------------------------------------- BORDER POINTS LOOP

        // Read data from the input file and create an array of points
        for (i = 0; i < N; i++) {
            // Find k nearest neighbors for each point and the mean point
            getNeighbors(ptrPoints, ptrPoints[i].x, ptrPoints[i].y, ptrKnnPoint, ptrMeanPoint);

            // Find directional angles between the center, its k nearest neighbors, and the mean point
            for (j = 0; j < K; j++) {
                ptrDirectionalAnglesPoint[j] = getDirectionalAngle(ptrPoints[i], ptrMeanPoint, ptrKnnPoint[j]);
            }

            // Find the border points with enclosing angle for each point and border degree
            if (isBorderPoint(getEnclosingAngle(ptrDirectionalAnglesPoint)) == 1) {
                ptrBorderPointsAll[counter].x = ptrPoints[i].x;
                ptrBorderPointsAll[counter].y = ptrPoints[i].y;
                ptrBorderPointsAll[counter].z = getBorderDegree(ptrDirectionalAnglesPoint);
                ++counter;
            }
        }

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

        // Scatter the data to all processes
        int *sendcounts = calloc(size, sizeof(int));
        int *displs = calloc(size, sizeof(int));
        int remainder = factor % size;
        for (int i = 0; i < size; i++) {
            sendcounts[i] = factor / size + (i < remainder ? 1 : 0);
            displs[i] = (i > 0) ? displs[i - 1] + sendcounts[i - 1] : 0;
        }

        MPI_Scatterv(borderPointsAndLabels, sendcounts, displs, mpi_point_label_type, MPI_IN_PLACE, 0, mpi_point_label_type, 0, MPI_COMM_WORLD);

        free(sendcounts);
        free(displs);
    }

    MPI_Bcast(&factor, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int local_factor = factor / size + (rank < factor % size ? 1 : 0);
    printf("Rank %d has %d points\n", rank, local_factor);
    
    struct point_label *local_borderPointsAndLabels = calloc(local_factor, sizeof(struct point_label));
    if (local_borderPointsAndLabels == NULL) {
        printErrorAllocation();
    }

    MPI_Scatterv(NULL, NULL, NULL, mpi_point_label_type, local_borderPointsAndLabels, local_factor, mpi_point_label_type, 0, MPI_COMM_WORLD);

    for (d = 0; d < local_factor; d++) {
        printf("Rank %d: %f, %f, %d\n", rank, local_borderPointsAndLabels[d].x, local_borderPointsAndLabels[d].y, local_borderPointsAndLabels[d].label);
    }

    // Call the original getLabelsBorderPoints function on the local data
    getLabelsBorderPoints(local_factor, 19000, 3, local_borderPointsAndLabels);

    // Send the computed results back to the root process
    MPI_Gatherv(local_borderPointsAndLabels, local_factor, mpi_point_label_type, borderPointsAndLabels, sendcounts, displs, mpi_point_label_type, 0, MPI_COMM_WORLD);
    free(local_borderPointsAndLabels);

    if (rank == 0) {
        free(sendcounts);
        free(displs);
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
