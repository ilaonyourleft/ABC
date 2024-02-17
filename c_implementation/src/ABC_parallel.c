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

int main(int argc, char *argv[])
{
    // Open input and output files
    FILE *file = fopen("../data/dataset1.csv", "r");
    FILE *output = fopen("../results/results_parallel.txt", "w");

    // Variable declarations
    int i, j, g, d, e, counter = 0, factor = BETA * N, otherFactor = 0, localFactor = 0, localOtherFactor = 0;
    float x, y = 0.0;
    int *ptrLabels, *ptrNonBorderLabels, *localPtrLabels, *localPtrNonBorderLabels;
    float *ptrMeanPoint, *ptrDirectionalAnglesPoint, *ptrEnclosingAnglesPoint, *ptrBorderDegreesPoint, **ptrBorderPointsAll, **ptrPoints, **ptrKnnPoint, **ptrBorderPoints, **ptrNonBorderPoints, **localPtrBorderPoints, **localPtrNonBorderPoints;

    // Variable declarations for parallelization
    int rank, size;

    // Error handling for file opening
    if (file == NULL || output == NULL)
    {
        printf("Could not open file.\n");
        exit(OPEN_FILE_ERROR);
    }

    // MPI initialization
    MPI_Init(&argc, &argv);
    // Process rank in a communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Number of processes in a communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // --------------------------------------------------- MEMORY ALLOCATIONS

    ptrPoints = calloc(N, sizeof(float *));
    if (ptrPoints == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            ptrPoints[i] = calloc(2, sizeof(float));
            if (ptrPoints[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    ptrEnclosingAnglesPoint = calloc(K, sizeof(float));
    if (ptrEnclosingAnglesPoint == NULL)
    {
        printErrorAllocation();
    }

    ptrBorderDegreesPoint = calloc(K, sizeof(float));
    if (ptrBorderDegreesPoint == NULL)
    {
        printErrorAllocation();
    }

    ptrKnnPoint = calloc(K, sizeof(float *));
    if (ptrKnnPoint == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < K; i++)
        {
            ptrKnnPoint[i] = calloc(2, sizeof(float));
            if (ptrKnnPoint[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    ptrMeanPoint = calloc(2, sizeof(float));
    if (ptrMeanPoint == NULL)
    {
        printErrorAllocation();
    }

    ptrDirectionalAnglesPoint = calloc(K, sizeof(float));
    if (ptrDirectionalAnglesPoint == NULL)
    {
        printErrorAllocation();
    }

    ptrBorderPointsAll = calloc(N, sizeof(float *));
    if (ptrBorderPointsAll == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            ptrBorderPointsAll[i] = calloc(3, sizeof(float));
            if (ptrBorderPointsAll[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    ptrBorderPoints = calloc(factor, sizeof(float *));
    if (ptrBorderPoints == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < factor; i++)
        {
            ptrBorderPoints[i] = calloc(2, sizeof(float));
            if (ptrBorderPoints[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    ptrLabels = calloc(factor, sizeof(int));
    if (ptrLabels == NULL)
    {
        printErrorAllocation();
    }

    // Create an array of points
    for (i = 0; i < N; i++)
    {
        fscanf(file, "%f,%f", &x, &y);
        ptrPoints[i][0] = x;
        ptrPoints[i][1] = y;
    }

    // --------------------------------------------------- BORDER POINTS LOOP

    // Read data from the input file and create an array of points
    for (i = 0; i < N; i++)
    {
        // Find k nearest neighbors for each point and the mean point
        getNeighbors(ptrPoints, ptrPoints[i][0], ptrPoints[i][1], ptrKnnPoint, ptrMeanPoint);

        // Find directional angles between the center, its k nearest neighbors, and the mean point
        for (j = 0; j < K; j++)
        {
            ptrDirectionalAnglesPoint[j] = getDirectionalAngle(ptrPoints[i], ptrMeanPoint, ptrKnnPoint[j]);
        }

        // Find the border points with enclosing angle for each point and border degree
        if (isBorderPoint(getEnclosingAngle(ptrDirectionalAnglesPoint)) == 1)
        {
            ptrBorderPointsAll[counter][0] = ptrPoints[i][0];
            ptrBorderPointsAll[counter][1] = ptrPoints[i][1];
            ptrBorderPointsAll[counter][2] = getBorderDegree(ptrDirectionalAnglesPoint);
            ++counter;
        }
    }

    // Free allocated memory
    free(ptrEnclosingAnglesPoint);
    free(ptrBorderDegreesPoint);
    free(ptrKnnPoint);
    free(ptrMeanPoint);
    free(ptrDirectionalAnglesPoint);

    // Get factor border points
    getBorderPoints(ptrBorderPointsAll, counter, ptrBorderPoints);
    free(ptrBorderPointsAll);

    // Check if counter is less than factor
    if (counter < factor)
    {
        factor = counter;
    }

    otherFactor = N - factor;
    ptrNonBorderPoints = calloc(otherFactor, sizeof(float *));
    if (ptrNonBorderPoints == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < otherFactor; i++)
        {
            ptrNonBorderPoints[i] = calloc(2, sizeof(float));
            if (ptrNonBorderPoints[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    ptrNonBorderLabels = calloc(otherFactor, sizeof(int));
    if (ptrNonBorderLabels == NULL)
    {
        printErrorAllocation();
    }

    // Get non-border points
    getNonBorderPoints(ptrPoints, ptrBorderPoints, factor, ptrNonBorderPoints);
    free(ptrPoints);

    if (rank == 0)
    {
        localFactor = factor % size;
        localOtherFactor = otherFactor % size;
    }
    else
    {
        localFactor = factor / size;
        localOtherFactor = otherFactor / size;
    }

    // Memory allocation for local pointers, used for parallelization
    localPtrBorderPoints = calloc(localFactor, sizeof(float *));
    if (localPtrBorderPoints == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < localFactor; i++)
        {
            localPtrBorderPoints[i] = calloc(2, sizeof(float));
            if (localPtrBorderPoints[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    localPtrLabels = calloc(localFactor, sizeof(int));
    if (localPtrLabels == NULL)
    {
        printErrorAllocation();
    }

    localPtrNonBorderPoints = calloc(localOtherFactor, sizeof(float *));
    if (localPtrNonBorderPoints == NULL)
    {
        printErrorAllocation();
    }
    else
    {
        for (i = 0; i < localOtherFactor; i++)
        {
            localPtrNonBorderPoints[i] = calloc(2, sizeof(float));
            if (localPtrNonBorderPoints[i] == NULL)
            {
                printErrorAllocation();
            }
        }
    }

    localPtrNonBorderLabels = calloc(localOtherFactor, sizeof(int));
    if (localPtrNonBorderLabels == NULL)
    {
        printErrorAllocation();
    }

    printf("Process %d, localFactor = %d, localOtherFactor = %d\n", rank, localFactor, localOtherFactor);

    // Scatter the border points data and labels
    MPI_Scatter(ptrBorderPoints, localFactor, MPI_FLOAT, localPtrBorderPoints, localFactor, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(ptrLabels, localFactor, MPI_INT, localPtrLabels, localFactor, MPI_INT, 0, MPI_COMM_WORLD);

    // Scatter the non-border points data and labels
    MPI_Scatter(ptrNonBorderPoints, localOtherFactor, MPI_FLOAT, localPtrNonBorderPoints, localOtherFactor, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(ptrNonBorderLabels, localOtherFactor, MPI_INT, localPtrNonBorderLabels, localOtherFactor, MPI_INT, 0, MPI_COMM_WORLD);

    // Each process computes labels for its portion of the data
    getLabelsBorderPoints(localPtrBorderPoints, localFactor, 19000, 3, localPtrLabels);
    getLabelsNonBorderPoints(localPtrBorderPoints, localFactor, localPtrLabels, localPtrNonBorderPoints, localOtherFactor, localPtrNonBorderLabels);

    // Gather the results back to the process with rank 0
    MPI_Gather(localPtrBorderPoints, localFactor, MPI_FLOAT, ptrBorderPoints, localFactor, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(localPtrLabels, localFactor, MPI_INT, ptrLabels, localFactor, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(localPtrNonBorderPoints, localOtherFactor, MPI_FLOAT, ptrNonBorderPoints, localOtherFactor, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(localPtrNonBorderLabels, localOtherFactor, MPI_INT, ptrNonBorderLabels, localOtherFactor, MPI_INT, 0, MPI_COMM_WORLD);

    // process the gathered results
    if (rank == 0)
    {
        // Write border points and labels to the output file
        for (d = 0; d < factor; d++)
        {
            if (ptrBorderPoints[d][0] != 0 && ptrBorderPoints[d][1] != 0)
            {
                fprintf(output, "%f,%f,%d\n", ptrBorderPoints[d][0], ptrBorderPoints[d][1], ptrLabels[d]);
            }
        }

        // Write non border points and labels to the output file
        for (e = 0; e < otherFactor; e++)
        {
            if (ptrNonBorderPoints[e][0] != 0 && ptrNonBorderPoints[e][1] != 0)
            {
                fprintf(output, "%f,%f,%d\n", ptrNonBorderPoints[e][0], ptrNonBorderPoints[e][1], ptrNonBorderLabels[e]);
            }
        }

        // Free remaining allocated memory
        free(ptrBorderPoints);
        free(ptrLabels);
        free(ptrNonBorderPoints);
        free(ptrNonBorderLabels);
        free(localPtrBorderPoints);
        free(localPtrNonBorderPoints);
        free(localPtrLabels);
        free(localPtrNonBorderLabels);
    }

    // Close output file
    fclose(output);

    // MPI finalization
    MPI_Finalize();

    return 0;
}
