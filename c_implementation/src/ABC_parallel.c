#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include "ABC_sequential_lib.h"

int main(int argc, char *argv[]) {
    // Open input and output files
	FILE *file = fopen("../data/dataset1.csv", "r");
	FILE *output = fopen("../results/results_sequential.txt", "w");

    int rank, size;
    int start, end, local_counter = 0, global_counter;

    // Variable declarations
	int i, j, g, d, e, counter = 0, factor = BETA * N, otherFactor;
	float x, y = 0.0;
	int *ptrLabels, *ptrNonBorderLabels;
	float *ptrMeanPoint, *ptrDirectionalAnglesPoint, *local_ptrDirectionalAnglesPoint, *ptrEnclosingAnglesPoint, *ptrBorderDegreesPoint, **ptrBorderPointsAll, **ptrPoints, **local_ptrPoints, **ptrKnnPoint, **local_ptrKnnPoint, **ptrBorderPoints, **ptrNonBorderPoints;

    // inizializzazione per gestire le chiamate alla libreria
    MPI_Init(&argc, &argv);
    // determina il rank del processo in un comunicatore
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // determina il numero di processi in un comunicatore
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate the range of loop iterations for each process
    int chunk_size = N / size;
    start = rank * chunk_size;
    end = (rank == size - 1) ? N : start + chunk_size;

    ptrPoints = calloc(N, sizeof(float *));
	if (ptrPoints == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < N; i++) {
			ptrPoints[i] = calloc(2, sizeof(float));
			if (ptrPoints[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

    local_ptrPoints = calloc(N, sizeof(float *));
	if (local_ptrPoints == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < N; i++) {
			local_ptrPoints[i] = calloc(2, sizeof(float));
			if (local_ptrPoints[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	ptrEnclosingAnglesPoint = calloc(K, sizeof(float));
	if (ptrEnclosingAnglesPoint == NULL) {
		printErrorAllocation();
	}

	ptrBorderDegreesPoint = calloc(K, sizeof(float));
	if (ptrBorderDegreesPoint == NULL) {
		printErrorAllocation();
	}

	ptrKnnPoint = calloc(K, sizeof(float *));
	if (ptrKnnPoint == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < K; i++) {
			ptrKnnPoint[i] = calloc(2, sizeof(float));
			if (ptrKnnPoint[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	local_ptrKnnPoint = calloc(K, sizeof(float *));
	if (local_ptrKnnPoint == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < K; i++) {
			local_ptrKnnPoint[i] = calloc(2, sizeof(float));
			if (local_ptrKnnPoint[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	ptrMeanPoint = calloc(2, sizeof(float));
	if (ptrMeanPoint == NULL) {
		printErrorAllocation();
	}

	ptrDirectionalAnglesPoint = calloc(K, sizeof(float));
	if (ptrDirectionalAnglesPoint == NULL) {
		printErrorAllocation();
	}

	local_ptrDirectionalAnglesPoint = calloc(K, sizeof(float));
	if (local_ptrDirectionalAnglesPoint == NULL) {
		printErrorAllocation();
	}

	ptrBorderPointsAll = calloc(N, sizeof(float *));
	if (ptrBorderPointsAll == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < N; i++) {
			ptrBorderPointsAll[i] = calloc(3, sizeof(float));
			if (ptrBorderPointsAll[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	// Create an array of points
	for (i = 0; i < N; i++) {
		fscanf(file, "%f,%f", &x, &y);
		ptrPoints[i][0] = x;
		ptrPoints[i][1] = y;
	}

    // Scatter the ptrPoints array among processes
    MPI_Scatter(ptrPoints, chunk_size * 2, MPI_FLOAT, &local_ptrPoints[0][0], chunk_size * 2, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Scatter the ptrKnnPoint array among processes
    MPI_Scatter(ptrKnnPoint, K * 2, MPI_FLOAT, &local_ptrKnnPoint[0][0], K * 2, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Parallelize the loop for border points
    for (int i = start; i < end; i++) {
        // Find k nearest neighbors for each point and the mean point
        getNeighbors(local_ptrPoints, local_ptrPoints[i][0], local_ptrPoints[i][1], local_ptrKnnPoint, ptrMeanPoint);

        // Find directional angles between the center, its k nearest neighbors, and the mean point
        for (int j = 0; j < K; j++) {
            local_ptrDirectionalAnglesPoint[j] = getDirectionalAngle(local_ptrPoints[i], ptrMeanPoint, local_ptrKnnPoint[j]);
        }

        // Find the border points with enclosing angle for each point and border degree
        if (isBorderPoint(getEnclosingAngle(local_ptrDirectionalAnglesPoint)) == 1) {
            ptrBorderPointsAll[counter][0] = local_ptrPoints[i][0];
            ptrBorderPointsAll[counter][1] = local_ptrPoints[i][1];
            ptrBorderPointsAll[counter][2] = getBorderDegree(local_ptrDirectionalAnglesPoint);
            ++local_counter;
        }
    }

    // operazione di riduzione fra tutti i processi in un comunicatore e ottiene il risultato sul processo root
    // Use MPI_Allreduce to get the sum of local_counter from all processes
    MPI_Allreduce(&local_counter, &global_counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    if (rank == 0) {
        printf("Total border points: %d\n", global_counter);
    }


    // Free memory and finalize MPI
    free(ptrPoints);
    free(local_ptrPoints);
	free(ptrKnnPoint);
    free(local_ptrKnnPoint);
	free(ptrDirectionalAnglesPoint);
    free(local_ptrDirectionalAnglesPoint);
	free(ptrEnclosingAnglesPoint);
	free(ptrBorderDegreesPoint);
	free(ptrMeanPoint);
    free(ptrBorderPointsAll);

    MPI_Finalize();

    return 0;
}