/*
 ============================================================================
 Name        : ABC_sequential.c
 Author      : Ilaria Malinconico
 Version     :
 Copyright   : Your copyright notice
 Description : ABC sequential implementation
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ABC_sequential_lib.h"

int main(int argc, char *argv[]) {
	// puts("Angle Based Clustering Program - Sequential");

	// Open input and output files
	FILE *file = fopen("../data/dataset1.csv", "r");
	FILE *output = fopen("../results/results_sequential.txt", "w");

	// Variable declarations
	int i, j, g, d, e, counter = 0, factor = BETA * N, otherFactor;
	float x, y = 0.0;
	int *ptrLabels, *ptrNonBorderLabels;
	float *ptrMeanPoint, *ptrDirectionalAnglesPoint, *ptrEnclosingAnglesPoint, *ptrBorderDegreesPoint, **ptrBorderPointsAll, **ptrPoints, **ptrKnnPoint, **ptrBorderPoints, **ptrNonBorderPoints;

	// Error handling for file opening
	if (file == NULL || output == NULL) {
		printf("Could not open file.\n");
		exit(OPEN_FILE_ERROR);
	}

	// --------------------------------------------------- MEMORY ALLOCATIONS

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

	ptrMeanPoint = calloc(2, sizeof(float));
	if (ptrMeanPoint == NULL) {
		printErrorAllocation();
	}

	ptrDirectionalAnglesPoint = calloc(K, sizeof(float));
	if (ptrDirectionalAnglesPoint == NULL) {
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

	ptrBorderPoints = calloc(factor, sizeof(float *));
	if (ptrBorderPoints == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < factor; i++) {
			ptrBorderPoints[i] = calloc(2, sizeof(float));
			if (ptrBorderPoints[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	ptrLabels = calloc(factor, sizeof(int));
	if (ptrLabels == NULL) {
		printErrorAllocation();
	}

	// Create an array of points
	for (i = 0; i < N; i++) {
		fscanf(file, "%f,%f", &x, &y);
		ptrPoints[i][0] = x;
		ptrPoints[i][1] = y;
	}

	// --------------------------------------------------- BORDER POINTS LOOP
	
	// Read data from the input file and create an array of points
	for (i = 0; i < N; i++) {
		// Find k nearest neighbors for each point and the mean point
		getNeighbors(ptrPoints, ptrPoints[i][0], ptrPoints[i][1], ptrKnnPoint, ptrMeanPoint);

		// Write the k nearest neighbors to the output file
		//for (g = 0; g < K; g++) {
			//if (ptrKnnPoint[g][0] != 0 && ptrKnnPoint[g][1] != 0) {
				//fprintf (output, "%d, %d\n", ptrKnnPoint[g][0], ptrKnnPoint[g][1]);
			//}
		//}
		//fprintf(output, "\n");

		// Find directional angles between the center, its k nearest neighbors, and the mean point
		for (j = 0; j < K; j++) {
			ptrDirectionalAnglesPoint[j] = getDirectionalAngle(ptrPoints[i], ptrMeanPoint, ptrKnnPoint[j]);
		}

		// Find the border points with enclosing angle for each point and border degree
		if (isBorderPoint(getEnclosingAngle(ptrDirectionalAnglesPoint)) == 1) {
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
	if (counter < factor) {
		factor = counter;
	}

	// Get label for each border point
	getLabelsBorderPoints(ptrBorderPoints, factor, 19000, 3, ptrLabels);

	// Write border points and labels to the output file
	for (d = 0; d < factor; d++) {
		if (ptrBorderPoints[d][0] != 0 && ptrBorderPoints[d][1] != 0) {
			fprintf (output, "%f,%f,%d\n", ptrBorderPoints[d][0], ptrBorderPoints[d][1], ptrLabels[d]);
		}
	}

	otherFactor = N - factor;
	ptrNonBorderPoints = calloc(otherFactor, sizeof(float *));
	if (ptrNonBorderPoints == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < otherFactor; i++) {
			ptrNonBorderPoints[i] = calloc(2, sizeof(float));
			if (ptrNonBorderPoints[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	ptrNonBorderLabels = calloc(otherFactor, sizeof(int));
	if (ptrNonBorderLabels == NULL) {
		printErrorAllocation();
	}

	// Get non-border points
	getNonBorderPoints(ptrPoints, ptrBorderPoints, factor, ptrNonBorderPoints);
	free(ptrPoints);

	// Get labels for non-border points
	getLabelsNonBorderPoints(ptrBorderPoints, factor, ptrLabels, ptrNonBorderPoints, otherFactor, ptrNonBorderLabels);

	// Write non border points and labels to the output file
	for (e = 0; e < otherFactor; e++) {
		if (ptrNonBorderPoints[e][0] != 0 && ptrNonBorderPoints[e][1] != 0) {
			fprintf (output, "%f,%f,%d\n", ptrNonBorderPoints[e][0], ptrNonBorderPoints[e][1], ptrNonBorderLabels[e]);
		}
	}
	
	// Free remaining allocated memory
	free(ptrBorderPoints);
	free(ptrLabels);
	free(ptrNonBorderPoints);
	free(ptrNonBorderLabels);

	// Close output file
	fclose(output);
	
	return 0;
}
