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
	// Open input and output files
	FILE *file = fopen("../data/dataset1.csv", "r");
	FILE *output = fopen("../results/results_sequential.txt", "w");

	// Variable declarations
	int i = 0, j = 0, g = 0, d = 0, e = 0, counter = 0, factor = BETA * N, otherFactor = 0;
	float x, y = 0.0;

	// Error handling for file opening
	if (file == NULL || output == NULL) {
		printf("Could not open file.\n");
		exit(OPEN_FILE_ERROR);
	}

	// --------------------------------------------------- MEMORY ALLOCATIONS

	struct double_float *ptrPoints = calloc(N, sizeof(struct double_float));
	if (ptrPoints == NULL) {
		printErrorAllocation();
	}

	struct double_float *ptrKnnPoint = calloc(K, sizeof(struct double_float));
	if (ptrKnnPoint == NULL) {
		printErrorAllocation();
	}

	float *ptrMeanPoint = calloc(2, sizeof(float));
	if (ptrMeanPoint == NULL) {
		printErrorAllocation();
	}

	float *ptrDirectionalAnglesPoint = calloc(K, sizeof(float));
	if (ptrDirectionalAnglesPoint == NULL) {
		printErrorAllocation();
	}

	struct triple_float *ptrBorderPointsAll = calloc(N, sizeof(struct triple_float));
	if (ptrBorderPointsAll == NULL) {
		printErrorAllocation();
	}

	// Create an array of points
	for (i = 0; i < N; i++) {
		fscanf(file, "%f,%f", &x, &y);
		ptrPoints[i].x = x;
		ptrPoints[i].y = y;
	}

	struct point_label *borderPointsAndLabels = calloc(factor, sizeof(struct point_label));
	if (borderPointsAndLabels == NULL) {
		printErrorAllocation();
	}

	// --------------------------------------------------- BORDER POINTS LOOP
	
	// Read data from the input file and create an array of points
	for (i = 0; i < N; i++) {
		// Find k nearest neighbors for each point and the mean point
		getNeighbors(ptrPoints, ptrPoints[i].x, ptrPoints[i].y, ptrKnnPoint, ptrMeanPoint);
		for (j = 0; j < K; j++) {
			printf("Rank 0: %f, %f, %f, %f, %f\n", ptrPoints[i].x, ptrPoints[i].y, ptrKnnPoint[j].x, ptrKnnPoint[j].y, ptrMeanPoint[j]);
		}
		

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

	// Close output file
	fclose(output);
	
	return 0;
}
