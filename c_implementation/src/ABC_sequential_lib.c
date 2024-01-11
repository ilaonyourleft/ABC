/*
 * ABC_sequential.h
 *
 *  Created on: 25 giu 2022
 *      Author: ilaria
 */

#ifndef ABC_SEQUENTIAL_H_
#define ABC_SEQUENTIAL_H_
#endif /* ABC_SEQUENTIAL_H_ */

#define K 12  // 12
#define BETA 0.3  // 0.2
#define N 500  // 500

#define OPEN_FILE_ERROR -1
#define MEMORY_ALLOCATION_ERROR -2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

// MISCELLANEOUS

/**
 * Print an error message and exit the program if memory allocation fails.
 */
void printErrorAllocation() {
	printf("Could not allocate memory to pointer.\n");
	exit(MEMORY_ALLOCATION_ERROR);
}

// KNN methods

/**
 * Calculate the Euclidean distance between two points in a two-dimensional space.
 *
 * @param x1 The x-coordinate of the first point.
 * @param y1 The y-coordinate of the first point.
 * @param x2 The x-coordinate of the second point.
 * @param y2 The y-coordinate of the second point.
 * @return The Euclidean distance between the two points.
 */
float euclideanDistance(float x1, float y1, float x2, float y2) {
	return sqrtf(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

/**
 * Sorts an array of distances in ascending order based on the third element. (Bubble sort)
 *
 * @param distancesPoints The array of distances to be sorted.
 */
void sortArrayDistances(float **distancesPoints) {
	float tmp[3];

	for (int i = 0; i < N-1; i++) {
		for (int j = i + 1; j < N; j++) {
			if (distancesPoints[i][2] > distancesPoints[j][2]) {
				tmp[0] = distancesPoints[i][0];
				tmp[1] = distancesPoints[i][1];
				tmp[2] = distancesPoints[i][2];
				distancesPoints[i][0] = distancesPoints[j][0];
				distancesPoints[j][0] = tmp[0];
				distancesPoints[i][1] = distancesPoints[j][1];
				distancesPoints[j][1] = tmp[1];
				distancesPoints[i][2] = distancesPoints[j][2];
				distancesPoints[j][2] = tmp[2];
			}
		}
	}
}

/**
 * Get the nearest neighbors of a given point based on Euclidean distances.
 *
 * @param points     The array of points.
 * @param x          The x-coordinate of the reference point.
 * @param y          The y-coordinate of the reference point.
 * @param knn        The array to store the nearest neighbors.
 * @param meanPoint  The array to store the mean point of the nearest neighbors.
 */
// non permette di ritornare un array, ma si può ritornare il puntatore all'array specificandone il nome senza indice
void getNeighbors(float **points, float x, float y, float **knn, float *meanPoint) {
	// non ritorna l'indirizzo di una variabile locale all'esterno della funzione, quindi serve static nella definizione della variabile locale
	float **distances;
	float tmp_x = 0.00, tmp_y = 0.00;
	int i, j;

	distances = calloc(N, sizeof(float *));
	if (distances == NULL) {
		printErrorAllocation();
	} else {
		for (int i = 0; i < N; i++) {
			distances[i] = calloc(3, sizeof(float));
			if (distances[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	for (i = 0; i < N; i++) {
		float distance = euclideanDistance(x, y, points[i][0], points[i][1]);
		for (j = 0; j < 3; j++) {
			if (j != 2) {
				distances[i][j] = points[i][j];
			} else {
				distances[i][j] = distance;
			}
		}
	}

	sortArrayDistances(distances);

	for (i = 0; i < K; i++) {
		knn[i][0] = distances[i+1][0];
		knn[i][1] = distances[i+1][1];
	}

	free(distances);

	for (i = 0; i < K; i++) {
		tmp_x += knn[i][0];
		tmp_y += knn[i][1];
	}

	meanPoint[0] = tmp_x / K;
	meanPoint[1] = tmp_y / K;
}

// ENCLOSING ANGLES

/**
 * Calculate the directional angle between the line segments formed by the center, mean point, and neighbor.
 *
 * @param center     The coordinates of the center point.
 * @param meanPoint  The coordinates of the mean point.
 * @param neighbor   The coordinates of the neighbor point.
 * @return The directional angle in degrees.
 */
float getDirectionalAngle(float *center, float *meanPoint, float *neighbor) {
	float uX = 0.00, uY = 0.00, vX = 0.00, vY = 0.00, directionalAngle = 0.00, directionalAngleDegree = 0.00;

	// mean point - center
	uX = meanPoint[0] - center[0];
	uY = meanPoint[1] - center[1];

	// neighbor - center
	vX = neighbor[0] - center[0];
	vY = neighbor[1] - center[1];

	directionalAngle = atan2(vY, vX) - atan2(uY, uX);
	directionalAngleDegree = directionalAngle * (180 / M_PI);

	if (directionalAngleDegree < 0) {
		directionalAngleDegree += 360;
	}

	return directionalAngleDegree;
}

/**
 * Find the number of directional angles that are greater than or equal to 180 degrees.
 *
 * @param directionalAngles The array of directional angles.
 * @return The number of directional angles greater than or equal to 180 degrees.
 */
int findSize(float *directionalAngles) {
	int counter = 0;
	for (int i = 0; i < K; i++) {
		if (directionalAngles[i] >= 180) {
			++counter;
		}
	}
	return counter;
}

/**
 * Calculate the enclosing angle based on the directional angles greater than or equal to 180 degrees.
 *
 * @param directionalAngles The array of directional angles.
 * @return The enclosing angle in degrees.
 */
float getEnclosingAngle(float *directionalAngles) {
	int sizeTmp = findSize(directionalAngles), i;
	float tmpDirectionalAngles[sizeTmp], enclosingAngle = 0.00;
	int counter = 0;
	for (i = 0; i < K; i++) {
		if (directionalAngles[i] >= 180) {
			tmpDirectionalAngles[counter] = directionalAngles[i];
			++counter;
		}
	}

	float minimumAngle = tmpDirectionalAngles[0];
	for (i = 0; i < sizeTmp; i++) {
		if (minimumAngle > tmpDirectionalAngles[i]) {
			minimumAngle = tmpDirectionalAngles[i];
		}
	}

	enclosingAngle = 360 - minimumAngle;
	return enclosingAngle;
}

/**
 * Calculate the border degree based on the minimum directional angle.
 *
 * @param directionalAngles The array of directional angles.
 * @return The border degree in degrees.
 */
float getBorderDegree(float *directionalAngles) {
	float minimumAngle = directionalAngles[0], borderDegree = 0.00;
	for (int i = 0; i < K; i++) {
		if (minimumAngle > directionalAngles[i]) {
			minimumAngle = directionalAngles[i];
		}
	}

	borderDegree = 360 - minimumAngle;
	return borderDegree;
}

/**
 * Check if a point is a border point based on the enclosing angle.
 *
 * @param enclosingAngle The enclosing angle in degrees.
 * @return 1 if the point is a border point, 0 otherwise.
 */
int isBorderPoint(float enclosingAngle) {
	if (enclosingAngle < 60) {
		return 1;
	} else {
		return 0;
	}
}

/**
 * Sort the array of border degrees in descending order based on the third element of each row. (Bubble sort)
 *
 * @param borderDegrees The 2D array of border degrees.
 */
void sortArrayBorderDegrees(float **borderDegrees) {
	float tmp[3];

	for (int i = 0; i < N-1; i++) {
		for (int j = i + 1; j < N; j++) {
			if (borderDegrees[i][2] < borderDegrees[j][2]) {
				tmp[0] = borderDegrees[i][0];
				tmp[1] = borderDegrees[i][1];
				tmp[2] = borderDegrees[i][2];
				borderDegrees[i][0] = borderDegrees[j][0];
				borderDegrees[j][0] = tmp[0];
				borderDegrees[i][1] = borderDegrees[j][1];
				borderDegrees[j][1] = tmp[1];
				borderDegrees[i][2] = borderDegrees[j][2];
				borderDegrees[j][2] = tmp[2];
			}
		}
	}
}

/**
 * Get the border points from the sorted array of border points.
 *
 * @param borderPointsAll The 2D array of border points.
 * @param sizeArray       The size of the array.
 * @param borderPoints    The 2D array to store the border points.
 */
void getBorderPoints(float **borderPointsAll, int sizeArray, float **borderPoints) {
	sortArrayBorderDegrees(borderPointsAll);
	for (int i = 0; i < sizeArray; i++) {
		borderPoints[i][0] = borderPointsAll[i][0];
		borderPoints[i][1] = borderPointsAll[i][1];
	}
}

// DBSCAN

/**
 * Calculate the module (magnitude) of a vector given its x and y components.
 *
 * @param x The x component of the vector.
 * @param y The y component of the vector.
 * @return The module of the vector.
 */
float moduleVector(float x, float y) {
	return sqrt(pow(x, 2) + pow(y, 2));
}

/**
 * Calculate the modified distance between two points based on the direction angle between their vectors.
 *
 * @param aX The x component of the first point.
 * @param aY The y component of the first point.
 * @param bX The x component of the second point.
 * @param bY The y component of the second point.
 * @return The modified distance between the two points.
 */
float directionAngleModifiedDistanceFunction(float aX, float aY, float bX, float bY) {
	double product = 0.00;
	float angleBetweenVectors = 0.00;
	product = (double) aX * (double) bX + (double) aY * (double) bY;
	angleBetweenVectors = product / (moduleVector(aX, aY) * moduleVector(bX, bY));
	return euclideanDistance(aX, aY, bX, bY) * (1 + ((0.5 - 1) / M_PI) * angleBetweenVectors);
}

/**
 * Perform a region query to find the neighboring points within a specified distance threshold.
 *
 * @param borderPoints The 2D array of border points.
 * @param neighbors    The 2D array to store the neighboring points.
 * @param factor       The number of border points.
 * @param x            The x component of the point to query.
 * @param y            The y component of the point to query.
 * @param epsilon      The distance threshold.
 * @return The number of neighboring points found.
 */
int regionQuery(float **borderPoints, float **neighbors, int factor, float x, float y, int epsilon) {
	int counter = 0;
	for (int i = 0; i < factor; i++) {
		float disComputed = directionAngleModifiedDistanceFunction(x, y, borderPoints[i][0], borderPoints[i][1]);
		if (disComputed < epsilon) {
			neighbors[counter][0] = i;
			neighbors[counter][1] = borderPoints[i][0];
			neighbors[counter][2] = borderPoints[i][1];
			++counter;
		}
	}
	return counter;
}

/**
 * Check if a point is already a neighbor in the list of neighbors.
 *
 * @param neighbors    The 2D array of neighbors.
 * @param lenNeighbors The length of the neighbors array.
 * @param index        The index of the point to check.
 * @return 1 if the point is already a neighbor, 0 otherwise.
 */
int checkIfAlreadyNeighbor(float **neighbors, int lenNeighbors, float *index) {
	int flag = 0;
	for (int i = 0; i < lenNeighbors; i++) {
		if (neighbors[i][0] >= 0) {
			if (neighbors[i][1] == index[1] && neighbors[i][2] == index[2]) {
				flag = 1;
				break;
			}
		}
	}
	return flag;
}

/**
 * Grow a cluster by expanding it with neighboring points.
 *
 * @param borderPoints      The 2D array of border points.
 * @param factor            The number of border points.
 * @param labels            The array to store the cluster labels.
 * @param index             The index of the initial point.
 * @param x                 The x component of the initial point.
 * @param y                 The y component of the initial point.
 * @param neighbors         The 2D array of neighboring points.
 * @param lenNeighbors      The length of the neighbors array.
 * @param clusterId         The ID of the cluster.
 * @param epsilon           The distance threshold for neighboring points.
 * @param minNumberPoints   The minimum number of points required to form a cluster.
 */
void growCluster(float **borderPoints, int factor, int *labels, int index, float x, float y, float **neighbors, int lenNeighbors, int clusterId, int epsilon, int minNumberPoints) {
	labels[index] = clusterId;
	int counter = 0, i, j, neighborsIncrement;
	float **ptrNextNeighbors;

	ptrNextNeighbors = calloc(factor, sizeof(float *));
	if (ptrNextNeighbors == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < factor; i++) {
			ptrNextNeighbors[i] = calloc(3, sizeof(float));
			if (ptrNextNeighbors[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	while (counter < lenNeighbors) {
		int lenNextNeighbors = 0;
		if (counter != 0) {
			for (j = 0; j < factor; j++) {
				ptrNextNeighbors[j][0] = 0;
				ptrNextNeighbors[j][1] = 0;
				ptrNextNeighbors[j][2] = 0;
			}
		}
		int next_index = neighbors[counter][0];
		if (labels[next_index] == -1) {
			labels[next_index] = clusterId;
		} else if (labels[next_index] == 0) {
			labels[next_index] = clusterId;
			lenNextNeighbors = regionQuery(borderPoints, ptrNextNeighbors, factor, borderPoints[next_index][0], borderPoints[next_index][1], epsilon);
			if (lenNextNeighbors >= minNumberPoints) {
				neighborsIncrement = 0;
				for (i = 0; i < lenNextNeighbors; i++) {
					if (checkIfAlreadyNeighbor(neighbors, lenNeighbors, ptrNextNeighbors[i]) == 0) {
						neighbors[lenNeighbors + neighborsIncrement][0] = ptrNextNeighbors[i][0];
						neighbors[lenNeighbors + neighborsIncrement][1] = ptrNextNeighbors[i][1];
						neighbors[lenNeighbors + neighborsIncrement][2] = ptrNextNeighbors[i][2];
						++lenNeighbors;
						++neighborsIncrement;
					}
				}
			}
		}
		++counter;
	}
	free(ptrNextNeighbors);
}

/**
 * Assign cluster labels to the border points based on density connectivity.
 *
 * @param borderPoints      The 2D array of border points.
 * @param factor            The number of border points.
 * @param epsilon           The distance threshold for neighboring points.
 * @param minNumberPoints   The minimum number of points required to form a cluster.
 * @param labels            The array to store the cluster labels.
 */
void getLabelsBorderPoints(float **borderPoints, int factor, int epsilon, int minNumberPoints, int *labels) {
	int clusterId = 0, i, j;
	float **ptrNeighbors;

	ptrNeighbors = calloc(factor, sizeof(float *));
	if (ptrNeighbors == NULL) {
		printErrorAllocation();
	} else {
		for (i = 0; i < factor; i++) {
			ptrNeighbors[i] = calloc(3, sizeof(float));
			if (ptrNeighbors[i] == NULL) {
				printErrorAllocation();
			}
		}
	}

	for (i = 0; i < factor; i++) {
		int lenNeighbors = 0;
		if (i != 0) {
			for (j = 0; j < factor; j++) {
				ptrNeighbors[j][0] = 0;
				ptrNeighbors[j][1] = 0;
				ptrNeighbors[j][2] = 0;
			}
		}
		if (labels[i] == 0) {
			lenNeighbors = regionQuery(borderPoints, ptrNeighbors, factor, borderPoints[i][0], borderPoints[i][1], epsilon);
			if (lenNeighbors < minNumberPoints) {
				labels[i] = -1;
			} else {
				++clusterId;
				growCluster(borderPoints, factor, labels, i, borderPoints[i][0], borderPoints[i][1], ptrNeighbors, lenNeighbors, clusterId, epsilon, minNumberPoints);
			}
		}
	}

	free(ptrNeighbors);
}

// CLUSTER

/**
 * Check if a point is a border point.
 *
 * @param borderPoints    The 2D array of border points.
 * @param factor          The number of border points.
 * @param x               The x component of the point to check.
 * @param y               The y component of the point to check.
 * @return 1 if the point is a border point, 0 otherwise.
 */
int checkIfBorderPoint(float **borderPoints, int factor, float x, float y) {
	for (int i = 0; i < factor; i++) {
		if (borderPoints[i][0] == x && borderPoints[i][1] == y) {
			return 1;
		}
	}
	return 0;
}

/**
 * Get the non-border points from the given points array.
 *
 * @param points            The 2D array of points.
 * @param borderPoints      The 2D array of border points.
 * @param factor            The number of border points.
 * @param nonBorderPoints   The 2D array to store the non-border points.
 */
void getNonBorderPoints(float **points, float **borderPoints, int factor, float **nonBorderPoints) {
	int counter = 0;
	for (int i = 0; i < N; i++) {
		if (checkIfBorderPoint(borderPoints, factor, points[i][0], points[i][1]) == 0) {
			nonBorderPoints[counter][0] = points[i][0];
			nonBorderPoints[counter][1] = points[i][1];
			++counter;
		}
	}
}

/**
 * Find the minimum distance from the given distances array.
 *
 * @param distances   The 2D array of distances.
 * @param factor      The number of distances.
 * @return The minimum distance found.
 */
float findMinimumDistance(float **distances, int factor) {
	float minimumDistance = 0.0;
	for (int i = 0; i < factor; i++) {
		if (i == 0) {
			minimumDistance = distances[i][0];
		} else {
			if (minimumDistance > distances[i][0] && distances[i][0] != 0.00) {
				minimumDistance = distances[i][0];
			}
		}
	}
	return minimumDistance;
}

/**
 * Assign cluster labels to the non-border points based on distance and labels of border points.
 *
 * @param borderPoints        The 2D array of border points.
 * @param factor              The number of border points.
 * @param labels              The array of cluster labels for border points.
 * @param nonBorderPoints     The 2D array of non-border points.
 * @param otherFactor         The number of non-border points.
 * @param nonBorderLabels     The array to store the cluster labels for non-border points.
 */
void getLabelsNonBorderPoints(float **borderPoints, int factor, int *labels, float **nonBorderPoints, int otherFactor, int *nonBorderLabels) {
	float **distancesAndLabels, minDistance;
	int labelMin;

	distancesAndLabels = calloc(factor, sizeof(int *));
	if (distancesAndLabels == NULL) {
		printErrorAllocation();
	} else {
		for (int k = 0; k < factor; k++) {
			distancesAndLabels[k] = calloc(2, sizeof(float));
			if (distancesAndLabels[k] == NULL) {
				printErrorAllocation();
			}
		}
	}

	for (int i = 0; i < otherFactor; i++) {
		minDistance = 0.0;
		labelMin = 0;
		if (i != 0) {
			for (int h = 0; h < factor; h++) {
				distancesAndLabels[h][0] = 0;
				distancesAndLabels[h][1] = 0;
			}
		}
		for (int j = 0; j < factor; j++) {
			if (labels[j] != -1) {
				distancesAndLabels[j][0] = directionAngleModifiedDistanceFunction(nonBorderPoints[i][0], nonBorderPoints[i][1], borderPoints[j][0], borderPoints[j][1]);
				distancesAndLabels[j][1] = labels[j];
			}
		}
		minDistance = findMinimumDistance(distancesAndLabels, factor);
		for (int k = 0; k < factor; k++) {
			if ((minDistance == distancesAndLabels[k][0]) && (distancesAndLabels[k][1] != -1)) {
				labelMin = distancesAndLabels[k][1];
				break;
			}
		}
		nonBorderLabels[i] = labelMin;
	}

	free(distancesAndLabels);
}
