/*
 * ABC_sequential.h
 *
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

struct point_label {
    float x;
    float y;
    int label;
};

struct double_float {
    float x;
    float y;
};

struct triple_float {
    float x;
    float y;
    float z;
};

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
 * @param distancesPoints The array of distances to be sorted. Struct of three float elements: x, y, z.
 */
void sortArrayDistances(struct triple_float *distancesPoints) {
	float tmp[3];

	for (int i = 0; i < N-1; i++) {
		for (int j = i + 1; j < N; j++) {
			if (distancesPoints[i].z > distancesPoints[j].z) {
				tmp[0] = distancesPoints[i].x;
				tmp[1] = distancesPoints[i].y;
				tmp[2] = distancesPoints[i].z;
				distancesPoints[i].x = distancesPoints[j].x;
				distancesPoints[j].x = tmp[0];
				distancesPoints[i].y = distancesPoints[j].y;
				distancesPoints[j].y = tmp[1];
				distancesPoints[i].z = distancesPoints[j].z;
				distancesPoints[j].z = tmp[2];
			}
		}
	}
}

/**
 * Get the nearest neighbors of a given point based on Euclidean distances.
 *
 * @param points     The array of points. Struct of two float elements: x, y.
 * @param x          The x-coordinate of the reference point.
 * @param y          The y-coordinate of the reference point.
 * @param knn        The array to store the nearest neighbors. Struct of two float elements: x, y.
 * @param meanPoint  The array to store the mean point of the nearest neighbors.
 */
// non permette di ritornare un array, ma si può ritornare il puntatore all'array specificandone il nome senza indice
void getNeighbors(struct double_float *points, float x, float y, struct double_float *knn, float *meanPoint) {
	// non ritorna l'indirizzo di una variabile locale all'esterno della funzione, quindi serve static nella definizione della variabile locale
	float tmp_x = 0.00, tmp_y = 0.00;
	int i;

	struct triple_float *distances = calloc(N, sizeof(struct triple_float));
	if (distances == NULL) {
		printErrorAllocation();
	}

	for (i = 0; i < N; i++) {
		float distance = euclideanDistance(x, y, points[i].x, points[i].y);
		distances[i].x = points[i].x;
		distances[i].y = points[i].y;
		distances[i].z = distance;
	}

	sortArrayDistances(distances);

	for (i = 0; i < K; i++) {
		knn[i].x = distances[i+1].x;
		knn[i].y = distances[i+1].y;
	}

	free(distances);

	for (i = 0; i < K; i++) {
		tmp_x += knn[i].x;
		tmp_y += knn[i].y;
	}

	meanPoint[0] = tmp_x / K;
	meanPoint[1] = tmp_y / K;
}

// ENCLOSING ANGLES

/**
 * Calculate the directional angle between the line segments formed by the center, mean point, and neighbor.
 *
 * @param center     The coordinates of the center point. Struct of two float elements: x, y.
 * @param meanPoint  The coordinates of the mean point.
 * @param neighbor   The coordinates of the neighbor point. Struct of two float elements: x, y.
 * @return The directional angle in degrees.
 */
float getDirectionalAngle(struct double_float center, float *meanPoint, struct double_float neighbor) {
	float uX = 0.00, uY = 0.00, vX = 0.00, vY = 0.00, directionalAngle = 0.00, directionalAngleDegree = 0.00;

	// mean point - center
	uX = meanPoint[0] - center.x;
	uY = meanPoint[1] - center.y;

	// neighbor - center
	vX = neighbor.x - center.x;
	vY = neighbor.y - center.y;

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
 * @param borderDegrees The 2D array of border degrees. Struct of three float elements: x, y, z.
 */
void sortArrayBorderDegrees(struct triple_float *borderDegrees) {
	float tmp[3];

	for (int i = 0; i < N-1; i++) {
		for (int j = i + 1; j < N; j++) {
			if (borderDegrees[i].z < borderDegrees[j].z) {
				tmp[0] = borderDegrees[i].x;
				tmp[1] = borderDegrees[i].y;
				tmp[2] = borderDegrees[i].z;
				borderDegrees[i].x = borderDegrees[j].x;
				borderDegrees[j].x = tmp[0];
				borderDegrees[i].y = borderDegrees[j].y;
				borderDegrees[j].y = tmp[1];
				borderDegrees[i].z = borderDegrees[j].z;
				borderDegrees[j].z = tmp[2];
			}
		}
	}
}

/**
 * Get the border points from the sorted array of border points.
 *
 * @param borderPointsAll The 2D array of border points. Struct of three float elements: x, y, z.
 * @param sizeArray       The size of the array.
 * @param borderPointsAndLabels    The 2D array to store the border points. Struct of three elements: x (float), y (float), label (int).
 */
void getBorderPoints(struct triple_float *borderPointsAll, int sizeArray, struct point_label *borderPointsAndLabels) {
	sortArrayBorderDegrees(borderPointsAll);
	for (int i = 0; i < sizeArray; i++) {
		borderPointsAndLabels[i].x = borderPointsAll[i].x;
		borderPointsAndLabels[i].y = borderPointsAll[i].y;
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
 * @param borderPointsAndLabels The 2D array of border points. Struct of three elements: x (float), y (float), label (int).
 * @param neighbors    			The 2D array to store the neighboring points. Struct of three float elements: x, y, z.
 * @param factor       			The number of border points.
 * @param x            			The x component of the point to query.
 * @param y            			The y component of the point to query.
 * @param epsilon      			The distance threshold.
 * @return The number of neighboring points found.
 */
int regionQuery(struct point_label *borderPointsAndLabels, struct triple_float *neighbors, int factor, float x, float y, int epsilon) {
	int counter = 0;
	for (int i = 0; i < factor; i++) {
		float disComputed = directionAngleModifiedDistanceFunction(x, y, borderPointsAndLabels[i].x, borderPointsAndLabels[i].y);
		if (disComputed < epsilon) {
			neighbors[counter].x = i;
			neighbors[counter].y = borderPointsAndLabels[i].x;
			neighbors[counter].z = borderPointsAndLabels[i].y;
			++counter;
		}
	}
	return counter;
}

/**
 * Check if a point is already a neighbor in the list of neighbors.
 *
 * @param neighbors    The 2D array of neighbors. Struct of three float elements: x, y, z.
 * @param lenNeighbors The length of the neighbors array.
 * @param index        The index of the point to check. Struct of three float elements: x, y, z.
 * @return 1 if the point is already a neighbor, 0 otherwise.
 */
int checkIfAlreadyNeighbor(struct triple_float *neighbors, int lenNeighbors, struct triple_float index) {
	int flag = 0;
	for (int i = 0; i < lenNeighbors; i++) {
		if (neighbors[i].x >= 0) {
			if (neighbors[i].y == index.y && neighbors[i].z == index.z) {
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
 * @param factor            	The number of border points.
 * @param borderPointsAndLabels	The array to store the cluster labels. Struct of three elements: x (float), y (float), label (int).
 * @param index             	The index of the initial point.
 * @param x                 	The x component of the initial point.
 * @param y                 	The y component of the initial point.
 * @param neighbors         	The 2D array of neighboring points. Struct of three float elements: x, y, z.
 * @param lenNeighbors      	The length of the neighbors array.
 * @param clusterId         	The ID of the cluster.
 * @param epsilon           	The distance threshold for neighboring points.
 * @param minNumberPoints   	The minimum number of points required to form a cluster.
 */
void growCluster(int factor, struct point_label *borderPointsAndLabels, int index, float x, float y, struct triple_float *neighbors, int lenNeighbors, int clusterId, int epsilon, int minNumberPoints) {
	borderPointsAndLabels[index].label = clusterId;
	int counter = 0, i, j, neighborsIncrement;

	struct triple_float *ptrNextNeighbors = calloc(factor, sizeof(struct triple_float));
	if (ptrNextNeighbors == NULL) {
		printErrorAllocation();
	}

	while (counter < lenNeighbors) {
		int lenNextNeighbors = 0;
		if (counter != 0) {
			for (j = 0; j < factor; j++) {
				ptrNextNeighbors[j].x = 0;
				ptrNextNeighbors[j].y = 0;
				ptrNextNeighbors[j].z = 0;
			}
		}
		int next_index = neighbors[counter].x;
		if (borderPointsAndLabels[next_index].label == -1) {
			borderPointsAndLabels[next_index].label = clusterId;
		} else if (borderPointsAndLabels[next_index].label == 0) {
			borderPointsAndLabels[next_index].label = clusterId;
			lenNextNeighbors = regionQuery(borderPointsAndLabels, ptrNextNeighbors, factor, borderPointsAndLabels[next_index].x, borderPointsAndLabels[next_index].y, epsilon);
			if (lenNextNeighbors >= minNumberPoints) {
				neighborsIncrement = 0;
				for (i = 0; i < lenNextNeighbors; i++) {
					if (checkIfAlreadyNeighbor(neighbors, lenNeighbors, ptrNextNeighbors[i]) == 0) {
						neighbors[lenNeighbors + neighborsIncrement].x = ptrNextNeighbors[i].x;
						neighbors[lenNeighbors + neighborsIncrement].y = ptrNextNeighbors[i].y;
						neighbors[lenNeighbors + neighborsIncrement].z = ptrNextNeighbors[i].z;
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
 * @param factor            	The number of border points.
 * @param epsilon           	The distance threshold for neighboring points.
 * @param minNumberPoints   	The minimum number of points required to form a cluster.
 * @param borderPointsAndLabels	The array to store the cluster labels. Struct of three elements: x (float), y (float), label (int).
 */
void getLabelsBorderPoints(int factor, int epsilon, int minNumberPoints, struct point_label *borderPointsAndLabels) {
	int clusterId = 0, i, j;

	struct triple_float *ptrNeighbors = calloc(factor, sizeof(struct triple_float));
	if (ptrNeighbors == NULL) {
		printErrorAllocation();
	}

	for (i = 0; i < factor; i++) {
		int lenNeighbors = 0;
		if (i != 0) {
			for (j = 0; j < factor; j++) {
				ptrNeighbors[j].x = 0;
				ptrNeighbors[j].y = 0;
				ptrNeighbors[j].z = 0;
			}
		}
		if (borderPointsAndLabels[i].label == 0) {
			lenNeighbors = regionQuery(borderPointsAndLabels, ptrNeighbors, factor, borderPointsAndLabels[i].x, borderPointsAndLabels[i].y, epsilon);
			if (lenNeighbors < minNumberPoints) {
				borderPointsAndLabels[i].label = -1;
			} else {
				++clusterId;
				growCluster(factor, borderPointsAndLabels, i, borderPointsAndLabels[i].x, borderPointsAndLabels[i].y, ptrNeighbors, lenNeighbors, clusterId, epsilon, minNumberPoints);
			}
		}
	}

	free(ptrNeighbors);
}

// CLUSTER

/**
 * Check if a point is a border point.
 *
 * @param borderPointsAndLabels	The 2D array of border points. Struct of three elements: x (float), y (float), label (int).
 * @param factor          		The number of border points.
 * @param x               		The x component of the point to check.
 * @param y               		The y component of the point to check.
 * @return 1 if the point is a border point, 0 otherwise.
 */
int checkIfBorderPoint(struct point_label *borderPointsAndLabels, int factor, float x, float y) {
	for (int i = 0; i < factor; i++) {
		if (borderPointsAndLabels[i].x == x && borderPointsAndLabels[i].y == y) {
			return 1;
		}
	}
	return 0;
}

/**
 * Get the non-border points from the given points array.
 *
 * @param points            		The 2D array of points. Struct of two float elements: x, y.
 * @param borderPointsAndLabels		The 2D array of border points. Struct of three elements: x (float), y (float), label (int).
 * @param factor            		The number of border points.
 * @param internalPointsAndLabels	The 2D array to store the non-border points. Struct of three elements: x (float), y (float), label (int).
 */
void getNonBorderPoints(struct double_float *points, struct point_label *borderPointsAndLabels, int factor, struct point_label *internalPointsAndLabels) {
	int counter = 0;
	for (int i = 0; i < N; i++) {
		if (checkIfBorderPoint(borderPointsAndLabels, factor, points[i].x, points[i].y) == 0) {
			internalPointsAndLabels[counter].x = points[i].x;
			internalPointsAndLabels[counter].y = points[i].y;
			++counter;
		}
	}
}

/**
 * Find the minimum distance from the given distances array.
 *
 * @param distances   The 2D array of distances. Struct of two float elements: x, y.
 * @param factor      The number of distances.
 * @return The minimum distance found.
 */
float findMinimumDistance(struct double_float *distances, int factor) {
	float minimumDistance = 0.0;
	for (int i = 0; i < factor; i++) {
		if (i == 0) {
			minimumDistance = distances[i].x;
		} else {
			if (minimumDistance > distances[i].x && distances[i].x != 0.00) {
				minimumDistance = distances[i].x;
			}
		}
	}
	return minimumDistance;
}

/**
 * Assign cluster labels to the non-border points based on distance and labels of border points.
 *
 * @param factor              		The number of border points.
 * @param borderPointsAndLabels		The array of cluster labels for border points. Struct of three elements: x (float), y (float), label (int).
 * @param internalPointsAndLabels	The 2D array of non-border points. Struct of three elements: x (float), y (float), label (int).
 * @param otherFactor         		The number of non-border points.
 */
void getLabelsNonBorderPoints(int factor, struct point_label *borderPointsAndLabels, struct point_label *internalPointsAndLabels, int otherFactor) {
	float minDistance;
	int labelMin;

	struct double_float *distancesAndLabels = calloc(factor, sizeof(struct double_float));
	if (distancesAndLabels == NULL) {
		printErrorAllocation();
	}

	for (int i = 0; i < otherFactor; i++) {
		minDistance = 0.0;
		labelMin = 0;
		if (i != 0) {
			for (int h = 0; h < factor; h++) {
				distancesAndLabels[h].x = 0;
				distancesAndLabels[h].y = 0;
			}
		}
		for (int j = 0; j < factor; j++) {
			if (borderPointsAndLabels[j].label != -1) {
				distancesAndLabels[j].x = directionAngleModifiedDistanceFunction(internalPointsAndLabels[i].x, internalPointsAndLabels[i].y, borderPointsAndLabels[j].x, borderPointsAndLabels[j].y);
				distancesAndLabels[j].y = borderPointsAndLabels[j].label;
			}
		}
		minDistance = findMinimumDistance(distancesAndLabels, factor);
		for (int k = 0; k < factor; k++) {
			if ((minDistance == distancesAndLabels[k].x) && (distancesAndLabels[k].y != -1)) {
				labelMin = distancesAndLabels[k].y;
				break;
			}
		}
		internalPointsAndLabels[i].label = labelMin;
	}

	free(distancesAndLabels);
}
