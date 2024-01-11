/*
 * ABC_sequential.h
 *
 *  Created on: 25 giu 2022
 *      Author: ilaria
 */

#define K 12  // 12
#define BETA 0.3  // 0.2
#define N 500  // 500

#define OPEN_FILE_ERROR -1
#define MEMORY_ALLOCATION_ERROR -2

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

// MISCELLANEOUS
void printErrorAllocation();

// KNN methods
float euclideanDistance(float x1, float y1, float x2, float y2);
void sortArrayDistances(float **distancesPoints);
void getNeighbors(float **points, float x, float y, float **knn, float *meanPoint);

// ENCLOSING ANGLES
float getDirectionalAngle(float *center, float *meanPoint, float *neighbor);
int findSize(float *directionalAngles);
float getEnclosingAngle(float *directionalAngles);
float getBorderDegree(float *directionalAngles);
int isBorderPoint(float enclosingAngle);
void sortArrayBorderDegrees(float *borderDegrees);
void getBorderPoints(float **borderPointsAll, int sizeArray, float **borderPoints);

// DBSCAN
float moduleVector(float x, float y);
float directionAngleModifiedDistanceFunction(float aX, float aY, float bX, float bY);
int regionQuery(float **borderPoints, float **neighbors, int factor, float x, float y, int epsilon);
void growCluster(float **borderPoints, int factor, int *labels, int index, float x, float y, float **neighbors, int lenNeighbors, int clusterId, int epsilon, int minNumberPoints);
void getLabelsBorderPoints(float **borderPoints, int factor, int epsilon, int minNumberPoints, int *labels);

// CLUSTER
int checkIfBorderPoint(float **borderPoints, int factor, float x, float y);
void getNonBorderPoints(float **points, float **borderPoints, int factor, float **nonBorderPoints);
float findMinimumDistance(float **distances, int factor);
void getLabelsNonBorderPoints(float **borderPoints, int factor, int *labels, float **nonBorderPoints, int otherFactor, int *nonBorderLabels);
