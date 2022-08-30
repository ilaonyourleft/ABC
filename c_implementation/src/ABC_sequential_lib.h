/*
 * ABC_sequential.h
 *
 *  Created on: 25 giu 2022
 *      Author: ilaria
 */

#define K 12
#define BETA 0.2
#define N 500

#define OPEN_FILE_ERROR -1
#define MEMORY_ALLOCATION_ERROR -2

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

// KNN methods
float euclideanDistance(int x1, int y1, int x2, int y2);
void sortArrayDistances(float **distancesPoints);
void getNeighbors(int **points, int x, int y, int **knn, float *meanPoint);

// ENCLOSING ANGLES
float getDirectionalAngle(int *center, float *meanPoint, int *neighbor);
int findSize(float *directionalAngles);
float getEnclosingAngle(float *directionalAngles);
float getBorderDegree(float *directionalAngles);
int isBorderPoint(float enclosingAngle);
void sortArrayBorderDegrees(float *borderDegrees, int sizeArray);
void getBorderPoints(float **borderPointsAll, int sizeArray, int **borderPoints, int factor);

// DBSCAN
int scalarProduct(int aX, int aY, int bX, int bY);
float moduleVector(int x, int y);
float angleBetweenVectors(int aX, int aY, int bX, int bY);
float directionAngleModifiedDistanceFunction(int aX, int aY, int bX, int bY);
int regionQuery(int **borderPoints, int **neighbors, int factor, int x, int y, int epsilon);
void growCluster(int **borderPoints, int factor, int *labels, int index, int x, int y, int **neighbors, int lenNeighbors, int clusterId, int epsilon, int minNumberPoints);
void getLabelsBorderPoints(int **borderPoints, int factor, int epsilon, int minNumberPoints, int *labels);
