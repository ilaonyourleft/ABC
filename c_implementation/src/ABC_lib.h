/*
 * ABC_lib.h
 *
 *      Author: ilaria
 */

#define K 12  // 12
#define BETA 0.2  // 0.2
#define N 10000  // 500

#define OPEN_FILE_ERROR -1
#define MEMORY_ALLOCATION_ERROR -2

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
void printErrorAllocation();
void printExitError();

// KNN methods
float euclideanDistance(float x1, float y1, float x2, float y2);
void swap(struct triple_float *a, struct triple_float *b);
int partitionDistances(struct triple_float *distancesPoints, int low, int high);
int partitionDegrees(struct triple_float *borderDegrees, int low, int high);
void quicksort(struct triple_float *arrays, int low, int high, int flag);
void sortArrayDistances(struct triple_float *distancesPoints);
void getNeighbors(struct double_float *points, float x, float y, struct double_float *knn, float *meanPoint);

// ENCLOSING ANGLES
float getDirectionalAngle(struct double_float center, float *meanPoint, struct double_float neighbor);
int findSize(float *directionalAngles);
float getEnclosingAngle(float *directionalAngles);
float getBorderDegree(float *directionalAngles);
int isBorderPoint(float enclosingAngle);
void sortArrayBorderDegrees(struct triple_float *borderDegrees);
void getBorderPoints(struct triple_float *borderPointsAll, int sizeArray, struct point_label *borderPointsAndLabels);

// DBSCAN
float moduleVector(float x, float y);
float directionAngleModifiedDistanceFunction(float aX, float aY, float bX, float bY);
int regionQuery(struct point_label *borderPointsAndLabels, struct triple_float *neighbors, int factor, float x, float y, int epsilon);
int checkIfAlreadyNeighbor(struct triple_float *neighbors, int lenNeighbors, struct triple_float index);
void growCluster(int factor, struct point_label *borderPointsAndLabels, int index, float x, float y, struct triple_float *neighbors, int lenNeighbors, int clusterId, int epsilon, int minNumberPoints);
void getLabelsBorderPoints(int factor, float epsilon, int minNumberPoints, struct point_label *borderPointsAndLabels);

// CLUSTER
int checkIfBorderPoint(struct point_label *borderPointsAndLabels, int factor, float x, float y);
void getNonBorderPoints(struct double_float *points, struct point_label *borderPointsAndLabels, int factor, struct point_label *internalPointsAndLabels);
float findMinimumDistance(struct double_float *distances, int factor);
void getLabelsNonBorderPoints(int factor, struct point_label *borderPointsAndLabels, struct point_label *internalPointsAndLabels, int otherFactor);
