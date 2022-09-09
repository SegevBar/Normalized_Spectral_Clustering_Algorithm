#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Clusters, Points matrix and its' dimentions
* Action: Executes the kmeans algorithm
* Return: Print output centroids
*/
void kmeansmain(CLUSTER *clusters, double **dataPoints, int numOfFeatures,
                int numOfVectors) {
    int i;
    int j;
    int k;
    double minval;
    int mincluster;
    double arg;
    int stop;

    for (i = 0; i < MAX_ITER_KMEANS; i++) {
        for (j = 0; j < numOfVectors; j++) {
            minval = euclideanNorm(dataPoints[j], clusters[0].centroid,
                                   numOfFeatures);
            mincluster = 0;
            for (k = 1; k < numOfFeatures; k++) {
                arg = euclideanNorm(dataPoints[j], clusters[k].centroid,
                                    numOfFeatures);
                if (arg < minval) {
                    minval = arg;
                    mincluster = k;
                }
            }
            updateClosest(&clusters[mincluster], dataPoints[j], numOfFeatures);
        }
        stop = 1;
        for (k = 0; k < numOfFeatures; k++) {
            updateCentroid(&clusters[k], numOfFeatures);
            if (!clusters[k].equalTolLast) {
                stop = 0;
            }
        }
        if (stop) {
            break;
        }
    }
    printCentroids(clusters, numOfFeatures);
    freeMatrix(dataPoints);
    freeClusters(clusters, numOfFeatures);
}



/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Points matrix, Point size
* Action: Create k clusters 
* Return: Array of clusters
*/
CLUSTER *initializeClusters(double **dataPoints, int K) {

    CLUSTER *clusters;
    int i;

    clusters = (CLUSTER *) calloc(K, sizeof(CLUSTER));
    ourAssert(clusters != NULL);

    for (i = 0; i < K; i++) {
        initCluster(&clusters[i], dataPoints[i], K);
    }

    return clusters;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Cluster, a Point, Point size
* Action: Init cluster
* Return: None
*/
void initCluster(CLUSTER *curCluster, double *dataPoint, int numOfFeatures) {

    curCluster->centroid = (double *) calloc(numOfFeatures, sizeof(double));
    ourAssert(curCluster->centroid != NULL);
    memcpy(curCluster->centroid, dataPoint, numOfFeatures * sizeof(double));
    curCluster->centroid_closest =
            (double *) calloc(numOfFeatures, sizeof(double));
    ourAssert(curCluster->centroid_closest != NULL);

    curCluster->size = 0;
    curCluster->equalTolLast = 0; /* 0 is False */
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Cluster, a Point, Point size
* Action: Updates centroid_closest of cluster to next centroid
* Return: None
*/
void
updateClosest(CLUSTER *curCluster, const double *datapoint, int numOfFeatures) {

    int i;
    for (i = 0; i < numOfFeatures; i++) {
        curCluster->centroid_closest[i] =
                curCluster->centroid_closest[i] + datapoint[i];
    }
    curCluster->size++;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Cluster, Point size
* Action: Updates cluster centroid and zeroing the centroid_closest
* Return: None
*/
void updateCentroid(CLUSTER *curCluster, int numOfFeatures) {

    int i;
    curCluster->equalTolLast = 1;
    for (i = 0; i < numOfFeatures; i++) {
        curCluster->centroid_closest[i] =
                curCluster->centroid_closest[i] / curCluster->size;
        if (curCluster->centroid_closest[i] != curCluster->centroid[i]) {
            curCluster->equalTolLast = 0;
        }
    }
    free(curCluster->centroid);
    curCluster->centroid = curCluster->centroid_closest;
    curCluster->centroid_closest = (double *) calloc(numOfFeatures,
                                                     sizeof(double));
    ourAssert(curCluster->centroid_closest != NULL);
    curCluster->size = 0;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 2 points, Point size
* Action: Computes euclidean norm between 2 points
* Return: square norm
*/
double euclideanNorm(const double *datapoint1, const double *datapoint2,
                     int numOfFeatures) {

    int i;
    double sum = 0;
    double temp;

    for (i = 0; i < numOfFeatures; i++) {
        temp = datapoint1[i] - datapoint2[i];
        temp *= temp;
        sum += temp;
    }

    return sum;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Clusters, point size
* Action: Prints final centroids from the K-means algorithm
* Return: None
*/
void printCentroids(CLUSTER *clusters, int numOfFeatures) {

    int i;
    int j;

    for (i = 0; i < numOfFeatures; i++) {
        for (j = 0; j < numOfFeatures; j++) {
            if (j != numOfFeatures - 1) {
                printf("%.4f,", round(clusters[i].centroid[j]));
            } else {
                if (i != numOfFeatures - 1) {
                    printf("%.4f\n", round(clusters[i].centroid[j]));
                } else {
                    printf("%.4f", round(clusters[i].centroid[j]));
                }
            }
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Clusters, k
* Action: Frees allocated memory of CLUSTER
* Return: None
*/
void freeClusters(CLUSTER *clusters, int K) {

    int i;
    for (i = 0; i < K; i++) {
        free(clusters[i].centroid);
        free(clusters[i].centroid_closest);
    }
    free(clusters);
}
