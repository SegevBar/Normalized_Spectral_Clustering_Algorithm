#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input file
* Action: Calculate and output the Weighted Adjacency Matrix
* Return: Prints Weighted Adjacency Matrix 
*/
void wam(char* filename) {
    int numOfVectors, numOfFeatures;
    double **dataPoints, **weightedAdjacencyMatrix;
    
    dataPoints = getDataPoints(&numOfVectors, &numOfFeatures, filename);
    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix(dataPoints,         
                              numOfVectors, numOfFeatures);
    freeMatrix(dataPoints);
    printSymmetricMatrix(weightedAdjacencyMatrix, numOfVectors);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input points matrix, and it's dimentions
* Action: Create Weighted Adjacency Matrix
* Return: Weighted Adjacency Matrix
*/
double **createWeightedAdjacencyMatrix(double **dataPoints, int numOfVectors,
                                       int numOfFeatures) {

    double **weightedAdjacencyMatrix;
    int i, j;
    weightedAdjacencyMatrix = createSymmetricMatrix(numOfVectors);
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < i; j++) {
            weightedAdjacencyMatrix[i][j] = exp(
                    -sqrt(euclideanNorm(dataPoints[i], dataPoints[j],
                                        numOfFeatures)) / 2);
        }
    }
    return weightedAdjacencyMatrix;
}

