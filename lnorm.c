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
* Action: Calculate and output the Normalized Graph Laplacian
* Return: Prints Normalized Graph Laplacian
*/
void lnorm(char* filename) {
    int numOfVectors, numOfFeatures;
    double **dataPoints, **weightedAdjacencyMatrix, *diagonalDegreeArray;
    double **lnorm;

    dataPoints = getDataPoints(&numOfVectors, &numOfFeatures, filename);
    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix(dataPoints,         
                              numOfVectors, numOfFeatures);
    freeMatrix(dataPoints);              
    diagonalDegreeArray = calculateDiagonalDegreeMatrix(
                                    weightedAdjacencyMatrix, numOfVectors);
    lnorm = createLnorm(diagonalDegreeArray, 
                                    weightedAdjacencyMatrix, numOfVectors);
    printSymmetricMatrix(lnorm, numOfVectors);

    freeMatrix(weightedAdjacencyMatrix);
    free(diagonalDegreeArray);
    freeMatrix(lnorm);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Array of values in diagonal, Weighted Adjacency Matrix, Diagonal 
*         Degree Matrix, Points amount
* Action: Creates Lnorm
* Return: Lnorm
*/
double **
createLnorm(double *diagonalDegreeArray, double **weightedAdjacencyMatrix, 
            int numOfVectors) {
    double **lnormMatrix;
    int i, j;

    /* Sets values x in diagonal to (x)^(-1/2) */
    for (i = 0; i < numOfVectors; i++) {
        diagonalDegreeArray[i] = 1 / sqrt(diagonalDegreeArray[i]);
    }

    lnormMatrix = createSymmetricMatrix(numOfVectors);
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j <= i; j++) {
            if (i != j) {
                lnormMatrix[i][j] = -diagonalDegreeArray[i] *
                                    weightedAdjacencyMatrix[i][j] *
                                    diagonalDegreeArray[j];
            } else {
                lnormMatrix[i][j] =
                        1 -
                        diagonalDegreeArray[i] * weightedAdjacencyMatrix[i][j] *
                        diagonalDegreeArray[j];
            }
        }
    }
    return lnormMatrix;
}
