#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Points matrix, Vector count, Vector dimension
* Action: Calculate and output the Normalized Graph Laplacian
* Return: Prints Normalized Graph Laplacian
*/
void lnorm(double** vectorsMatrix, int N, int vectorDim) {
    double **weightedAdjacencyMatrix, *diagonalDegreeArray;
    double **lnorm;

    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix
                                            (vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    diagonalDegreeArray = calculateDiagonalDegreeMatrix
                                            (weightedAdjacencyMatrix, N);
    lnorm = createLnorm(diagonalDegreeArray, 
                                    weightedAdjacencyMatrix, N);
    printSymmetricMatrix(lnorm, N);

    freeMatrix(weightedAdjacencyMatrix, N);
    free(diagonalDegreeArray);
    freeMatrix(lnorm, N);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Array of values in diagonal, Weighted Adjacency Matrix, Diagonal 
*         Degree Matrix, Points amount
* Action: Creates Lnorm
* Return: Lnorm
*/
double** createLnorm(double *diagonalDegreeArray, 
                     double **weightedAdjacencyMatrix, int N) {
    double **lnormMatrix;
    int i, j;

    /* Sets values x in diagonal to (x)^(-1/2) */
    for (i = 0; i < N; i++) {
        diagonalDegreeArray[i] = 1 / sqrt(diagonalDegreeArray[i]);
    }

    lnormMatrix = createSymmetricMatrix(N);
    for (i = 0; i < N; i++) {
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
