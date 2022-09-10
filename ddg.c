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
* Action: Calculate and output the Diagonal Degree Matrix
* Return: Prints Diagonal Degree Matrix
*/
void ddg(double** vectorsMatrix, int N, int vectorDim) {
    double **weightedAdjacencyMatrix, **diagonalDegreeMatrix;

    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix
                                            (vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    diagonalDegreeMatrix = createDDGMatrixforDDG(weightedAdjacencyMatrix,   
                           N);
    printMatrix(diagonalDegreeMatrix, N, vectorDim);
    
    freeMatrix(weightedAdjacencyMatrix, N);
    freeMatrix(diagonalDegreeMatrix, N);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Points amount
* Action: Creates Diagonal Degree Matrix
* Return: Diagonal Degree Matrix
*/
double** createDDGMatrixforDDG(double **weightedAdjacencyMatrix, 
                               int N) {
    double *DiagonalDegreeArray;
    double **DiagonalDegreeMatrix;
    int i;

    DiagonalDegreeArray = calculateDiagonalDegreeMatrix
                                    (weightedAdjacencyMatrix, N);
    DiagonalDegreeMatrix = createRegularSquareMatrix(N);
    for (i = 0; i < N; i++) {
        DiagonalDegreeMatrix[i][i] = DiagonalDegreeArray[i];
    }
    free(DiagonalDegreeArray);
    return DiagonalDegreeMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Points amount
* Action: Calculates the diagonal of the diagonal degree matrix
* Return: Array of values in diagonal
*/
double *calculateDiagonalDegreeMatrix(double **weightedAdjacencyMatrix,
                                      int N) {

    double *diagonalDegreeArray, sum;
    int i, j;
    diagonalDegreeArray = (double *) calloc(N, sizeof(double));
    validateAction(diagonalDegreeArray != NULL);
    for (i = 0; i < N; i++) {
        sum = 0;
        for (j = 0; j < N; j++) {
            if (j <= i) {
                sum += weightedAdjacencyMatrix[i][j];
            } else {
                sum += weightedAdjacencyMatrix[j][i];
            }
        }
        diagonalDegreeArray[i] = sum;
    }
    return diagonalDegreeArray;
}
