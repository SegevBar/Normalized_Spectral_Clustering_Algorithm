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
* Action: Calculate and output the Diagonal Degree Matrix
* Return: Prints Diagonal Degree Matrix
*/
void ddg(char* filename) {
    int numOfVectors, numOfFeatures;
    double **dataPoints, **weightedAdjacencyMatrix, **diagonalDegreeMatrix;

    dataPoints = getDataPoints(&numOfVectors, &numOfFeatures, filename);
    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix(dataPoints,         
                              numOfVectors, numOfFeatures);
    freeMatrix(dataPoints);              
    diagonalDegreeMatrix = createDDGMatrixforDDG(weightedAdjacencyMatrix,   
                           numOfVectors);
    printMatrix(diagonalDegreeMatrix, numOfVectors, numOfVectors);
    
    freeMatrix(weightedAdjacencyMatrix);
    freeMatrix(diagonalDegreeMatrix);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Points amount
* Action: Creates Diagonal Degree Matrix
* Return: Diagonal Degree Matrix
*/
double **
createDDGMatrixforDDG(double **weightedAdjacencyMatrix, int numOfVectors) {
    double *DiagonalDegreeArray;
    double **DiagonalDegreeMatrix;
    int i;

    DiagonalDegreeArray = calculateDiagonalDegreeMatrix
                                    (weightedAdjacencyMatrix, numOfVectors);
    DiagonalDegreeMatrix = createRegularSquareMatrix(numOfVectors);
    for (i = 0; i < numOfVectors; i++) {
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
                                      int numOfVectors) {

    double *diagonalDegreeArray, sum;
    int i, j;
    diagonalDegreeArray = (double *) calloc(numOfVectors, sizeof(double));
    ourAssert(diagonalDegreeArray != NULL);
    for (i = 0; i < numOfVectors; i++) {
        sum = 0;
        for (j = 0; j < numOfVectors; j++) {
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
