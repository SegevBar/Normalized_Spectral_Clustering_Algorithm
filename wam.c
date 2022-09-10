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
* Action: Calculate and output the Weighted Adjacency Matrix
* Return: Prints Weighted Adjacency Matrix 
*/
void wam(double** vectorsMatrix, int N, int vectorDim) {
    double **weightedAdjacencyMatrix;

    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix
                                    (vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);
    printMatrix(weightedAdjacencyMatrix, N, N);
    freeMatrix(weightedAdjacencyMatrix, N);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input points matrix, and it's dimentions
* Action: Create Weighted Adjacency Matrix
* Return: Weighted Adjacency Matrix
*/
double **createWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
                                       int vectorDim) {
    double **weightedAdjacencyMatrix;
    int i, j;
    double wij;
    
    weightedAdjacencyMatrix = createSquareMatrix(N);
    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            wij = exp(-sqrt
            (euclideanNorm(vectorsMatrix[i], vectorsMatrix[j],vectorDim)) / 2);
            weightedAdjacencyMatrix[i][j] = wij;
            weightedAdjacencyMatrix[j][i] = wij;
        }
    }
    return weightedAdjacencyMatrix;
}

