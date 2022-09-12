#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: void wam(double** vectorsMatrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: Vectors matrix, Vector count, Vector dimension
* Action: Calculate and output the Weighted Adjacency Matrix
* Return: Prints Weighted Adjacency Matrix 
*/
void wam(double** vectorsMatrix, int N, int vectorDim) {
    double **wam;

    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);
    printMatrix(wam, N, N);
    freeMatrix(wam, N);
}

/*
* Funcion: double **getWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
*                                       int vectorDim)
* -----------------------------------------------------------------------------
* Params: Input vectors matrix, and it's dimentions
* Action: Create Weighted Adjacency Matrix
* Return: Weighted Adjacency Matrix
*/
double **getWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
                                       int vectorDim) {
    double **wam, wij, norm;
    int i, j;
    
    wam = createSquareMatrix(N); /* allocate memory */
    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            norm = euclideanNorm(vectorsMatrix[i], vectorsMatrix[j], vectorDim);
            wij = exp(-norm/2);
            wam[i][j] = wij;
            wam[j][i] = wij;
        }
    }
    return wam;
}
