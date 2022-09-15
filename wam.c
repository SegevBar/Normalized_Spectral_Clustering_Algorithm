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

    /* Form The Weighted Adjacency Matrix W from X */
    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);

    /* print the weighted adjacency */
    printMatrix(wam, N, N);
    freeMatrix(wam, N);
}

/*
* Funcion: double **getWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
*                                       int vectorDim)
* -----------------------------------------------------------------------------
* Params: Input vectors matrix, and it's dimentions N*(vector dimension)
* Action: Create Weighted Adjacency Matrix
* Return: Weighted Adjacency Matrix
*/
double **getWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
                                       int vectorDim) {
    double **wam, wij, norm;
    int i, j;
    
    wam = createSquareMatrix(N); /* allocate memory */
    /* calculate wij = exp(-(1/2)*(euclideanNorm(xi-xj))) */
    /* run on lower triangle of matrix and update it symetricly */
    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            /* calculate norm = sqrt(euclideanNorm(xi-xj)) */
            norm = sqrt(euclideanNorm(
                        vectorsMatrix[i], vectorsMatrix[j], vectorDim));
            /* calculate wij = exp(-(1/2)*norm) */
            wij = exp(-norm/2);
            /* update it symetricly */
            wam[i][j] = wij;
            wam[j][i] = wij;
        }
    }
    return wam;
}
