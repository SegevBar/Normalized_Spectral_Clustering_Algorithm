#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: void lnorm(double** vectorsMatrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: Vectors matrix, Vector amount N, Vector dimension
* Action: Calculate and output the Normalized Graph Laplacian
* Return: Prints Normalized Graph Laplacian
*/
void lnorm(double** vectorsMatrix, int N, int vectorDim) {
    double **wam, *ddgDiagonal, **lnorm;

    /* Form The Weighted Adjacency Matrix W from X */
    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    
    /* Form the diagonal of The Diagonal Degree Matrix D from W */
    ddgDiagonal = getDdgDiagonal(wam, N);
    /* Form The Normalized Graph Laplacian */
    lnorm = getLnorm(ddgDiagonal, wam, N);
    freeMatrix(wam, N);
    free(ddgDiagonal);

    /* print The Normalized Graph Laplacian */
    printMatrix(lnorm, N, N);
    freeMatrix(lnorm, N);
}

/*
* Funcion: double** getLnorm(double *ddgDiagonal, double **wam, int N)
* -----------------------------------------------------------------------------
* Params: Array of values in diagonal, Weighted Adjacency Matrix, Vectors 
*         amount N
* Action: Creates Lnorm
* Return: Lnorm
*/
double** getLnorm(double *ddgDiagonal, double **wam, int N) {
    double **lnormMatrix;
    int i, j;

    /* Sets values x in diagonal to (x)^(-1/2) */
    for (i = 0; i < N; i++) {
        ddgDiagonal[i] = 1 / sqrt(ddgDiagonal[i]);
    }

    lnormMatrix = createSquareMatrix(N); /* allocate memory */
    /* calculate I-(D^-1/2)(W)(D^-1/2) */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            lnormMatrix[i][j] = -ddgDiagonal[i]*wam[i][j]*ddgDiagonal[j];
            if (i == j) {  /* add 1 to diagonal */
                lnormMatrix[i][j] = lnormMatrix[i][j] + 1;         
            }
        }
    }
    return lnormMatrix;
}
