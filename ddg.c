#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: void ddg(double** vectorsMatrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: Vectors matrix, Vector amount N, Vector dimension
* Action: Calculate and output the Diagonal Degree Matrix
* Return: Prints Diagonal Degree Matrix
*/
void ddg(double** vectorsMatrix, int N, int vectorDim) {
    double **wam, **ddg;

    /* Form The Weighted Adjacency Matrix W from X */
    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);

    /* Form The Diagonal Degree Matrix D from W */         
    ddg = getDiagonalDegreeMatrix(wam, N);
    freeMatrix(wam, N);
    
    /* print The Diagonal Degree Matrix */
    printMatrix(ddg, N, N);
    freeMatrix(ddg, N);
}

/*
* Funcion: double** getDiagonalDegreeMatrix(double **wam, int N)
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Vectors amount N
* Action: Creates Diagonal Degree Matrix
* Return: Diagonal Degree Matrix
*/
double** getDiagonalDegreeMatrix(double **wam, int N) {
    double *diagonal, **ddg;
    int i;

    /* calculate the diagonal of DDG matrix */
    diagonal = getDdgDiagonal(wam, N);
    
    ddg = createSquareMatrix(N); /* allocate memory */
    /* update values of diagonal in matrix */
    for (i = 0; i < N; i++) {
        ddg[i][i] = diagonal[i];
    }
    free(diagonal);
    return ddg;
}

/*
* Funcion: double* getDdgDiagonal(double **wam, int N)
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Vectors amount N
* Action: Calculates the diagonal of the diagonal degree matrix
* Return: Array of values in diagonal
*/
double* getDdgDiagonal(double **wam, int N) {
    double *diagonal, sum;
    int i, j;
    
    diagonal = (double*) calloc(N, sizeof(double));
    validateAction(diagonal != NULL);

    /* sum wij at each row of wam */
    /* dij = sum(wiz) for z = 1,..,n if i = j */
    for (i = 0; i < N; i++) {
        sum = 0;
        for (j = 0; j < N; j++) {
            sum += wam[i][j];
        }
        diagonal[i] = sum;
    }
    return diagonal;
}
