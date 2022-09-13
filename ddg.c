#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: void ddg(double** vectorsMatrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: Vectors matrix, Vector count, Vector dimension
* Action: Calculate and output the Diagonal Degree Matrix
* Return: Prints Diagonal Degree Matrix
*/
void ddg(double** vectorsMatrix, int N, int vectorDim) {
    double **wam, **ddg;

    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    ddg = getDiagonalDegreeMatrix(wam, N);
    freeMatrix(wam, N);
    
    printMatrix(ddg, N, vectorDim);
    freeMatrix(ddg, N);
}

/*
* Funcion: double** getDiagonalDegreeMatrix(double **wam, int N)
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Vectors amount
* Action: Creates Diagonal Degree Matrix
* Return: Diagonal Degree Matrix
*/
double** getDiagonalDegreeMatrix(double **wam, int N) {
    double *diagonal, **ddg;
    int i;

    diagonal = getDdgDiagonal(wam, N);
    
    ddg = createSquareMatrix(N); /* allocate memory */
    for (i = 0; i < N; i++) {
        ddg[i][i] = diagonal[i];
    }
    free(diagonal);
    return ddg;
}

/*
* Funcion: double *getDdgDiagonal(double **wam, int N)
* -----------------------------------------------------------------------------
* Params: Weighted Adjacency Matrix, Vectors amount
* Action: Calculates the diagonal of the diagonal degree matrix
* Return: Array of values in diagonal
*/
double* getDdgDiagonal(double **wam, int N) {
    double *diagonal, sum;
    int i, j;
    
    diagonal = (double*) calloc(N, sizeof(double));
    validateAction(diagonal != NULL);

    /* sum wij at each row of wam */
    for (i = 0; i < N; i++) {
        sum = 0;
        for (j = 0; j < N; j++) {
            sum += wam[i][j];
        }
        diagonal[i] = sum;
    }
    return diagonal;
}
