#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: double** getNormalizedKEigenvectorsMatrix(int *kp, 
*          double** vectorsMatrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: k pointer, Vectors matrix, Vector count, Vector dimension
* Action: performs The Normalized Spectral Clustering Algorithm *WITHOUT*
*         kmeans (first step executes at python program)
* Return: Normalized K Eigenvectors Matrix (T)
*/
double** getNormalizedKEigenvectorsMatrix(int *kp, double** vectorsMatrix, 
                                          int N, int vectorDim) {
    double **wam, *ddgDiagonal, **lnorm, **eigenVecMatrix, **T;
    EIGEN* eigens;

    int k = *kp;

    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    
    ddgDiagonal = getDdgDiagonal(wam, N);
    lnorm = getLnorm(ddgDiagonal, wam, N);

    freeMatrix(wam, N);
    free(ddgDiagonal);

    eigenVecMatrix = jacobiAlgorithm(lnorm, N);
    eigens = createEigensArr(lnorm, N); 
    freeMatrix(lnorm, N);
    
    descendingSort(eigens, N); /*sort eigans from largest to smallest*/
    *kp = (k == 0) ? eigengapHeuristic(eigens, N) : k;

    T = createT(eigens, eigenVecMatrix, k, N);
    freeMatrix(eigenVecMatrix, N);

    return T;
}

/*
* Funcion: EIGEN *createEigensArr(double **eigenVectors, double **eiganVals, 
*          int n)
* -----------------------------------------------------------------------------
* Params: Eigenvectors Matrix, Eigenvalues Matrix, Matrix size (1D)
* Action: Create array of EIGENS representing eigenvalue and its' eigenvector
* Return: Array of EIGENS representing eigenvalue and its' eigenvector
*/
EIGEN *createEigensArr(double **eiganVals, int n) {
    EIGEN *eigenArray;
    int i;

    eigenArray = (EIGEN*) calloc(n, sizeof(EIGEN));
    validateAction(eigenArray != NULL);

    for (i = 0; i < n; i++) {
        eigenArray[i].eigenValue = eiganVals[i][i];
        eigenArray[i].eiganIndex = i;
    }
    return eigenArray;
}

/*
* Funcion: int eigengapHeuristic(EIGEN *eigenArray, int n)
* -----------------------------------------------------------------------------
* Params: EIGENS array and it's size
* Action: Execute the eigengap heuristic
* Return: k - number of clusters
*/
int eigengapHeuristic(EIGEN *eigenArray, int n) {
    double max, curMax;
    int maxIndex, i;
    
    maxIndex = 0;
    max = fabs(eigenArray[1].eigenValue - eigenArray[0].eigenValue);
    for (i = 1; i < n / 2; i++) {
        curMax = fabs(eigenArray[i+1].eigenValue - eigenArray[i].eigenValue);
        if (curMax > max) {
            maxIndex = i;
            max = curMax;
        }
    }
    return maxIndex + 1;
}

/*
* Funcion: double** createT(EIGEN* eigens, double** eigenVecMatrix, int k, 
*          int N)
* -----------------------------------------------------------------------------
* Params: Descending sorted eigens array, eigan vectors, k (col), N (rows)
* Action: Create Normalized K Eigenvectors Matrix (T) from k largest eigans
* Return: Matrix With Eigenvectors As Columns
*/
double** createT(EIGEN* eigens, double** eigenVecMatrix, int k, int N) {
    double **T;
    int i, j;
    
    T = createMatrix(N, k);  /* allocate memory */
    /* copy eiganvectors of k largest eigan values */
    for (j = 0; j < k; j++) {        
        for (i = 0; i < N; i++) {
            T[i][j] = eigenVecMatrix[i][eigens[j].eiganIndex];
        }
    }
    free(eigens);
    normalizeMatrixByRows(T, N, k);

    return T;
}

/*
* Funcion: void normalizeMatrixByRows(double **matrix, int row, int col)
* -----------------------------------------------------------------------------
* Params: a Matrix and its' dimensions
* Action: Normalizes values in matrix
* Return: None
*/
void normalizeMatrixByRows(double **matrix, int row, int col) {
    double *denominatorsArr;
    int i, j;

    denominatorsArr = getNormalizeDenominators(matrix, row, col);
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            if (denominatorsArr[i] != 0) {
                matrix[i][j] = matrix[i][j] / denominatorsArr[i];
            }
        }
    }
    free(denominatorsArr);
}

/*
* Funcion: double* getNormalizeDenominators(double **matrix, int row, int col)
* -----------------------------------------------------------------------------
* Params: a Matrix and its' dimensions
* Action: Calaulate sqrt of sum of (u_ij)^2 for each row for normalize algo
* Return: Array with calculated sum of each row
*/
double* getNormalizeDenominators(double **matrix, int row, int col) {
    double *denominatorsArr, sum;
    int j, i;
    
    denominatorsArr = (double *) calloc(row, sizeof(double));
    validateAction(denominatorsArr != NULL);

    sum = 0;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            sum += pow(matrix[i][j], 2);
        }
        denominatorsArr[i] = sqrt(sum);
        sum = 0;
    }
    return denominatorsArr;
}

/*
* Funcion: void descendingSort(EIGEN* eigens, int n)
* -----------------------------------------------------------------------------
* Params: EIGEN array, EIGEN array length
* Action: Sorts EIGEN array in descending order
* Return: None
*/
void descendingSort(EIGEN* eigens, int n) {
    int i, j, tmpEigenIndex;
    double tmpEigenVal;
 
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            if (eigens[i].eigenValue < eigens[j].eigenValue) {
                tmpEigenVal = eigens[i].eigenValue;
                eigens[i].eigenValue = eigens[j].eigenValue;
                eigens[j].eigenValue = tmpEigenVal;

                tmpEigenIndex = (eigens[i].eiganIndex);
                eigens[i].eiganIndex = (eigens[j].eiganIndex);
                eigens[j].eiganIndex = tmpEigenIndex;
            }
        }
    }
}
