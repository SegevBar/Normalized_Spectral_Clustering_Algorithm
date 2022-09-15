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
* Params: k pointer, Vectors matrix, Vector amount N, Vector dimension
* Action: performs The Normalized Spectral Clustering Algorithm *WITHOUT*
*         kmeans (the step after - kmeans++ - executes at python program)
* Return: Normalized K Eigenvectors Matrix (T)
*/
double** getNormalizedKEigenvectorsMatrix(int *kp, double** vectorsMatrix, 
                                          int N, int vectorDim) {
    double **wam, *ddgDiagonal, **lnorm, **eigenVecMatrix, **T;
    EIGEN* eigens;

    /* Form The Weighted Adjacency Matrix W from X */
    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    
    /* Form the diagonal of The Diagonal Degree Matrix D from W */
    ddgDiagonal = getDdgDiagonal(wam, N);
    /* Form The Normalized Graph Laplacian */
    lnorm = getLnorm(ddgDiagonal, wam, N);

    freeMatrix(wam, N);
    free(ddgDiagonal);

    /* Find Eigenvalues and Eigenvectors using Jacobi algorithm */
    eigenVecMatrix = jacobiAlgorithm(lnorm, N);  /* diagonalized lnorm */
    eigens = createEigensArr(lnorm, N); 
    freeMatrix(lnorm, N);
    
    /* Sort Eigenvalues from largest to smallest*/
    descendingSort(eigens, N); 
    /* if k = 0 find k using the Eigengap Heuristic */
    *kp = (*kp == 0) ? eigengapHeuristic(eigens, N) : *kp;

    /* Form Normalized K Eigenvectors Matrix T */
    T = createT(eigens, eigenVecMatrix, *kp, N);
    freeMatrix(eigenVecMatrix, N);

    return T;
}

/*
* Funcion: EIGEN *createEigensArr(double **eiganVals, int n)
* -----------------------------------------------------------------------------
* Params: Eigenvalues Matrix, Matrix size n (size is n*n) 
* Action: Create array of EIGENS representing eigenvalue and its' matching 
*         eigenvectors' index
* Return: Array of EIGENS representing eigenvalue and its' eigenvector index
*/
EIGEN *createEigensArr(double **eiganVals, int n) {
    EIGEN *eigenArray;
    int i;

    /* create EIGEN array */
    eigenArray = (EIGEN*) calloc(n, sizeof(EIGEN));
    validateAction(eigenArray != NULL);

    for (i = 0; i < n; i++) {
        /* get Eigenvalue from Eigenvalues Matrix diagonal */
        eigenArray[i].eigenValue = eiganVals[i][i];
        /* get index of Eigenvalues and Eigenvector (same index) */
        eigenArray[i].eiganIndex = i;
    }
    return eigenArray;
}

/*
* Funcion: int eigengapHeuristic(EIGEN *eigenArray, int n)
* -----------------------------------------------------------------------------
* Params: EIGENS array and it's size n
* Action: Execute the eigengap heuristic
* Return: k - number of clusters
*/
int eigengapHeuristic(EIGEN *eigenArray, int n) {
    double max, curMax;
    int maxIndex, i;
    
    maxIndex = 0;
    /* find argmaxi(deltai) for i=1,..,n/2 */
    max = fabs(eigenArray[0].eigenValue - eigenArray[1].eigenValue);
    for (i = 1; i < n / 2; i++) {
        /* calculate dletai = |lambda(i)-lambda(i+1)| */
        curMax = fabs(eigenArray[i].eigenValue - eigenArray[i+1].eigenValue);
        if (curMax > max) {  /* update max if needed */
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
* Params: Descending sorted eigens array, eigan vectors, k (colunms), N (rows)
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
            /* copy eiganvector (colunm) according to eiganvalue index */
            T[i][j] = eigenVecMatrix[i][eigens[j].eiganIndex];
        }
    }
    free(eigens);
    /* normalize matrix */
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

    /* get array of calculated denominators = sqrt(sum((uij)^2)) for each row */
    denominatorsArr = getNormalizeDenominators(matrix, row, col);
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            if (denominatorsArr[i] != 0) {
                /* update matrix to normalized values */
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
    
    /* allocate memory for array rows size */
    denominatorsArr = (double *) calloc(row, sizeof(double));
    validateAction(denominatorsArr != NULL);

    /* for each row - sum (u_ij)^2 and save in denominators array */
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
* Params: EIGEN array, EIGEN array length n
* Action: Sorts EIGEN array in descending order
* Return: None
*/
void descendingSort(EIGEN* eigens, int n) {
    int i, j, tmpEigenIndex;
    double tmpEigenVal;
    
    /* kind of selection Sort */
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            /* if current eiganvalue smaller then - swap */
            if (eigens[i].eigenValue < eigens[j].eigenValue) {
                /* swap eiganvalues */
                tmpEigenVal = eigens[i].eigenValue;
                eigens[i].eigenValue = eigens[j].eigenValue;
                eigens[j].eigenValue = tmpEigenVal;

                /* swap eiganvalues index in original diagonal matrix */
                tmpEigenIndex = (eigens[i].eiganIndex);
                eigens[i].eiganIndex = (eigens[j].eiganIndex);
                eigens[j].eiganIndex = tmpEigenIndex;
            }
        }
    }
}
