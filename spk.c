#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: double** spk(int k, double** vectorsMatrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: k, Vectors matrix, Vector count, Vector dimension
* Action: performs The Normalized Spectral Clustering Algorithm *WITHOUT*
*         kmeans (first step executes at python program)
* Return: runKMeans ? kkmeansmain : Eigenvectors matrix
*/
double** spk(int k, double** vectorsMatrix, int N, int vectorDim) {
    double **wam, *ddgDiagonal;
    double **lnorm, **matrix;

    wam = getWeightedAdjacencyMatrix(vectorsMatrix, N, vectorDim);
    freeMatrix(vectorsMatrix, N);              
    
    ddgDiagonal = getDdgDiagonal(wam, N);
    lnorm = getLnorm(ddgDiagonal, wam, N);
    freeMatrix(wam, N);
    free(ddgDiagonal);

    printMatrix(lnorm, N, N);
    freeMatrix(lnorm, N);

    matrix = getEigenvectorsMatrix(lnorm, &k, N);
    freeMatrix(lnorm, N);
    
    normalizeMatrix(matrix, N, k);
    
    return matrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Symmetric Matrix, k pointer, Matrix size (1D)
* Action: Execute jacobi algorithm, eigengap heuristic and create returned 
*         matrix
* Return: Matrix With Eigenvectors As Columns
*/
double **getEigenvectorsMatrix(double **matrix, int *kp,
                                                  int n) {

    EIGEN *eigenArray; /* array of EIGENS -
     * structs that contain an eigenvalue and a pointer to its eigenvector. */
    double **eigenvectorsMatrix;
    double **newMatrix;
    int i;
    int j;
    eigenvectorsMatrix = jacobiAlgorithm(matrix, n);
    eigenArray = getEigensArray(eigenvectorsMatrix, matrix, n);
    freeMatrix(matrix, n);
    
    descendingSort(eigenArray, n);
    /* get k with eigengap heuristic if k == 0 */
    if (*kp == 0) {
        *kp = eigengapHeuristic(eigenArray, n);
    }
    newMatrix = createRegularMatrix(n, *kp);
    for (j = 0; j < *kp; j++) {
        for (i = 0; i < n; i++) {
            newMatrix[i][j] = *(eigenArray[j].eigenVector + n * i);
        }
    }
    freeMatrix(eigenvectorsMatrix, n);
    free(eigenArray);
    return newMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Eigenvectors Matrix, Eigenvalues Matrix, Matrix size (1D)
* Action: Create array of EIGENS representing eigenvalue and its' eigenvector
* Return: Array of EIGENS representing eigenvalue and its' eigenvector
*/
EIGEN *getEigensArray(double **vectorsMatrix, double **valuesMatrix,
                           int n) {
    EIGEN *eigenArray;
    int i;
    eigenArray = (EIGEN *) calloc(n, sizeof(EIGEN));
    validateAction(eigenArray != NULL);
    for (i = 0; i < n; i++) {
        eigenArray[i].eigenValue = valuesMatrix[i][i];
        eigenArray[i].eigenVector = &vectorsMatrix[0][i];
    }
    return eigenArray;
}


/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: EIGENS array and it's size
* Action: Execute the eigengap heuristic
* Return: k - number of clusters
*/
int eigengapHeuristic(EIGEN *eigenArray, int arrLength) {
    double max;
    double cur;
    int max_index;
    int i;
    max_index = 0;
    max = eigenArray[1].eigenValue - eigenArray[0].eigenValue;
    for (i = 1; i < arrLength / 2; i++) {
        cur = eigenArray[i + 1].eigenValue - eigenArray[i].eigenValue;
        if (cur > max) {
            max_index = i;
            max = cur;
        }
    }
    return max_index + 1;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Matrix and its' dimensions
* Action: Normalizes values in matrix
* Return: None
*/
void normalizeMatrix(double **matrix, int rows, int columns) {
 
    double *arraySumOfSquaresRows;
    int i;
    int j;
    arraySumOfSquaresRows = calculateRootOfSumOfSquaresRows(matrix, rows,
                                                            columns);
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            matrix[i][j] = matrix[i][j] / arraySumOfSquaresRows[i];
        }
    }
    free(arraySumOfSquaresRows);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Matrix and its' dimensions
* Action: Calaulate sum of (u_ij)^2 for each row for normalize algo
* Return: Array with calculated sum of each row
*/
double* calculateRootOfSumOfSquaresRows(double **matrix, int rows, 
                                        int columns) {
    double *arraySumOfSquaresColumns;
    double sum;
    int j;
    int i;
    sum = 0;
    arraySumOfSquaresColumns = (double *) calloc(rows, sizeof(double));
    validateAction(arraySumOfSquaresColumns != NULL);
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            sum += pow(matrix[i][j], 2);
        }
        arraySumOfSquaresColumns[i] = sqrt(sum);
        sum = 0;
    }
    return arraySumOfSquaresColumns;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: EIGEN array, EIGEN array length
* Action: Sorts EIGEN array in descending order
* Return: None
*/
void descendingSort(EIGEN* eigenArray, int arrLength) {
    int i, j;
    double tmpEigenVal, *tmpEigenVector;
 
    for (i = 0; i < arrLength; ++i) {
        for (j = i + 1; j < arrLength; ++j) {
            if (eigenArray[i].eigenValue < eigenArray[j].eigenValue) {
                tmpEigenVal = eigenArray[i].eigenValue;
                eigenArray[i].eigenValue = eigenArray[j].eigenValue;
                eigenArray[j].eigenValue = tmpEigenVal;

                tmpEigenVector = (eigenArray[i].eigenVector);
                eigenArray[i].eigenVector = (eigenArray[j].eigenVector);
                eigenArray[j].eigenVector = tmpEigenVector;
            }
        }
    }
}
