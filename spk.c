#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/* ***unite with normalizedSpectralClustering
* Funcion: 
* -----------------------------------------------------------------------------
* Params: k, file name
* Action: performs The Normalized Spectral Clustering Algorithm
* Return: None
*/
void spk(int k, char* filename) {
    normalizedSpectralClustering(k, filename, 1);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: k, file name
* Action: performs The Normalized Spectral Clustering Algorithm
* Return: runKMeans ? kkmeansmain : Eigenvectors matrix
*/
double** normalizedSpectralClustering(int k, char *filename, int runKMeans) {
    int numOfVectors, numOfFeatures;
    double **dataPoints, **weightedAdjacencyMatrix, *diagonalDegreeArray; 
    double **lnorm, **matrix;
    CLUSTER *clusters;

    dataPoints = getDataPoints(&numOfVectors, &numOfFeatures, filename);
    weightedAdjacencyMatrix = createWeightedAdjacencyMatrix(dataPoints,         
                              numOfVectors, numOfFeatures);
    freeMatrix(dataPoints);              
    diagonalDegreeArray = calculateDiagonalDegreeMatrix(
                                    weightedAdjacencyMatrix, numOfVectors);
    lnorm = createLnorm(diagonalDegreeArray, 
                                    weightedAdjacencyMatrix, numOfVectors);
    freeMatrix(weightedAdjacencyMatrix);
    free(diagonalDegreeArray);
   
    matrix = calculateMatrixWithEigenvectorsAsColumns(lnorm, &k, numOfVectors);
    freeMatrix(lnorm);
    
    normalizeMatrix(matrix, numOfVectors, k);
    clusters = initializeClusters(matrix, k);
    if (runKMeans == 1) {
        kmeansmain(clusters, matrix, k, numOfVectors);
    }
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
double **calculateMatrixWithEigenvectorsAsColumns(double **matrix, int *p_k,
                                                  int lenMatrix) {

    EIGEN *eigenArray; /* array of EIGENS -
     * structs that contain an eigenvalue and a pointer to its eigenvector. */
    double **eigenvectorsMatrix;
    double **newMatrix;
    int i;
    int j;
    eigenvectorsMatrix = jacobiAlgorithm(matrix, lenMatrix);
    eigenArray = createArrayOfEigens(eigenvectorsMatrix, matrix, lenMatrix);
    freeMatrix(matrix);
    
    descendingSort(eigenArray, lenMatrix);
    if (*p_k == 0) {
        *p_k = eigengapHeuristic(eigenArray, lenMatrix);
    }
    newMatrix = createRegularMatrix(lenMatrix, *p_k);
    for (j = 0; j < *p_k; j++) {
        for (i = 0; i < lenMatrix; i++) {
            newMatrix[i][j] = *(eigenArray[j].eigenVector + lenMatrix * i);
        }
    }
    freeMatrix(eigenvectorsMatrix);
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
EIGEN *createArrayOfEigens(double **vectorsMatrix, double **valuesMatrix,
                           int lenMatrix) {
    EIGEN *eigenArray;
    int i;
    eigenArray = (EIGEN *) calloc(lenMatrix, sizeof(EIGEN));
    ourAssert(eigenArray != NULL);
    for (i = 0; i < lenMatrix; i++) {
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
int eigengapHeuristic(EIGEN *eigenArray, int lenArray) {

    double max;
    double cur;
    int max_index;
    int i;
    max_index = 0;
    max = eigenArray[1].eigenValue - eigenArray[0].eigenValue;
    for (i = 1; i < lenArray / 2; i++) {
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
double *
calculateRootOfSumOfSquaresRows(double **matrix, int rows, int columns) {

    double *arraySumOfSquaresColumns;
    double sum;
    int j;
    int i;
    sum = 0;
    arraySumOfSquaresColumns = (double *) calloc(rows, sizeof(double));
    ourAssert(arraySumOfSquaresColumns != NULL);
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
void descendingSort(EIGEN* eigenArray, int lenArray) {
    int i, j;
    double tmpEigenVal, *tmpEigenVector;
 
    for (i = 0; i < lenArray; ++i) {
        for (j = i + 1; j < lenArray; ++j) {
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


/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: EIGEN array, EIGEN array length
* Action: Performs merge sort algorithm
* Return: None
*/
void mergeSort(EIGEN eigenArray[], int lenArray) {

    int currSize; /* cur size of the subarrays that we merge */
    int leftStart;
    for (currSize = 1; currSize <= lenArray - 1; currSize = 2 * currSize) {
        for (leftStart = 0;
             leftStart < lenArray - 1; leftStart += 2 * currSize) {
            int leftEnd = min(leftStart + currSize - 1, lenArray - 1);
            int rightEnd = min(leftStart + 2 * currSize - 1, lenArray - 1);
            merge(eigenArray, leftStart, leftEnd, rightEnd);
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: EIGEN array, pointers to indexes in arrays
* Action: Merges 2 arrays sorted arrays (merge sort)
* Return: None
*/
void merge(EIGEN eigenArray[], int leftStart, int leftEnd, int rightEnd) {

    int i, j, k;
    EIGEN *leftArray, *rightArray;
    int lenLeftArray = leftEnd - leftStart + 1;
    int lenRightArray = rightEnd - leftEnd; /* rightStart = leftEnd + 1 */

    leftArray = (EIGEN *) calloc(lenLeftArray, sizeof(EIGEN));
    ourAssert(leftArray != NULL);
    rightArray = (EIGEN *) calloc(lenRightArray, sizeof(EIGEN));
    ourAssert(rightArray != NULL);

    /* copy values to temp leftArray and rightArray */
    for (i = 0; i < lenLeftArray; i++)
        leftArray[i] = eigenArray[leftStart + i];
    for (j = 0; j < lenRightArray; j++)
        rightArray[j] = eigenArray[leftEnd + 1 + j];

    i = 0;
    j = 0;
    k = leftStart;
    /* merge the two subarrays */
    while (i < lenLeftArray && j < lenRightArray) {
        if (compareEIGEN(leftArray[i], rightArray[j]) <= 0) {
            eigenArray[k] = leftArray[i];
            i++;
        } else {
            eigenArray[k] = rightArray[j];
            j++;
        }
        k++;
    }

    /* copy the reminded values */
    while (i < lenLeftArray) {
        eigenArray[k] = leftArray[i];
        i++;
        k++;
    }
    while (j < lenRightArray) {
        eigenArray[k] = rightArray[j];
        j++;
        k++;
    }
    free(leftArray);
    free(rightArray);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 2 EIGENs
* Action: Compares 2 EIGENs
* Return: if EIGENs1 == EIGENs2 -> 0
*         if EIGENs1 > EIGENs2 -> 1
*         if EIGENs1 <>> EIGENs2 -> -1
*/
int compareEIGEN(EIGEN e1, EIGEN e2) {

    if (e1.eigenValue > e2.eigenValue) {
        return 1;
    } else if (e1.eigenValue < e2.eigenValue) {
        return -1;
    }
    return 0;
}
