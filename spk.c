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
* Funcion TO DELETE
*/
double **stepsOneToFive(double **dataPoints, int *p_k, int numOfVectors,
                        int numOfFeatures) {
    
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
    
    mergeSort(eigenArray, lenMatrix);
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
