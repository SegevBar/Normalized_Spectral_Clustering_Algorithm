#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input file, pointers to points amount and point size
* Action: Reads points from inpuf file and save in metrix
* Return: Points matrix
*/
double **
getDataPoints(int *p_numOfVectors, int *p_numOfFeatures, char *filename) {
    
    int lineCounter, j;
    double *dataPointsArray, **dataPoints;
    char line[MAX_CHARS_LINE];
    char *left;
    FILE *f;

    f = fopen(filename, "r");
    ourAssert(f != NULL);

    fscanf(f, "%s", line);
    *p_numOfFeatures = featuresCount(line);
    dataPoints = (double **) calloc(MAX_LINES, sizeof(double *));
    ourAssert(dataPoints != NULL);

    dataPointsArray = (double *) calloc(MAX_LINES * (*p_numOfFeatures),
                                        sizeof(double));
    ourAssert(dataPointsArray != NULL);
    
    lineCounter = 0;
    do {
        dataPoints[lineCounter] =
                dataPointsArray + lineCounter * (*p_numOfFeatures);
        left = line;
        for (j = 0; j < *p_numOfFeatures; j++) {
            dataPoints[lineCounter][j] = strtod(left, &left);
            left++;
        }
        lineCounter++;
    } while (fscanf(f, "%s", line) != EOF);
    *p_numOfVectors = lineCounter;

    fclose(f);
    return dataPoints;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input file, pointers to points amount
* Action: Counts number of points in file
* Return: Points count
*/
int getVectorCount(char *filename) {

    int lineCounter;
    char line[MAX_CHARS_LINE];
    FILE *f;

    f = fopen(filename, "r");
    ourAssert(f != NULL);
    fscanf(f, "%s", line);

    lineCounter = 0;
    do {
        lineCounter++;
    } while (fscanf(f, "%s", line) != EOF);
    
    fclose(f);
    
    return lineCounter;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Point
* Action: Counts features in point
* Return: Point size
*/
int featuresCount(const char *line) {
    int i;
    int count = 1;
    for (i = 0; line[i] != '\0'; i++) {
        if (line[i] == ',') {
            count++;
        }
    }
    return count;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Symmetric Matrix size (1D)
* Action: Creates symetric matrix - values in the diagonal or bottom triangle
* Return: Symmetric Matrix
*/
double **createSymmetricMatrix(int n) {

    int i;
    double *array;
    double *cur;
    double **symMatrix;
    double d_n = (double) n;
    /* calculate how many values we need to save for symmetric matrix - how many
     * values are in the diagonal or the bottom triangle */
    int lenSymMatrix = (int) ((d_n * d_n) / 2 + d_n / 2);
    array = (double *) calloc(lenSymMatrix, sizeof(double));
    ourAssert(array != NULL);

    symMatrix = (double **) calloc(n, sizeof(double *));
    ourAssert(symMatrix != NULL);

    cur = array;
    for (i = 0; i < n; i++) {
        symMatrix[i] = cur + i;
        cur = symMatrix[i];
    }
    return symMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Length of squared identity matrix
* Action: Creates regular squared matrix
* Return: Squared matrix
*/
double **createRegularSquareMatrix(int n) {

    return createRegularMatrix(n, n);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: rows and columns number 
* Action: Create empty matrix
* Return: Matrix
*/
double **createRegularMatrix(int rows, int columns) {

    int i;
    double *array;
    double **matrix;

    array = (double *) calloc(rows * columns, sizeof(double));
    ourAssert(array != NULL);
    matrix = (double **) calloc(rows, sizeof(double *));
    ourAssert(matrix != NULL);

    for (i = 0; i < rows; i++) {
        matrix[i] = array + i * columns;
    }
    return matrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Length of squared identity matrix
* Action: Create identity matrix
* Return: Identity matrix
*/
double **identityMatrix(int n) {

    int i;
    double **matrix;
    matrix = createRegularSquareMatrix(n);
    for (i = 0; i < n; i++) {
        matrix[i][i] = 1;
    }
    return matrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Symmetric Matrix, Matrix size (1D)
* Action: Prints matrix
* Return: None
*/
void printSymmetricMatrix(double **matrix, int lenMatrix) {
    int i;
    int j;

    for (i = 0; i < lenMatrix; i++) {
        for (j = 0; j < lenMatrix; j++) {
            if (i < j) {
                printf("%.4f", round(matrix[j][i]));
            } else {
                printf("%.4f", round(matrix[i][j]));
            }
            if (j != lenMatrix - 1) {
                printf(",");
            } else {
                printf("\n");
            }
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Matrix, Matrix size (2D)
* Action: Prints matrix
* Return: None
*/
void printMatrix(double **matrix, int rows, int columns) {

    int i;
    int j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            if (j != columns - 1) {
                printf("%.4f,", round(matrix[i][j]));
            } else {
                if (i != rows - 1) {
                    printf("%.4f\n", round(matrix[i][j]));
                } else {
                    printf("%.4f", round(matrix[i][j]));
                }
            }
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Matrix, Matrix size (2D)
* Action: Prints transposed matrix
* Return: None
*/
void printTransposedMatrix(double **matrix, int rows, int columns) {

    int i;
    int j;
    for (j = 0; j < columns; j++) {
        for (i = 0; i < rows; i++) {
            if (i != rows - 1) {
                printf("%.4f,", round(matrix[i][j]));
            } else {
                if (j != columns - 1) {
                    printf("%.4f\n", round(matrix[i][j]));
                } else {
                    printf("%.4f", round(matrix[i][j]));
                }
            }
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Matrix, Matrix size (1D)
* Action: Prints the transposed matrix
* Return: None
*/
void printDiagonal(double **matrix, int lenMatrix) {

    int i;
    for (i = 0; i < lenMatrix; i++) {
        if (i != lenMatrix - 1) {
            printf("%.4f,", round(matrix[i][i]));
        } else {
            printf("%.4f\n", round(matrix[i][i]));
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Matrix
* Action: Frees matrix memory
* Return: None
*/
void freeMatrix(double **matrix) {

    free(matrix[0]);
    free(matrix);
}

/*
* Funcion TO DELETE
*/
void ourAssert(int trueOrFalse) {

    if (trueOrFalse == 0) {
        printf("An Error Has Occured");
        abort();
    }
}