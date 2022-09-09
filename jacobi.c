#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"


/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input file
* Action: Calculate and output the eigenvalues and eigenvectors
* Return: Prints eigenvalues and eigenvectors
*/
void jacobi(char* filename) {
    int numOfVectors;
    double **matrix, **eigenvectorsMatrix;

    matrix = readSymatricMatrixFromFile(filename, &numOfVectors);
    eigenvectorsMatrix = jacobiAlgorithm(matrix, numOfVectors);
    
    printDiagonal(matrix, numOfVectors);
    printTransposedMatrix(eigenvectorsMatrix, numOfVectors, numOfVectors);
    
    freeMatrix(matrix);
    freeMatrix(eigenvectorsMatrix);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Input file, pointer to points amount
* Action: Reads points from inpuf file and save in metrix
* Return: Points matrix
*/
double **readSymatricMatrixFromFile(char *filename, int *p_lenMatrix) {

    FILE *f;
    int i;
    int j;
    char line[MAX_CHARS_LINE_MATRIX];
    char *left;
    double *array;
    double *cur;
    double **symMatrix;
    double d_n;
    int lenArray;

    f = fopen(filename, "r");
    ourAssert(f != NULL);
    fscanf(f, "%s", line);
    *p_lenMatrix = featuresCount(line);

    d_n = (double) *p_lenMatrix;
    /* calculate how many values we need to save for symmetric matrix - how many
     * values are in the diagonal or the bottom triangle */
    lenArray = (int) ((d_n * d_n) / 2 + d_n / 2);

    array = (double *) calloc(lenArray, sizeof(double));
    ourAssert(array != NULL);
    symMatrix = (double **) calloc(*p_lenMatrix, sizeof(double *));
    ourAssert(symMatrix != NULL);

    cur = array;
    for (i = 0; i < *p_lenMatrix; i++) {
        symMatrix[i] = cur + i;
        cur = symMatrix[i];
        left = line;
        for (j = 0; j <= i; j++) {
            symMatrix[i][j] = strtod(left, &left);
            left++;
        }
        fscanf(f, "%s", line);
    }
    fclose(f);
    return symMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix and it's dimension (1D)
* Action: Execute jacobi algorithm
* Return: Eigenvectors matrix
*/
double **jacobiAlgorithm(double **matrix, int lenMatrix) {

    int i;
    int j;
    int q;
    double c;
    double s;
    double offMatrixPre;
    double offMatrixPost;
    double **eigenvectorsMatrix;
    offMatrixPre = off(matrix, lenMatrix);
    eigenvectorsMatrix = identityMatrix(lenMatrix);

    for (q = 0; q < MAX_ITER_JACOBI; q++) {
        if (checkIfDiagonalMatrix(matrix, lenMatrix)) {
            break;
        }
        calculateMax(matrix, lenMatrix, &i, &j);
        rotateMatrix(matrix, i, j, &s, &c);
        transformMatrix(matrix, lenMatrix, i, j, s, c);
        updateEigenvectors(eigenvectorsMatrix, lenMatrix, i, j, s, c);
        offMatrixPost = off(matrix, lenMatrix);
        if (offMatrixPre - offMatrixPost <= EPSLION) {
            break;
        }
        offMatrixPre = offMatrixPost;
    }
    return eigenvectorsMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix and it's dimension (1D)
* Action: Checks if the matrix is a diagonal matrix
* Return: diagonal matrix ? 1 : 0
*/
int checkIfDiagonalMatrix(double **matrix, int lenMatrix) {

    int i;
    int j;
    for (i = 0; i < lenMatrix; i++) {
        for (j = 0; j < i; j++) {
            if (matrix[i][j] != 0) {
                return 0; /* False */
            }
        }
    }
    return 1; /* True */
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix and it's dimension (1D), 
*         pointers to the indexes ij of maximum off-diagonal value
* Action: Finds maximum absolut off-diagonal value and update ij
* Return: None
*/
void calculateMax(double **matrix, int lenMatrix, int *p_i, int *p_j) {

    int t;
    int q;
    double max;
    *p_i = 1;
    *p_j = 0;
    max = fabs(matrix[1][0]);
    for (t = 0; t < lenMatrix; t++) {
        for (q = 0; q < t; q++) {
            if (fabs(matrix[t][q]) > max) {
                max = fabs(matrix[t][q]);
                *p_i = t;
                *p_j = q;
            }
        }
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix, indexes of max off-diagonal value, 
*         pointers to s and c
* Action: calculates and updates s and c
* Return: None
*/
void rotateMatrix(double **matrix, int i, int j, double *p_s, double *p_c) {

    double theta;
    double t;
    /* matrix[i][j] is in the bottom triangle so the formulas have been
     * updated accordingly. */
    theta = (matrix[i][i] - matrix[j][j]) / (2 * matrix[i][j]);
    if (theta >= 0) {
        t = 1 / (theta + sqrt(pow(theta, 2) + 1));
    } else {
        t = -1 / (-theta + sqrt(pow(theta, 2) + 1));
    }
    *p_c = 1 / sqrt(pow(t, 2) + 1);
    *p_s = t * *p_c;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix and it's dimension (1D), indexes of max 
*         off-diagonal value, s, c
* Action: Transforms the matrix according to the rotation matrix
* Return: None
*/
void transformMatrix(double **matrix, int lenMatrix, int i, int j, double s,
                     double c) {

    int r;
    double *p_r_i;
    double *p_r_j;
    double temp_r_i;
    double temp_r_j;
    double temp_i_i;
    double temp_j_j;
    double temp_i_j;
    /* matrix[i][j] is in the bottom triangle, so the formulas have been
     * updated accordingly. */
    for (r = 0; r < lenMatrix; r++) {
        if (r != i && r != j) {
            /* we save only values where rows >= columns because the matrix is
             * symmetric so we want to update only the values that we save */
            if (r < i) {
                p_r_i = &matrix[i][r];
            } else {
                p_r_i = &matrix[r][i];
            }
            if (r < j) {
                p_r_j = &matrix[j][r];
            } else {
                p_r_j = &matrix[r][j];
            }
            temp_r_i = *p_r_i;
            temp_r_j = *p_r_j;
            *p_r_i = c * temp_r_i + s * temp_r_j;
            *p_r_j = c * temp_r_j - s * temp_r_i;
        }
    }
    temp_i_i = matrix[i][i];
    temp_j_j = matrix[j][j];
    temp_i_j = matrix[i][j];
    matrix[i][i] =
            pow(s, 2) * temp_j_j + pow(c, 2) * temp_i_i + 2 * s * c * temp_i_j;
    matrix[j][j] =
            pow(c, 2) * temp_j_j + pow(s, 2) * temp_i_i - 2 * s * c * temp_i_j;
    matrix[i][j] = 0;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix and it's dimension (1D), indexes of max 
*         off-diagonal value, s, c
* Action: Multiplies matrix with the rotation matrix
* Return: None
*/
void updateEigenvectors(double **matrix, int lenMatrix, int i, int j, double s,
                        double c) {

    int q;
    double temp_i;
    double temp_j;
    for (q = 0; q < lenMatrix; q++) {
        temp_i = matrix[q][i];
        temp_j = matrix[q][j];
        /* other values of the matrix are without change, because except of the
         * 4 s-c numbers, the rotation matrix is like the identity matrix. */
        matrix[q][j] = temp_j * c + temp_i * -s;
        matrix[q][i] = temp_j * s + temp_i * c;
    }
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: a Symmetric Matrix and it's dimension (1D)
* Action: Executes off function
* Return: off funtion result
*/
double off(double **matrix, int lenMatrix) {

    double sum;
    int i;
    int j;
    sum = 0;
    for (i = 0; i < lenMatrix; i++) {
        for (j = 0; j < i; j++) {
            sum += pow(matrix[i][j], 2);
        }
    }
    return 2 * sum;
}
