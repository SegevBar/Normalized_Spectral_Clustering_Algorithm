#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"


/*
* Funcion: void jacobi(double** matrix, int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: Vectors matrix, Vector count, Vector dimension
* Action: Calculate and output the eigenvalues and eigenvectors
* Return: Prints eigenvalues and eigenvectors
*/
void jacobi(double** matrix, int N, int vectorDim) {
    double **eigenvectorsMatrix;

    eigenvectorsMatrix = jacobiAlgorithm(matrix, N);
    
    printMatrixDiagonal(matrix, N);
    printMatrix(eigenvectorsMatrix, N, vectorDim);
    
    freeMatrix(matrix, N);
    freeMatrix(eigenvectorsMatrix, N);
}

/*
* Funcion: double **jacobiAlgorithm(double **A, int n)
* -----------------------------------------------------------------------------
* Params: A matrix and it's dimension (1D)
* Action: Execute jacobi algorithm
* Return: Eigenvectors matrix
*/
double **jacobiAlgorithm(double **A, int n) {
    int i, j, k;
    double c, s, offA, offAt;
    double **V; /* Eigenvectors matrix */
    
    offA = off(A, n);
    V = createIdentityMatrix(n);

    for (k = 0; k < MAX_ITER_JACOBI; k++) {
        if (matrixIsDiagonal(A, n)) {  /*break loop if A is diagonal*/
            break;
        }
        getIJOfLargestOffDiag(A, n, &i, &j); /*updates i j via pointers*/
        getCSOfP(A, i, j, &c, &s); /*updates s c via pointers*/
        transformA(A, n, i, j, s, c); 
        getCurrentV(V, n, i, j, s, c);
        
        /* check convergence */
        offAt = off(A, n);
        if (offA - offAt <= EPSLION) {
            break;
        }
        offA = offAt;
    }
    return V;
}

/*
* Funcion: double off(double **A, int n)
* -----------------------------------------------------------------------------
* Params: A matrix and it's dimension (1D)
* Action: Executes off function
* Return: off funtion result
*/
double off(double **A, int n) {
    int i, j;
    double sum = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                sum += pow(A[i][j], 2);
            }
        }
    }
    return sum; 
}

/*
* Funcion: double **createIdentityMatrix(int n)
* -----------------------------------------------------------------------------
* Params: Length of squared identity matrix
* Action: Create identity matrix
* Return: Identity matrix
*/
double **createIdentityMatrix(int n) {
    int i;
    double **matrix;

    matrix = createSquareMatrix(n); /* allocate memory */
    for (i = 0; i < n; i++) {
        matrix[i][i] = 1;
    }
    return matrix;
}

/*
* Funcion: int matrixIsDiagonal(double **A, int n)
* -----------------------------------------------------------------------------
* Params: A matrix and it's dimension (1D)
* Action: Checks if the matrix is a diagonal matrix
* Return: diagonal matrix ? 1 : 0
*/
int matrixIsDiagonal(double **A, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j && A[i][j] != 0.0) {
                return 0; 
            }
        }
    }
    return 1; 
}

/*
* Funcion: void getIJOfLargestOffDiag(double **A, int n, int* pi, int* pj)
* -----------------------------------------------------------------------------
* Params: A matrix and it's dimension (1D), pointers to indexes ij of maximum
*         off-diagonal value
* Action: Finds maximum absolut off-diagonal value and update ij
* Return: None
*/
void getIJOfLargestOffDiag(double **A, int n, int* pi, int* pj) {
    int i, j;
    double currMax = 0.0;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            if (i != j && fabs(A[i][j]) > currMax) {
                currMax = fabs(A[i][j]);
                *pi = i;
                *pj = j;
            }
        }
    }
}

/*
* Funcion: void getCSOfP(double **A, int i, int j, double *cp, double *sp)
* -----------------------------------------------------------------------------
* Params: A matrix, indexes of max off-diagonal value, 
*         pointers to s and c
* Action: calculates and updates s and c
* Return: None
*/
void getCSOfP(double **A, int i, int j, double *cp, double *sp) {
    double theta, t;
    int sign;

    /* theta = (Ajj-Aii)/(2Ajj) */
    theta = (A[i][i] - A[j][j]) / (2 * A[i][j]);
    /* t = sign(theta)/(|thetha|+sqrt((theta)^2+1) */
    sign = (theta >= 0) ? 1 : -1;
    t = sign / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    /* c = 1/sqrt(t^2+1) */
    *cp = 1 / sqrt(pow(t, 2) + 1);
    /* s = t*c */
    *sp = (t)*(*cp);
}

/*
* Funcion: void transformA(double **A, int n, int i, int j, double s, double c)
* -----------------------------------------------------------------------------
* Params: A matrix and it's dimension (1D), indexes of max 
*         off-diagonal value, s, c
* Action: Transforms the matrix formulas
* Return: None
*/
void transformA(double **A, int n, int i, int j, double s, double c) {
    int r;
    double Aii, Ajj, Aij, Ari, Arj;

    Aii = A[i][i];
    Ajj = A[j][j];
    Aij = A[i][j]; 
    
    /* a'ij = 0 */  
    A[i][j] = 0;
    /* a'ii = c^2*aii + s^2*ajj - 2*s*c*aij */
    A[i][i] = pow(c, 2) * Aii + pow(s, 2) * Ajj - 2 * s * c * Aij;
    /* a'jj = s^2*aii + c^2*ajj + 2*s*c*aij */
    A[j][j] = pow(s, 2) * Aii + pow(c, 2) * Ajj + 2 * s * c * Aij;
    
    /* if r != i,j */
    for (r = 0; r < n; r++) {
        if (r != i && r != j) {
            Ari = A[r][i];
            Arj = A[r][j];
            /* a'ri = c*ari - s*arj */
            A[r][i] = c * Ari + s * Arj;
            /* a'rj = c*arj + s*ari */
            A[r][j] = c * Arj - s * Ari;
        }
    }
    
}

/*
* Funcion: void getCurrentV(double **V, int n, int i, int j, double s, double c)
* -----------------------------------------------------------------------------
* Params: Eigenvectors matrix (V) and it's dimension (1D), indexes of max 
*         off-diagonal value, s, c
* Action: Multiplies current V matrix with the new P
* Return: None
*/
void getCurrentV(double **V, int n, int i, int j, double s, double c) {
    int k;
    double Vki, Vkj;

    for (k = 0; k < n; k++) {
        Vki = V[k][i];
        Vkj = V[k][j];
        V[k][i] = Vkj * s + Vki * c;
        V[k][j] = Vkj * c + Vki * -s;
    }
}
