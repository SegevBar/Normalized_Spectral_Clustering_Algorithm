#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"


/*
* Funcion: int getVectorDim(char *filename)
* -----------------------------------------------------------------------------
* Params: Input file
* Action: Counts features in point
* Return: Point dimension
*/
int getVectorDim(char *filename) {
    int vectorDim = 1;
    int c;
    FILE *ifp;

    ifp = fopen(filename, "r");
    validateAction(ifp != NULL);

    /*find vector dimensions*/
    while ((c = fgetc(ifp)) != 10) {  /*run until end of line*/
        if (c == 44) {  /*if c == "," increment dimension*/
            vectorDim++; 
        }
    }
    fclose(ifp);
    return vectorDim;
} 

/*
* Funcion: int getVectorCount(char *filename)
* -----------------------------------------------------------------------------
* Params: Input file
* Action: Counts number of vectors in file
* Return: vectors count
*/
int getVectorCount(char *filename) {
    int N = 0;
    int c;
    FILE *ifp;

    ifp = fopen(filename, "r");
    validateAction(ifp != NULL);

    /*find N*/
    while ((c = fgetc(ifp)) != EOF) {  /*run until end of file*/
        
        if(c == 10) { /*if (c == "\n") increment N*/
            N++;
        }
    }
    fclose(ifp);
    return N;
}

/*
* Funcion: double** getVectorsMatrix(char *filename, int N, int dim)
* -----------------------------------------------------------------------------
* Params: Input file, pointers to vectors amount and point size
* Action: Reads vectors from inpuf file and save in metrix
* Return: Vectors matrix
*/
double** getVectorsMatrix(char *filename, int N, int dim) {
    int i, j;
    double** vectorsMatrix;
    FILE *ifp;
    char ch;
    double value;

    ifp = fopen(filename, "r");
    validateAction(ifp != NULL);

    /* allocate memory of vector matrix */ 
    vectorsMatrix = (double**) calloc(N, sizeof(*vectorsMatrix));
    validateAction(vectorsMatrix != NULL);
    for(i = 0; i < N; i++){
        vectorsMatrix[i] = (double*) calloc(dim ,sizeof(*vectorsMatrix[i]));
        validateAction(vectorsMatrix[i] != NULL);
    }
    /* fill matrix with data from file */
    for(i = 0; i < N; i++){
        for(j = 0; j < dim; j++){
            fscanf(ifp, "%lf%c", &value, &ch);
            vectorsMatrix[i][j] = value;
        }
    }
    fclose(ifp);

    return vectorsMatrix;
}

/*
* Funcion: double** createMatrix(int n, int m)
* -----------------------------------------------------------------------------
* Params: rows and columns number 
* Action: Allocate memory to zeros matrix n*m
* Return: Matrix
*/
double** createMatrix(int n, int m) {
    int i;
    double** matrix;

    /* allocate memory of matrix */ 
    matrix = (double**) calloc(n, sizeof(*matrix));
    validateAction(matrix != NULL);
    for(i = 0; i < n; i++){
        matrix[i] = (double*) calloc(m ,sizeof(*matrix[i]));
        validateAction(matrix[i] != NULL);
    }
    return matrix;
}

/*
* Funcion: double** createSquareMatrix(int n)
* -----------------------------------------------------------------------------
* Params: Matrix dimension
* Action: Allocate memory to zeros square matrix n*n
* Return: Matrix
*/
double** createSquareMatrix(int n) {
    return createMatrix(n, n);
}

/*
* Funcion: void printMatrix(double **matrix, int n, int m)
* -----------------------------------------------------------------------------
* Params: Matrix, Matrix size (2D)
* Action: Prints matrix
* Return: None
*/
void printMatrix(double **matrix, int n, int m) {
    int i;
    int j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%.4f", matrix[i][j]);
            if (j != m-1) {
                printf(",");
            } else {
                printf("\n");
            }
        }
    }
}

/*
* Funcion: void printMatrixDiagonal(double **matrix, int n)
* -----------------------------------------------------------------------------
* Params: Matrix, Matrix size (1D)
* Action: Prints the transposed matrix
* Return: None
*/
void printMatrixDiagonal(double **matrix, int n) {
    int i;
    
    for (i = 0; i < n; i++) {
        printf("%.4f", matrix[i][i]);
        if (i != n-1) {
            printf(",");
        } else {
            printf("\n");
        }
    }
}

/*
* Funcion: double euclideanNorm(double* vector1, double* vector2, int dim)
* -----------------------------------------------------------------------------
* Params: 2 vectors of the same dimension and dimension
* Action: Calculates res = sum(xi-v)^2 i=0,..,k-1
* Return: res
*/
double euclideanNorm(double* vector1, double* vector2, int dim) {
    double sum = 0.0;
    int j = 0;

    for (j = 0; j < dim; j++) {
        sum += (vector1[j]-vector2[j])*(vector1[j]-vector2[j]);
    } 
    return sum;
}

/*
* Funcion: void freeMatrix(double **matrix, int n)
* -----------------------------------------------------------------------------
* Params: Matrix, Matrix lines amount
* Action: Frees matrix memory
* Return: None
*/
void freeMatrix(double **matrix, int n) {
    int i;
    for(i = 0; i < n; i++){
        free(matrix[i]);
    }
    free(matrix);
}

/*
* Funcion: void validateAction(int bool)
* -----------------------------------------------------------------------------
* Params: boolean
* Action: Abort program if an error has occured
* Return: None
*/
void validateAction(int bool) {
    if (bool == 0) {
        printf("An Error Has Occurred\n");
        abort();
    }
}

/*
* Funcion: void validateInput(int bool)
* -----------------------------------------------------------------------------
* Params: boolean
* Action: Abort program if the input isn't valid
* Return: None
*/
void validateInput(int bool) {
    if (bool == 0) {
        printf("Invalid Input!\n");
        abort();
    }
    
}
