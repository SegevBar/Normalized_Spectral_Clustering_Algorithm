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
* Funcion: 
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
* Funcion: 
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

double** createSquareMatrix(int n) {
    int i;
    double** matrix;

    /* allocate memory of matrix */ 
    matrix = (double**) calloc(n, sizeof(*matrix));
    validateAction(matrix != NULL);
    for(i = 0; i < n; i++){
        matrix[i] = (double*) calloc(n ,sizeof(*matrix[i]));
        validateAction(matrix[i] != NULL);
    }
    return matrix;
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
    validateAction(array != NULL);
    matrix = (double **) calloc(rows, sizeof(double *));
    validateAction(matrix != NULL);

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
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Symmetric Matrix, Matrix size (1D)
* Action: Prints matrix
* Return: None
*/
void printSymmetricMatrix(double **matrix, int n) {
    int i;
    int j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i < j) {
                printf("%.4f", round(matrix[j][i]));
            } else {
                printf("%.4f", round(matrix[i][j]));
            }
            if (j != n - 1) {
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
void printMatrix(double **matrix, int n, int m) {
    int i;
    int j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%.4f", round(matrix[i][j]));
            if (j != m-1) {
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
void printDiagonal(double **matrix, int n) {

    int i;
    for (i = 0; i < n; i++) {
        if (i != n - 1) {
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
