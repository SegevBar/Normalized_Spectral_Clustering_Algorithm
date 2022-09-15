#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"


/*
* Funcion: int main(int argc, char *argv[])
* -----------------------------------------------------------------------------
* Params: input from CMD
* Action: runs C program according to CMD arguments
* Return: Output according to goal
*/
int main(int argc, char *argv[]) {
    char *goal, *filename;

    /* parse CMD parameters and validate amount */
    validateInput(argc == 3);
    goal = argv[1];
    filename = argv[2];
    /* run one of the goals : wam, ddg, lnorm, jacobi */
    runGoal(goal, filename); 

    return 0;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: k (required clusters), goal, file name
* Action: runs program according to CMD input. prints output according goal
* Return: None
*/
void runGoal(char* goal, char* filename) {
    int N, vectorDim;
    double **vectorsMatrix;

    /* create matrix of vectors from file data points */
    vectorDim = getVectorDim(filename);  /* get vector dimension */
    N = getVectorCount(filename);  /* get vectors amount N */
    vectorsMatrix = getVectorsMatrix(filename, N, vectorDim);

    /* run wam */
    if (strcmp(goal, "wam") == 0) {
        wam(vectorsMatrix, N, vectorDim);
    }
    /* run ddg */
    else if (strcmp(goal, "ddg") == 0) {
        ddg(vectorsMatrix, N, vectorDim);
    }
    /* run lnorm */
    else if (strcmp(goal, "lnorm") == 0) {
        lnorm(vectorsMatrix, N, vectorDim);
    }
    /* run jacobi */
    else if (strcmp(goal, "jacobi") == 0) {
        jacobi(vectorsMatrix, N, vectorDim);
    }
    /* goal is not one of the above options */
    else {
        validateInput(0);
    }
}
