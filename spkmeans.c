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

    validateInput(argc == 3);
    goal = argv[1];
    filename = argv[2];
    runGoal(0, goal, filename);

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

    vectorDim = getVectorDim(filename);
    N = getVectorCount(filename);
    vectorsMatrix = getVectorsMatrix(filename, N, vectorDim);

    if (strcmp(goal, "wam") == 0) {
        wam(vectorsMatrix, N, vectorDim);
    }
    else if (strcmp(goal, "ddg") == 0) {
        ddg(vectorsMatrix, N, vectorDim);
    }
    else if (strcmp(goal, "lnorm") == 0) {
        lnorm(vectorsMatrix, N, vectorDim);
    }
    else if (strcmp(goal, "jacobi") == 0) {
        jacobi(vectorsMatrix, N, vectorDim);
    }
    else {
        validateInput(0);
    }
}
