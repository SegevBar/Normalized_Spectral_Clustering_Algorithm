#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"


/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: input from CMD
* Action: runs C program according to CMD arguments
* Return: 0 (None)
*/
int main(int argc, char *argv[]) {
   
    int k;
    char *goal;
    char *filename;
    if (argc == 4) {
        k = atoi(argv[1]);
        goal = argv[2];
        filename = argv[3];
        goalFunc(k, goal, filename);
    }
    return 0;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: k (required clusters), goal, file name
* Action: runs program according to CMD input. prints output according goal
* Return: None
*/
void goalFunc(int k, char* goal, char* filename) {

    if (strcmp(goal, "spk") == 0) {
        spk(k, filename);
    }
    if (strcmp(goal, "wam") == 0) {
        wam(filename);
    }
    if (strcmp(goal, "ddg") == 0) {
        ddg(filename);
    }
    if (strcmp(goal, "lnorm") == 0) {
        lnorm();
    }
    if (strcmp(goal, "jacobi") == 0) {
        jacobi(filename);
    }
}
