#!/bin/bash
# Script to compile and execute a c program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c wam.c ddg.c lnorm.c jacobi.c utils.c normalizedKEigenvectorsMatrix.c kmeans.c -lm -o spkmeans
