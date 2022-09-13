#ifndef KMEANS_H_
#define KMEANS_H_
#define MAX_ITER_JACOBI 100
#define MAX_ITER_KMEANS 300
#define EPSLION 0.00001

/* Structs */ 
typedef struct 
{
    double eigenValue;
    int eiganIndex;
} EIGEN;

typedef struct 
{
    int vectors_count;
    double* vectors_sum;
    double* centroid;
} CLUSTER;

/* spkmeans.c */
void runGoal(char *goal, char *filename);

/* Utils.c */
int getVectorDim(char *filename);
int getVectorCount(char *filename);
double** getVectorsMatrix(char *filename, int N, int dim);
double** createMatrix(int n, int m);
double** createSquareMatrix(int n);
void printMatrix(double **matrix, int n, int m);
void printMatrixDiagonal(double **matrix, int n);
double euclideanNorm(double* vector1, double* vector2, int dim);
void freeMatrix(double **matrix, int n);
void validateAction(int bool);
void validateInput(int bool);

/* wam.c */
void wam(double** vectorsMatrix, int N, int vectorDim);
double **getWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
                                       int vectorDim);

/* ddg.c */
void ddg(double** vectorsMatrix, int N, int vectorDim);
double** getDiagonalDegreeMatrix(double **wam, int N);
double* getDdgDiagonal(double **wam, int N);

/* lnorm.c */
void lnorm(double** vectorsMatrix, int N, int vectorDim);
double** getLnorm(double *ddgDiagonal, double **wam, int N);

/* jacobi */
void jacobi(double** matrix, int N, int vectorDim);
double **jacobiAlgorithm(double **A, int n);
double off(double **A, int n);
double **createIdentityMatrix(int n);
int matrixIsDiagonal(double **A, int n);
void getIJOfLargestOffDiag(double **A, int n, int* pi, int* pj);
void getCSOfP(double **A, int i, int j, double *cp, double *sp);
void transformA(double **A, int n, int i, int j, double s, double c);
void getCurrentV(double **V, int n, int i, int j, double s, double c);

/* spk.c */
double** getNormalizedKEigenvectorsMatrix(int k, double** vectorsMatrix, 
                                          int N, int vectorDim);
EIGEN *createEigensArr(double **eiganVals, int n);
int eigengapHeuristic(EIGEN *eigenArray, int n);
double** createT(EIGEN* eigens, double** eigenVecMatrix, int k, int N);
void normalizeMatrixByRows(double **matrix, int row, int col);
double* getNormalizeDenominators(double **matrix, int row, int col);
void descendingSort(EIGEN* eigens, int n);

#endif
