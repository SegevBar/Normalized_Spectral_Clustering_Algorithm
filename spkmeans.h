#ifndef KMEANS_H_
#define KMEANS_H_
#define MAX_ITER_JACOBI 100
#define MAX_ITER_KMEANS 300
#define EPSLION 0.00001
#define round(x)((((x)>-0.00005)&&((x)<=0))?(0):(x))

/* Structs */ 
typedef struct 
{
    double *eigenVector;
    double eigenValue;
} EIGEN;

typedef struct 
{
    double *centroid;
    double *centroid_closest;
    int size;
    int equalTolLast;
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
EIGEN *createEigensArr(double **eigenVectors, double **eiganVals, int n);
int eigengapHeuristic(EIGEN *eigenArray, int n);
double** createT(EIGEN* eigens, int k, int N);
void normalizeMatrixByRows(double **matrix, int row, int col);
double* getNormalizeDenominators(double **matrix, int row, int col);
void descendingSort(EIGEN* eigens, int n);

/* kmeans */
void kmeansmain(CLUSTER *clusters, double **vectorsMatrix, int vectorDim,
                int N);
CLUSTER *initializeClusters(double **vectorsMatrix, int K);
void initCluster(CLUSTER *curCluster, double *dataPoint, int vectorDim);
void updateClosest(CLUSTER *curCluster, const double *datapoint, int vectorDim);
void updateCentroid(CLUSTER *curCluster, int vectorDim);
double euclideanNorm(const double *datapoint1, const double *datapoint2,
                     int vectorDim);
void printCentroids(CLUSTER *clusters, int vectorDim);
void freeClusters(CLUSTER *clusters, int K);

#endif
