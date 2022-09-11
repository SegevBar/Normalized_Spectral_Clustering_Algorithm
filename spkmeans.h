#ifndef KMEANS_H_
#define KMEANS_H_
#define MAX_ITER_JACOBI 100
#define MAX_ITER_KMEANS 300
#define EPSLION 0.00001
#define round(x)((((x)>-0.00005)&&((x)<=0))?(0):(x))

/* Structs */ 
struct cluster {
    double *centroid;
    double *centroid_closest;
    int size;
    int equalTolLast;
};
typedef struct cluster CLUSTER;

struct eigen {
    double *eigenVector;
    double eigenValue;
};
typedef struct eigen EIGEN;

/* spkmeans.c */
void runGoal(int k, char *goal, char *filename);

/* Utils.c */
double** getVectorsMatrix(char *filename, int N, int dim);
int getVectorCount(char *filename);
int getVectorDim(char *filename);
double** createSquareMatrix(int n);
void freeMatrix(double **matrix, int n);
void validateAction(int bool);
void validateInput(int bool);



double **createRegularMatrix(int rows, int columns);
double **identityMatrix(int n);
void printSymmetricMatrix(double **matrix, int n);
void printMatrix(double **matrix, int rows, int columns);
void printTransposedMatrix(double **matrix, int rows, int columns);
void printDiagonal(double **matrix, int n);


/* spk.c */
double** spk(int k, double** vectorsMatrix, int N, int vectorDim);
double **getEigenvectorsMatrix(double **matrix, int *kp, int n);
EIGEN *getEigensArray(double **vectorsMatrix, double **valuesMatrix,
                           int n);
int eigengapHeuristic(EIGEN *eigenArray, int arrLength);
void normalizeMatrix(double **matrix, int rows, int columns);
double* getRootOfSumOfSquares(double **matrix, int rows, 
                                        int columns);
void descendingSort(EIGEN* eigenArray, int arrLength);


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
int matrixIsDiagonal(double **A, int n);
void getIJOfLargestOffDiag(double **A, int n, int* pi, int* pj);
void getCSOfP(double **A, int i, int j, double *cp, double *sp);
void transformA(double **A, int n, int i, int j, double s, double c);
void getCurrentEigenvectors(double **V, int n, int i, int j, double s,
                           double c);

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
