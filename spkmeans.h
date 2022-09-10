#ifndef KMEANS_H_
#define KMEANS_H_

#define MAX_ITER_KMEANS 300
#define MAX_ITER_JACOBI 100
#define EPSLION 0.000000000000001
#define MAX_CHARS_LINE 111 /* at most 10 features '-XXXX.XXXX,' +\0  */
#define MAX_CHARS_LINE_MATRIX 551 /* at most 50 features '-XXXX.XXXX,' +\0  */
#define MAX_LINES 1000
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

double **createSymmetricMatrix(int n);
double **createRegularSquareMatrix(int n);
double **createRegularMatrix(int rows, int columns);
double **identityMatrix(int n);
void printSymmetricMatrix(double **matrix, int lenMatrix);
void printMatrix(double **matrix, int rows, int columns);
void printTransposedMatrix(double **matrix, int rows, int columns);
void printDiagonal(double **matrix, int lenMatrix);
void freeMatrix(double **matrix, int n);
void validateAction(int trueOrFalse);

/* spk.c */
double** spk(int k, double** vectorsMatrix, int N, int vectorDim);
double **calculateMatrixWithEigenvectorsAsColumns(double **matrix, int *kp,
                                                  int lenMatrix);
EIGEN *createArrayOfEigens(double **vectorsMatrix, double **valuesMatrix,
                           int lenMatrix);
int eigengapHeuristic(EIGEN *eigenArray, int arrLength);
void normalizeMatrix(double **matrix, int rows, int columns);
double* calculateRootOfSumOfSquaresRows(double **matrix, int rows, 
                                        int columns);
void descendingSort(EIGEN* eigenArray, int arrLength);


/* wam.c */
void wam(double** vectorsMatrix, int N, int vectorDim);
double **createWeightedAdjacencyMatrix(double **vectorsMatrix, int N,
                                       int vectorDim);

/* ddg.c */
void ddg(double** vectorsMatrix, int N, int vectorDim);
double** createDDGMatrixforDDG(double **weightedAdjacencyMatrix, 
                               int N);
double *calculateDiagonalDegreeMatrix(double **weightedAdjacencyMatrix,
                                      int N);

/* lnorm.c */
void lnorm(double** vectorsMatrix, int N, int vectorDim);
double** createLnorm(double *diagonalDegreeArray, 
                     double **weightedAdjacencyMatrix, int N);

/* jacobi */
void jacobi(double** vectorsMatrix, int N, int vectorDim);
double **jacobiAlgorithm(double **matrix, int lenMatrix);
int checkIfDiagonalMatrix(double **matrix, int lenMatrix);
void calculateMax(double **matrix, int lenMatrix, int *p_i, int *p_j);
void rotateMatrix(double **matrix, int i, int j, double *p_s, double *p_c);
void transformMatrix(double **matrix, int lenMatrix, int i, int j, double s,
                     double c);
void updateEigenvectors(double **matrix, int lenMatrix, int i, int j, double s,
                        double c);
double off(double **matrix, int lenMatrix);

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
