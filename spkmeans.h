#ifndef KMEANS_H_
#define KMEANS_H_

#define MAX_ITER_KMEANS 300
#define MAX_ITER_JACOBI 100
#define EPSLION 0.000000000000001
#define MAX_CHARS_LINE 111 /* at most 10 features '-XXXX.XXXX,' +\0  */
#define MAX_CHARS_LINE_MATRIX 551 /* at most 50 features '-XXXX.XXXX,' +\0  */
#define MAX_LINES 50
#define round(x)((((x)>-0.00005)&&((x)<=0))?(0):(x))
#define min(a,b) (((a) < (b)) ? (a) : (b))

/* ------------------------spk - step 6, Kmeans ----------------------------- */

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

void spk(int k, char* filename);

void wam(char* filename);

void ddg(char* filename);

void lnorm(char* filename);

void jacobi(char* filename);

void initCluster(CLUSTER *curCluster, double *dataPoint, int numOfFeatures);

CLUSTER *initializeClusters(double **dataPoints, int K);

void
updateClosest(CLUSTER *curCluster, const double *datapoint, int numOfFeatures);

void updateCentroid(CLUSTER *curCluster, int numOfFeatures);

double euclideanNorm(const double *datapoint1, const double *datapoint2,
                     int numOfFeatures);

void kmeansmain(CLUSTER *clusters, double **dataPoints, int numOfFeatures,
                int numOfVectors);

void printCentroids(CLUSTER *clusters, int numOfFeatures);

void freeClusters(CLUSTER *clusters, int K);

/* --------------------------- spk - steps 1+2 ------------------------------ */

double **createWeightedAdjacencyMatrix(double **dataPoints, int numOfVectors,
                                       int numOfFeaturs);

double *calculateDiagonalDegreeMatrix(double **weightedAdjacencyMatrix,
                                      int numOfVectors);

double **createDDGMatrixforDDG(double **weightedAdjacencyMatrix, int numOfVectors);

/*void computePowOfMinusHalf(double *diagonalDegreeMatrix, int numOfVectors);*/

double **
createLnorm(double *diagonalDegreeArray, double **weightedAdjacencyMatrix, 
            int numOfVectors);

/* ----------------------- spk - steps 3-5, Jacobi -------------------------- */


double **calculateMatrixWithEigenvectorsAsColumns(double **matrix, int *p_k,
                                                  int lenMatrix);

double **jacobiAlgorithm(double **matrix, int lenMatrix);

void updateEigenvectors(double **matrix, int lenMatrix, int i, int j, double s,
                        double c);

void rotateMatrix(double **matrix, int i, int j, double *p_s, double *p_c);

void transformMatrix(double **matrix, int lenMatrix, int i, int j, double s,
                     double c);

void calculateMax(double **matrix, int lenMatrix, int *p_i, int *p_j);

double off(double **matrix, int lenMatrix);

int checkIfDiagonalMatrix(double **matrix, int lenMatrix);

EIGEN *createArrayOfEigens(double **vectorsMatrix, double **valuesMatrix,
                           int lenMatrix);

int eigengapHeuristic(EIGEN *eigenArray, int lenArray);

void normalizeMatrix(double **matrix, int rows, int columns);

double *
calculateRootOfSumOfSquaresRows(double **matrix, int rows, int columns);

/* --------------------------------- utils ---------------------------------- */

double **createSymmetricMatrix(int n);

double **createRegularSquareMatrix(int n);

double **createRegularMatrix(int rows, int columns);

double **identityMatrix(int n);

void printSymmetricMatrix(double **matrix, int lenMatrix);

void printMatrix(double **matrix, int rows, int columns);

void printTransposedMatrix(double **matrix, int rows, int columns);

void printDiagonal(double **matrix, int lenMatrix);

void freeMatrix(double **matrix);

void mergeSort(EIGEN eigenArray[], int lenArray);

void merge(EIGEN eigenArray[], int leftStart, int leftEnd, int rightEnd);

int compareEIGEN(EIGEN e1, EIGEN e2);

void ourAssert(int trueOrFalse);

/* ------------------------------ input data -------------------------------- */

double **getDataPoints(int *numOfVectors, int *numOfFeatures, char *filename);

double **readSymatricMatrixFromFile(char *filename, int *p_lenMatrix);

int featuresCount(const char *line);

/* --------------------------------- main ----------------------------------- */

double **normalizedSpectralClustering(int k, char *filename, int runKMeans);

void goalFunc(int k, char *goal, char *filename);

#endif
