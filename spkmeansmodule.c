#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"

/* spkmeansmodule.c */
static PyObject *spkWithoutKmeans(PyObject *self, PyObject *args);
static PyObject *nkMatrixToPython(double **nkMatrix, int N, int k);
static PyObject *goalsOtherThenSpk(PyObject *self, PyObject *args);
static PyObject *kmeans(PyObject *self, PyObject *args);
static double **getVectorsFromPython(PyObject *pythonVectorsMatrix, int N, 
                                    int vectorDim);
static CLUSTER *pythonInitializeClusters(PyObject *pythonClusters, int N);
static void pythonInitCluster(CLUSTER *curCluster, int vectorDim);

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static PyObject *spkWithoutKmeans(PyObject *self, PyObject *args) {
    /* This function get the argument k and the file name from python.
    The function then returns a PyObject containing the vectors which are
    calculated in steps 1 - 5 of the algorithm. */

    int k, N, vectorDim;
    char *filename;
    double** T, vectorsMatrix;
    
    if (!PyArg_ParseTuple(args, "is", &k, &filename)) {
        return NULL;
    }

    N = getVectorCount(filename);
    vectorDim = getVectorDim(filename);
    vectorsMatrix = getVectorsMatrix(filename, N, vectorDim);

    T = getNormalizedKEigenvectorsMatrix(k, vectorsMatrix, N, vectorDim);

    return sendMatrixToPython(T, N, k);
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: Matrix and it's dimensions n*m
* Action: Creates a PyObject matrix from input matrix
* Return: PyObject matrix
*/
static PyObject *sendMatrixToPython(double **matrix, int n, int m) {
    PyObject *pyMatrix;
    PyObject *pyList;
    PyObject *pyVal;
    int i; int j;

    pyMatrix = PyList_New(N);
    for (i = 0; i < n; ++i) {
        pyList = PyList_New(k);
        for (j = 0; j < m; ++j) {
            pyVal = Py_BuildValue("d", matrix[i][j]);
            PyList_SetItem(pyList, j, pyVal);
        }
        PyList_SetItem(pyMatrix, i, pyList);
    }
    return pyMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static PyObject *goalsOtherThenSpk(PyObject *self, PyObject *args) {
    /* This function calls the appropriate function, per requested by the goal variable:
    wam, ddg, lnorm or jacobi. */

    char *goal;
    char *filename;

    if (!PyArg_ParseTuple(args, "ss", &goal, &filename)) {
        return NULL;
    }

    runGoal(goal, filename);
    Py_RETURN_NONE;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static PyObject *kmeans(PyObject *self, PyObject *args) {
    /* this function will be called by python.
    * It gets a PyObjects containing the data points and a PyObject containing
    * the clusters and number of features.
    * The function runs the kmeans algorithm */
    PyObject *pythonVectors;
    PyObject *pythonClusters;
    CLUSTER *clusters;
    int k, N;
    double **vectorsMatrix;

    if (!PyArg_ParseTuple(args, "00i", &pythonVectors, &pythonClusters, &k)) {
        return NULL;
    }
    N = PyList_Size(pythonVectors);
    vectorsMatrix = getVectorsFromPython(pythonVectors, N, k);
    clusters = pythonInitializeClusters(pythonClusters, k);
    kmeansmain(clusters, vectorsMatrix, k, N);

    Py_RETURN_NONE;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static double **getVectorsFromPython(PyObject *pythonVectorsMatrix, int N, 
                                    int vectorDim) {
    /* This function recieves a python object containing points, the number of data points
    and the number of features.
    This function copies the python object data points into a dynamic array and returns it. */

    double **vectorsMatrix, *vectorsArray;
    int i; int j;
    PyObject *point;

    vectorsMatrix = (double **) calloc(N, sizeof(double *));
    validateAction(vectorsMatrix != NULL);
    vectorsArray = (double *) calloc(vectorDim * N, sizeof(double));
    validateAction(vectorsArray != NULL);
    
    for (i = 0; i < N; i++) {
        vectorsMatrix[i] = vectorsArray + i * (vectorDim);
        point = PyList_GetItem(pythonVectorsMatrix, i);

        for (j = 0; j < vectorDim; j++) {
            vectorsMatrix[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    return vectorsMatrix;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static CLUSTER *pythonInitializeClusters(PyObject *pythonClusters, int vectorDim) {
    /* This function gets a python object containing the clusters and the number
    of features. The function then copies the clusters into a struct array,
    CLUSTERS, initializes it and returns the array. */

    CLUSTER *clusters;
    PyObject *centroid;
    int i; int j;
    
    clusters = calloc(vectorDim, sizeof(CLUSTER));
    validateAction(clusters != NULL);

    for (i = 0; i < vectorDim; i++) {
        pythonInitCluster(&clusters[i], vectorDim);
        centroid = PyList_GetItem(pythonClusters, i);

        for (j = 0; j < vectorDim; j++) {
            *(clusters[i].centroid + j) = PyFloat_AsDouble(PyList_GetItem(centroid, j));
        }
    }

    return clusters;
}

/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static void pythonInitCluster(CLUSTER *curCluster, int vectorDim) {
    /* This function receives a specific cluster struct pointer and the number
    of features. The function then initializes the values that are related to the
    cluster struct */

    curCluster -> centroid_closest = calloc(vectorDim, sizeof(double));
    validateAction(curCluster -> centroid_closest != NULL);

    curCluster -> size = 0;
    curCluster -> equalTolLast = 0; // 0 == False
}

/* API */ 
static PyMethodDef _methods[] = {
    {"kmeans",
            (PyCFunction) kmeans,
                METH_VARARGS,
                    PyDoc_STR("find a partition of N unlabeled observations in k distinct "
                              "clusters using the k-means algorith")},
    {"spkWithoutKmeans",
            (PyCFunction) spkWithoutKmeans,
                METH_VARARGS,
                    PyDoc_STR("perform normalized spectral clustering without  "
                              "k-means level")},
    {"goalsOtherThenSpk",
            (PyCFunction) goalsOtherThenSpk,
                METH_VARARGS,
                    PyDoc_STR("perform wam, ddg, lnorm or jacobi")},
    {NULL, NULL, 0, NULL} 
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void) {
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
