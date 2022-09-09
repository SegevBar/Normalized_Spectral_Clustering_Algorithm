#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"

static PyObject *kmeans(PyObject *self, PyObject *args);

static double **pythonDataPointsMatrixToDataPoints(PyObject *pythonDataPoints, int numOfDataPoints, int numOfFeatures);

static CLUSTER *pythonInitializeClusters(PyObject *pythonClusters, int numOfVectors);

static void pythonInitCluster(CLUSTER *curCluster, int numOfFeatures);

static PyObject *spkWithoutKmeans(PyObject *self, PyObject *args);

static PyObject *nkMatrixToPython(double **nkMatrix, int numOfVectors, int k);

static PyObject *goalsOtherThenSpk(PyObject *self, PyObject *args);

static PyObject *kmeans(PyObject *self, PyObject *args) {
    /* this function will be called by python.
    * It gets a PyObjects containing the data points and a PyObject containing
    * the clusters and number of features.
    * The function runs the kmeans algorithm */
    PyObject *pythonDataPoints;
    PyObject *pythonClusters;
    CLUSTER *clusters;
    int k, numOfVectors;
    double **dataPoints;

    if (!PyArg_ParseTuple(args, "00i", &pythonDataPoints, &pythonClusters, &k)) {
        return NULL;
    }
    numOfVectors = PyList_Size(pythonDataPoints);
    dataPoints = pythonDataPointsMatrixToDataPoints(pythonDataPoints, numOfVectors, k);
    clusters = pythonInitializeClusters(pythonClusters, k);
    kmeansmain(clusters, dataPoints, k, numOfVectors);

    Py_RETURN_NONE;
}

static double **pythonDataPointsMatrixToDataPoints(PyObject *pythonDataPoints, int numOfVectors, int numOfFeatures) {
    /* This function recieves a python object containing dataPoints, the number of data points
    and the number of features.
    This function copies the python object data points into a dynamic array and returns it. */

    double **dataPoints, *dataPointsArray;
    int i; int j;
    PyObject *dataPoint;

    dataPoints = (double **) calloc(numOfVectors, sizeof(double *));
    ourAssert(dataPoints != NULL);
    dataPointsArray = (double *) calloc(numOfFeatures * numOfVectors, sizeof(double));
    ourAssert(dataPointsArray != NULL);
    
    for (i = 0; i < numOfVectors; i++) {
        dataPoints[i] = dataPointsArray + i * (numOfFeatures);
        dataPoint = PyList_GetItem(pythonDataPoints, i);

        for (j = 0; j < numOfFeatures; j++) {
            dataPoints[i][j] = PyFloat_AsDouble(PyList_GetItem(dataPoint, j));
        }
    }
    
    return dataPoints;
}

static CLUSTER *pythonInitializeClusters(PyObject *pythonClusters, int numOfFeatures) {
    /* This function gets a python object containing the clusters and the number
    of features. The function then copies the clusters into a struct array,
    CLUSTERS, initializes it and returns the array. */

    CLUSTER *clusters;
    PyObject *centroid;
    int i; int j;
    
    clusters = calloc(numOfFeatures, sizeof(CLUSTER));
    ourAssert(clusters != NULL);

    for (i = 0; i < numOfFeatures; i++) {
        pythonInitCluster(&clusters[i], numOfFeatures);
        centroid = PyList_GetItem(pythonClusters, i);

        for (j = 0; j < numOfFeatures; j++) {
            *(clusters[i].centroid + j) = PyFloat_AsDouble(PyList_GetItem(centroid, j));
        }
    }

    return clusters;
}

static void pythonInitCluster(CLUSTER *curCluster, int numOfFeatures) {
    /* This function receives a specific cluster struct pointer and the number
    of features. The function then initializes the values that are related to the
    cluster struct */

    curCluster -> centroid_closest = calloc(numOfFeatures, sizeof(double));
    ourAssert(curCluster -> centroid_closest != NULL);

    curCluster -> size = 0;
    curCluster -> equalTolLast = 0; // 0 == False
}

static PyObject *spkWithoutKmeans(PyObject *self, PyObject *args) {
    /* This function get the argument k and the file name from python.
    The function then returns a PyObject containing the dataPoints which are
    calculated in steps 1 - 5 of the algorithm. */

    int k, numOfVectors;
    char *filename;
    double **nkMatrix;
    
    if (!PyArg_ParseTuple(args, "is", &k, &filename)) {
        return NULL;
    }
    numOfVectors = getVectorCount(filename);
    nkMatrix = normalizedSpectralClustering(k, filename, 0);

    return nkMatrixToPython(nkMatrix, numOfVectors, k);
}

static PyObject *nkMatrixToPython(double **nkMatrix, int numOfVectors, int k) {
    /* This functions receives a dynamic matrix with data points, the number of vectors
    and the number of features. The function then returns a PyObject containing the
    data points. */

    PyObject *python_listOflist;
    PyObject *python_list;
    PyObject *python_val;
    int i; int j;
    python_listOflist = PyList_New(numOfVectors);

    for (i = 0; i < numOfVectors; ++i) {
        python_list = PyList_New(k);
        for (j = 0; j < k; ++j) {
            python_val = Py_BuildValue("d", nkMatrix[i][j]);
            PyList_SetItem(python_list, j, python_val);
        }
        PyList_SetItem(python_listOflist, i, python_list);
    }

    return python_listOflist;
}

static PyObject *goalsOtherThenSpk(PyObject *self, PyObject *args) {
    /* This function calls the appropriate function, per requested by the goal variable:
    wam, ddg, lnorm or jacobi. */

    int k;
    char *goal;
    char *filename;

    if (!PyArg_ParseTuple(args, "iss", &k, &goal, &filename)) {
        return NULL;
    }

    goalFunc(k, goal, filename);
    Py_RETURN_NONE;
}

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
