#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

/* spkmeansmodule.c */
static PyObject *getPythonNormalizedKEigenvectorsMatrix(PyObject *self, 
                                                        PyObject *args);
static PyObject *runGoalOfCProgram(PyObject *self, PyObject *args);
static PyObject *sendMatrixToPython(double **matrix, int N, int k);
static PyObject *runKmeansFromCProgram(PyObject *self, PyObject *args);
static CLUSTER *initPyClusters(PyObject *pyCentroids, int k);
static PyObject *kmeans(PyObject *vectors_py, CLUSTER *clusters, int k, int N);
void calcCluster(double* vector, CLUSTER* clusters, int k, int dim);
int updateCentroids(CLUSTER* clusters, int k, int dim, double epsilon);
static PyObject *sendCentroidsToPython(CLUSTER *clusters, int k);


/*
* Funcion: static PyObject *getPythonNormalizedKEigenvectorsMatrix(
           PyObject *self, PyObject *args)
* -----------------------------------------------------------------------------
* Params: k and filename from python program
* Action: Parses arguments from python and runs C func to find T matrix
* Return: PyObject Normalized K Eigenvectors Matrix (T)
*/
static PyObject *getPythonNormalizedKEigenvectorsMatrix(PyObject *self, 
                                                        PyObject *args) {
    int k, N, vectorDim;
    char *filename;
    double **T, **vectorsMatrix;

    /* Parses arguments from python */
    if (!PyArg_ParseTuple(args, "is", &k, &filename)) {
        return NULL;
    }
    /* create matrix of vectors from file data points */
    N = getVectorCount(filename);
    vectorDim = getVectorDim(filename);
    vectorsMatrix = getVectorsMatrix(filename, N, vectorDim);

    T = getNormalizedKEigenvectorsMatrix(&k, vectorsMatrix, N, vectorDim);

    return sendMatrixToPython(T, N, k);
}

/*
* Funcion: static PyObject *runGoalofCProgram(PyObject *self, PyObject *args)
* -----------------------------------------------------------------------------
* Params: goal and filename from python program
* Action: Parses arguments from python and runs goal from c program
*         Goals: wam, ddg, lnorm, jacobi
* Return: None
*/
static PyObject *runGoalOfCProgram(PyObject *self, PyObject *args) {
    char *goal, *filename;

    /* Parses arguments from python */
    if (!PyArg_ParseTuple(args, "ss", &goal, &filename)) {
        return NULL;
    }
    runGoal(goal, filename); /* runs goal from C program*/
    Py_RETURN_NONE;
}

/*
* Funcion: static PyObject *sendMatrixToPython(double **matrix, int N, int k)
* -----------------------------------------------------------------------------
* Params: Matrix and it's dimensions N*k
* Action: Creates a PyObject matrix from input matrix
* Return: PyObject matrix
*/
static PyObject *sendMatrixToPython(double **matrix, int N, int k) {
    PyObject *pyMatrix, *pyList, *pyVal;
    int i; int j;

    /* parse C matrix to python matrix */
    pyMatrix = PyList_New(N);
    for (i = 0; i < N; ++i) {
        pyList = PyList_New(k);
        for (j = 0; j < k; ++j) {
            pyVal = Py_BuildValue("d", matrix[i][j]);
            PyList_SetItem(pyList, j, pyVal);
        }
        PyList_SetItem(pyMatrix, i, pyList);
    }
    freeMatrix(matrix, N);
    return pyMatrix;
}

/*
* Funcion: static PyObject *runKmeansFromCProgram(PyObject *self, 
*          PyObject *args)
* -----------------------------------------------------------------------------
* Params: Vectors, Centroids, N and k from python kmeans++ program
* Action: Parses arguments from python and run kmeans from C program
* Return: None
*/
static PyObject *runKmeansFromCProgram(PyObject *self, PyObject *args) {
    PyObject *pyVectors, *pyCentroids;
    CLUSTER *clusters;
    int k, N;

    /* Parses arguments from python */
    if (!(PyArg_ParseTuple(args, "00ii", &pyVectors, &pyCentroids, &N, &k))) {
        return NULL;
    }

    clusters = initPyClusters(pyCentroids, k);    
    kmeans(pyVectors, clusters, k, N);
    Py_RETURN_NONE;
}

/*
* Funcion: static CLUSTER *initPyClusters(PyObject *pyCentroids, int k)
* -----------------------------------------------------------------------------
* Params: pyObject of list of centroids, k
* Action: Init clusters with centroids from kmeans++ 
* Return: pointer to array of clusters
*/
static CLUSTER *initPyClusters(PyObject *pyCentroids, int k) {
    CLUSTER *clusters;
    int i, j, curr;
      
    clusters = (CLUSTER *)calloc(k, sizeof(CLUSTER));
    validateAction(clusters == NULL);

    curr = 0;
    for (i = 0; i < k; i++) {
        clusters[i].centroid = (double *)calloc(k, sizeof(double));
        validateAction(clusters[i].centroid == NULL);

        for (j = 0; j < k; j++) {
            clusters[i].centroid[j] = PyFloat_AsDouble(PyList_GetItem(pyCentroids, curr));
            curr++;
        }

        clusters[i].vectors_count = 0;
        clusters[i].vectors_sum = (double *)calloc(k, sizeof(double));
        validateAction(clusters[i].vectors_sum == NULL);
    }
    return clusters;
}

/*
* Funcion: static PyObject *kmeans(PyObject vectors_py, CLUSTER *clusters, 
*          int k, int N)
* -----------------------------------------------------------------------------
* Params: pyObject of list of data points, clusters array, k, N
* Action: Performs kmeans main loop
* Return: pyObject matrix of calculated centroids
*/
static PyObject *kmeans(PyObject *vectors_py, CLUSTER *clusters, int k, int N) {
    double epsilon = 0;
    double *curr_vector;
    int cnt, i, j, curr;
    int has_converged = 0;

    /* kmeans main loop */
    cnt = 0;
    while ((cnt < MAX_ITER_KMEANS) && (!has_converged)) {
        curr = 0;
        /*find current vector cluster*/
        for (i = 0; i < N; i++) {
            curr_vector = (double*)calloc(k, sizeof(double));
            validateAction(curr_vector == NULL);

            for (j = 0; j < k; j++) {
                curr_vector[j] = PyFloat_AsDouble(
                                        PyList_GetItem(vectors_py, curr));
                curr++;
            }
            calcCluster(curr_vector, clusters, k, k);
            free(curr_vector);
        }
        /*update centroids*/
        has_converged = updateCentroids(clusters, k, k, epsilon);
        
        /*reset*/
        for(i = 0; i < k; i++) {
            clusters[i].vectors_count = 0;
            for(j = 0; j < k; j++) {
                clusters[i].vectors_sum[j] = 0;
            }
        }
        cnt++;
    }
    return sendCentroidsToPython(clusters, k);
}

/*
* Funcion: void calcCluster(double* vector, CLUSTER* clusters, int k, int dim)
* -----------------------------------------------------------------------------
* Params: a vector, clusters array, k and vector dimension
* Action: Finds closest cluster to current vector and update closest cluster
* Return: None
*/
void calcCluster(double* vector, CLUSTER* clusters, int k, int dim) {
    double min_distance = -1.0;
    int closest_cluster = -1;
    double distance;
    int i, j;

    /*find closest cluster to current vector*/
    for (i = 0; i < k; i++) {
        distance = euclideanNorm(vector, clusters[i].centroid, dim);
        if ((distance < min_distance) || (min_distance < 0)) {
            min_distance = distance; 
            closest_cluster = i;
        }
    }
    /*update closest cluster*/
    clusters[closest_cluster].vectors_count++; 
    for (j = 0; j < dim; j++) {
        clusters[closest_cluster].vectors_sum[j] += vector[j];
    }
}

/*
* Funcion: updateCentroids(CLUSTER* clusters, int k, int dim, double epsilon)
* -----------------------------------------------------------------------------
* Params: clusters array, k, vector dimension and epsilon
* Action: Calculates new centroid from clusters (vectors sum)/(vectors count)
*         and checks convergence
* Return: boolean representing convergence (1-yes, 0-no)
*/
int updateCentroids(CLUSTER* clusters, int k, int dim, double epsilon) {
    int i, j;
    int has_converged = 1;
    double* new_centroid = NULL;
    double dist = 0;

    /*calculate new centroid*/
    for (i = 0; i < k; i++) {
        new_centroid = (double*)calloc(dim, sizeof(double));
        validateAction(new_centroid == NULL);

        for (j = 0; j < dim; j++) {
            new_centroid[j] = (clusters[i].vectors_sum[j]/
                               clusters[i].vectors_count);
        }
        dist = sqrt(euclideanNorm(clusters[i].centroid, new_centroid, dim));

        /*check if convergence did not accured*/
        has_converged = (dist >= epsilon) ? 0 : 1;

        /*update centroid*/
        memcpy(clusters[i].centroid, new_centroid, sizeof(double)*dim);

        free(new_centroid);
    }
    return has_converged;
}

/*
* Funcion: static PyObject *sendCentroidsToPython(CLUSTER *clusters, int k)
* -----------------------------------------------------------------------------
* Params: clusters array, k
* Action: Converts centroids from C to python
* Return: pyObject matrix of calculated centroids
*/
static PyObject *sendCentroidsToPython(CLUSTER *clusters, int k) {
    PyObject *clusters_py, *value, *curr_vector;
    int i, j;

    clusters_py = PyList_New(k);
    for (i = 0; i < k; i++) {
        curr_vector = PyList_New(k);
        for (j = 0; j < k; j++) {
            value = Py_BuildValue("d", clusters[i].centroid[j]);
            PyList_SetItem(curr_vector, j, value);
        }
        /*add PyObject centroid to PyList clusters*/
        PyList_SetItem(clusters_py, i, curr_vector);
    }
    /*free clusters memory*/
    for (i = 0; i < k; i++) {
        free(clusters[i].centroid);
        free(clusters[i].vectors_sum);
    }
    free(clusters);
    
    return clusters_py;
}


/* API */ 
static PyMethodDef _methods[] = {
    {"runKmeansFromCProgram",
            (PyCFunction)runKmeansFromCProgram,
            METH_VARARGS,
            PyDoc_STR("Performs k-means algorithm")},
    {"getPythonNormalizedKEigenvectorsMatrix",
            (PyCFunction)getPythonNormalizedKEigenvectorsMatrix,
            METH_VARARGS,
            PyDoc_STR("Calculates Normalized K Eigenvectors Matrix")},
    {"runGoalOfCProgram",
            (PyCFunction)runGoalOfCProgram,
            METH_VARARGS,
            PyDoc_STR("Performs wam, ddg, lnorm or jacobi steps")},
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
