#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"

/* spkmeansmodule.c */
static PyObject *getPythonNormalizedKEigenvectorsMatrix(PyObject *self, 
                                                        PyObject *args);
static PyObject *sendMatrixToPython(double **matrix, int n, int m);
static PyObject *runGoalOfCProgram(PyObject *self, PyObject *args);
static PyObject *runKmeansFromCProgram(PyObject *self, PyObject *args);
static double **getVectorsFromPython(PyObject *pythonVectorsMatrix, int N, 
                                    int vectorDim);
static CLUSTER *pythonInitializeClusters(PyObject *pythonClusters, int N);
static void pythonInitCluster(CLUSTER *curCluster, int vectorDim);

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
    double** T, vectorsMatrix;
    
    /* Parses arguments from python */
    if (!PyArg_ParseTuple(args, "is", &k, &filename)) {
        return NULL;
    }
    /* create matrix of vectors from file data points */
    N = getVectorCount(filename);
    vectorDim = getVectorDim(filename);
    vectorsMatrix = getVectorsMatrix(filename, N, vectorDim);

    T = getNormalizedKEigenvectorsMatrix(k, vectorsMatrix, N, vectorDim);

    return sendMatrixToPython(T, N, k);
}

/*
* Funcion: static PyObject *sendMatrixToPython(double **matrix, int n, int m)
* -----------------------------------------------------------------------------
* Params: Matrix and it's dimensions n*m
* Action: Creates a PyObject matrix from input matrix
* Return: PyObject matrix
*/
static PyObject *sendMatrixToPython(double **matrix, int n, int m) {
    PyObject *pyMatrix, *pyList, *pyVal;
    int i; int j;

    /* parse C matrix to python matrix */
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
* Funcion: static PyObject *runKmeansFromCProgram(PyObject *self, 
*          PyObject *args)
* -----------------------------------------------------------------------------
* Params: Vectors, Centroids, N and k from python kmeans++ program
* Action: Parses arguments from python and run kmeans from C program
* Return: None
*/
static PyObject *runKmeansFromCProgram(PyObject *self, PyObject *args) {
    PyObject *pythonVectors, *pythonCentroids;
    CLUSTER *clusters;
    int k, N;
    double **vectorsMatrix;

    /* Parses arguments from python */
    if (!PyArg_ParseTuple(args, "00ii", &pythonVectors, &pythonClusters, &N, 
    &k)) {
        return NULL;
    }
    /* parse N, data points matrix and cluster from python to c */
    vectorsMatrix = getVectorsFromPython(pythonVectors, N, k);
    clusters = pythonInitializeClusters(pythonClusters, k);
    
    kmeans(clusters, vectorsMatrix, k, N);

    Py_RETURN_NONE;
}

/*
* Funcion: static double **getVectorsFromPython(PyObject *pythonVectorsMatrix, 
*          int N, int vectorDim)
* -----------------------------------------------------------------------------
* Params: Python matrix and it's dimension
* Action: Creates C vectors matrix based on Python matrix
* Return: C vectors matrix
*/
static double **getVectorsFromPython(PyObject *pythonVectorsMatrix, int N, 
                                    int vectorDim) {
    double **vectorsMatrix;
    int i, j;
    PyObject *pyVector;

    vectorsMatrix = createMatrix(N, vectorDim);
    /* fill matrix with python matrix values */
    for (i = 0; i < N; i++) {
        pyVector = PyList_GetItem(pythonVectorsMatrix, i);
        for (j = 0; j < vectorDim; j++) {
            vectorsMatrix[i][j] = PyFloat_AsDouble(PyList_GetItem(pyVector, j));
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
static CLUSTER *pythonInitializeClusters(PyObject *pyClusters, int k) {
    /* This function gets a python object containing the clusters and the number
    of features. The function then copies the clusters into a struct array,
    CLUSTERS, initializes it and returns the array. */

    CLUSTER *clusters;
    double *curr_vector;
    int has_converged = 0;
    int cnt = 0;
    int i = 0;
    int j = 0;
    int curr = 0;
      
    clusters = (CLUSTER *)calloc(k, sizeof(CLUSTER));
    validateAction(clusters == NULL);

    curr = 0;
    for (i = 0; i < k; i++) {
        clusters[i].centroid = (double *)calloc(k, sizeof(double));
        validateAction(clusters[i].centroid == NULL);

        for (j = 0; j < k; j++) {
            clusters[i].centroid[j] = PyFloat_AsDouble(PyList_GetItem(centroids_py, curr));
            curr++;
        }

        clusters[i].vectors_count = 0;
        clusters[i].vectors_sum = (double *)calloc(vectorDim, sizeof(double));
        validateAction(clusters[i].vectors_sum == NULL);

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


/*
* Funcion: 
* -----------------------------------------------------------------------------
* Params: 
* Action: 
* Return: 
*/
static PyObject *kmeans(int k, int max_iter_py, double epsilon_py, int dim_py, int N_py, PyObject *centroids_py, PyObject *vectors_py) {
    int N = N_py;
    int max_iter = max_iter_py;
    double epsilon = epsilon_py;
    int dim = dim_py;
    Cluster *clusters;
    double *curr_vector;
    int has_converged = 0;
    int cnt = 0;
    int i = 0;
    int j = 0;
    int curr = 0;

    /*convert k centroids from python to C*/
    

    /*main loop*/
    cnt = 0;
    while ((cnt < max_iter) && (!has_converged))
    {
        curr = 0;

        /*find current vector cluster*/
        for (i = 0; i < N; i++)
        {
            curr_vector = (double*)calloc(dim, sizeof(double));
            if (curr_vector == NULL)
            {
                printf("An Error Has Occurred\n");
                exit(1);
            }
            for (j = 0; j < dim; j++)
            {
                curr_vector[j] = PyFloat_AsDouble(PyList_GetItem(vectors_py, curr));
                curr++;
            }
            calcCluster(curr_vector, clusters, k, dim);
            free(curr_vector);
        }
        
        /*update centroids*/
        has_converged = updateCentroids(clusters, k, dim, epsilon);
        
        /*reset*/
        for(i = 0; i < k; i++)
        {
            clusters[i].vectors_count = 0;
            for(j = 0; j < dim; j++)
            {
                clusters[i].vectors_sum[j] = 0;
            }
        }
        cnt++;
    }

    return cToPyObject(clusters, k, dim, N);
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
