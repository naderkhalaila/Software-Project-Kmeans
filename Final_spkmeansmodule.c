#define PY_SSIZE_T_CLEANS
#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include <malloc.h>
#include "spkmeans.h"

static PyObject* kmeans_capi(PyObject *self, PyObject *args);
static double **parse_arrays(PyObject* _list, int num_row, int num_col);
static PyObject *pseudo_main(PyObject* Py_data, PyObject* Py_centroids, int N, int d, int K, int MAX_ITER, int goal , int Part); //static?


/* the wrapping function for the pseudo_main - parses PyObjects */
static PyObject* kmeans_capi(PyObject *self, PyObject *args){

    PyObject *data, *centroids;
    int rows, dim, K, MAX_ITER;
    int goal , Part;
    if(!PyArg_ParseTuple(args, "OOiiiiii", &data, &centroids, &rows, &dim, &K, &MAX_ITER , &goal , &Part)){
        return NULL;
    }
    return Py_BuildValue("O", pseudo_main(data, centroids, rows, dim, K, MAX_ITER, goal , Part));
}

/* functino that parses the data and puts them in arrays */
static double **parse_arrays(PyObject* _list, int num_row, int num_col) {

    int i, j;
    Py_ssize_t Py_i, Py_j;
    double **parsed_data;
    parsed_data = malloc(num_row * sizeof(double*));
    assert(parsed_data!=NULL);
    PyObject* item; PyObject* num;
    for (i = 0; i < num_row; i++) {
        Py_i = (Py_ssize_t)i;
        parsed_data[Py_i] = malloc(num_col * sizeof(double));
        assert(parsed_data[Py_i]!=NULL);
        item = PyList_GetItem(_list, Py_i);
        if (!PyList_Check(item)){ /* Skips non-lists */
            continue;
        }
        for (j = 0; j < num_col; j++) {
            Py_j = (Py_ssize_t)j;
            num = PyList_GetItem(item, Py_j);
            if (!PyFloat_Check(num)) continue; /* Skips non-floats */
            parsed_data[Py_i][Py_j] = PyFloat_AsDouble(num);
        }
    }return parsed_data;
}

/* this array tells python what methods this module has */
static PyMethodDef capiMethods[] = {

        {"fit",
                (PyCFunction) kmeans_capi,
                     METH_VARARGS,
                        PyDoc_STR("calculates the centroids using kmeans algorithm")},
        {NULL, NULL, 0, NULL}
};

/* This struct initiates the module using the above definition. */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",           //name of module
        NULL,                   //module documentation
        -1,                     //size of per-interpreter
        capiMethods
};

/* Module Creation */
PyMODINIT_FUNC
PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


static PyObject *pseudo_main(PyObject* Py_data, PyObject* Py_centroids, int num_rows, int dim, int K, int MAX_ITER,
                             int goal , int Part){

    PyObject *array, *vec, *num;
    PyObject *lst_centroids ;
    int dimension = dim;
    int rows = num_rows;
    int Goal = goal;
    int i, j;
    double ** centroids_list, **sum_array, **DataPoints , **t , **kmatrix;
    int *count_array;
    int delta = 0;
    int iter;
    int Temp = 1;
    double** WeightedAdjacencyMatrix, **DiagonalDegreeMatrix, **NormalizedGraphLaplacian ,**jacobi ,**Vectors  , **matrix;
    int k = K;

    /* initialising centroids, data, cSum, cNum, data_index*/
    DataPoints = parse_arrays(Py_data, rows, dimension);

    if (Part == 0){
        if (Goal == 1) {

            t = malloc(sizeof(double *) * rows);

            if (K == 0) {
                k = NormalizedSpectralClustering(rows, K, dimension, DataPoints, t);
                kmatrix = malloc(sizeof(double *) * 1);
                for (i = 0; i < 1; i++) {
                    kmatrix[i] = (double *) malloc(sizeof(double *) * 1);
                }
                kmatrix[0][0] = k;

                array = PyList_New(1);
                if (!array){
                    return NULL;
                }
                for(i=0; i<1; i++){
                    vec = PyList_New(1);
                    if (!vec){
                        return NULL;
                    }
                    for (j = 0; j < 1; j++) {
                        num = PyFloat_FromDouble(kmatrix[i][j]);
                        if (!num) {
                            Py_DECREF(vec);
                            return NULL;
                        }PyList_SET_ITEM(vec, j, num);
                    }PyList_SET_ITEM(array, i, vec);
                }

                return array;
            }

            if (K != 0) {
                for (i = 0; i < rows; i++) {
                    t[i] = (double *) malloc(sizeof(double *) * k);
                }

                k = NormalizedSpectralClustering(rows, k, dimension, DataPoints, t);

                array = PyList_New(rows);
                if (!array){
                    return NULL;
                }
                for(i=0; i<rows; i++){
                    vec = PyList_New(k);
                    if (!vec){
                        return NULL;
                    }
                    for (j = 0; j < k; j++) {
                        num = PyFloat_FromDouble(t[i][j]);
                        if (!num) {
                            Py_DECREF(vec);
                            return NULL;
                        }PyList_SET_ITEM(vec, j, num);
                    }PyList_SET_ITEM(array, i, vec);
                }
                return array;
            }
        }
        if (Temp != Goal){
            Temp ++;
            WeightedAdjacencyMatrix = malloc(sizeof(double *) * rows);
            for (i = 0; i<rows; i++) {
                WeightedAdjacencyMatrix[i] = (double *) malloc(sizeof(double *) * rows);
            }
            TheWeightedAdjacencyMatrix(rows , dimension , WeightedAdjacencyMatrix , DataPoints);
            if(Temp == Goal){

                array = PyList_New(rows);
                if (!array){
                    return NULL;
                }
                for(i=0; i<rows; i++){
                    vec = PyList_New(rows);
                    if (!vec){
                        return NULL;
                    }
                    for (j = 0; j < rows; j++) {
                        num = PyFloat_FromDouble(WeightedAdjacencyMatrix[i][j]);
                        if (!num) {
                            Py_DECREF(vec);
                            return NULL;
                        }PyList_SET_ITEM(vec, j, num);
                    }PyList_SET_ITEM(array, i, vec);
                }

                return array;
            }
        }
        if (Temp != Goal){
            Temp ++;
            DiagonalDegreeMatrix = malloc(sizeof(double *) * rows);
            for (i = 0; i<rows; i++) {
                DiagonalDegreeMatrix[i] = (double *) malloc(sizeof(double *) * rows);
            }
            TheDiagonalDegreeMatrix(rows, DiagonalDegreeMatrix , WeightedAdjacencyMatrix);
            if(Temp == Goal){

                array = PyList_New(rows);
                if (!array){
                    return NULL;
                }
                for(i=0; i<rows; i++){
                    vec = PyList_New(rows);
                    if (!vec){
                        return NULL;
                    }
                    for (j = 0; j < rows; j++) {
                        num = PyFloat_FromDouble(DiagonalDegreeMatrix[i][j]);
                        if (!num) {
                            Py_DECREF(vec);
                            return NULL;
                        }PyList_SET_ITEM(vec, j, num);
                    }PyList_SET_ITEM(array, i, vec);
                }

                return array;
            }
        }
        if (Temp != Goal){
            Temp ++;
            NormalizedGraphLaplacian = malloc(sizeof(double *) * rows);
            for (i = 0; i<rows; i++) {
                NormalizedGraphLaplacian[i] = (double *) malloc(sizeof(double *) * rows);
            }
            TheNormalizedGraphLaplacian(rows , NormalizedGraphLaplacian ,DiagonalDegreeMatrix , WeightedAdjacencyMatrix);
            if(Temp == Goal){

                array = PyList_New(rows);
                if (!array){
                    return NULL;
                }
                for(i=0; i<rows; i++){
                    vec = PyList_New(rows);
                    if (!vec){
                        return NULL;
                    }
                    for (j = 0; j < rows; j++) {
                        num = PyFloat_FromDouble(NormalizedGraphLaplacian[i][j]);
                        if (!num) {
                            Py_DECREF(vec);
                            return NULL;
                        }PyList_SET_ITEM(vec, j, num);
                    }PyList_SET_ITEM(array, i, vec);
                }

                return array;
            }
        }
        if (Temp != Goal){
            Temp ++;

            jacobi = malloc(sizeof(double *) * rows);
            Vectors = malloc(sizeof(double *) * rows);
            for (i = 0; i<rows; i++) {
                jacobi[i] = (double *) malloc(sizeof(double *) * rows);
                Vectors[i] = (double *) malloc(sizeof(double *) * rows);
            }

            matrix = malloc(sizeof(double *) * rows+1);
            for (i = 0; i<rows+1; i++) {
                matrix[i] = (double *) malloc(sizeof(double *) * rows);
            }

            Jacobi(rows, jacobi , Vectors ,NormalizedGraphLaplacian);

            for (i = 0 ; i < rows+1 ; i++){
                matrix[0][i] = jacobi[i][i];
            }
            for (i =1 ; i < rows ; i++){
                for(j = 0 ; j < rows ; j++){
                    matrix[i][j] = Vectors[i-1][j];
                }
            }

            array = PyList_New(rows+1);
            if (!array){
                return NULL;
            }
            for(i=0; i<rows+1; i++){
                vec = PyList_New(rows);
                if (!vec){
                    return NULL;
                }
                for (j = 0; j < rows; j++) {
                    num = PyFloat_FromDouble(matrix[i][j]);
                    if (!num) {
                        Py_DECREF(vec);
                        return NULL;
                    }PyList_SET_ITEM(vec, j, num);
                }PyList_SET_ITEM(lst_centroids, i, vec);
            }

            return array;
        }
    }
    else{
        centroids_list = parse_arrays(Py_centroids, k, dimension);
        sum_array = malloc(sizeof(double *) * k);
        for (i = 0; i<k; i++) {
            sum_array[i] = (double *) malloc((dimension) * sizeof(double));
        }

        count_array = (int *) malloc(sizeof(int) * K);

        iter = 0;
        while(delta==0 && iter < MAX_ITER) {
            Init(k, dimension, count_array, sum_array);
            clustering(k, rows, dimension, count_array, sum_array, DataPoints, centroids_list);
            delta = calculate_delta(k, dimension, centroids_list, sum_array);
            iter++;
        }

        lst_centroids = PyList_New(k);
        if (!lst_centroids){
            return NULL;
        }
        for(i=0; i<k; i++){
            vec = PyList_New(dimension);
            if (!vec){
                return NULL;
            }
            for (j = 0; j < dimension; j++) {
                num = PyFloat_FromDouble(centroids_list[i][j]);
                if (!num) {
                    Py_DECREF(vec);
                    return NULL;
                }PyList_SET_ITEM(vec, j, num);
            }PyList_SET_ITEM(lst_centroids, i, vec);
        }

        for(i=0 ; i<K ; i++){
            free(centroids_list[i]);
            free((sum_array[i]));
        }
        free(centroids_list);
        for(i=0 ; i<rows ; i++){
            free(DataPoints[i]);
        }
        free(DataPoints);
        free(sum_array);
        free(count_array);

        return lst_centroids;
    }
}

