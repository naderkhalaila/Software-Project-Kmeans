//
// Created by nadsa on 12/04/2022.
//

#ifndef FINALPROJECT_SPKMEANS_H
#define FINALPROJECT_SPKMEANS_H

void print(double **matrix, int rows, int col);

void check(double ** Data);

void checkSub(double * Data);

void TheWeightedAdjacencyMatrix(int N ,int dimension ,double **matrix , double **DataPoints);

void TheDiagonalDegreeMatrix(int N , double **matrix , double **WeightedAdjacencyMatrix);

void MatrixMultiplication (int N ,double **matrix , double **matrix1 , double **matrix2);

void TheNormalizedGraphLaplacian (int N , double **matrix ,double **DiagonalDegreeMatrix ,
                                  double **WeightedAdjacencyMatrix);

void CreatePmatrix(int N , double **matrix, double **Amatrix);

double absDouble(double d);

void Ptrans(int N, double **matrix , double **P);

double convergence(int N, double **matrix1 , double **matrix2);

void Jacobi(int N, double **matrix , double **Vectors ,double **Lmatrix);

void LnormSort(int n, double * arr , int *index);

void insertionSort(int n, double * arr);

int Eigengap(int N ,double *eigenvalues);

int NormalizedSpectralClustering(int N ,int K , int dimension , double**DataPoints ,  double** t);

double dist(double *point1, double *point2, int dimension);

int mindist(int k, int dimension, double **centroids_list, double *point);

void Init(int k, int dimension, int *count_array, double **sum_array);

void clustering(int k, int rows, int dimension, int *count_array,
                double **sum_array, double **DataPoints, double **centroids_list);

int calculate_delta(int k, int dimension, double **centroids_list, double **sum_array);

void Sizefile(char *filename, int *dimension, int *rows);

void Getpoints(char filename[] , int rows , int dimension , double **DataPoints);

int kmeans(char filename[], int Goal);

int main(int argc, char** argv);



#endif //FINALPROJECT_SPKMEANS_H
