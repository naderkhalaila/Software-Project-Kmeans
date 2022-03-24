#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>


void TheWeightedAdjacencyMatrix(int N ,int dimension ,double **matrix , double **DataPoints){
    int i,j,k;
    double norm;
    for (i = 0; i < N; i++) {
        matrix[i][i] = 0;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++){
            if(j == i){
                continue;
            }

            norm = 0;
            for(k=0 ; k< dimension ; k++){
                norm += pow((DataPoints[i][k] - DataPoints[j][k]), 2);
            }
            matrix[i][j] = exp(-(sqrt(norm)/2));
        }
    }
}


void TheDiagonalDegreeMatrix(int N , double **matrix , double **WeightedAdjacencyMatrix){

    /* Matrix must be all zeros */

    int i , j ;
    double value;
    for (i = 0; i < N; i++) {
        value =0;
        for (j = 0; j < N; j++){
            value += WeightedAdjacencyMatrix[i][j];
            value = 1/ sqrt(value);
        }
        matrix[i][i] = value;
    }
}


void MatrixMultiplication (int N ,double **matrix , double **matrix1 , double **matrix2){
    int i , j , k ;
    double value;
    for (i=0 ; i < N ; i++){
        for (j=0 ; j<N ; j++){
            value = 0;
            for (k=0; k<N ; k++){
                value+= matrix1[i][k]  * matrix2[k][j];
            }
            matrix[i][j] = value;
        }
    }
}


void TheNormalizedGraphLaplacian (int N , double **matrix ,double **DiagonalDegreeMatrix ,
                                  double **WeightedAdjacencyMatrix){
    int i , j ;
    double Identity[N][N];
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++){
            if(i==j){
                matrix[i][j] = 1;
            }
            else{
                matrix[i][j]=0;
            }
        }
    }

    double **matrix1;
    matrix1 = malloc(sizeof(double *) * N);
    for (i = 0; i<N; i++) {
        matrix1[i] = (double *) malloc(sizeof(double *) * N);
    }
    MatrixMultiplication(N , matrix1, DiagonalDegreeMatrix , WeightedAdjacencyMatrix);
    MatrixMultiplication(N , matrix, matrix1 , DiagonalDegreeMatrix);

    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            matrix[i][j] = Identity[i][j] - matrix[i][j];
        }
    }

    for (i=0 ; i<N ; i++){
        free(matrix1[i]);
    }
    free(matrix1);
}



void CreatePmatrix(int N , double **matrix, double **Amatrix) {
    double max ;
    int i,j;
    int ii , jj;
    int sign;
    double tita, s, c, t;
    max =Amatrix[0][0];
    ii =0;
    jj =0;
    for ( i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if(i!=j && Amatrix[i][j]>max){
                max= Amatrix[i][j];
                ii = i;
                jj = j;
            }
        }
    }

    tita = (Amatrix[jj][jj] - Amatrix[ii][ii])/(2*max);
    if(tita>=0){
        sign = 1;
    }else{
        sign = -1;
    }


    t = (sign)/((sign*(tita)) + (sqrt(pow(tita,2) +1)));
    c = 1/sqrt(pow(t,2) +1);
    s = t*c;

    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++){
            if(i==j){
                matrix[i][j] = 1;
            }
            else{
                matrix[i][j]=0;
            }
        }
    }

    matrix[ii][ii] = c;
    matrix[jj][jj] = c;

    /*check*/

    matrix[ii][jj] = s;
    matrix[jj][ii] = -s;
}

void Ptrans(int N, double **matrix , double **P){

    int i ,j;
    for ( i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {

            matrix[i][j]= P[j][i];
        }
    }

}

double convergence(int N, double **matrix1 , double **matrix2){
    int i , j;
    double offMatrix1 , offMatrix2;
    offMatrix1 = 0;
    offMatrix2 = 0;
    for (i =0 ; i<N ; i++){
        for (j=0 ; j<N ; j++){
            if (i!=j){
                offMatrix1 += pow(matrix1[i][j],2);
                offMatrix2 += pow(matrix2[i][j],2);
            }
        }
    }

    return offMatrix1 - offMatrix2;
}


void Jacobi(int N, double **matrix , double **Vectors ,double **Lmatrix){
    double eps = 1.0 * exp(-15);
    int Max_Iter = 100;
    double conv ;
    int iter =0;
    int i , j;

    double **matrixP ,**TmatrixP , **A , **A_tag , **temp;

    matrixP = malloc(sizeof(double *) * N);
    TmatrixP = malloc(sizeof(double *) * N);
    for (i = 0; i<N; i++) {
        matrixP[i] = (double *) malloc((N) * sizeof(double));
        TmatrixP[i] = (double *) malloc((N) * sizeof(double));
    }

    CreatePmatrix(N ,matrixP , Lmatrix);
    Ptrans(N , TmatrixP , matrixP);

    A = malloc(sizeof(double *) * N);
    A_tag = malloc(sizeof(double *) * N);
    temp = malloc(sizeof(double *) * N);
    for (i = 0; i<N; i++) {
        A[i] = (double *) malloc((N) * sizeof(double));
        A_tag[i] = (double *) malloc((N) * sizeof(double));
        temp[i] = (double *) malloc((N) * sizeof(double));
    }

    for( i=0 ; i<N ; i++) {
        for ( j = 0; j < N; j++) {
            A[i][j] = Lmatrix[i][j];
        }
    }

    MatrixMultiplication(N,temp , TmatrixP , A);
    MatrixMultiplication(N ,A_tag , temp , matrixP);

    for( i=0 ; i<N ; i++) {
        for ( j = 0; j < N; j++) {
            Vectors[i][j] = matrixP[i][j];
        }
    }

    conv = convergence(N ,A, A_tag);
    if(conv <= eps){
        for( i=0 ; i<N ; i++) {
            for ( j = 0; j < N; j++) {
                matrix[i][j] = A_tag[i][j];
            }
        }
    }

    while(conv > eps && iter < Max_Iter){

        double **V;
        V = malloc(sizeof(double *) * N);
        for (i = 0; i<N; i++) {
            V[i] = (double *) malloc((N) * sizeof(double));
        }

        for( i=0 ; i<N ; i++) {
            for ( j = 0; j < N; j++) {
                A[i][j] = A_tag[i][j];
            }
        }

        CreatePmatrix(N,matrixP, A);
        Ptrans(N ,TmatrixP , matrixP);

        MatrixMultiplication(N ,V , Vectors , matrixP);
        for( i=0 ; i<N ; i++) {
            for ( j = 0; j < N; j++) {
                Vectors[i][j] = V[i][j];
            }
        }

        MatrixMultiplication(N, temp , TmatrixP , A);
        MatrixMultiplication(N, A_tag , temp , matrixP);

        conv = convergence(N, A, A_tag);
        iter++;
    }

    for( i=0 ; i<N ; i++) {
        for ( j = 0; j < N; j++) {
            matrix[i][j] = A_tag[i][j];
        }
    }
}

void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void selectionSort(int N , double *arr)
{
    int i, j, min_idx;

    // One by one move boundary of unsorted subarray
    for (i = 0; i < N-1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i+1; j < N; j++)
            if (arr[j] < arr[min_idx])
                min_idx = j;

        // Swap the found minimum element with the first element
        swap(&arr[min_idx], &arr[i]);
    }
}

int Eigengap(int N ,double *eigenvalues){
    int k =0;
    int i;
    double arr[N];
    double max = fabs(eigenvalues[0] - eigenvalues[1]);
    arr[0]= max;

    selectionSort(N, eigenvalues);

    for( i=1 ; i< floor(N/2) ; i++){
        arr[i] = (eigenvalues[i] - eigenvalues[i+1]);
        if(arr[i]>max){
            max=arr[i];
            k=i;
        }
    }

    return k;
}

int NormalizedSpectralClustering(int N ,int K , int dimension , double**DataPoints ,  double** t){
    double ** WeightedAdjacencyMatrix ,**DiagonalDegreeMatrix , **NormalizedGraphLaplacian , **eigenvectors  ,** eigenvalues , **U;
    double *eign;
    int i , j ,  k;

    eign = malloc(sizeof(double *) * N);

    eigenvalues =  malloc(sizeof(double *) * N);
    eigenvectors = malloc(sizeof(double *) * N);
    NormalizedGraphLaplacian = malloc(sizeof(double *) * N);
    DiagonalDegreeMatrix = malloc(sizeof(double *) * N);
    WeightedAdjacencyMatrix = malloc(sizeof(double *) * N);
    for (i = 0; i<N; i++) {
        NormalizedGraphLaplacian[i] = (double *)malloc(sizeof(double *) * N);
        eigenvectors[i] = (double *)malloc(sizeof(double *) * N);
        eigenvalues[i] = (double *)malloc(sizeof(double *) * N);
        DiagonalDegreeMatrix[i] = (double *) malloc((N) * sizeof(double));
        WeightedAdjacencyMatrix[i] = (double *) malloc((N) * sizeof(double));
    }

    TheWeightedAdjacencyMatrix(N, dimension , WeightedAdjacencyMatrix , DataPoints);
    TheDiagonalDegreeMatrix(N , DiagonalDegreeMatrix , WeightedAdjacencyMatrix);
    TheNormalizedGraphLaplacian(N , NormalizedGraphLaplacian ,DiagonalDegreeMatrix , WeightedAdjacencyMatrix );
    Jacobi(N , eigenvalues , eigenvectors , NormalizedGraphLaplacian);

    if(K != 0){
        k=K;
    }

    if(K == 0) {
        for (i = 0; i < N; i++) {
            eign[i] = eigenvalues[i][i];
        }

        k = Eigengap(N, eign);
    }

    U = malloc(sizeof(double *) * k);
    for (i = 0; i<k; i++) {
        U[i] = malloc(sizeof(double *) * N);
    }



    for (i = 0; i<k; i++) {
        for (j = 0 ; j<N ; j++){
            U[i][j] = eigenvectors[i][j];
        }
    }

    int sum;
    for(i = 0 ; i<N ; i++){
        sum = 0;
        for(j = 0 ; j<k ; j++){
            sum += pow(U[i][j],2);
        }
        sum = pow(sum,0.5);
        for(j = 0 ; j<k ; j++){
            t[i][j] = U[i][j]/sum;
        }
    }
    return k;

}

double dist(double *point1, double *point2, int dimension) {

    int i;
    double sum = 0;
    for (i = 0; i < dimension; i++) {
        sum += pow((point1[i]-point2[i]), 2);
    }

    sum = pow(sum, 0.5);
    return sum;
}

int mindist(int k, int dimension, double **centroids_list, double *point) {

    double min = dist(centroids_list[0], point, dimension);
    int index = 0;
    int j;
    for (j = 1; j < k; j++) {
        double distance = dist(centroids_list[j], point, dimension);
        if (distance < min) {
            min = distance;
            index = j;
        }
    }
    return index;
}

void Init(int k, int dimension, int *count_array, double **sum_array) {
    int i, j;
    for (i = 0; i < k; i++) {
        count_array[i] = 0;
        for (j = 0; j < dimension; j++) {
            sum_array[i][j] = 0;
        }
    }
}

void clustering(int k, int rows, int dimension, int *count_array,
                double **sum_array, double **DataPoints, double **centroids_list) {
    int i;
    int j;
    int index;
    for (i = 0; i < rows; i++) {
        index = mindist(k, dimension, centroids_list, DataPoints[i]);
        count_array[index]++;
        for (j = 0; j < dimension; j++) {
            sum_array[index][j] += DataPoints[i][j];
        }
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < dimension; j++) {
            sum_array[i][j] = sum_array[i][j] / count_array[i];
        }
    }
}

int calculate_delta(int k, int dimension, double **centroids_list, double **sum_array) {
    double distance;
    int i, j;
    int bol = 1;
    for (i = 0; i < k; i++) {
        distance = dist(centroids_list[i], sum_array[i], dimension);
        if (distance > 0.01) {
            bol = 0;
            break;
        }
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < dimension; j++) {
            centroids_list[i][j] = sum_array[i][j];
        }
    }
    return bol;
}


void Sizefile(char *filename, int *dimension, int *rows) {
    int tmp_dimension;
    int tmp_rows;
    FILE *f;
    char c;
    f = fopen(filename, "rt");
    tmp_dimension=0;
    tmp_rows=0;


    while ((c = fgetc(f)) != EOF) {
        if (c == ',') {
            tmp_dimension++;
        }
        if (c == '\n') {
            tmp_rows++;
        }
    }

    tmp_dimension = tmp_dimension / tmp_rows + 1;
    fclose(f);
    *dimension = tmp_dimension;
    *rows = tmp_rows;
}


void Getpoints(char filename[] , int rows , int dimension , double **DataPoints){
    int negative;
    int p = 1;
    int i, j;
    double num, fr ,tmp;
    FILE *file;
    char C;

    file = fopen(filename, "rt");
    for (i = 0; i < rows; i++) {
        for (j = 0; j < dimension; j++) {
            tmp = 0;
            negative = 0;
            while ((C = fgetc(file)) != '.' && C != EOF) {
                num = 1;
                if (C == '-') {
                    negative = 1;
                } else {
                    tmp = tmp*10;
                    num = num * ((double) C - 48);
                    tmp += num;
                    num = tmp;
                }

            }
            fr = 0;
            p = 1;
            while ((C = fgetc(file)) != EOF && C != ',' && C != '\n') {
                fr = fr + ((double) C - 48) * pow(10, -1 * p);
                p++;
            }
            num = num + fr;
            if (negative == 1) {
                num = -1 * num;
            }
            DataPoints[i][j] = num;
        }
    }
    fclose(file);

}


void print(double **matrix, int rows, int col) {
    int i;
    int j;
    for(i=0;i<rows;i++){
        for ( j = 0; j < col; j++) {
            printf("%.3f",matrix[i][j]);
            if(j<col-1){
                printf(",");
            }
        }
        printf("\n");
    }
}

int kmeans(char filename[], int Goal) {
    int dimension = 0;
    int rows = 0;
    int Temp =0;
    int i, j ,column;
    double **DataPoints;
    double** WeightedAdjacencyMatrix, **DiagonalDegreeMatrix, **NormalizedGraphLaplacian ,**jacobi ,**Vectors  , **matrix;


    Sizefile(filename, &dimension, &rows);
    DataPoints = malloc(sizeof(double *) * rows);
    for (i = 0; i < rows; i++) {
        DataPoints[i] = (double *) malloc((dimension) * sizeof(double));
    }
    Getpoints(filename , rows , dimension , DataPoints);


    if (Temp != Goal){
        Temp ++;
        WeightedAdjacencyMatrix = malloc(sizeof(double *) * rows);
        for (i = 0; i<rows; i++) {
            WeightedAdjacencyMatrix[i] = (double *) malloc(sizeof(double *) * rows);
        }
        TheWeightedAdjacencyMatrix(rows , dimension , WeightedAdjacencyMatrix , DataPoints);

        if(Temp == Goal){
            column = rows;
            print(WeightedAdjacencyMatrix , rows , column );
            return 1;
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
            column = rows;
            print(DiagonalDegreeMatrix , rows , column);
            return 1;
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
            column = rows;
            print(NormalizedGraphLaplacian , rows , column);
            return 1;
        }
    }

    if (Temp != Goal) {
        Temp++;
        jacobi = malloc(sizeof(double *) * rows);
        Vectors = malloc(sizeof(double *) * rows);
        for (i = 0; i < rows; i++) {
            jacobi[i] = (double *) malloc(sizeof(double *) * rows);
            Vectors[i] = (double *) malloc(sizeof(double *) * rows);
        }

        matrix = malloc(sizeof(double *) * rows + 1);
        for (i = 0; i < rows + 1; i++) {
            matrix[i] = (double *) malloc(sizeof(double *) * rows);
        }

        Jacobi(rows, jacobi, Vectors, NormalizedGraphLaplacian);

        for (i = 0; i < rows; i++) {
            matrix[0][i] = jacobi[i][i];
        }
        for (i = 1; i < rows; i++) {
            for (j = 0; j < rows; j++) {
                matrix[i][j] = Vectors[i - 1][j];
            }
        }
        column = rows;
        print(matrix , rows+1 , column);
        return 1;
    }

    return 0;

}
int main() {

    char path[] = "C:\\Users\\weamm\\Documents\\example.txt";
    kmeans(path, 1);
    return 0;
}
