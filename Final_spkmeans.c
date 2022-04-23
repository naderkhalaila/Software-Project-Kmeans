#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PY_SSIZE_T_CLEANS
#include <malloc.h>
#include <string.h>

double absDouble(double d);

void check(double ** Data){
    if (Data == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
}

void checkSub(double * Data){
    if (Data == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
}

void print(double **matrix, int rows, int col) {
    int i;
    int j;
    for(i=0;i<rows;i++){
        for ( j = 0; j < col; j++) {
            printf("%.4f",matrix[i][j]);
            if(j<col-1){
                printf(",");
            }
        }
        printf("\n");
    }
}

void TheWeightedAdjacencyMatrix(int N ,int dimension ,double **matrix , double **DataPoints){
    int i,j,k;
    double norm;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++){
            if(j == i){
                matrix[i][i] = 0;
                continue;
            }

            norm = 0;
            for(k=0 ; k< dimension ; k++){
                norm += pow((DataPoints[i][k] - DataPoints[j][k]), 2);
            }
            matrix[i][j] =exp(-(sqrt(norm)/2));
        }
    }
}

void TheDiagonalDegreeMatrix(int N , double **matrix , double **WeightedAdjacencyMatrix){
    int i , j ;
    double value , d;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            matrix[i][j] = 0;
        }
    }
    for (i = 0; i < N; i++) {
        value =0;
        for (j = 0; j < N; j++){
            value += WeightedAdjacencyMatrix[i][j];
        }
        d = 1/ sqrt(value);
        matrix[i][i] = d;
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
    double **matrix1 , **Identity;

    Identity = malloc(sizeof(double *) * N);
    check(Identity);
    for (i = 0; i<N; i++) {
        Identity[i] = (double *) malloc(sizeof(double *) * N);
        checkSub(Identity[i]);
    }
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++){
            if(i==j){
                Identity[i][j] = 1;
            }
            else{
                Identity[i][j]=0;
            }
        }
    }


    matrix1 = malloc(sizeof(double *) * N);
    check(matrix1);
    for (i = 0; i<N; i++) {
        matrix1[i] = (double *) malloc(sizeof(double *) * N);
        checkSub(matrix1[i]);
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
        free(Identity[i]);
    }
    free(matrix1);
    free(Identity);
}

void CreatePmatrix(int N , double **matrix, double **Amatrix) {
    double max ;
    int i,j, ii , jj ,sign;
    double tita, s, c, t;
    max = absDouble(Amatrix[0][1]);
    ii =0;
    jj =1;
    for ( i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if(i!=j && absDouble(Amatrix[i][j])>max){
                max= absDouble(Amatrix[i][j]);
                ii = i;
                jj = j;

            }
        }
    }
    if(Amatrix[ii][jj]< 0){
        max = -max;
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

    matrix[ii][jj] = s;
    matrix[jj][ii] = -s;
}

double absDouble(double d) {
    if(d >= 0 ){
        return d;
    }
    else{
        return -d;
    }
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


void Jacobi(int N, double **matrix , double **Vectors ,double **Lmatrix) {
    double eps = 1.0 * pow(10,-5);
    int Max_Iter = 100;
    double conv  , tempdouble;
    int iter = 0;
    int i, j , tempint;

    double **matrixP, **TmatrixP, **A, **A_tag, **temp;

    matrixP = malloc(sizeof(double *) * N);
    TmatrixP = malloc(sizeof(double *) * N);
    check(matrixP);
    check(TmatrixP);
    for (i = 0; i < N; i++) {
        matrixP[i] = (double *) malloc((N) * sizeof(double));
        TmatrixP[i] = (double *) malloc((N) * sizeof(double));
        checkSub(matrixP[i]);
        checkSub(TmatrixP[i]);
    }

    CreatePmatrix(N, matrixP, Lmatrix);
    Ptrans(N, TmatrixP, matrixP);

    A = malloc(sizeof(double *) * N);
    A_tag = malloc(sizeof(double *) * N);
    temp = malloc(sizeof(double *) * N);
    check(A);
    check(A_tag);
    check(temp);
    for (i = 0; i < N; i++) {
        A[i] = (double *) malloc((N) * sizeof(double));
        A_tag[i] = (double *) malloc((N) * sizeof(double));
        temp[i] = (double *) malloc((N) * sizeof(double));
        checkSub(A[i]);
        checkSub(A_tag[i]);
        checkSub(temp[i]);
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = Lmatrix[i][j];
        }
    }

    MatrixMultiplication(N, temp, TmatrixP, A);
    MatrixMultiplication(N, A_tag, temp, matrixP);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Vectors[i][j] = matrixP[i][j];
        }
    }

    conv = convergence(N, A, A_tag);
    if (conv <= eps) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                matrix[i][j] = A_tag[i][j];
            }
        }
    }

    while (conv > eps && iter < Max_Iter) {

        double **V;
        V = malloc(sizeof(double *) * N);
        check(V);
        for (i = 0; i < N; i++) {
            V[i] = (double *) malloc((N) * sizeof(double));
            checkSub(V[i]);
        }

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                A[i][j] = A_tag[i][j];
            }
        }

        CreatePmatrix(N, matrixP, A);
        Ptrans(N, TmatrixP, matrixP);

        MatrixMultiplication(N, V, Vectors, matrixP);
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                Vectors[i][j] = V[i][j];
            }
        }

        for(i=0 ; i <N ; i++){
            free(V[i]);
        }
        free(V);

        MatrixMultiplication(N, temp, TmatrixP, A);
        MatrixMultiplication(N, A_tag, temp, matrixP);

        conv = convergence(N, A, A_tag);
        iter++;
    }

    for(i=0 ; i<N ; i++){
        for(j=0 ; j<N ; j++){
            tempdouble = Vectors[i][j] * 10000;
            tempint = (int)tempdouble;
            tempdouble= ((double)tempint) * 10000 ;
            if(tempdouble == -0.0000)
                Vectors[i][j] = 0.0000;
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            matrix[i][j] = A_tag[i][j];
        }
    }

    for (i = 0; i < N; i++) {
        free(matrixP[i]);
        free(A[i]);
        free(temp[i]);
        free(A_tag[i]);
        free(TmatrixP[i]);
    }

    free(matrixP);
    free(A);
    free(temp);
    free(A_tag);
    free(TmatrixP);
}

void insertionSort(int n, double * arr)
{
    int i, j;
    double key;

    for (i = 1; i < n; i++) {

        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void LnormSort(int n, double * arr , int *index)
{
    int i, j;
    double key;

    for (i = 1; i < n; i++) {

        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            index[j + 1] = index[j];
            j = j - 1;
        }
        arr[j + 1] = key;
        index[j + 1] = i;
    }
}

int Eigengap(int N ,double *eigenvalues){
    int k =0;
    int i;
    double *arr;
    double max = fabs(eigenvalues[0] - eigenvalues[1]);
    arr = malloc(sizeof(double *) * N);
    check(arr);
    arr[0]= max;

    insertionSort(N, eigenvalues);
    for( i=1 ; i< floor(N/2) ; i++){
        arr[i] = fabs(eigenvalues[i] - eigenvalues[i+1]);
        if(arr[i]>=max){
            max=arr[i];
            k=i;
        }
    }
    free(arr);
    return k+1;
}

int NormalizedSpectralClustering(int N ,int K , int dimension , double**DataPoints ,  double** t){
    double ** WeightedAdjacencyMatrix ,**DiagonalDegreeMatrix , **NormalizedGraphLaplacian , **eigenvectors  ,**sortedL,** eigenvalues , **U;
    double *eign;
    int i , j ,  k;
    double sum;
    int * index;

    eign = malloc(sizeof(double *) * N);
    eigenvalues =  malloc(sizeof(double *) * N);
    eigenvectors = malloc(sizeof(double *) * N);
    NormalizedGraphLaplacian = malloc(sizeof(double *) * N);
    DiagonalDegreeMatrix = malloc(sizeof(double *) * N);
    WeightedAdjacencyMatrix = malloc(sizeof(double *) * N);
    check(eign);
    check(eigenvalues);
    check(eigenvectors);
    check(NormalizedGraphLaplacian);
    check(DiagonalDegreeMatrix);
    check(WeightedAdjacencyMatrix);
    for (i = 0; i<N; i++) {
        NormalizedGraphLaplacian[i] = (double *)malloc(sizeof(double *) * N);
        eigenvectors[i] = (double *)malloc(sizeof(double *) * N);
        eigenvalues[i] = (double *)malloc(sizeof(double *) * N);
        DiagonalDegreeMatrix[i] = (double *) malloc((N) * sizeof(double));
        WeightedAdjacencyMatrix[i] = (double *) malloc((N) * sizeof(double));
        checkSub(eigenvalues[i]);
        checkSub(eigenvectors[i]);
        checkSub(NormalizedGraphLaplacian[i]);
        checkSub(DiagonalDegreeMatrix[i]);
        checkSub(WeightedAdjacencyMatrix[i]);
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
        for (i = 0; i < N; i++) {
            t[i] = (double *) malloc(sizeof(double *) * k);
            checkSub(t[i]);
        }
    }

    index =  malloc(sizeof(int *) * N);
    check(index);
    for (i = 0; i < N; i++) {
        eign[i] = eigenvalues[i][i];
        index[i] = i;
    }

    LnormSort(N, eign, index);
    sortedL = malloc(sizeof(double *) * N);
    check(sortedL);
    for (i = 0; i<N; i++) {
        sortedL[i] = malloc(sizeof(double *) * N);
        checkSub(sortedL[i]);
    }
    for (i = 0 ; i < N ; i++){
        for (j = 0 ; j<N ; j++){
            sortedL[i][j] = eigenvectors[i][index[j]] ;
        }
    }

    U = malloc(sizeof(double *) * N);
    check(U);
    for (i = 0; i<N; i++) {
        U[i] = malloc(sizeof(double *) * k);
        checkSub(U[i]);
    }

    for (i = 0; i<N; i++) {
        for (j = 0 ; j<k ; j++){
            U[i][j] = sortedL[i][j];
        }
    }

    for(i = 0 ; i<N ; i++) {
        sum = 0.0;
        for (j = 0; j < k; j++) {
            sum += pow(U[i][j], 2);
        }
        for (j = 0; j < k; j++) {

            if (sum != 0) {
                t[i][j] = U[i][j] / sqrt(sum);
            }
            if (sum == 0) {
                t[i][j] = U[i][j];
            }


        }

    }

    for (i=0; i < N; i++){
        free(eigenvalues[i]);
        free(eigenvectors[i]);
        free(NormalizedGraphLaplacian[i]);
        free(DiagonalDegreeMatrix[i]);
        free(WeightedAdjacencyMatrix[i]);
        free(U[i]);
        free(sortedL[i]);
    }
    free(sortedL);
    free(index);
    free(eigenvalues);
    free(eigenvectors);
    free(NormalizedGraphLaplacian);
    free(DiagonalDegreeMatrix);
    free(WeightedAdjacencyMatrix);
    free(U);
    free(eign);

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

int symetric(double ** Datapoints ,int dimension , int rows){
    int i ,j ;
    if (dimension != rows){
        return 0;
    }
    else{
        for(i = 0 ; i<rows ; i++)
            for (j = 0 ; j<rows ; j++)
                if(Datapoints[i][j] != Datapoints[j][i])
                    return 0;
    }
    return 1;
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
    double  fr ,tmp;
    FILE *file;
    char C;
    double num =0;

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

int kmeans(char filename[], int Goal) {
    int dimension = 0;
    int rows = 0;
    int Temp = 0;
    int i, j, column;
    double **DataPoints;
    double **WeightedAdjacencyMatrix, **DiagonalDegreeMatrix, **NormalizedGraphLaplacian, **jacobi, **Vectors, **matrix;


    Sizefile(filename, &dimension, &rows);
    DataPoints = malloc(sizeof(double *) * rows);
    check(DataPoints);
    for (i = 0; i < rows; i++) {
        DataPoints[i] = (double *) malloc((dimension) * sizeof(double));
        checkSub(DataPoints[i]);
    }
    Getpoints(filename, rows, dimension, DataPoints);

    if(Goal == 4){

        if(symetric(DataPoints ,dimension , rows) == 0){
            printf("Invalid Input!");
            return 0;
        }
        jacobi = malloc(sizeof(double *) * rows);
        Vectors = malloc(sizeof(double *) * rows);
        check(jacobi);
        check(Vectors);
        for (i = 0; i < rows; i++) {
            jacobi[i] = (double *) malloc(sizeof(double *) * rows);
            Vectors[i] = (double *) malloc(sizeof(double *) * rows);
            checkSub(jacobi[i]);
            checkSub(Vectors[i]);
        }

        matrix = malloc(sizeof(double *) * rows + 1);
        check(matrix);
        for (i = 0; i < rows + 1; i++) {
            matrix[i] = (double *) malloc(sizeof(double *) * rows);
            checkSub(matrix[i]);
        }

        Jacobi(rows, jacobi, Vectors, DataPoints);

        for (i = 0; i < rows; i++) {
            matrix[0][i] = jacobi[i][i];
        }
        for (i = 0; i < rows; i++) {
            for (j = 0; j < rows; j++) {
                matrix[i + 1][j] = Vectors[i][j];
            }
        }
        column = rows;
        print(matrix, rows + 1, column);

        for (i=0; i < rows; i++){
            free(jacobi[i]);
            free(Vectors[i]);
            free(matrix[i]);
        }
        free(matrix[rows]);
        free(matrix);
        free(jacobi);
        free(Vectors);

        return 1;
    }
    if (Temp != Goal) {
        Temp++;
        WeightedAdjacencyMatrix = malloc(sizeof(double *) * rows);
        check(WeightedAdjacencyMatrix);
        for (i = 0; i < rows; i++) {
            WeightedAdjacencyMatrix[i] = (double *) malloc(sizeof(double *) * rows);
            checkSub(WeightedAdjacencyMatrix[i]);
        }
        TheWeightedAdjacencyMatrix(rows, dimension, WeightedAdjacencyMatrix, DataPoints);

        for (i=0; i < rows; i++){
            free(DataPoints[i]);
        }
        free(DataPoints);

        if (Temp == Goal) {
            column = rows;
            print(WeightedAdjacencyMatrix, rows, column);

            for (i=0; i < rows; i++){
                free(WeightedAdjacencyMatrix[i]);
            }
            free(WeightedAdjacencyMatrix);

            return 1;
        }
    }
    if (Temp != Goal) {
        Temp++;
        DiagonalDegreeMatrix = malloc(sizeof(double *) * rows);
        check(DiagonalDegreeMatrix);
        for (i = 0; i < rows; i++) {
            DiagonalDegreeMatrix[i] = (double *) malloc(sizeof(double *) * rows);
            checkSub(DiagonalDegreeMatrix[i]);
        }
        TheDiagonalDegreeMatrix(rows, DiagonalDegreeMatrix, WeightedAdjacencyMatrix);
        if (Temp == Goal) {
            column = rows;
            for (i = 0 ; i < rows ; i++){
                DiagonalDegreeMatrix[i][i] = 1/pow(DiagonalDegreeMatrix[i][i],2);
            }
            print(DiagonalDegreeMatrix, rows, column);

            for (i=0; i < rows; i++){
                free(DiagonalDegreeMatrix[i]);
                free(WeightedAdjacencyMatrix[i]);
            }
            free(DiagonalDegreeMatrix);
            free(WeightedAdjacencyMatrix);

            return 1;
        }
    }
    if (Temp != Goal) {
        Temp++;
        NormalizedGraphLaplacian = malloc(sizeof(double *) * rows);
        check(NormalizedGraphLaplacian);
        for (i = 0; i < rows; i++) {
            NormalizedGraphLaplacian[i] = (double *) malloc(sizeof(double *) * rows);
            checkSub(NormalizedGraphLaplacian[i]);
        }
        TheNormalizedGraphLaplacian(rows, NormalizedGraphLaplacian, DiagonalDegreeMatrix, WeightedAdjacencyMatrix);

        for (i = 0; i < rows; i++) {
            free(DiagonalDegreeMatrix[i]);
            free(WeightedAdjacencyMatrix[i]);
        }
        free(DiagonalDegreeMatrix);
        free(WeightedAdjacencyMatrix);

        if (Temp == Goal) {
            column = rows;
            print(NormalizedGraphLaplacian, rows, column);

            for (i = 0; i < rows; i++) {
                free(NormalizedGraphLaplacian[i]);
            }
            free(NormalizedGraphLaplacian);

            return 1;
        }
    }

    return 0;
}

int main(int argc, char** argv) {
    if(argc != 3){
        printf("Invalid Input!");
        return 0;
    }
    if(strcmp(argv[1],"wam")==0){
        return kmeans(argv[2] , 1);
    }
    else{
        if (strcmp(argv[1],"ddg")==0){
            return kmeans(argv[2] , 2);
        }
        else{
            if(strcmp(argv[1],"lnorm")==0){
                return kmeans(argv[2] , 3);
            }
            else{
                if(strcmp(argv[1],"jacobi")==0){
                    return kmeans(argv[2] , 4);
                }
                else{
                    printf("Invalid Input!");
                    return 0;
                }
            }
        }
    }

    return 0;
}


