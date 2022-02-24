#include <stdio.h>
#include <math.h>
#include <stdlib.h>

 void TheWeightedAdjacencyMatrix(int N ,int dimension ,double **matrix , double **DataPoints){
    int i,j,k;
    double norm;
    for (i = 0; i < N; i++) {
        matrix[i][i] = 0;
    }

     for (i = 0; i < N; i++) {
         for (j = i+1; j < N; j++){
             norm = 0;
             for(k=0 ; k< dimension ; k++){
                 norm += pow((DataPoints[i][k] - DataPoints[j][k]),2);
             }
             matrix[i][j] = exp(-(sqrt(norm)/2));
        }
     }

     for (i = 1; i < N; i++) {
         for (j = 0; j < i; j++){
             matrix[i][j] = matrix[j][i];
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
    int i;

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
        for (int j = 0; j < N; j++) {
            A[i][j] = Lmatrix[i][j];
        }
    }

    MatrixMultiplication(N,temp , TmatrixP , A);
    MatrixMultiplication(N ,A_tag , temp , matrixP);

    for( i=0 ; i<N ; i++) {
        for (int j = 0; j < N; j++) {
            Vectors[i][j] = matrixP[i][j];
        }
    }

    conv = convergence(N ,A, A_tag);
    if(conv <= eps){
        for( i=0 ; i<N ; i++) {
            for (int j = 0; j < N; j++) {
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
            for (int j = 0; j < N; j++) {
                A[i][j] = A_tag[i][j];
            }
        }

        CreatePmatrix(N,matrixP, A);
        Ptrans(N ,TmatrixP , matrixP);

        MatrixMultiplication(N ,V , Vectors , matrixP);
        for( i=0 ; i<N ; i++) {
            for (int j = 0; j < N; j++) {
                Vectors[i][j] = V[i][j];
            }
        }

        MatrixMultiplication(N, temp , TmatrixP , A);
        MatrixMultiplication(N, A_tag , temp , matrixP);

        conv = convergence(N, A, A_tag);
        iter++;
    }

    for( i=0 ; i<N ; i++) {
        for (int j = 0; j < N; j++) {
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
    double arr[N];
    double max = fabs(eigenvalues[0] - eigenvalues[1]);
    arr[0]= max;

    selectionSort(N, eigenvalues);

    for(int i=1 ; i< floor(N/2) ; i++){
        arr[i] = (eigenvalues[i] - eigenvalues[i+1]);
        if(arr[i]>max){
            max=arr[i];
            k=i;
        }
    }

    return k;
}

void NormalizedSpectralClustering(int N ,int K , int dimension , double**DataPoints){
    double ** WeightedAdjacencyMatrix ,**DiagonalDegreeMatrix , **NormalizedGraphLaplacian , **eigenvectors  ,** eigenvalues , **U , **t;
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

    if(K == -1) {
        for (i = 0; i < N; i++) {
            eign[i] = eigenvalues[i][i];
        }

        k = Eigengap(N, eign);
    }

    U = malloc(sizeof(double *) * k);
    t = malloc(sizeof(double *) * k);
    for (i = 0; i<K; i++) {
        U[i] = (double *) malloc(sizeof(double *) * N);
        t[i] = (double *) malloc(sizeof(double *) * N);
    }

    for (i = 0; i<K; i++) {
        for (j = 0 ; j<N ; j++){
            U[i][j] = NormalizedGraphLaplacian[i][j];
        }
    }

    
    for (i = 0; i<K; i++) {
        for (j = 0 ; j<N ; j++){

        }
        for (j = 0 ; j<N ; j++){

        }
    }


};


int main(){
    double A[3][3] = {{3, 2, 4},
                      {2, 0, 2},
                      {4, 2, 3}};

    double B[3][3] = {{1,0,0},
                      {0,4/(3* sqrt(2)),1/3},
                      {0,-1/3,4/(3* sqrt(2))}};

    double**AP;
    double**BP;
    double **matrix;
    double **matrixP;
    double **TmatrixP;
    int i ,j;

    AP = malloc(sizeof(double *) * 3);
    BP = malloc(sizeof(double *) * 3);
    matrix = malloc(sizeof(double *) * 3);
    matrixP = malloc(sizeof(double *) * 3);
    TmatrixP = malloc(sizeof(double *) * 3);
    for (i = 0; i<3; i++) {
        AP[i] = (double *) malloc((3) * sizeof(double));
        BP[i] = (double *) malloc((3) * sizeof(double));
        matrix[i] = (double *) malloc((3) * sizeof(double));
        matrixP[i] = (double *) malloc((3) * sizeof(double));
        TmatrixP[i] = (double *) malloc((3) * sizeof(double));
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            AP[i][j] = A[i][j];
            BP[i][j] = B[i][j];
        }
    }

    Jacobi(3 ,matrix , matrixP, AP);


    for( i=0 ; i<3 ; i++) {
        printf("\n");
        for (j = 0; j < 3; j++) {
            printf("%f ,", matrix[i][j]);
        }
    }
    return 1;
}
