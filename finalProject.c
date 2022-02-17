#include <stdio.h>
#include <math.h>


 void TheWeightedAdjacencyMatrix(int N ,int dimension ,double matrix[N][N] , double DataPoints[N][dimension]){
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


void TheDiagonalDegreeMatrix(int N , int dimension , double matrix[N][N] , double WeightedAdjacencyMatrix[N][N]){

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


void MatrixMultiplication (int N ,double matrix[N][N] , double matrix1[N][N] , double matrix2[N][N]){
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


void TheNormalizedGraphLaplacian (int N , double matrix[N][N] ,double DiagonalDegreeMatrix[N][N] ,
                                  double WeightedAdjacencyMatrix[N][N]){
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
    double matrix1[N][N];
    MatrixMultiplication(N , matrix1, DiagonalDegreeMatrix , WeightedAdjacencyMatrix);
    MatrixMultiplication(N , matrix, matrix1 , DiagonalDegreeMatrix);

    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            matrix[i][j] = Identity[i][j] - matrix[i][j];
        }
    }
}



void CreatePmatrix(int N , double matrix[N][N], double Amatrix[N][N]) {
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

    tita = (Amatrix[jj][jj] - Amatrix[ii][ii])/2*max;
    if(tita>=0){
        sign = 1;
    }else{
        sign = -1;
    }

    t = sign/(abs(tita) + sqrt(pow(tita,2) +1));
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

void Ptrans(int N, double matrix[N][N] , double P[N][N]){
    int i ,j;
    for ( i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            matrix[i][j]= P[j][i];
        }
    }
}

double convergence(int N, double matrix1[N][N] , double matrix2[N][N]){
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


void Jacobi(int N, double matrix[N][N] , double Vectors[N][N] ,double Lmatrix[N][N]){
    double eps = 1.0 * exp(-15);
    int Max_Iter = 100;
    double conv ;
    int iter =0;

    double matrixP[N][N];
    double TmatrixP[N][N];
    CreatePmatrix(N ,matrixP , Lmatrix);
    Ptrans(N , TmatrixP , matrixP);

    double **A = Lmatrix;
    double A_tag[N][N];
    double temp[N][N];

    MatrixMultiplication(N,temp , TmatrixP , A);
    MatrixMultiplication(N ,A_tag , temp , matrixP);

    Vectors =  matrixP;
    conv = convergence(N ,A, A_tag);
    if(conv <= eps){
        matrix = A_tag;
    }


    while(conv > eps && iter < 100){
        double **V;
        A = A_tag;

        CreatePmatrix(N,matrixP, A);
        Ptrans(N ,TmatrixP , matrixP);

        MatrixMultiplication(N ,V , Vectors , matrixP);
        Vectors=V;

        MatrixMultiplication(N, temp , TmatrixP , A);
        MatrixMultiplication(N, A_tag , temp , matrixP);

        conv = convergence(N, A, A_tag);
        iter++;

    }

    matrix = A_tag;
}


void Eigengap(double *eigenvalues, double** A_tag , double ** eigenvectors){

}

void main(){
    double A[3][3] = {{3,2,4},
                  {2,0,2},
                  {4,2,3}};

    double B[3][3] = {{2,3,2},
                  {1,2,3},
                  {0,5,6}};

    double matrixP[3][3];
    double TmatrixP[3][3];

    /*CreatePmatrix(matrixP , A , 3);

    printf("r\n");

    for(int i=0 ; i<3 ; i++){
        printf("\n");
        for(int j=0 ; j<3 ; j++) {
            printf("%f ," , matrixP[i][j]);
        }
    }*/

    MatrixMultiplication(matrixP , A, B , 3);

    for(int i=0 ; i<3 ; i++) {
        printf("\n");
        for (int j = 0; j < 3; j++) {
            printf("%f ,", matrixP[i][j]);
        }
    }
}
