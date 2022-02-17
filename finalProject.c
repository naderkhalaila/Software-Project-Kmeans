#include <stdio.h>
#include <math.h>


 void TheWeightedAdjacencyMatrix(double **matrix , double **DataPoints , int dimension , int N){
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


void TheDiagonalDegreeMatrix(double **matrix , double **WeightedAdjacencyMatrix, int dimension , int N){

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


void MatrixMultiplication (double ** matrix , double ** matrix1 , double ** matrix2  , int N){
    int i , j , k ;
    double value;
    for (i=0 ; i < N ; i++){
        printf("%d\n" , i);

        for (j=0 ; j<N ; j++){
            printf("%d\n" , j);

            value = 0;
            for (k=0; k<N ; k++){
                printf("%f\n" , matrix1[i][k]);
                printf("%f\n" , matrix2[k][j]);

                value+= matrix1[i][k]  * matrix2[k][j];
            }
            matrix[i][j] = value;
        }
    }
}


void TheNormalizedGraphLaplacian (double **matrix ,double **DiagonalDegreeMatrix ,
                                  double **WeightedAdjacencyMatrix, int dimension , int N){
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
    MatrixMultiplication(matrix1, DiagonalDegreeMatrix , WeightedAdjacencyMatrix, N);
    MatrixMultiplication(matrix, matrix1 , DiagonalDegreeMatrix, N);

    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            matrix[i][j] = Identity[i][j] - matrix[i][j];
        }
    }
}



void CreatePmatrix(double **matrix, double ** Amatrix, int N) {
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

void Ptrans(double **matrix , double** P, int N){
    int i ,j;
    for ( i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            matrix[i][j]= P[j][i];
        }
    }
}

double convergence(double ** matrix1 , double ** matrix2 , int N){
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


void Jacobi(double ** matrix , double **Vectors ,double** Lmatrix , int N){
    double eps = 1.0 * exp(-15);
    int Max_Iter = 100;
    double conv ;
    int iter =0;

    double matrixP[N][N];
    double TmatrixP[N][N];
    CreatePmatrix(matrixP , Lmatrix , N);
    Ptrans(TmatrixP , matrixP, N);

    double **A = Lmatrix;
    double **A_tag[N][N];
    double **temp[N][N];

    MatrixMultiplication(temp , TmatrixP , A, N);
    MatrixMultiplication(A_tag , temp , matrixP , N);

    Vectors =  matrixP;
    conv = convergence(A, A_tag, N);
    if(conv <= eps){
        matrix = A_tag;
    }


    while(conv > eps && iter < 100){
        double **V;
        A = A_tag;

        CreatePmatrix(matrixP, A , N);
        Ptrans(TmatrixP , matrixP , N);

        MatrixMultiplication(V , Vectors , matrixP ,N);
        Vectors=V;

        MatrixMultiplication(temp , TmatrixP , A, N);
        MatrixMultiplication(A_tag , temp , matrixP , N);

        conv = convergence(A, A_tag , N);
        iter++;

    }

    matrix = A_tag;
}


void Eigengap(double *eigenvalues, double** A_tag , double ** eigenvectors){

}

void main(){
    double **A = {{3,2,4},
                  {2,0,2},
                  {4,2,3}};

    double **B = {{2,3,2},
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
