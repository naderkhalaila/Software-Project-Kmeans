#include <stdio.h>
#include <malloc.h>
#include <math.h>

void Sizefile(char *filename, int *dimension, int *rows) {

    FILE *f;
    char c;
    f = fopen(filename, "rt");

    int tmp_dimension = 0;
    int tmp_rows = 0;

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

double dist(double *point1, double *point2, int dimension) {
    int i;
    double sum = 0;
    for (i = 0; i < dimension; i++) {
        sum += pow((point1[i] - point2[i]), 2);
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
        for (j = 0; j < dimension; j++) {
            sum_array[i][j] = 0;

        }
    }

    for (i = 0; i < k; i++) {
        count_array[i] = 0;
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
            //printf( "%d , %d ,%f\n" , i, j,sum_array[i][j]);
        }
    }
}

int calculate_delta(int k, int dimension, double **centroids_list, double **sum_array) {
    double distance;
    int i, j;
    int bol = 1;
    for (i = 0; i < k; i++) {
        distance = dist(centroids_list[i], sum_array[i], dimension);
        if (distance > 0.001) {
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


void kmeans(int k, int max_iter, char filename[]) {
    int dimension = 0;
    int rows = 0;
    int negative;
    int p = 1;
    int i, j;
    double num, fr;
    double **centroids_list;
    double **sum_array;
    double **DataPoints;
    int delta = 0;
    int iter;

    Sizefile(filename, &dimension, &rows);
    DataPoints = malloc(sizeof(double *) * rows);
    for (i = 0; i < rows; i++) {
        DataPoints[i] = (double *) malloc((dimension) * sizeof(double));
    }
    FILE *file;
    char C;
    file = fopen(filename, "rt");
    for (i = 0; i < rows; i++) {
        for (j = 0; j < dimension; j++) {
            num = 1;
            negative = 0;
            while ((C = fgetc(file)) != '.' && C != EOF) {
                if (C == '-') {
                    negative = 1;
                } else {
                    num = num * ((double) C - 48);
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
    //printf("1\n");


    centroids_list = malloc(sizeof(double *) * k);
    for (i = 0; i < k; i++) {
        centroids_list[i] = (double *) malloc((dimension) * sizeof(double));
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < dimension; j++) {
            centroids_list[i][j] = DataPoints[i][j];
        }
    }

    sum_array = malloc(sizeof(double *) * k);
    for (i = 0; i<k; i++) {
        sum_array[i] = (double *) malloc((dimension) * sizeof(double));
    }
    //printf("2\n");


    int *count_array;
    count_array = (int *) malloc(sizeof(int) * k);

    iter = 0;
    while(delta==0 && iter < max_iter) {

        iter++;
        Init(k, dimension, count_array, sum_array);
        clustering(k, rows, dimension, count_array, sum_array, DataPoints, centroids_list);
        delta = calculate_delta(k, dimension, centroids_list, sum_array);


    }
    //printf("3\n");


    for (i = 0; i < k; i++) {
        for (j = 0; j < dimension; j++) {
            printf("%f", centroids_list[i][j]);
            printf(",");
        }
        printf("\n");
    }
    //printf("4\n");

}

/*************************/
int main() {

    char filename[] = "C:\\Users\\nadsa\\Documents\\input_1.txt";
    kmeans(3, 600, filename);
}



/*
    // FILE *in_file = fopen("C:\\Users\\nadsa\\Documents\\input_3.txt", "r");
    FILE *f;
    char c;
    f=fopen("C:\\Users\\nadsa\\Documents\\input_3.txt", "rt");

    int cnt1 = 0;
    int cnt2 = 0;
    while((c=fgetc(f))!=EOF){
        if(c==','){
            cnt1++;
        }
        if(c=='\n'){
            cnt2++;
        }
        //printf("%c",c);
    }
    fclose(f);
    printf("%d, %d",cnt1,cnt2);
    double *arr[cnt2];
    for(int i=0 ; i<cnt2 ; i++){
        arr[i] = (double*)malloc((cnt1/cnt2 +1)* sizeof(double));
    }
    FILE *file;
    char C;
    file=fopen("C:\\Users\\nadsa\\Documents\\input_3.txt", "rt");
    for(int i=0 ; i<cnt2 ; i++){
        for(int j=0 ; j<cnt1/cnt2 +1 ; j++){
            double num = 1;
            int negative = 0;
            while((C=fgetc(file))!='.' && C!=EOF){
                if(C=='-'){
                    negative = 1;
                }
                else{
                    num = num * ((double)C - 48);
                }

            }
            double fr = 0;
            int p = 1;
            while((C=fgetc(file))!=EOF && C!=',' && C!='\n'){
                fr = fr + ((double)C - 48) * pow(10 , -1*p);
                p++;
            }
            num = num + fr;
            if(negative == 1){
                num = -1*num;
            }
            arr[i][j] = num;
        }
    }
    fclose(file);



    for(int i=0 ; i<cnt2 ; i++) {
        for (int j = 0; j < cnt1/cnt2 + 1; j++) {
            printf("%f", arr[i][j]);
            printf(",");
        }
        printf("\n");
    }
    return 0;
      */
