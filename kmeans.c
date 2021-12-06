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


void kmeans(int k, int max_iter, char filename[], char *output_filename) {
    int dimension = 0;
    int rows = 0;
    int negative;
    int p = 1;
    int i, j;
    double num, fr;
    double **centroids_list;
    double **sum_array;
    double **DataPoints;
    int *count_array;
    FILE *output;
    int delta = 0;
    int iter;
    double tmp;

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
            tmp = 0;
            negative = 0;
            while ((C = fgetc(file)) != '.' && C != EOF) {
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



    count_array = (int *) malloc(sizeof(int) * k);

    iter = 0;
    while(delta==0 && iter < max_iter) {
        Init(k, dimension, count_array, sum_array);
        clustering(k, rows, dimension, count_array, sum_array, DataPoints, centroids_list);
        delta = calculate_delta(k, dimension, centroids_list, sum_array);
        iter++;
    }
    output = fopen(output_filename,"w");

    printf("\n");
    for (i = 0; i < k; i++) {
        for (j = 0; j < dimension; j++) {
            if(j!= (dimension-1)){
                fprintf(output,"%.4f,",centroids_list[i][j]);
            }
            else{
                fprintf(output,"%.4f",centroids_list[i][j]);
            }
        }
        fprintf(output,"\n");
    }
    fclose(output);
    free(centroids_list);
    free(DataPoints);
    free(sum_array);
    free(count_array);

}

/*************************/
int main() {

    char filename[] = "C:\\Users\\nadsa\\Documents\\SW-project\\hw1\\tests\\input_20.txt";
    char output[] = "C:\\Users\\nadsa\\Documents\\SW-project\\result_1.txt";
    kmeans(7, 200, filename,output);

}

