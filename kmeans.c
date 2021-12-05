#include <stdio.h>
#include <malloc.h>
#include <math.h>

void Sizefile(char filename[],int *dimension ,int *rows){

    FILE *f;
    char c;
    f=fopen(filename, "rt");

    int tmp_dimension = 0;
    int tmp_rows = 0;

    while((c=fgetc(f))!=EOF){
        if(c==','){
            tmp_dimension++;
        }
        if(c=='\n') {
            tmp_rows++;
        }
    }

    tmp_dimension = tmp_dimension/tmp_rows +1;
    fclose(f);
    *dimension = tmp_dimension;
    *rows = tmp_rows;
}

double dist(double point1[],double point2[],int dimension){
    double sum = 0;
    for(int i=0 ; i<dimension ; i++){
        sum += pow((point1[i] - point2[i]),2);
    }
    sum = pow(sum,0.5);
    return sum;
}

int mindist(int k ,int dimension, double centroids_list[][dimension], double point[]){
    double min = dist(centroids_list[0] , point,dimension);
    int index = 0;
    for (int j =1; j<k; j++){
        double distance = dist(centroids_list[j] , point, dimension);
        if( distance < min){
            min =distance;
            index = j;
        }
    }
    return index;
}
void clustering(int dimension,double centroids_list[][dimension], double DataPoints[][dimension], double indexingList[])

void kmeans(int k, int max_iter, char filename[]){
    int dimension=0;
    int rows=0;
    Sizefile(filename,&dimension,&rows);
    double *DataPoints[rows];
    for(int i=0 ; i<rows ; i++){
        DataPoints[i] = (double*)malloc((dimension)* sizeof(double));
    }
    FILE *file;
    char C;
    file=fopen(filename, "rt");
    for(int i=0 ; i<rows ; i++){
        for(int j=0 ; j<dimension ; j++){
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
            DataPoints[i][j] = num;
        }
    }
    fclose(file);


    double *centroids_list[k];
    for(int i=0 ; i<k ; i++){
        centroids_list[i] = (double*)malloc((dimension)* sizeof(double));
    }
    for(int i=0 ; i<k ; i++){
        for(int j=0 ; j<dimension ; j++){
            centroids_list[i][j] = DataPoints[i][j];
        }
    }
    //int indexingList[k][rows]
    int size = 2;
    double **sumarrays;
    sumarrays=malloc(sizeof(int*)*k);
    for(int i=0 ; i<k ; i++){
        sumarrays[i] = (int*)malloc((dimension)* sizeof(int));
    }

    for(int i=0 ; i<k ; i++){
        for(int j=0 ; j<dimension ; j++){
            printf("%f",centroids_list[i][j]);
            printf(",");
        }
        printf("\n");
    }

    int** a = malloc(2*sizeof (int*));
    int* b = malloc(sizeof(int) * 2);
    b[0]=1;b[1]=2;
    int* c = malloc(sizeof(int) * 3);
    c[0]=1;c[1]=2;c[2]=3;
    a[0]=b;
    a[1]=c;


}



int main() {

    char filename[] = "C:\\Users\\nadsa\\Documents\\input_1.txt";
    kmeans(3,200,filename);



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
