#include <stdio.h>
#include <malloc.h>
#include <math.h>

int main() {



    // FILE *in_file = fopen("C:\\Users\\nadsa\\Documents\\input_3.txt", "r");
    FILE *f;
    char c;
    f=fopen("C:\\Users\\nadsa\\Documents\\input_1.txt", "rt");

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
    //printf("%d, %d",cnt1,cnt2);
    double *arr[cnt1];
    for(int i=0 ; i<cnt2 ; i++){
        arr[i] = (double*)malloc((cnt1/cnt2 +1)* sizeof(double));
    }
    FILE *file;
    char C;
    file=fopen("C:\\Users\\nadsa\\Documents\\input_1.txt", "rt");
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
}

--------------------------------------------------------
    

#include <stdio.h>
#include <malloc.h>
#include <math.h>

/*double Pow (double x, double y){
    double result = x ;
    for (int i =0 ; i<y ; i++){
        result = result * x;
    }
}

double *DataPoints[rows];
    for(int i=0 ; i<rows ; i++){
        DataPoints[i] = (double*)malloc((dimension)*sizeof(double));
    }

 */


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


void Init(char filename[] ,int rows ,int dimension ,double DataPoints[][dimension]){

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

    //print arr
    /*for(int i=0 ; i<rows ; i++) {
        for (int j = 0; j < dimension; j++) {
            printf("%f", DataPoints[i][j]);
            printf(", ");
        }
        printf("\n");
    }*/

}

int main() {

    int dimension =0 ;
    int rows =0;

    char filename[] = "C:\\Users\\weamm\\Downloads\\input_3.txt";
    Sizefile(filename , &dimension , &rows);

    double DataPoints[rows][dimension];
    Init(filename , rows, dimension, &DataPoints[0]);

    //print
    for(int i=0 ; i<rows ; i++) {
        //printf("%d, " , i);
        for (int j = 0; j < dimension; j++) {
            //printf("%d, " , j);
            printf("%f" , DataPoints[i][j]);
            printf(",");
        }
        printf("\n");
    }

    return 0;
}
