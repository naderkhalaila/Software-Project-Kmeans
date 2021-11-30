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
