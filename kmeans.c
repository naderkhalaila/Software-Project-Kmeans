#include <stdio.h>

int main() {

   int array_length = 100;
    int *array = (int*) malloc(array_length * sizeof(int));;
    printf("%d\n", sizeof(array));
    for(int i=0 ; i< sizeof(array) ; i++){
        printf("%f\n",array[i] );
    }
   
   // FILE *in_file = fopen("C:\\Users\\nadsa\\Documents\\input_3.txt", "r");
    FILE *f;
    char c;
    f=fopen("C:\\Users\\nadsa\\Documents\\input_3.txt", "rt");


    while((c=fgetc(f))!=EOF){
        printf("%c",c);
    }

    fclose(f);
    return 0;
}
