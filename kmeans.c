#include <stdio.h>

int main() {

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
