//matrix multiplication
//Ellie Sceeles
#include <stdio.h>
#include <stdlib.h>

int main(){
FILE *f1;
f1 = fopen("Data4.txt", "r");

int n;

fscanf(f1, "%d",&n);

//creates 1D dynamic arrays
int *A = malloc(n*n*sizeof(int));
int *B = malloc(n*n*sizeof(int));
int *C = malloc(n*n*sizeof(int));

//scans file into two arrays for multiplication
for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
      fscanf(f1, "%d", &A[(i*n) + j]);
   }
}

for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
      fscanf(f1, "%d", &B[(i*n) + j]);
   }
}

int sum = 0;

//does the multiplication
for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
      for (int k = 0; k < n; k++){
         sum += (A[(i*n)+k] * B[(k*n)+j]);
      }
   C[(i*n)+j] = sum;
   sum = 0;
   }
}

//print final array
for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
      printf("%d ", C[(i*n)+j]);
   }
printf("\n");

}

return 0;
}
