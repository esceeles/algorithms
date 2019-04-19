#include <stdio.h>
#include <stdlib.h>

int main(){
FILE *f1;
f1 = fopen("Data1.txt", "r");

int n;

fscanf(f1, "%d", &n);

int A[n];

for (int i = 0; i < n; i++){
   fscanf(f1, "%d", &A[i]);
}

for (int i = 0; i < n-1; i++){
   for (int j = 0; j < n-1; j++){
      if (A[j+1] < A[j]){
         int temp = A[j];
         A[j] = A[j+1];
         A[j+1] = temp;
      }
   }
}


for (int i = 0; i < n; i++){
   printf("%d ", A[i]);
}

return 0;
}
