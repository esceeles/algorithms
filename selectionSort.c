//selectionSort
//Ellie Sceeles

#include <stdlib.h>
#include <stdio.h>

int main(){

FILE *f1;
f1 = fopen("Data1.txt", "r");

int n;
fscanf(f1, "%d", &n);

int *A = malloc(n *sizeof(int));

for (int i = 0; i < n; i++){
   fscanf(f1, "%d", &A[i]);
}

int i, j, min;
for (i = 0; i < (n-1); i++){
   min = i;
   for (j = i+1; j < n; j++){
      if (A[j] < A[min]) min = j;}
   int temp = A[i];
   A[i] = A[min];
   A[min] = temp;
}

printf("Selection Sorted Array:\n");

for (int i = 0; i < n; i++){
   printf("%d ", A[i]);
}

return 0;
}
