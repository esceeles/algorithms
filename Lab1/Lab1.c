#include <stdio.h>
#include <stdlib.h>

//looking for the max element in a list of positive integers
int maxInt(int n, int *A){
   int i;
   int m = 0;
   for (i = 0; i < n; i++){
      if (A[i] > m)
          m = A[i];
   }
return m;
}

//this is a sequential search for an integer in a list
int seqSearch(int start, int elem, int n, int *A){
   //start is pos within array to begin
   int i;
   if (start >= n)
      return -1;
   for (i = start; i < n; i++){
      if (A[i] == elem)
         return i;
   }
   return -1;
}

//this is a search for repeated elements
int repeatedEl(int num, int *A){
int re = 0;
   for (int i = 0; i < num; i++){
      for (int j = 0; j < num; j++){
         if (i != j){
            if (A[i] == A[j]){
               re = 1;
               printf("index: %d, value: %d\n", i, A[i]);
            }
         }
      }
   }
return re;
}

//finding the dot product
int dotProd(int num, int *A, int *B){
int dp = 0;
   for (int i = 0; i < num; i++){
      dp += (A[i]*B[i]);
   }

return dp;
}


int main(){ 

   //saves index of max int in list
   int A[1000];
   int max = A[0];

   FILE *f1;
   f1 = fopen("data1.txt", "r");
   int num;
   //read into variable num, the number of numbers in list
   fscanf(f1, "%d", &num);

   //reads file stuff into the array of ints
   for (int i = 0; i < num; i++){
      fscanf(f1, "%d", &A[i]);
   }
   printf("\nData1:\nnumInFile = %d\n", num);

   //find maxElement
   max = maxInt(num, A);
   printf("maxElement = %d\n", max);
   int x;

   //find index of number in dataset
   printf("Enter number to search: ");
   scanf("%d", &x);
   int s = seqSearch(0, x, num, A);
   if (s != -1)
      printf("index = %d\n\n", s); 
   else
      printf("NOT FOUND\n\n");

   //rewrites array A and num values for data2.txt
   FILE *f2;
   f2 = fopen("data2.txt", "r");
   fscanf(f1, "%d", &num);
   for (int i = 0; i < num; i++){
      fscanf(f2, "%d", &A[i]);
   }
   printf("Data2:\nnumInFile = %d\n", num);

   //search for repeated elements
   printf("Repeated Elements: \n");
   if (repeatedEl(num, A) == 0)
      printf("NO REPEATED ELEMENTS\n");

   //dot product
   FILE *f3;
   f3 = fopen("data3.txt", "r");
   //rewrites array A and num values for first vector from data3.txt
   fscanf(f3, "%d", &num);
   for (int i = 0; i < num; i++){
      fscanf(f3, "%d", &A[i]);
   }
   //writes second vector into Array B   
   int B[100];
   for (int i = 0; i < num; i++){
      fscanf(f3, "%d", &B[i]);
   }

   //displays arrays
   printf("\nData3, Array 1: \n");
   for (int i = 0; i < num; i++){
      printf("%d ", A[i]);
   }
   printf("\n\nData3, Array 2: \n");
   for (int i = 0; i < num; i++){
      printf("%d ", B[i]);
   }
   printf("\n\n");
   
   //calculate the dot product
   printf("Dot product: %d\n\n", dotProd(num, A, B));


return 0;
}
