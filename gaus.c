#include <stdlib.h>
#include <stdio.h>

double* gaus(double*a, double* b, int n){
   double* A = malloc(n*(n+1)*sizeof(double));
   for (int i = 0; i < n; i++){
      for (int j = 0; j< n; j++){
         A[i*(n+1)+j] = a[i*n+j];
      }
   }

for (int i = 0; i < n; i++)
      A[i*(n+1)+(n)] = b[i];
  
  
for (int i = 0; i < n-1; i++){
      for (int j = i+1; j < n; j++){
         double temp = A[j*(n+1)+i];
         for (int k = i; k < n+1; k++){
            A[j*(n+1)+k] = A[j*(n+1)+k] - A[i*(n+1)+k] * temp / A[i*(n+1)+i];
         }
      }
      
   for (int i = 0 ; i < n; i++){
      for (int j = 0; j < n+1; j++){
         printf("%lf ", A[i*(n+1)+j]);
      }
      printf("\n");
   }
   printf("\n");

}

return A;
}



int main(){
   double A[9] = {2, -1, 1, 4, 1, -1, 1, 1, 1};
   double b[3] = {1, 5, 0};
   
   gaus(A, b, 3);

}
