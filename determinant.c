#include <stdio.h>
#include <stdlib.h>

double det2(double* B){
   int n = 2;
   double ans = (B[0 * n + 0] * B[1*n+1]) - (B[1*n+0] * B[0*n+1]);
return ans;
};

void getCofactor(double* A, double* temp, int p, int q, int n){
   int i = 0, j = 0;

   for (int row = 0; row < n; row++){
      for (int col = 0; col < n; col++){
         if (row != p && col != q){
            temp[i*(n)+(j++)] = A[row * (n) + col];
            if (j == n-1){
               j = 0;
               i++;
            }
         }
      }
   }
}

double detn(double* A, int n){
   double ans = 0;
   if (n ==1)
      return A[0];
   double* temp = malloc(n*n*sizeof(double));
   int sign = 1;

   for (int i = 0; i < n; i++){
      getCofactor(A, temp, 0, i, n);
      ans += sign * A[i] * detn(temp, n-1);
      sign = -sign;
   }
return ans;
};

int main(){

double B[4] = {2, 3, 5, 8};
double A[9] = {2, 1, 1, 4, 1, 1, 1, 1, 1};


printf("%lf\n\n", det2(B));

printf("%lf\n\n", detn(A, 3));

}
