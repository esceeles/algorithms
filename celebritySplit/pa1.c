#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//creates all possible combinations of sets, length n as a bitfield and fills in 2d array with arrangments. 0 is empty, all other numbers mean present

void makeSets(int** M, int n){
   for(int i=1;i<(1<<n);i++)
   {
      for(int j=0;j<n;j++)
           if ((1<<j)&i){
               M[i][j] =1;
           }
           else {M[i][j] = 0;}
   }

return;
};

//take the combinations of sets and sums the prices assoriated with them in col-l
void sumSets(int**M, int* prices, int cols, int rows){
   for (int i =0; i < rows; i++){
      int sum = 0;
      for (int k = 0; k < cols; k++)
         if (M[i][k] != 0)
            sum += prices[k];
      M[i][cols-1] = sum;
   }

return;
}

//make sure that if jack gets zero, so does jill, and if they sell nothing, there must be houses to jack and jill instead
int checkValid(int x, int y, int z){
   if (x == 0 && y != 0) return 0;
   if (z == 0 && x == 0) return 0;
   if (z == 0 && y == 0) return 0;

//if valid, return 1
return 1;
}

//makes sure that arrangements for allocating houses to jack, jill, and liquidation do not allocate the same house to dif. people
int checkValid2(int **M, int cols, int total, int x, int y){
   for (int i = 0; i < cols-1; i++){
      if (M[x][i] != 0 && M[y][i] != 0) return 0;
   }
return 1;
}

//looks for ways to arrange the different sums in last column (col-1) so that jack and jill get the same ammount and it is maximized
int findArrangement(int **M, int cols, int rows, int total){
//only look at M[i][cols-1] for sums
int maxSum = 0;
int sumSell = 0;
int sell = 0;
for (int x = 0; x < rows; x++){
   for (int y = 0; y< rows; y++){
         if (x != y){
            int jack = M[x][cols-1];
            int jill = M[y][cols-1];
            if ((jack == jill)){
               sell = total - (2*jack);
               if (checkValid(jack, jill, sell) == 1){
                     if (checkValid2(M, cols, total, x, y) == 1){
                        if (jack > maxSum){
                              maxSum = jack;
                              sumSell = total - (jack + jill);
                        }
                     }
                  }
               }
            }
                           //}
         }}               

return sumSell;
};

int main(){
   char fileName[20];
   printf("Enter a file name: ");
   scanf("%s", fileName);
   FILE *f;
   f = fopen(fileName, "r");
   int n;
   fscanf(f, "%d", &n);
   int rows = pow(2, n);          //formula for sum of nth row of pascals triangle+1
   int *M[rows];
   int * prices = malloc (n * sizeof(int));
   int cols = n+1;
   
   for (int i = 0; i < rows; i++){        //creates rows/possible combinations
      M[i] = (int*)calloc((n+1),sizeof(int));
   }

   int total = 0;
   int x;
   for (int i = 0; i < n; i++){         //prices are read into prices array. 0 is left blank
      fscanf(f, "%d", &x);
      prices[i] = x;
      total += x;
   }
   //printf("\nTotal houses to divide: %d", n);
   //printf("\nTotal cost of those houses: $%d\n", total);
   makeSets(M, n);  
   sumSets(M, prices, cols, rows);
   int ammtToSell = findArrangement(M, cols, rows, total);
   
   printf("Ammount to liquidate: $%d\n", ammtToSell);

return 0;
}
