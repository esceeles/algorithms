//Lab 4, DFS and BFS
//#include <queue>
#include <stdlib.h>
#include <stdio.h>

void DFS(int vert, int i, int G[][15], int* visited){
   int j;
   printf("\n%d", i);
   visited[i] = 1;

   for (j = 0; j < vert; j++){
      if (!visited[j] && G[i][j]==1)
         DFS(vert, j, G, visited);
   };
};

int main(){

   FILE *f;
   f = fopen("data7.txt", "r");

   int vert;

   fscanf(f, "%d", &vert);

   int matrix[15][15];

   for (int i = 0; i < vert; i++){
      for (int j = 0; j < vert; j++){
         fscanf(f, "%d", &matrix[i][j]);
      }
   }

   int visited[15] ={0};

//begins recurstion
   DFS(vert, 0, matrix, visited);
   printf("\n");




return 0;
}
