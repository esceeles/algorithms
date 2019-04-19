#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void paren(int pos, int n, int open, int closed){

   
   static char str[100];
   
   if (closed == n){
      for (int i = 0; i < n*2; i++){
         printf("%c", str[i]);}
      return;
   }
   else{
      if (open > closed){
         str[pos] = ')';
      paren(pos+1, n, open, closed+1);
      }
      if (open < n){
         str[pos] = '(';
         paren(pos+1, n, open+1, closed);
      }
   }

}


int main(){

int n = 3;
paren(0, n, 0, 0 );

return 0;
}
