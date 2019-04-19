#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef unsigned int *word_pointer;

float get_final_frac(word_pointer start, int n){
   float frac = 0.0;
   unsigned int i, j;
   unsigned int mask = 0x800000000;

   i = start[0];
   mask>>=(32-n);
   for (j = 0; j < n; ++j){
      if (i & mask)
         frac += 1/pow(2, (j+1));
      mask >>= 1;
   }

return frac;
}

float get_frac(int frac_shift, int n){
   float frac = get_final_frac((word_pointer) & frac_shift, n);
   return frac;
}

int calculate_sign_bit(word_pointer start, int n, int k){
   unsigned int i = start[0];
   unsigned int mask = 0x80000000;

   i =  start[0];
   mask >>=31;
   if (i & mask)
      return -1;
   return 0;
}

int get_sign_bit(int sign_shift, int n, int k){
   signed int sign = calculate_sign_bit((word_pointer) & sign_shift, n, k);
   return sign;
}

int check_hex(int hex, int n, int k){
   unsigned int value = hex >> (n+k+1);
   if (value == 0)
      return 1;
   return 0;
}

int main(){




return 0;
}
