#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

unsigned int expB = 4, fracB = 4, hex = 0x1af;

unsigned int frac_shift = hex <<(32 -fracB) >>(32-fracB);

printf("expB: %d, fracB: %d, hex: %x\n", expB, fracB, hex);
printf("frac shift: %d\n", frac_shift);

unsigned int exp = hex << (32-(fracB + expB)) >> (32 -expB);

printf("exp: %d\n", exp);

//get frac
float fracU = 0;
unsigned int i, j;
//unsigned int mask = 0x80000000;
//printf("mask: %x\n", mask);

unsigned int *start = (unsigned int *) &frac_shift;
i = start[0];

//mask >>=(32- fracB);

//printf("mask after 32b shift: %x\n", mask);

for (j = 0; j < fracB; ++j)
{
   if (i & frac_shift)
      fracU += 1/(pow(2, (j+1)));
   frac_shift >>= 1;
}
printf("frac: %f\n", fracU);


      int bias = (pow(2,(expB-1)) -1);
      int E;
      float M;
      int S = 0;
E = exp-bias;
         M = 1 + fracU;
printf("bias: %d, S: %d, M: %f, E: %d\n", bias, S, M, E);
      float V = pow((-1),S) * M * pow(2,E);
      printf("output: %f\n", V);


//find sign value
printf("hex: %x\n", hex);

return 0;
}
