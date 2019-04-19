//some code snippets attributed to randyrollofson on github

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

//check number of arguments
if (argc < 4)
   printf("Invalid number of command line arguments\n");

else{
//translate command arguements to something you can use
   unsigned int fracB, expB, hex;
   sscanf(argv[1], "%d", &fracB);
   sscanf(argv[2], "%d", &expB);
   sscanf(argv[3], "%x", &hex);

//check that values from command line are within appropriate range   
   if (fracB < 2 || fracB > 10)
      printf("Illegal number of fraction bits (%d). Should be between 2 and 10\n", fracB);
   else if (expB < 3 || expB > 5)
      printf("Illegal number of exponent bits (%d). Should be between 3 and 5\n", expB);
   else{

      //get frac bitfield out of hex bitfield
      unsigned int frac_shifted = hex <<(32-fracB) >>(32-fracB);
      float frac = 0;
      unsigned int mask = frac_shifted;

      //calculate frac value out of frac bitfield
      for (int j = 0; j < fracB; ++j)
      {
         if (mask & frac_shifted)
            frac += 1/(pow(2, (j+1)));
         frac_shifted >>= 1;
      }

      //get exp value from hex bitfield
      unsigned int exp = hex <<(32-(fracB + expB)) >> (32 - expB);

      //get the sign value from bitfield
      int S = hex << (32- (fracB + expB +1)) >> 31;
      if (S == 1)
         S = -1;
      else
         S = 0;
      int bias = (pow(2,(expB-1)) -1);
      int E;
      float M;

      //calculate odd cases
      if(exp == 7 || exp == 15 || exp == 31){
         if (frac == 0){
            if (S == -1)
               printf("-inf\n\n");
            else
               printf("+inf\n\n");
         }
         else
            printf("NaN\n\n");
      return 0;
      }
      //calculate E and M for denormalized cases
      else if (exp == 0){
         E = 1-bias;
         M = frac;   
      }
      //calculate E and M for normalized cases
      else{
         E = exp-bias;
         M = 1 + frac;
      }

      float V = pow((-1),S) * M * pow(2,E);
      printf("%f\n\n", V);
      }
}
return 0;
}
