//submission.c
#ifndef _submission_c_
#define _submission_c_

#include <stdbool.h>
#include <stdlib.h>
#include "vectors.h"
#include <stdio.h>

bool vector_fma(struct doubleVector * a, const struct doubleVector * b, const struct doubleVector * c)
{
//length = length / 4;
double * ad = a->data;
double * bd = b->data;
double * cd = c->data;
int length = a->length;
if (length != b->length || length != c->length) return false;
int limit = length/4;
int i;

if (length< 20){
   __asm__ volatile (
      "1:\n"
      "movsd (%0), %%xmm0\n"
      "movsd (%1), %%xmm1\n"
      "movsd (%2), %%xmm2\n"
      "vfmadd231sd %%xmm1, %%xmm2, %%xmm0\n"
      "movsd %%xmm0, (%0)\n"
      "addq $8, %0\n"
      "addq $8, %1\n"
      "addq $8, %2\n"
      "loop 1b\n"
      :"+r" (ad)
      :"r" (bd), "r" (cd), "c" (length)
      :"xmm0", "xmm1", "xmm2"
      );
}

else{

   for (i = 0; i < limit-1; i++){
      __asm__ volatile (
         "1:\n"
         "vmovupd (%0), %%ymm0\n"
         "vmovupd (%1), %%ymm1\n"
         "vmovupd (%2), %%ymm2\n"
         "vfmadd231pd %%ymm1, %%ymm2, %%ymm0\n"
         "vmovupd %%ymm0, (%0)\n"
         "addq $32, %0\n"
         "addq $32, %1\n"
         "addq $32, %2\n"
      : "+r" (ad)
      : "r" (bd), "r" (cd)
      : "ymm0", "ymm1", "ymm2"
      );
   }
   for(i = i*4 ; i < length; i++){
      ad[i] += bd[i]*cd[i];  
   }
};
return true;
}

#endif
