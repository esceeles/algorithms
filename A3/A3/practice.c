#include <stdlib.h>
#include <stdio.h>
//a[i] += b[i] * c[i]
void func(int *a, int *b, int *c)
{
//rdi = adr a; 
   asm ("movl (%0), %%eax; movl (%1), %%ecx; movl (%2), %%edx; imul %%ecx, %%edx; addl %%edx, %%eax; movl %%eax, (%0);"
   : "+r" (a) 
   : "r" (b), "r" (c)   
   : "eax", "ecx", "edx", "ebx", "esi", "edi"
   );

/*for (int i = 1; i <= 2; i++){
   asm ("addl $32, %%ebx; addl $32, %%ecx; addl $32, %%edx; imull %%ecx, %%edx; addl %%edx, %%ebx; movl %%ebx, (%0)"
   : "+r" (a) 
   : "r" (b), "r" (c)   
   : "eax", "ecx", "edx"
   );
};
*/}

int main(){
int a[3] = {2, 3, 4};
int b[3] = {3, 4, 5};
int c[3] = {4, 5, 6};
int *pta = a;
int *ptb = b;
int *ptc = c;

for (int i = 0; i < 3; i++){
   func(pta, ptb, ptc);
   pta++;
   ptb++;
   ptc++;
}

for (int i = 0; i < 3; i++){
   printf("%d ", a[i]);
}

return 0;
}
