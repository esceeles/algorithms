#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

int n;

printf("\nEnter an integer: ");
scanf("%d", &n);

printf("Your number in decimal: %d\n", n);

int count = 0;

while(n > 0){
   n = n/2;
   count += 1;
}

printf("Digits needed for binary: %d\n\n", count);

return 0;
}
