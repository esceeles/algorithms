//#define _BSD_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

int get_user_input(int r, char* arr[], int l, int i){
   char input[50];
   int match = 1;
//get input
   scanf("%49s", input);
   int il = strlen(input);
//check user input with word in list
   match = strncmp(input, arr[r], l);
   if (match != 0 || l != il){
      printf("Incorrect. Try again. \n");
      printf("word #%d is %s: ", i+1, arr[r]);
      match = 1;
   }
   else
      match = 0;

   return match;
}

int main(){
   char *arr[] = {
         "the", 
         "quick", 
         "brown", 
         "fox", 
         "jumps", 
         "over", 
         "the", 
         "lazy", 
         "dog"
   };
   int used[9] = {0};   
   struct timeval now, then, dif;
   gettimeofday(&now, NULL);
   int s = now.tv_sec;
   srand(s);

   printf("This is a game that tests typing speed.\n\nType the following words:\n");
   for (int i = 0; i < 9; i++){
//generate words and get length
      int r = rand() % 9;
      while(used[r] == 1){
            r = rand() % 9;
      }
            printf("word #%d is %s: ", i+1, arr[r]);
            used[r] = 1;
      int l = strlen(arr[r]);
//get user input
      while (get_user_input(r, arr, l, i) != 0){}
   }
   gettimeofday(&then, NULL);
   timersub(&now, &then, &dif);
   printf("\nCorrect! You're time is: %ld sec %6ld usec", (-1* dif.tv_sec), dif.tv_usec);
   printf("\n\n");
   return 0;
}
