// tools.c
//
// Ellie Sceeles & Bryant W. York
// November 2019
#include "math.h"
#include "tools.h"


int minDistance(int* dist, int* sptSet, int V) 
{ 
   int min = 10000, min_index; 
   
   for (int v = 0; v < V; v++) 
     if (sptSet[v] == 0 && dist[v] <= min) 
         min = dist[v], min_index = v; 
   
   return min_index; 
} 
   
int printSolution(int* dist, int n, int V) 
{ 
   printf("Vertex   Distance from Source\n"); 
   for (int i = 0; i < V; i++) 
      printf("%d tt %d\n", i, dist[i]); 
} 
   
void dijkstra(int* graph, int src, int V) 
{ 
     int* dist = malloc(V * sizeof(int));     //shortest distance from src to i 
   
     int* sptSet= malloc(V * sizeof(int)); //true if vertex is included in shortest path tree 
   
     //all distances begin as MAX and all vertexes not in SPT 
     for (int i = 0; i < V; i++) 
        dist[i] = 10000, sptSet[i] = 0; 
   
     dist[src] = 0; 
   
     for (int count = 0; count < V-1; count++) 
     { 
       // Pick next min dist vertex. u == src in first iteration
       int u = minDistance(dist, sptSet, V); 
   
       // Mark the vertex as processed 
       sptSet[u] = 1; 
   
       // Update dist value of the adjacent vertices to the picked vertex. 
       for (int v = 0; v < V; v++) 
   
         // Update dist[v] iff it is not in sptSet, there is an edge from u to v, 
         //and total path from src to v through u is smaller than current value of dist[v] 
         if (!sptSet[v] && graph[u* V + v] && dist[u] != 10000 && dist[u]+graph[u* V + v] < dist[v]) 
            dist[v] = dist[u] + graph[u* V + v]; 
     } 
   
     // print the constructed distance array 
     printSolution(dist, V, V); 
}


int printMST(int* final, int V, graph* G) 
{ 
printf("Edge \tWeight\n"); 
for (int i = 1; i < V; i++) 
    printf("%d - %d \t%d \n", final[i], i, G->adjmat[i* V + final[i]]); 
}; 

int minWeight(int* weight, int *set, int V){
   int min = 10000;
   int minIndex = 0;
   for (int v = 0; v < V; v++){
      if (set[v] == 0 && weight[v] < min){
         min = weight[v];
         minIndex = v;
      }
   }
return minIndex;      
}
  
void Prim(graph* G, int V){
   int* VT = calloc(V, sizeof(int));         //final tree
   int* weight = calloc(V, sizeof(int));     //holds weights
   int* set = calloc(V, sizeof(int));        //0 if not in final tree yet, 1 if in final tree

   for (int i = 0; i < V; i++){
      weight[i] = 10000;            //set weights to MAX
      set[i] = 0;                   //says no verts are in final tree yet
   }
   weight[0] = 0;                   //the first vertex will be weight of 0
   VT[0] = -1;                      //the first vertex will be in the final tree

   for (int i = 0; i < V-1; i++){
      int nextMin = minWeight(weight, set, V);           //get the next smallest weight for something that isn't already in VT
      set[nextMin] = 1;                   //set it as visited and in final tree
      for (int v = 0; v < V; v++){        //go through all the verts, check if theres an edge there and if it hasn't been already put in VT, and that the weight of that edge is less than its weight in our set
         if (G->adjmat[nextMin*V + v] && set[v] == 0 && G->adjmat[nextMin *V + v] < weight[v]){    
            VT[v] = nextMin;           //add the vertex into the final tree
            weight[v] = G->adjmat[nextMin * V + v];      //update the weight of that vertex
         }
      }
   }

printMST(VT, V, G);

}

void printDblArray(double *C, int n){

printf("\n");
for (int i = 1; i <= n+1; i++){
   for (int j = 0; j <= n; j++){
      printf("%lf ", C[i*(n+1)+j]);
   }
   printf("\n");
}
printf("\n\n");
}

//assuming p is size n + 1 and vals are between 1-n?
void avgCompOptBST(double* p, int n){
   double * C = calloc((n+2) * (n+1), sizeof(double));
   double * R = calloc((n+2)*(n+1), sizeof(double));
   for (int i = 1; i <= n; i++)         //1, 2, 3, 4
   {
      C[i * (n+1) + i] = p[i];
   }

   for (int i = 1; i <=n; i++) R[i*(n+1)+i] = i; 
   double minC = 1000.0;
   double temp = 0;
   int minR = 0;
   double sum = 0;
   printf("\n original:\n");
   printDblArray(C, n);
   
   for (int d = 1; d <= n-1; d++){          // k = 2, 3, 4
      for (int i = 1; i <= n - d; i++){
         int j = i + d;
	 minC = 1000;
         sum = 0;
         for(int z = i; z <=j; z++) {
		sum += C[z* (n+1) + z];
	 }
	 for(int k = i+1; k <= j; k++) {	
		if(C[i *(n+1) + k-1] + C[(k+1)* (n+1) + j] < minC) {
			minR = k;
			minC = C[i * (n+1) + (k-1)] + C[(k+1) *(n+1) + j]; 
		}
	 }
        C[i * (n+1) + j] = minC + (double) sum;
	R[i * (n+1) + j] = minR;
        printDblArray(C, n);
      }
   }
   
printf("\nFinal Array:\n");
printDblArray(C, n);
printf("\nRoot Table:\n");
printDblArray(R, n);
}

void floyds(int* p,int n) {
	int i,j,k;
	for (k=0;k<n;k++)
	  for (i=0;i<n;i++)
	   for (j=0;j<n;j++)
	    if(i==j)
	     p[i*n+j]=0; else
	     p[i*n+j]=min(p[i*n+j],p[i*n+k]+p[k*n+j]);
}


//access 2 day array as a 1d array
void warshal(int* p, int n)
{
 int i,j,k;
 for(k=0;k<n;k++)
  for(i=0;i<n;i++)
   for(j=0;j<n;j++)
    p[i*n+j]=max(p[i*n+j],p[i*n+k]&&p[k*n+j]);
}

int max(int a,int b){                                                       ;
if(a>b)
 return a;
else
 return b;
}

int min(int a, int b){
if (a<b)
   return a;
else
   return b;
}

void printCmpArray(comp* a, int n){
   printf("(real, imag)\n");
   for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
         printf("($%.3lf, %.3lf)", a[i*n+j].real, a[i*n+j].imag);
      }
      printf("\n");
   }
return;
}


comp * genCmpFn(int n){
   comp* f = malloc(n*sizeof(comp));
   comp* x;
   for (int i = 0; i < n; i++){
      x = randCmp(1000);
      f[i].real = x->real;
      f[i].imag = x->imag;
   }
for (int i = 0; i < n; i++)
   printf("%lf %lf, ", f[i].real, f[i].imag);
printf("\n");
return f;
}

comp* randCmp(int k){
   double a, b;
   a = rand() % k;
   b = rand() % k;
   comp* c = malloc(sizeof(comp));
   c->real = a;
   c->imag = b;

return c;
}

comp* powerC(comp* a, int n){
   polar* b = toPolar(a);
   double br = pow((b->radius), n);
   double bt = (b->theta) * n;
   polar * c = malloc(sizeof(polar));
   c->radius = br;
   c->theta = bt;
   comp *d = malloc(sizeof(comp));
   d = toComplex(c);
   free(b);
   free(c);

return d;
}

comp* make_dft(int n){
   comp* dft;
   comp* rootn;
   dft = malloc(n*n*sizeof(comp));
   rootn = malloc(1 * sizeof(comp));
   rootn = nthRoot(n);
   for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
         dft[(i*n)+j] = *powerC(rootn, i *j);
      }
   }

for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++)
      printf("(%.3f , %.3f) ", dft[i*n+j].real, dft[i*n+j].imag);
   printf("\n");
}

return dft;
}

comp* cvm(comp* a, comp* b, int n){
   comp* c = malloc(n*n*sizeof(comp));
   comp *sum = malloc(sizeof(comp)); 
   
   for (int i = 0; i < n; i++){
         c[i].real = 0.0; 
         c[i].imag = 0.0; 
         sum->real = 0.0;
         sum->imag = 0.0;
         for (int k = 0; k < n; k++){
            comp* mult = multc(&a[i+k], &b[k]);
            sum = addc(sum, mult);
         }    
         c[i].real = sum->real;
         c[i].imag = sum->imag;
   }

/*   // Different Version
for (int i = 0; i < n; i++)
         {
         int j = 0; // one column approach
             {
             for (int k = 0; k < n; k++)
                 {
                 // complex number style
                 comp* temp = multc(&a[i*n+k], &b[i*k +j]);
                 sum = addc(sum, temp);
                 }
             c[i*n+j] = *sum;
             }
         }
*/
for (int i = 0; i < n; i++){
   printf("(%lf, %lf) ", c[i].real, c[i].imag);
}
return c;
}

comp* cmm(comp* a, comp* b, int n){
   comp* c = malloc(n*n*sizeof(comp));
   comp *sum;
   for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
         c[i*n+j].real = 0.0;
         c[i*n+j].imag = 0.0;
         for (int k = 0; k < n; k++){
            comp* mult = multc(&a[i*n+k], &b[k*n+j]);
            sum = addc(sum, mult);
         }
      c[i *n +j].real = sum->real;
      c[i*n+j].imag = sum->imag;
      }
   }

return c;
}

polar* toPolar(comp* a){
   polar* c = malloc(sizeof(polar));
   c->radius = sqrt((a->real * a->real) + (a->imag * a->imag));
   c->theta = atan(a->imag/a->real);
return c;
}


comp* toComplex(polar* a){
   comp * c = malloc(sizeof(comp));
   c -> imag = a -> radius * sin(a -> theta);
   c -> real = a -> radius * cos(a -> theta);
   return c;
}

comp* multc(comp* a, comp * b){
   comp * c = malloc(sizeof(comp));
   c->real = (a->real * b->real) - (a->imag * b->imag);
   c->imag = (a->real * b->imag) + (a->imag * b->real);
return c;
}

comp* addc (comp* a, comp* b){
   comp * c = malloc(sizeof(comp));
   c->real = a->real + b->real;
   c->imag = a->imag + b->imag;
return c;
}

comp* conjugate(comp* a){
   comp * c = malloc(sizeof(comp));
   c->real = a->real;
   c->imag = -1 * a->imag;
return c;
}

comp* nthRoot(int n){
   comp * c = malloc(sizeof(comp));
   double theta = (2*M_PI)/(double)n;
   c->real = cos(theta);
   c->imag = sin(theta);
return c;
}

pt2df* jarvisMarchClass(pt2df* A, int n){
   double *next = malloc(n * sizeof(double));
   double *prev = malloc(n * sizeof(double));
   double * x = malloc(n * sizeof(double));
   double * y = malloc(n * sizeof(double));
   
   int q, p;
   int l = 1;
   for (int i = 0; i < n; i++){
      x[i] = A[i].x;
      y[i] = A[i].y;
   }
   for (int i = 2; i < n; i++){
      if (x[i] < x[l])
         l = i;
   }
   p = l;
   do{
   q = p + 1;
   for (int i = 0; i < n; i++){
     if ((p != i) && CCW(&A[p], &A[i], &A[q]))
      q = i;
   }
   next[p] = q;
   prev[q] = p;
   p = q;
   } while(p != l);

printf("\nnext order: \n");
for (int i = 0; i < n; i++){
   int x = next[i];
   if (x != 0) printf("%f, %f\n", A[x].x, A[x].y);
}

return A;
}


void CH3(pt2df * data, int n){
   pt2df *next = malloc(n * sizeof(pt2df));
   pt2df *prev = malloc(n * sizeof(pt2df));
   //p = 0, q = 1, r = 2
   pt2df p = data[0];
   pt2df q = data[1];
   pt2df r = data[2];
   double x = CCW(&p, &q, &r);
   if (x > 0){
      printf("CCW\n");
      next[0] = q;
      prev[1] = p;
      next[1] = r;
      prev[2] = q;
      next[2] = p;
      prev[0] = r;
   }
   else if (x == 0)
      printf("linear\n");
   else{
      printf("CW\n");
      next[0] = r;
      prev[2] = p;
      next[1] = p;
      prev[0] = q;
      next[2] = q;
      prev[1] = r;
   }

for (int i = 0 ; i < n; i++){
   printf("%f, %f ", next[i].x, next[i].y);
}
return;
}

double CCW(pt2df *p, pt2df *q, pt2df *r){
   double a, b, c, d, e, f;
   a = p->x;
   b = p-> y;
   c = q-> x;
   d = q-> y;
   e = r-> x;
   f = r->y;
   double ans;
   ans = ((d-b)*(e-a)) < ((f-b)*(c-a));

return ans;
}

//finds orientation of ordered triplet (p, q, r)
int orientation(pt2d p, pt2d q, pt2d r) {
    int val = ((q.y - p.y) * (r.x - q.x)) - ((q.x - p.x) * (r.y - q.y));

    if (val == 0) {
        return 0; //colinear
    }
    else if (val > 0) {
        return 1; //clockwise
    }
    else {
        return 2;  //counterclockwise
    }
}

void jarvisMarchMy(pt2d* points, int n) {
    if (n < 3) {
        return; //because there's no hull if there's not at least a triangle
    }

    //to store hull points in
    pt2d hull[n];
    int index = 0;

    //find the point closest to 0 on x axis (leftmost)
    int l = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].x < points[l].x) {
            l = i;
        }
    }

    //start at leftmost point, move counterclockwise until we're at start point again
    int p = l;
    int q;
    do
    {
        //add current point to result
        hull[index].x = points[p].x;
        hull[index].y = points[p].y;
			
        index += 1;

        //search for a point 'q' such that orientation(p, x, q) is counterclockwise for all points 'x' and to keep track of last visited most counterclockwise point in q
        q = (p+1) % n;
        for (int i = 0; i < n; i++) {
        //if i is more counterclockwise than current q, then update q
            if (orientation(points[p], points[i], points[q]) == 2) {
                q = i;
            }
        }

        //q is now most counterclockwise vs p
        //set p as q for next iteration, so that q is added to hull
        p = q;

    } while (p != l); //ends when we're back at start

    for (int i = 0; i < index; i++) {
        printf("(%d, %d)\n", hull[i].x, hull[i].y);
    }
}


int interpolationSearch(int* arr, int size, int key)
{
    int low = 0;
    int high = size - 1;
    int mid;

    while ((arr[high] != arr[low]) && (key >= arr[low]) && (key <= arr[high])) {
        mid = low + ((key - arr[low]) * (high - low) / (arr[high] - arr[low]));

        if (arr[mid] < key)
            low = mid + 1;
        else if (key < arr[mid])
            high = mid - 1;
        else
            return mid;
    }

    if (key == arr[low])
        return low ;
    else
        return -1;
}

int divAndConqExponentiation(int x, int n) {
     if (n < 0) {
         return divAndConqExponentiation(1/x, -n);
     }
     else if (n == 0) {
         return 1;
     }
     else if (n == 1) {
         return x;
      }
     else if (n % 2 == 0) {
         return  divAndConqExponentiation(x*x, n/2);
     }
     else if (n % 2 == 1) {
         return x* divAndConqExponentiation(x*x, (n-1)/2);
     }
}

int* matrixMultiply(int* A, int* B, int n){
int * C = malloc(n*n*sizeof(int));
int sum = 0;
for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
      for (int k = 0; k < n; k++){
         sum += (A[(i*n)+k] * B[(k*n)+j]);
      }   
   C[(i*n)+j] = sum;
   sum = 0;
   }   
}
return C;
}

//also CCW
double angleBtw(pt2d* v1, pt2d *v2){
   double n1, n2, dp, theta;
   n1 = f2norm(v1);
   n2 = f2norm(v2);
   dp = dotProd2D(v1, v2);
   theta = acos(dp/(n1*n2));
return theta;
}

double f2norm(pt2d* v){
   double res;
   int x1, y1;
   x1 = v->x;
   y1 = v->y;
   res = sqrt((x1*x1) + (y1*y1));
return res;
}

double dotProd2D(pt2d *A, pt2d *B){
int dp = 0;
int x1, y1, x2, y2;
x1 = A->x;
y1 = A->y;
x2 = B->x;
y2 = B->y;
dp = (x1*x2) + (y1*y2);

return dp;
}

int dotProd(int num, int *A, int *B){
int dp = 0;
   for (int i = 0; i < num; i++){
      dp += (A[i]*B[i]);
   }   

return dp; 
}

int gcd(int m, int n){
   if (n != 0){
      int r = m % n;
      gcd(n, r);
   }
   if (n == 0)
      return m;
}

int medianSelect(int *arr, int n){

   int pivot = arr[0];
   int *L = calloc(n, sizeof(int));
   int *R = calloc(n, sizeof(int));
   int indL = 0;
   int indR = 0;

   if (n == 1) return arr[0];

   for (int i = 0; i < n; i++){
      if (arr[i] <= pivot){
         L[indL] = arr[i];
         indL++;
      }
      else{
         R[indR] = arr[i];
         indR++;
      }
   }
   
   int *L2 = malloc((indL+1) * sizeof(int));
   int *R2 = malloc((indR+1) * sizeof(int));

   for (int i = 0; i < indL+1; i++)
      L2[i] = L[i];

   for (int i = 0; i < indR+1; i++)
      R2[i] = R[i];

   if (indL == indR && n%2!= 0) return pivot;
   
   else if(( indL == indR) && (n %2 == 0)){
      int min = maxInteger(indL, L2);
      min += minInteger(indR, R2);
      return (double) (min/2);
   }

   if (indR > indL)  medianSelect(R2, n-indL-1);
   if (indL > indR)  medianSelect(L2, n-indR-1);
}

int genRandInt(int n)
{
   int z;
   z = rand() % n;
   return z;
}

void genRandIntToFile(){
    int A[1000];
    char fileName[20];
    printf("enter file name to print to: ");
    scanf("%s", fileName);
    printf("enter number of ints to generate: ");
    int num;
    scanf("%d", &num);
    FILE *fip;
    fip = fopen(fileName, "w");
    
    fprintf(fip, "%d", num);
    fprintf(fip, "\n\n");
    for(int i =0; i < num; i++)
    {
        A[i] = genRandInt(num);
        fprintf(fip, "%d\n", A[i]);
    }
    
    fclose(fip);
return;
}


void swap(int* x, int* y){
   int temp = *x;
   *x = *y;
   *y = temp;
}

int partition(int * arr, int l, int r){
int x = arr[r], i = l;
for (int j = l; j <= r-1; j++){
   if (arr[j] <= x){
      swap(&arr[i], &arr[j]);
      i++;
   }
}
swap(&arr[i], &arr[r]);
return i;
}

int quickSelectMedian(int* arr, int l, int r, int k){          //need to pass in median index as k, in main... (n/2) +1;
   int pivot = partition(arr, l, r);
   if (pivot-l ==k -1)
      return arr[pivot];
   if (pivot -l > k-1)
      return quickSelectMedian(arr, l, pivot-1, k);
   return quickSelectMedian(arr, pivot+1, r, k-pivot +l-1);
}


void heapify(int * array, int n){ 
 n = n-1;
 int k;
    int v;
    int j;
    int heap_bool;
    int mid = n/2;

   for(int i = mid; i >= 1; i--){ //Important: notice start/end conditions/decrement format
        k = i;
        v = array[k];
        heap_bool = 0;
            while(heap_bool != 1 && 2*k <= n){ 
                j = 2 * k;
                if(j < n)  //there are two children
                {   
                    if (array[j] < array[j+1])
                        {j = j+1;}
                }   
                if (v >= array[j])
                    {heap_bool = 1;}
                else {
                    array[k] = array[j];
                    k = j;
                    }   
                array[k] = v;
        }   
    }   
}

void heapSort(int n, int* array)
{
    int temp;
    n = n-1;
    while(n > 0) {
        heapify(array, n);
        temp = array[1];
        array[1] = array[n];
        array[n] = temp;
        n -= 1;
      }

}


//DFA
void DFA(){
   int len;
   int state = 0;
   char inp[80];
   int i = 0;
   printf("DFA for even number of 1's\n");
   printf("Enter some input: ");
   scanf("%s", inp);
   len = strlen(inp);

   while (i < len){
      if (inp[i] == '1')
         state = !state;
      i++;
   }

   if (state == 0)
      printf("%d", 1);

   else
      printf("%d", 0);
}

int* readIntArray(char * fileName, int *num){
   int x;
   int numElements;
   //num of elements is first in file
   FILE *f1;
   f1 = fopen(fileName, "r");
   fscanf(f1, "%d", &numElements);
   int *A = malloc(numElements * sizeof(int));

   for (int i = 0; i < numElements; i++){
      fscanf(f1, "%d", &x);
      A[i] = x;
   }
   *num = numElements;

   fclose(f1);
return A;
}

void writeIntArray(char * fileName, int *A, int num){
   FILE *f1;
   f1 = fopen(fileName, "w");
   fprintf(f1, "%d\n\n", num);

   for (int i = 0; i < num; i++){
      fprintf(f1, "%d\n", A[i]);
   }

   fclose(f1);
return;
}

void writeDblArray(char * fileName, double *A, int num){
   FILE *f1;
   f1 = fopen(fileName, "w");
   fprintf(f1, "%d\n\n", num);

   for (int i = 0; i < num; i++){
      fprintf(f1, "%lf\n", A[i]);
   }

   fclose(f1);
return;
}

double* readDblArray(char* fileName, int * num){
   double x;
   int numElements;
   //num of elements is first in file
   FILE *f1;
   f1 = fopen(fileName, "r");
   fscanf(f1, "%d", &numElements);
   double *A = malloc(numElements * sizeof(double));

   for (int i = 0; i < numElements; i++){
      fscanf(f1, "%lf", &x);
      A[i] = x;
   }
   *num = numElements;

   fclose(f1);
return A;
}

char* readCharArray(char* fileName, int * num){
   FILE *f1;
   char *A;
   f1 = fopen(fileName, "r");
   int n;
   fscanf(f1, "%d", &n);
   A = malloc(n * sizeof(char));
   for (int i = 0; i<n; i++)
      fscanf(f1, "%c", &A[i]);
   *num = n;
   fclose(f1);
   return A;
}


//read string array
void writeStrArray(char* fileName, char ** A, int num){
        FILE *f1;
        int n;
        f1=fopen(fileName, "w");
        fprintf(f1, "%d\n\n", num);
         
        for (int i=0; i < num; i++)
        {   
         fprintf(f1, "%s\n", A[i]);
        } 
fclose(f1);
return;
}

char** readStrArray(char* fileName, int * num){
        FILE *f1;
        char** A;
        char xx[100];
        double avg;
        int n, i;
        f1 = fopen(fileName, "r");
        fscanf(f1,"%d",&n);
        A = malloc(n * sizeof(char*));
        for (i=0;i<n;i++)
        {   
                fscanf(f1,"%s",xx);
                n = strlen(xx)+1;
                A[i] = malloc((n) * sizeof(char));
                strcpy(A[i],xx);
        } 
*num = n;
fclose(f1);
return A;
}

void write2Dpts(char* fileName, pt2d *A, int n){
   FILE *f1;
   f1= fopen(fileName, "w");
   fprintf(f1, "%d", n);
   for (int i = 0; i < n; i++)
      fprintf(f1, "%d %d\n", A[i].x, A[i].y);

fclose(f1);
return;
}

pt2d* read2Dpts(char* fileName, int* points){
   FILE *f1;
   f1 = fopen(fileName, "r");
   double x, y;
   int p;
   fscanf(f1, "%d", &p);
   struct pt2d* A = malloc(p * sizeof(pt2d));
   for (int i = 0; i < p; i++){
      fscanf(f1, "%lf %lf", &x, &y);
      A[i].x =  x;
      A[i].y = y;
   }
   *points = p;
return A;
}


pt2df* read2Dptsf(char* fileName, int* points){
   FILE *f1;
   f1 = fopen(fileName, "r");
   float x, y;
   int p;
   fscanf(f1, "%d", &p);
   struct pt2df* A = malloc(p * sizeof(pt2df));
   for (int i = 0; i < p; i++){
      fscanf(f1, "%f %f", &x, &y);
      A[i].x =  x;
      A[i].y = y;
   }
   *points = p;
return A;
}

void closestPair(pt2d *X, int points, int* A){
   float minDist = 10000000;
   int j;
   for (int i = 0; i < points-1; i++){
      for (j = 0; j < points; j++){
         double sum = pow((X[i].x-X[j].x), 2) + pow((X[i].y-X[j].y), 2);
         double ans = sqrt(sum);
         if ((ans < minDist) && (i != j)){
            minDist = ans;
            A[0] = X[i].x;
            A[1] = X[i].y;
            A[2] = X[j].x;
            A[3] = X[j].y;}
      }   
   }   

printf("minDist: %lf\n", minDist);
return;
}

//Brute Force Char Matching
void bfCharMatch(char* text, int n){
   printf("Text: %s\n", text);
   int j;
   char pattern[100];
   printf("Enter patter to search for: ");
   scanf("%s", pattern);
   int lentext = n;
   int lenpattern = strlen(pattern);
      for (int i = 0; i < lentext; i++){
         for (j = 0; j < lenpattern; j++){
            if (pattern[j] != text[i+j])
               break;
         }
         if (j == lenpattern){
            printf("match at index: %d\n", i);
            return;
         }
      }
   printf("no match\n");
return;
}

//Brute Force String Matching
void bfStrMatch(char** A){
   int j;
   char *text = A[0];
   char pattern[100];
   printf("Enter pattern to search for: ");
   scanf("%s", pattern);
   int lentext = strlen(text);
   int lenpattern = strlen(pattern);
      for (int i = 0; i < lentext; i++){
         for (j = 0; j < lenpattern; j++){
            if (pattern[j] != text[i+j])
               break;
         }
         if (j == lenpattern){
            printf("match at index: %d\n", i);
            return;
         }
      }
   printf("no match\n");
return;
}

//elementUniqueness
int uniqueInt(int * A, int n){ 
   for (int i = 0; i < n-1; i++){
      for (int j = i+1; j < n; j++){
         if(A[i] == A[j]){
            return 0;
         }   
      }   
   }   
   return 1;
}


// Queue Functions
struct queue_node* make_queue() {
  struct queue_node* p;
  p = NULL;
  return p;
}

int isQueueEmpty(struct queue_node* head){
  if(head == NULL) {
    return 1;
  } else {
    return 0;
  }
}

int insert(struct queue_node** head, struct queue_node* element) {
  if(*head == NULL) {
    *head = element;
    return SUCCESS;
  }
  struct queue_node* end = *head;
  while(end->next != NULL) {
    end = end->next;
  }
  end->next = element;
  return SUCCESS;
}

struct queue_node* removeq(struct queue_node** head) {
  struct queue_node* temp = *head;
  *head = (*head)->next;
  return temp;
}

struct queue_node* make_queue_node() {
  struct queue_node* p = malloc(sizeof(struct queue_node));
  p->next = NULL;
  return p;
} 

void print_queue(struct queue_node * head){
	
	printf("\nprint_queue done\n");
}

// Stack Functions

struct stack_node* make_stack() {
  struct stack_node* p;
  p = NULL;
  return p;
}

int isStackEmpty(struct stack_node* head){
  if(head == NULL) {
    return 1;
  } else {
    return 0;
  }
}

int push(struct stack_node** head, struct stack_node* element) {
  if(*head == NULL) {
    *head = element;
    element->next = NULL;
    return SUCCESS;
  } else {
    struct stack_node* pushed = *head;
    *head = element;
    element->next = pushed;
    return SUCCESS;
  }
}

struct stack_node* pop(struct stack_node** head) {
  struct stack_node* temp = *head;
  *head = (*head)->next;
  return temp;
}

struct stack_node* make_stack_node() {
  struct stack_node* p = malloc(sizeof(struct stack_node));
  p->next = NULL;
  return p;
} 



// Search Functions
int maxInteger(int n, int *A){
	int m, i;
	m = A[0];
	for (i=0; i<n; i++)
		if (A[i] > m) m = A[i];
	return m;
}

int minInteger(int n, int *A){
	int m, i;
	m = A[0];
	for (i=0; i<n; i++)
		if (A[i] < m) m = A[i];
	return m;
}

int seqSearch(int start, int elem, int n, int *A){
	// start is the position within the array to begin
	// returns -1 if elem is not found
	int i;
	if (start >= n) return -1;
	
	for (i=start; i<n; i++)
		if (A[i] == elem) return i;
	
	return -1;
}

int recBS(int l,int u,int key,int n,int * A){

     int mid;

     if(l <= u){
          mid=(l+u)/2;
          if(key == A[mid]){
              return mid;
          }
          else if(key < A[mid]){
              return recBS(l,mid-1,key,n,A);
          }
          else
              return recBS(mid+1,u,key,n,A);
     }
     
	 return -1;
	 
}

int nonrecBS(int K, int n, int *A){
	int l, r, m;
	
	l = 0;
	r = n-1;
	
	while (l <= r){
		m = (l+r)/2;  // int div gives floor
		if (K == A[m]) return m;
		else if (K < A[m]) r = m-1;
		else l = m+1;
 	}
	
	return -1;
}

// Sort Functions
void swapint(int * x, int * y)
{
	int temp;
	temp = *x;
	*x = *y;
	*y = temp;	
}

void selectionSort(int n, int* A)
{
	int i, j, minx;
	int nm1 = n - 1;
	int nm2 = n - 2;
	for(i=0; i<=nm2; i++){
		minx = i;
		for(j=i+1; j<=nm1; j++){
		if(A[j]<A[minx]) minx = j;
		}
		swapint(&A[i],&A[minx]);	
	}
}

void bubbleSort(int n, int* A)
{
	int i, j;
	int nm2 = n - 2;
	for(i=0; i<=nm2; i++){
		for(j=0; j<=nm2-i; j++){
		if(A[j+1]<A[j]) swapint(&A[j],&A[j+1]);
		}
			
	}
}

void mergeSort(int n, int* A)
{}

void insertionSort(int n, int* A)
{}

void quickSort(int n, int* A)
{}

// Graph Functions
void make_adjmat(struct graph *g){
	int nv, ne;
	int i, e1, e2;
	struct edge *e;
	int *am;
	
	nv = g->num_verts;
	am = (int *)calloc((nv*nv),sizeof(int));
	ne = g->num_edges;
	e = g->edgelist;
	for(i=0; i<ne; i++){
		e1 = e->from; 
		e2 = e->to;
		am[e1*nv+e2] = 1;
		am[e2*nv+e1] = 1;
		e = e->next;
	}
	g->adjmat = am;
}

struct graph * readGraph(char * fname){
	
	FILE *f1;
	int numg, nv, ne;
	int i, e1, e2;
	char nm[10];
	int type;
	struct graph * g;
	struct vertex **verts, *v;
	struct edge *head, *e, *p;
	
	
	g = malloc(sizeof(struct graph));
	
	f1 = fopen(fname, "r");
	fscanf(f1,"%d", &numg);  // read the number of graphs
	fscanf(f1,"%s",nm);
	strcpy(g->name,nm);
	fscanf(f1, "%d" , &type);  //1=undirected, 2=directed
	g->gtype = type;
	fscanf(f1,"%d",&nv);
	g->num_verts = nv;
	verts = malloc(nv*sizeof(struct vertex *));
	for (i = 0; i < nv; i++) {
		v = malloc(sizeof(struct vertex));
		fscanf(f1,"%s",nm);
		strcpy(v->name,nm);
		v->num = i;
		verts[i] = v;
	}
	g->vertices = verts;
	
	fscanf(f1,"%d",&ne);
	g->num_edges = ne;
	fscanf(f1,"%d %d",&e1,&e2);
	head = malloc(sizeof(struct edge));
	head->from = e1;
	head->to = e2;
	p = head;
	
	for(i=1; i<ne; i++){
		fscanf(f1,"%d %d",&e1,&e2);
		e = malloc(sizeof(struct edge));
		e->from = e1;
		e->to = e2;
		p->next = e;
		p = e;	
	}

	g->edgelist = head;
	
	fclose(f1);
	return g;
	
}

void printGraph(struct graph * g){
	int nv, ne;
	int i, e1, e2;
	struct vertex *v;
	struct edge *e;
	
	printf("\nName = %s",g->name);
	printf("\nGraph Type = %d",g->gtype);
	printf("\nNumber of Vertices = %d",g->num_verts);
	nv = g->num_verts;
	for (i = 0; i < nv; i++) {
		v = g->vertices[i];
		printf("\nVertex number = %d",v->num);
		printf("\nVertex name = %s",v->name);
		printf("\n Vertex mark = %d",v->mark);
	}
	
	ne = g->num_edges;
	printf("\nNumber of edges = %d", g->num_edges);
	e = g->edgelist;
	for(i=0; i<ne; i++){
		e1 = e->from; 
		e2 = e->to;
		printf("\nEdge from = %d  to = %d", e1, e2);
		e = e->next;
	}

	printf("\nAdjacency Matrix\n");
	print_adjmat(nv,g->adjmat);
	return;
	
}

void print_adjmat(int n, int *mat)
{
	int i, j;
	for (i=0; i<n; i++){
		for(j=0;j<n;j++) {
			printf("%2d",mat[i*n+j]);
		}
		printf("\n");
	}
	
}


int adjacent(int v1, int v2, int n, int *mat){
	if(mat[v1*n+v2] == 1) return 1;
	else return 0;
}
	
//bfs

void bfs(struct graph * g, struct vertex *v, int *count){
	struct queue_node *nd, *q, *ndd;
	struct vertex *w;
	int i, wn, fn, val, ct;
	int nv, empty, ans;
	
//	printf("\nInside bfs\n");
//	scanf("%d",&ans);
	(*count)++;
	v->mark = *count;
// initialize a queue	
	q = make_queue_node();
	q->vert = v;
// get number of vertices
	nv = g->num_verts;
	empty = isQueueEmpty(q);
	
	while (!empty){
		// get the front of the queue
		fn = (q->vert)->num;
		for(i=0; i<nv; i++){
			w = g->vertices[i];
			wn = w->num;
			if (adjacent(fn,wn,nv,g->adjmat))
			{
			  if (w->mark == 0){
				(*count)++;
				w->mark = *count;
				nd = make_queue_node();
				nd->vert = w;
				nd->next = NULL;
				insert(&q,nd);
				//printf("\nInside set mark  count = %d w->num = %d\n",*count,w->num);
				//scanf("%d",&ans);
			  }
			}
		}
		// remove front of queue
		ndd = removeq(&q);
		val = (ndd->vert)->num;
		ct = (ndd->vert)->mark;
		printf("\nVertex = %d  count = %d\n",val,ct);
//		scanf("%d",&ans);
		free(ndd);
		empty = isQueueEmpty(q);
//		if (!empty) print_queue(q);
//		printf("\nEnd of while   empty = %d",empty);
	}
}


//BFS
void BFS(struct graph *g){
	int nv, count;
	int i;
	struct vertex *v;
	
	nv = g->num_verts;
	//mark all vertices 0
	for (i = 0; i<nv; i++) {
		v = g->vertices[i];
		v->mark = 0;
	}
	
	count = 0;
	
	for (i=0; i<nv; i++){
		v = g->vertices[i];
		if (v->mark == 0) bfs(g, v, &count);
	}
}


//dfs
void dfs(struct graph * g, struct vertex *v, int *count){
struct stack_node *nd, *st, *ndd;
	struct vertex *w;
	int i, wn, fn, val, ct;
	int nv, empty, ans;
	
//	printf("\nInside dfs\n");
//	scanf("%d",&ans);
	(*count)++;
	v->mark = *count;
// initialize a stack	
	st = make_stack_node();
	st->vert = v;
// get number of vertices
	nv = g->num_verts;
	empty = isStackEmpty(st);
	
	while (!empty){
		// get the top
		fn = (st->vert)->num;
		for(i=0; i<nv; i++){
			w = g->vertices[i];
			wn = w->num;
			if (adjacent(fn,wn,nv,g->adjmat))
			{
			  if (w->mark == 0){
				(*count)++;
				w->mark = *count;
				nd = make_stack_node();
				nd->vert = w;
				nd->next = NULL;
				push(&st,nd);
				//printf("\nInside set mark  count = %d w->num = %d\n",*count,w->num);
				//scanf("%d",&ans);
			  }
			}
		}
		// remove top of stack
		ndd = pop(&st);
		val = (ndd->vert)->num;
		ct = (ndd->vert)->mark;
		printf("\nVertex = %d  count = %d\n",val,ct);
//		scanf("%d",&ans);
		free(ndd);
		empty = isStackEmpty(st);
//		if (!empty) print_stack(q);
//		printf("\nEnd of while   empty = %d",empty);
	}
}	

//DFS
void DFS(struct graph *g){
	int nv, count;
	int i;
	struct vertex *v;
	
	nv = g->num_verts;
	//mark all vertices 0
	for (i = 0; i<nv; i++) {
		v = g->vertices[i];
		v->mark = 0;
	}
	
	count = 0;
	
	for (i=0; i<nv; i++){
		v = g->vertices[i];
		if (v->mark == 0) dfs(g, v, &count);
	}
}

// Utility functions
double avg_strlen(int n, char **A)
{		int i;
		double sum;
		char x[100];
		
		sum = 0.0;
		for (i=0; i<n; i++)
		{
			strcpy(x,A[i]);
			sum = sum + (double) strlen(x);
		}
		return sum /(double) n;
}


//funcs for effClosestPair
//needed to sort array of points according to X coord
int compare_x(const void* a, const void* b) {
    struct pt2d *p1 = (struct pt2d *)a;
    struct pt2d *p2 = (struct pt2d *)b;
    return (p1->x - p2->x);
};
//needed to sort array of points according to Y coord
int compare_y(const void* a, const void* b) {
    struct pt2d *p1 = (struct pt2d *)a;
    struct pt2d *p2 = (struct pt2d *)b;
    return (p1->y - p2->y);
};

float cart_dist(struct pt2d p1, struct pt2d p2) {
    return sqrt(((p1.x - p2.x)*(p1.x - p2.x)) + ((p1.y - p2.y)*(p1.y - p2.y)));
}

// A Brute Force method to return the smallest distance between two points
// in P[] of size n
float bruteForce(struct pt2d* arr, int n) {
    float min = 100000000000;
    float distance;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            distance = cart_dist(arr[i], arr[j]);
            if (min > distance) {
                min = distance;
            }
        }
    }
    return min;
};

// A utility function to find minimum of two float values
float minFloat(float x, float y) {
    return (x < y)? x : y;
};

//finds the distance beween the closest points of center area
float center_closest(struct pt2d* center, int size, float min_center) {
    float min = min_center;

    //sorting center[] according to y coord.
    qsort(center, size, sizeof(struct pt2d), compare_y);

    // Pick all points one by one and try the next points till the difference
    // between y coordinates is smaller than min_closest
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; ((j < size) && ((center[j].y - center[i].y) < min)); j++) {
            if (cart_dist(center[i], center[j]) < min) {
                min = cart_dist(center[i], center[j]);
            }
        }
    }
    return min;
};

//arr has been sorted according to distance of x coords
float effClosestPair(struct pt2d* arr, int size) {
    //use if <3 points
    if (size <= 3) {
        return bruteForce(arr, size);
    }

    //find smallest distance on l and r of midpoint on array
    int mid = size / 2;
    struct pt2d midpt = arr[mid];
    float dist_l = effClosestPair(arr, mid);
    float dist_r = effClosestPair(arr + mid, size-mid);
    float min_center = min(dist_l, dist_r);

    //center contains points closer to the midpoint than min_center
    struct pt2d center[size];
    int count = 0;
    for (int i = 0; i < size; i++) {
        if (min_center > abs(arr[i].x - midpt.x)) {     //abs() returns absolute value
            center[count] = arr[i];
            count++;
        }
    }

    // find the closest points in center area
    return min(min_center, center_closest(center, count, min_center) );
}

//funcs for BST search and insert
struct node* search(struct node* root, int key) 
{ 
    if (root == NULL || root->key == key) 
       return root; 
     
    if (root->key < key) 
       return search(root->right, key); 
  
    return search(root->left, key); 
}

struct node *newNode(int item) 
{ 
    struct node *temp =  (struct node *)malloc(sizeof(struct node)); 
    temp->key = item; 
    temp->left = temp->right = NULL; 
    return temp; 
}

void traverse(struct node *root) 
{ 
    if (root != NULL) 
    { 
        traverse(root->left); 
        printf("%d \n", root->key); 
        traverse(root->right); 
    } 
}

struct node* insertNode(struct node* node, int key) 
{ 
    if (node == NULL) return newNode(key); 
  
    if (key < node->key) 
        node->left  = insertNode(node->left, key); 
    else if (key > node->key) 
        node->right = insertNode(node->right, key);    
  
    return node; 
}

int binarySearch(int arr[], int l, int r, int x) 
{ 
    if (r >= l) { 
        int mid = l + (r - l) / 2; 
  
        if (arr[mid] == x) 
            return mid; 
  
        if (arr[mid] > x) 
            return binarySearch(arr, l, mid - 1, x); 
  
        return binarySearch(arr, mid + 1, r, x); 
    } 
  
    return -1; 
}
