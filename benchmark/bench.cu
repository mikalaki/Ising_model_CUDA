/*
*       Βenchmark code for Parallels and Distributed Systems Exercise 3
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "ising.h"
#include  <time.h>

struct timespec start, finish;
double elapsed;

static int n[]={100,250,500,750,1000,1500,2000,3000,5000};//9
static int k[]={1, 5, 10 , 25 ,40, 50 , 100 };//7
// static int n[]={1000};
// static int k[]={20};

int main(int argc, char const *argv[]) {

  printf("/----------------------CUDA BENCHMARK BEGINS----------------------\n");
  printf("/-----------------------------------------------------------------\n");
  double weights[25]={ 0.004,  0.016,  0.026,  0.016,   0.004,
                       0.016,  0.071,  0.117,  0.071,   0.016,
                       0.026,  0.117,    0  ,  0.117,   0.026,
                       0.016,  0.071,  0.117,  0.071,   0.016,
                       0.004,  0.016,  0.026,  0.016,   0.004};





  // clock_t t;
  // double time_taken;

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 7; j++) {
      FILE *pointerToFile;
      int * sample;

      sample=(int *)malloc((size_t)n[i]*n[i]*sizeof(int));
      for (size_t ii = 0; ii < n[i]*n[i]; ii++) {
        sample[ii]=(rand()%2 );
        if(sample[ii]==0){
          sample[ii]=-1;
        }
      }


      clock_gettime(CLOCK_MONOTONIC, &start);
      ising(sample, weights,k[j], n[i]);
      clock_gettime(CLOCK_MONOTONIC, &finish);
      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

      pointerToFile=fopen("datav3.csv","a");
      fprintf(pointerToFile,"%d,%d,%lf\n",n[i],k[j],elapsed);
      printf("Ising model evolution for n=%d, k=%d ,took %lf seconds! \n",n[i],k[j], elapsed );
      free(sample);
    }

  }
  printf("\n");
  printf("/-----------------------------------------------------------------\n");
  printf("/-----------------------CUDA BENCHMARK ENDS-----------------------\n");
  printf("/-----------------------------------------------------------------\n");


}
