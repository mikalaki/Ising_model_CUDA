/*
*       Tester for Parallels and Distributed Systems Exercise 3
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ising.h"
#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

int test(int *G1,int *G2, int n );

int main(int argc, char const *argv[]) {
  double weights[25]={ 0.004,  0.016,  0.026,  0.016,   0.004,
                       0.016,  0.071,  0.117,  0.071,   0.016,
                       0.026,  0.117,    0  ,  0.117,   0.026,
                       0.016,  0.071,  0.117,  0.071,   0.016,
                       0.004,  0.016,  0.026,  0.016,   0.004};
  // int n=atoi(arv[1]);
  // int k=atoi(arv[2]);
  int n=517;
  //int k=1;

  //Getting the initial situation of the lattice
  int *init_buffer =(int *)malloc((size_t)sizeof(int)*n*n);
  FILE * fp_init = fopen("conf-init.bin", "rb");
  fread(init_buffer,sizeof(int),n*n,fp_init); // read bytes to our buffer

  //Coppy the initial situation
  int *sit0 =(int *)malloc((size_t)sizeof(int)*n*n);
  memcpy (sit0, init_buffer, (size_t)sizeof(int)*n*n);

  //test for k=1 // In next version will implement it as a new function
  memcpy (sit0, init_buffer, (size_t)sizeof(int)*n*n);

  ising(sit0, weights,1, n);

  int *sit1 =(int *)malloc((size_t)sizeof(int)*n*n);
  FILE * fp_1 = fopen("conf-1.bin", "rb");

  fread(sit1,sizeof(int),n*n,fp_1); // read bytes to our buffer
  int test1=test(sit0,sit1,n);
  if(test1){
    printf("For k=1 correct \n");
  }
  else{
    printf("For k=1 wrong \n");
  }
  free(sit1);

  //Test for k=4 // In next version will implement it as a new function
  memcpy (sit0, init_buffer, (size_t)sizeof(int)*n*n);
  ising(sit0, weights,4, n);

  int *sit4 =(int *)malloc((size_t)sizeof(int)*n*n);
  FILE * fp_4 = fopen("conf-4.bin", "rb");
  fread(sit4,sizeof(int),n*n,fp_4); // read bytes to our buffer
  int test4=test(sit0,sit4,n);
  if(test4){
    printf("For k=4 correct \n");
  }
  else{
    printf("For k=4 wrong \n");
  }
  free(sit4);


  //Test for k=11 // In next version will implement it as a new function
  memcpy (sit0, init_buffer, (size_t)sizeof(int)*n*n);
  ising(sit0, weights,11, n);

  int *sit11 =(int *)malloc((size_t)sizeof(int)*n*n);
  FILE * fp_11 = fopen("conf-11.bin", "rb");
  fread(sit11,sizeof(int),n*n,fp_11); // read bytes to our buffer
  int test11=test(sit0,sit11,n);
  if(test11){
    printf("For k=11 correct \n");
  }
  else{
    printf("For k=11 wrong \n");
  }
  free(sit11);



}

int test(int *G1,int *G2, int n ){
  for(int i=0;i<n*n;i++){
    if(G1[i]!=G2[i]){
      return 0;
    }
  }

  return 1;
}
