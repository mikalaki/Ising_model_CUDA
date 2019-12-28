/*
*       Parallels and Distributed Systems Exercise 3
*       v1. CUDA modified ising model ,each block use block threads shared memory.
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ising.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
//The max threads per block for my gpu (gt 540m) is 1024 = 32*32 (1024 are run by a single processor)
#define BLOCK_DIM_X 32
#define BLOCK_DIM_Y 32
//In my gpu there 2(MPs)*48(SPs)=96 sqrt(96)>9 => grid dimensions:
#define GRID_DIM_X 9
#define GRID_DIM_Y 9
#define RADIUS 2

//Functions Declaration
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n);

__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int lRowIndex,
  int lColIndex,int gRowIndex,int gColIndex);

// __host__
// void checkKernelConfiguration

void pointer_swap(int **a , int **b);


///Functions Definition
void ising( int *G, double *w, int k, int n){

  int * d_G, *secondG, *d_secondG;
  double * d_w;

  //Allocate and Get the G Matrix in the Device
  cudaMalloc((void **)&d_G, (size_t)sizeof(int)*n*n);
  cudaMemcpy(d_G, G, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);

  //Allocate and Get the Weights Matrix in the Device
  cudaMalloc((void **)&d_w, (size_t)sizeof(double)*5*5);
  cudaMemcpy(d_w, w, (size_t)sizeof(double)*5*5, cudaMemcpyHostToDevice);

  //The second Matrix We use,allocation in CPU(host) and GPU(device)
  secondG= (int *)malloc((size_t)sizeof(int)*n*n);
  cudaMalloc((void **)&d_secondG, (size_t)sizeof(int)*n*n);

  //check for valid Kernel Configuration

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){


    //block and grid dimensions
    dim3 dimBlock(BLOCK_DIM_X,BLOCK_DIM_Y);
    dim3 dimGrid(GRID_DIM_X,GRID_DIM_Y);
    //Check for valid kernel configuration

    nextStateCalculation<<<dimGrid,dimBlock>>>(d_G,d_secondG,d_w,n);

    cudaMemcpy(secondG,d_secondG,(size_t)sizeof(int)*n*n,cudaMemcpyDeviceToHost);
    // //checkErrors
    // printf("%s\n", cudaGetErrorString(cudaGetLastError()));

    //Swapping the pointers between the two Matrices.
    pointer_swap(&G,&secondG);
    //Update data in device after the pointer swap.
    cudaMemcpy(d_G, G, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_secondG, secondG, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);




  }

  //Getting the right values to the initial Lattice matrix for odd number of spots
  if((k%2)!=0){
    memcpy (secondG, G, (size_t)sizeof(int)*n*n);
    //Freeing memory space I dont need from CPU and GPU to avoid memory leaks in both.
    free(G);
    cudaFree(d_G);
    cudaFree(d_secondG);
  }
  else{
    //Freeing memory space I dont need from CPU and GPU to avoid memory leaks in both.
    free(secondG);
    cudaFree(d_G);
    cudaFree(d_secondG);
  }
}
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n){

      //I will share the part of the G matrix which is used for the computations
      __shared__ int sharedGpart[(BLOCK_DIM_X+2*RADIUS)  *  (BLOCK_DIM_Y+2*RADIUS)];
      int sharedNcols=(BLOCK_DIM_X+2*RADIUS) ;
      //__shared__ double w_shared[25];
      //Get the index and the stride for each Thread
      int strideX = blockDim.x *gridDim.x ;
      int strideY = blockDim.y *gridDim.y ;
      int gIndex_X = threadIdx.x +blockDim.x*blockIdx.x;//global x index
      int gIndex_Y = threadIdx.y +blockDim.y*blockIdx.y;//global y index

      int lIndex_X=threadIdx.x+RADIUS;//local(in the block) x index
      int lIndex_Y=threadIdx.y+RADIUS;//local(in the block) y index

      for(int i=gIndex_Y; i<n +RADIUS ;i+=strideY){
        for(int j=gIndex_X; j<n +RADIUS;j+=strideX){

          //Everythread initialize its point
          sharedGpart[lIndex_Y*sharedNcols+lIndex_X]=Gptr[( (i + n)%n )*n + ( (j + n)%n )];


          if((threadIdx.x)<RADIUS){
            int sharedGAccessorX= (lIndex_Y)*sharedNcols+(lIndex_X -RADIUS);
            int GAccessorX=( (i + n)%n )*n+ ( ( (j-RADIUS)%n  + n) % n);
            sharedGpart[sharedGAccessorX]=Gptr[GAccessorX];

            sharedGAccessorX=(lIndex_Y)*sharedNcols+(lIndex_X+BLOCK_DIM_X);
            GAccessorX=( (i + n)%n )*n+( ( (j+BLOCK_DIM_X)%n  + n) % n);
            sharedGpart[sharedGAccessorX]=Gptr[GAccessorX];

            //FOR DIAGONIALS
            if((threadIdx.y)<RADIUS){
              //1st diagonial 4spots square
              int sharedDiagAccessorX= (lIndex_Y -RADIUS)*sharedNcols +(lIndex_X-RADIUS);
              int GDiagAccessorX=( ( (i-RADIUS)%n  + n) % n)*n+( ( (j-RADIUS)%n  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];

              //2nd diagonial 4spots square
              sharedDiagAccessorX= (lIndex_Y+BLOCK_DIM_Y)*sharedNcols +(lIndex_X-RADIUS);
              GDiagAccessorX=( ( (i+BLOCK_DIM_Y)%n  + n) % n)*n+( ( (j-RADIUS)%n  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];

              //3rd diagonial 4spots square
              sharedDiagAccessorX= (lIndex_Y+BLOCK_DIM_Y)*sharedNcols +(lIndex_X+BLOCK_DIM_X);
              GDiagAccessorX=( ( (i+BLOCK_DIM_Y)%n  + n) % n)*n+( ( (j+BLOCK_DIM_X)%n  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];

              //4rd diagonial 4spots square
              sharedDiagAccessorX= (lIndex_Y -RADIUS)*sharedNcols+(lIndex_X+BLOCK_DIM_X);
              GDiagAccessorX=( ( (i-RADIUS)%n  + n) % n)*n+( ( (j+BLOCK_DIM_X)%n  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];
            }
          }

          if((threadIdx.y)<RADIUS){
            int sharedGAccessorY= (lIndex_Y-RADIUS)*sharedNcols+lIndex_X;
            int GAccessorY=( ( (i-RADIUS)%n  + n) % n)*n+( (j + n)%n );
            sharedGpart[sharedGAccessorY]=Gptr[GAccessorY];

            sharedGAccessorY=(lIndex_Y+BLOCK_DIM_Y)*sharedNcols+lIndex_X;
            GAccessorY=( ( (i+BLOCK_DIM_Y)%n  + n) % n)*n+( (j + n)%n );
            sharedGpart[sharedGAccessorY]=Gptr[GAccessorY];
          }
          //Here we synchronize the block threads in order Shared G values are
          //updated for each thread
          __syncthreads();

          if((i<n)&&(j<n))
            getTheSpin(sharedGpart,newMat,w,n,lIndex_Y, lIndex_X,i,j);

          __syncthreads();
          //return;
        }
      }

}

void getTheSpin(int * Lat,int * newLat, double * weights , int n, int lRowIndex,int lColIndex,
int gRowIndex,int gColIndex ){

  double total=0;
  //Getting the total for a certain spin.
  for(int i=lRowIndex-2;i<lRowIndex+3;i++ ){
    for(int j=lColIndex-2;j<lColIndex+3;j++ ){
      if((i==lRowIndex) && (j==lColIndex))
        continue;

      total+=Lat[ i*(BLOCK_DIM_X+2*RADIUS) + j] *weights[(2+i-lRowIndex)*5 + (2+j-lColIndex)];

    }
  }

  //Checking the conditions
  //if (total ==0), with taking into account possible floating point errors
  if( (total<1e-6)  &&  (total>(-1e-6)) ){
    newLat[(gRowIndex)*n+(gColIndex)]=Lat[lRowIndex*(BLOCK_DIM_X+2*RADIUS)+lColIndex];
  }
  else if(total<0){
    newLat[(gRowIndex)*n+(gColIndex)]=-1;
  }
  else if(total>0){
    newLat[(gRowIndex)*n+(gColIndex)]=1;
  }

}

void pointer_swap(int **a , int **b){
  int * temp=*a;
  *a=*b;
  *b=temp;
}
