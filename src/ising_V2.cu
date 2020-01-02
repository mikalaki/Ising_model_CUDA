/*
*       Parallels and Distributed Systems Exercise 3
*       v2. CUDA modified ising model ,grid and block computes the magnetic moments.
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ising.h"
#include "essentials.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
//The max threads per block for my gpu (gt 540m) is 1024 = 32*32 (1024 are run by a single processor)
//(Preferably:set BLOCK_DIM_X and BLOCK_DIM_Y a multiple of 4)
#define BLOCK_DIM_X 32
#define BLOCK_DIM_Y 32
#define GRID_DIM_X 9
#define GRID_DIM_Y 9

//Functions Declaration
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n);

__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex);



///Functions Definition
void ising( int *G, double *w, int k, int n){

  int * d_G, *d_secondG;
  double * d_w;

  //Allocate and Get the G Matrix in the Device
  if(   cudaMalloc((void **)&d_G, (size_t)sizeof(int)*n*n)     != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }
  cudaMemcpy(d_G, G, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);

  //Allocate and Get the Weights Matrix in the Device
  if(  cudaMalloc((void **)&d_w, (size_t)sizeof(double)*5*5)   != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }
  cudaMemcpy(d_w, w, (size_t)sizeof(double)*5*5, cudaMemcpyHostToDevice);

  //The second Matrix We use,allocation only in GPU(device)
  if(cudaMalloc((void **)&d_secondG, (size_t)sizeof(int)*n*n) != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }


  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){


    //grid in order one thread to compute a block of moments.
    //dim3 BlockDim(BLOCK_SIZE_1D,BLOCK_SIZE_1D);

    dim3 dimBlock(BLOCK_DIM_X,BLOCK_DIM_Y);
    dim3 dimGrid(GRID_DIM_X,GRID_DIM_Y);


    nextStateCalculation<<<dimGrid,dimBlock>>>(d_G,d_secondG,d_w,n);
    cudaDeviceSynchronize();

    //Swapping the pointers between the two Matrices in device
    pointer_swap(&d_G,&d_secondG);

    //Passing updated values of G matrix in the CPU.
    cudaMemcpy(G,d_G,(size_t)sizeof(int)*n*n,cudaMemcpyDeviceToHost);


  }

  //Freeing memory space I dont need from GPU to avoid memory leaks.
  cudaFree(d_G);
  cudaFree(d_secondG);
  cudaFree(d_w);

}
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n){
      int strideX = blockDim.x *gridDim.x ;
      int strideY = blockDim.y *gridDim.y ;
      int index_X = threadIdx.x +blockDim.x*blockIdx.x;
      int index_Y = threadIdx.y +blockDim.y*blockIdx.y;


      for(int i=index_Y;i<n ;i+=strideY){
        for(int j=index_X; j<n;j+=strideX){
          getTheSpin(Gptr,newMat,w,n,i,j);
        }
      }
}
__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex){


  double total=0;
  int idxR,idxC;
  //Getting the total for a certain spin.
  for(int i=rowIndex-2;i<rowIndex+3;i++ ){
    for(int j=colIndex-2;j<colIndex+3;j++ ){
      if((i==rowIndex) && (j==colIndex))
        continue;

      //using modulus arithmetic for handle the edges
      //Getting the modulus from the remainder in negative values of Cmodulus operator
      idxR= (i + n) % n ;
      idxC= (j + n) % n ;

      total+=Lat[ idxR*n + idxC] *weights[(2+i-rowIndex)*5 + (2+j-colIndex)];
    }
  }

  //Checking the conditions
  //if (total ==0), with taking into account possible floating point errors
  if( (total<1e-6)  &&  (total>(-1e-6)) ){
    newLat[rowIndex*n+colIndex]=Lat[rowIndex*n+colIndex];
  }
  else if(total<0){
    newLat[rowIndex*n+colIndex]=-1;
  }
  else if(total>0){
    newLat[rowIndex*n+colIndex]=1;
  }

}
