/*
*       Parallels and Distributed Systems Exercise 3
*       v1. CUDA modified ising model evolution ,one thread computes a magnetic moment.
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "essentials.h"
#include "ising.h"
#include "cuda.h"
/*Setting the dimensions of the block,grid dimensions are computed int the programm
in order to have one thread per moment */
#define BLOCK_DIM_X 16
#define BLOCK_DIM_Y 16

//Functions'-Kernels' Declaration
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n, int * flag);
__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,
  int colIndex, int * flag);


///Functions'-kernels' Definitions
void ising( int *G, double *w, int k, int n){

  //Flag for indicate if there was no changes in the lattice during a step,in order to terminate the evolving.
  int no_changes_flag;

  int * d_G, *d_secondG, *d_no_changes_flag;
  double * d_w;


  //Allocate memory for the no changes flag in the Device
  if(   cudaMalloc(&d_no_changes_flag, (size_t)sizeof(int))    != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }


  //Allocate memory and "transfer" the given G Matrix in the Device
  cudaMalloc((void **)&d_G, (size_t)sizeof(int)*n*n);
  if(   cudaMalloc((void **)&d_G, (size_t)sizeof(int)*n*n)     != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }
  cudaMemcpy(d_G, G, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);

  //Allocate memory and "transfer" the Weights' Matrix in the Device
  if(  cudaMalloc((void **)&d_w, (size_t)sizeof(double)*5*5)   != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }
  cudaMemcpy(d_w, w, (size_t)sizeof(double)*5*5, cudaMemcpyHostToDevice);

  //Allocate memory for the second G matrix only in GPU(device)
  if(cudaMalloc((void **)&d_secondG, (size_t)sizeof(int)*n*n) != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }

  /*grid that matches the ising Model (gridDimX*BLOCK_DIM_X <=n,gridDimY*BLOCK_DIM_Y <=n )
  configuration for 1 thread per moment.*/

  //Calculating the grid dimensions, in order to match the ising model.
  int gridDimX= (n+BLOCK_DIM_X -1)/BLOCK_DIM_X;
  int gridDimY= (n+BLOCK_DIM_Y -1)/BLOCK_DIM_Y;

  //Setting grid and block dimensions
  dim3 dimGrid(gridDimX,gridDimY);
  dim3 dimBlock(BLOCK_DIM_X,BLOCK_DIM_Y);


  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){

    /*no_changes_flag=1, indicates no change in the lattice, if there are changes
    nextStateCalculation() kernel will update its value.*/
    no_changes_flag=1;
    cudaMemcpy(d_no_changes_flag, &no_changes_flag, (size_t)sizeof(int), cudaMemcpyHostToDevice);

    //calling the nextStateCalculation() kernel
    nextStateCalculation<<<dimGrid,dimBlock>>>(d_G,d_secondG,d_w,n,d_no_changes_flag);
    cudaDeviceSynchronize();

    //Swapping the pointers between the two Matrices in device
    pointer_swap(&d_G,&d_secondG);

    //The host get the value of the no changes flag as indication if no changes happened during the step.
    cudaMemcpy(&no_changes_flag, d_no_changes_flag,  (size_t)sizeof(int), cudaMemcpyDeviceToHost);
    //If there are no changes in the lattice we stop evolving the model
    if(no_changes_flag){
      break;
    }


  }
  //Passing updated values of G matrix in the host(CPU).
  cudaMemcpy(G,d_G,(size_t)sizeof(int)*n*n,cudaMemcpyDeviceToHost);


  //Freeing memory space I dont need from GPU to avoid memory leaks.
  cudaFree(d_G);
  cudaFree(d_secondG);
  cudaFree(d_w);

}

__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n, int * flag){
    //The unigue global indixes of the threads.
    int index_X = threadIdx.x +blockDim.x*blockIdx.x;
    int index_Y = threadIdx.y +blockDim.y*blockIdx.y;

    //getting rid of the odd trheads
    if((index_X < n) &&(index_Y < n) ){

      //Get each thread calcute spin of its spot
      getTheSpin(Gptr,newMat,w,n,index_X,index_Y, flag);
    }
}
__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex, int * flag){


  double total=0;
  int idxR,idxC;

  //Calculating the Total influence for a certain spot.
  for(int i=rowIndex-2;i<rowIndex+3;i++ ){
    for(int j=colIndex-2;j<colIndex+3;j++ ){
      if((i==rowIndex) && (j==colIndex))
        continue;

      //using modulus arithmetic for handle the boundaries' conditions
      //Getting the positive modulus
      idxR= (i + n) % n ;
      idxC= (j + n) % n ;

      //Total influence update
      total+=Lat[ idxR*n + idxC] *weights[(2+i-rowIndex)*5 + (2+j-colIndex)];
    }
  }

  //Checking the conditions in order to get the next state spin
  //if (total ==0), taking into account possible floating point errors
  if( (total<1e-6)  &&  (total>(-1e-6)) ){
    newLat[rowIndex*n+colIndex]=Lat[rowIndex*n+colIndex];
  }
  //if change in a certain spot happens we update no changes flag's value into 0.
  else if(total<0){
    //Checking if there is change in this certain spot
    if(Lat[rowIndex*n+colIndex]!=-1)
      *flag=0;
    newLat[rowIndex*n+colIndex]=-1;
  }
  else if(total>0){
    //Checking if there is change in this certain spot
    if(Lat[rowIndex*n+colIndex]!=1)
      *flag=0;
    newLat[rowIndex*n+colIndex]=1;
  }

}
