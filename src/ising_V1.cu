/*
*       Parallels and Distributed Systems Exercise 3
*       v1. CUDA modified ising model ,one thread computes a magnetic moment.
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "essentials.h"
#include "ising.h"
#include "cuda.h"

//Functions Declaration
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n);
__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex);

//! Ising model evolution
/*!

  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]

  NOTE: Both matrices G and w are stored in row-major format.
*/
void ising( int *G, double *w, int k, int n){

  int * d_G, *d_secondG;
  double * d_w;

  //Allocate and Get the G Matrix in the Device
  cudaMalloc((void **)&d_G, (size_t)sizeof(int)*n*n);
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
    //grid that matches the ising Model
    dim3 dimGrid(n,n);
    nextStateCalculation<<<dimGrid,1>>>(d_G,d_secondG,d_w,n);
    cudaDeviceSynchronize();

    //Swapping the pointers between the two Matrices in device
    pointer_swap(&d_G,&d_secondG);

    //Passing updated values of G matrix in the CPU
    cudaMemcpy(G,d_G,(size_t)sizeof(int)*n*n,cudaMemcpyDeviceToHost);


  }

  //Freeing memory space I dont need from GPU to avoid memory leaks.
  cudaFree(d_G);
  cudaFree(d_secondG);
  cudaFree(d_w);

}

__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n){
      getTheSpin(Gptr,newMat,w,n,blockIdx.y,blockIdx.x);
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
