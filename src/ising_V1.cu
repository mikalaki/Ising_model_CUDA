/*
*       Parallels and Distributed Systems Exercise 3
*       v0. Sequential version of Ising Model
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ising.h"
#include "cuda.h"

//Functions Declaration
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n);
__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex);

void pointer_swap(int **a , int **b);

//! Ising model evolution
/*!

  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]

  NOTE: Both matrices G and w are stored in row-major format.
*/
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

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){
    //grid that matches the ising Model
    dim3 dimGrid(n,n);
    nextStateCalculation<<<dimGrid,1>>>(d_G,d_secondG,d_w,n);
    cudaMemcpy(secondG,d_secondG,(size_t)sizeof(int)*n*n,cudaMemcpyDeviceToHost);

    //Swapping the pointers between the two Matrices
    pointer_swap(&G,&secondG);
    //Update data in device
    cudaMemcpy(d_G, G, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_secondG, secondG, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);
  }

  //Getting the right values to the initial Lattice matrix for odd number of spots
  if((k%2)!=0){
    memcpy (secondG, G, (size_t)sizeof(int)*n*n);
    //Freeing memory space I dont need from CPU and GPU to avoid memory leaks.
    free(G);
    cudaFree(d_G);
    cudaFree(d_secondG);
  }
  else{
    //Freeing memory space I dont need from CPU and GPU to avoid memory leaks.
    free(secondG);
    cudaFree(d_G);
    cudaFree(d_secondG);
  }




}

__global__ 
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n){
      getTheSpin(Gptr,newMat,w,n,blockIdx.x,blockIdx.y);
}

void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex){
  // int rowIndex= index/n;
  // int colIndex= index%n;

  double total=0;
  int idxR,idxC;
  //Getting the total for a certain spin.
  for(int i=rowIndex-2;i<rowIndex+3;i++ ){
    for(int j=colIndex-2;j<colIndex+3;j++ ){
      if((i==rowIndex) && (j==colIndex))
        continue;

      //using modulus arithmetic for handle the edges
      //Getting the modulus from the remainder in negative values of Cmodulus operator
      idxR= (i % n + n) % n ;
      idxC= (j % n + n) % n ;

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

void pointer_swap(int **a , int **b){
  int * temp=*a;
  *a=*b;
  *b=temp;
}
