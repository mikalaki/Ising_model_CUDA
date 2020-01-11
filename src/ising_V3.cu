/*
*       Parallels and Distributed Systems Exercise 3
*       v3. CUDA modified ising model evolution,each block use block threads' shared memory.
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
//The max threads per block for my gpu (gt 540m) is 1024 so it must be BLOCK_DIM_X* BLOCK_DIM_Y<=1024
//(Preferably:set BLOCK_DIM_X and BLOCK_DIM_Y a multiple of 4)
#define BLOCK_DIM_X 24
#define BLOCK_DIM_Y 24
#define GRID_DIM_X 4
#define GRID_DIM_Y 4
#define RADIUS 2

//Functions'-kernels' Declarations
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n, int * flag);

__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int lRowIndex,
  int lColIndex,int gRowIndex,int gColIndex, int * flag);



///Functions'-kernels' Definitions
void ising( int *G, double *w, int k, int n){

  //Flag for indicate if there was no changes in the lattice during a step,in order to terminate the evolving.
  int no_changes_flag;

  int * d_G,*d_secondG, *d_no_changes_flag;
  double * d_w;

  //Allocate memory for the no changes flag in the Device
  if(   cudaMalloc(&d_no_changes_flag, (size_t)sizeof(int))    != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }

  //Allocate memory and "transfer" the G Matrix in the Device
  if(   cudaMalloc((void **)&d_G, (size_t)sizeof(int)*n*n)     != cudaSuccess){
    printf("Couldn't allocate memory in device (GPU) !");
    exit(1);
  }
  cudaMemcpy(d_G, G, (size_t)sizeof(int)*n*n, cudaMemcpyHostToDevice);

  //Allocate memory and "transfer" the Weights Matrix in the Device
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

  //Setting block's and grid's dimensions
  dim3 dimBlock(BLOCK_DIM_X,BLOCK_DIM_Y);
  dim3 dimGrid(GRID_DIM_X,GRID_DIM_Y);


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

  //Freeing memory space I don't need from GPU to avoid memory leaks.
  cudaFree(d_G);
  cudaFree(d_secondG);
  cudaFree(d_w);

}
__global__
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n, int * flag){
      /* The part of the G matrix that is needed to be read in the block shared memory
      are the spots that their spin is going to get computed by the block
      and two "offset" spot around every edgy spot, that will be needed for the spin
      computation as they are neighbots of the edgy spots  */

      //The part of the G matrix that will pass in the shared memory.
      __shared__ int sharedGpart[(BLOCK_DIM_X+2*RADIUS)  *  (BLOCK_DIM_Y+2*RADIUS)];

      //The number of columns of the shared G part
      int sharedNcols=(BLOCK_DIM_X+2*RADIUS) ;

      //matrix to store the shared weight matrices
      __shared__ double w_shared[25];

      //The step of each thread
      int strideX = blockDim.x *gridDim.x ;
      int strideY = blockDim.y *gridDim.y ;

      //The unigue global indixes of the threads in the grid
      int gIndex_X = threadIdx.x +blockDim.x*blockIdx.x;//global x index
      int gIndex_Y = threadIdx.y +blockDim.y*blockIdx.y;//global y index

      //The local (in the block) Index
      int lIndex_X=threadIdx.x+RADIUS;//local(in the block) x index
      int lIndex_Y=threadIdx.y+RADIUS;//local(in the block) y index

      //Accessing the spins in the global lattice and pass them in the shared matrix.
      for(int i=gIndex_Y; i<n +RADIUS ;i+=strideY){
        for(int j=gIndex_X; j<n +RADIUS;j+=strideX){

          //Every thread read its own element in shared memory
          sharedGpart[lIndex_Y*sharedNcols+lIndex_X]=Gptr[( (i + n)%n )*n + ( (j + n)%n )];

          //Accessing and read read in shared memory the 2 left and 2 right "offset" elements on each row
          if((threadIdx.x)<RADIUS){
            int sharedGAccessorX= (lIndex_Y)*sharedNcols+(lIndex_X -RADIUS);
            int GAccessorX=( (i + n)%n )*n+ ( ( (j-RADIUS)  + n) % n);
            sharedGpart[sharedGAccessorX]=Gptr[GAccessorX];

            sharedGAccessorX=(lIndex_Y)*sharedNcols+(lIndex_X+BLOCK_DIM_X);
            GAccessorX=( (i + n)%n )*n+( ( (j+BLOCK_DIM_X)  + n) % n);
            sharedGpart[sharedGAccessorX]=Gptr[GAccessorX];

            //Accessing and read in shared memory "corner offset" elements(each corner has 4 elements)
            if((threadIdx.y)<RADIUS){
              //1st corner (4 points, up and left)
              int sharedDiagAccessorX= (lIndex_Y -RADIUS)*sharedNcols +(lIndex_X-RADIUS);
              int GDiagAccessorX=( ( (i-RADIUS)  + n) % n)*n+( ( (j-RADIUS)  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];

              //2nd diagonial (4 points, down and left)
              sharedDiagAccessorX= (lIndex_Y+BLOCK_DIM_Y)*sharedNcols +(lIndex_X-RADIUS);
              GDiagAccessorX=( ( (i+BLOCK_DIM_Y)  + n) % n)*n+( ( (j-RADIUS)  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];

              //3rd corner (4 points, down and right)
              sharedDiagAccessorX= (lIndex_Y+BLOCK_DIM_Y)*sharedNcols +(lIndex_X+BLOCK_DIM_X);
              GDiagAccessorX=( ( (i+BLOCK_DIM_Y)  + n) % n)*n+( ( (j+BLOCK_DIM_X)  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];

              //4rd diagonial (4 points, up and right)
              sharedDiagAccessorX= (lIndex_Y -RADIUS)*sharedNcols+(lIndex_X+BLOCK_DIM_X);
              GDiagAccessorX=( ( (i-RADIUS)  + n) % n)*n+( ( (j+BLOCK_DIM_X)  + n) % n);
              sharedGpart[sharedDiagAccessorX]=Gptr[GDiagAccessorX];
            }
          }

          //Accessing and read read in shared memory the 2 top and 2 bottom "offset" elements on each row
          if((threadIdx.y)<RADIUS){
            int sharedGAccessorY= (lIndex_Y-RADIUS)*sharedNcols+lIndex_X;
            int GAccessorY=( ( (i-RADIUS)  + n) % n)*n+( (j + n)%n );
            sharedGpart[sharedGAccessorY]=Gptr[GAccessorY];

            sharedGAccessorY=(lIndex_Y+BLOCK_DIM_Y)*sharedNcols+lIndex_X;
            GAccessorY=( ( (i+BLOCK_DIM_Y)  + n) % n)*n+( (j + n)%n );
            sharedGpart[sharedGAccessorY]=Gptr[GAccessorY];
          }

          /*Ιf (BLOCK_DIM_Y>=5) && (BLOCK_DIM_X>=5),we use shared memory  also for the weights matrix,
          I didn't implement it for smaller dimensions, because the benefit is very small anyway and it will
          make our code more complex .Also we choose BLOCK_DIM_X = BLOCK_DIM_Y =24 and BLOCK_DIM_X<5 or
          BLOCK_DIM_Y<5 aren't used in practice, so WE get the small benefit by transfering the weights' matrix.  */
          if((BLOCK_DIM_Y>=5) && (BLOCK_DIM_X>=5)){
            if(threadIdx.x<5 &&threadIdx.y<5)
              w_shared[threadIdx.x*5+ threadIdx.y]=w[threadIdx.x*5+ threadIdx.y];
          }



          //Here we synchronize the block threads in order Shared G values are
          //updated for each thread and w values are updated
          __syncthreads();

          if((i<n)&&(j<n)){
            if((BLOCK_DIM_Y>=5) && (BLOCK_DIM_X>=5))
              getTheSpin(sharedGpart,newMat,  w_shared,n,lIndex_Y, lIndex_X,i,j,flag);
            else //if((BLOCK_DIM_Y<5) && (BLOCK_DIM_X35))
              getTheSpin(sharedGpart,newMat,  w,n,lIndex_Y, lIndex_X,i,j,flag);

          }

          __syncthreads();

        }
      }

}
__device__ __forceinline__
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int lRowIndex,int lColIndex,
int gRowIndex,int gColIndex, int * flag ){

  double total=0;
  //Calculating the Total influence for a certain spot, by scanning the block shared part of G.
  for(int i=lRowIndex-2;i<lRowIndex+3;i++ ){
    for(int j=lColIndex-2;j<lColIndex+3;j++ ){
      if((i==lRowIndex) && (j==lColIndex))
        continue;

      //Total influence update
      total+=Lat[ i*(BLOCK_DIM_X+2*RADIUS) + j] *weights[(2+i-lRowIndex)*5 + (2+j-lColIndex)];

    }
  }

  //Checking the conditions in order to get the next state spin
  //  if (total ==0), with taking into account possible floating point errors
  if( (total<1e-6)  &&  (total>(-1e-6)) ){
    newLat[(gRowIndex)*n+(gColIndex)]=Lat[lRowIndex*(BLOCK_DIM_X+2*RADIUS)+lColIndex];
  }
  //if change in a certain spot happens we update no change flag's value into 0.
  else if(total<0){
    //Checking if there is change in this certain spot
    if(Lat[lRowIndex*(BLOCK_DIM_X+2*RADIUS)+lColIndex]!=1)
      *flag=0;
    newLat[(gRowIndex)*n+(gColIndex)]=-1;
  }
  else if(total>0){
    //Checking if there is change in this certain spot
    if(Lat[lRowIndex*(BLOCK_DIM_X+2*RADIUS)+lColIndex]!=1)
      *flag=0;
    newLat[(gRowIndex)*n+(gColIndex)]=1;
  }

}
