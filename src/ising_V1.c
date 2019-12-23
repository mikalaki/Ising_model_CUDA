/*
*       Parallels and Distributed Systems Exercise 3
*       v0. Sequential version of isinv_V1 Model
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include "ising.h"

//Functions Declaration
__global__
void step_V1(int **Gptr,int **newMat, double * w , int n){

  //scanning the Lattice
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if((threadIdx.x)==i*n+j)
        getTheSpin(*Gptr,*newMat,w,n,i,j);
    }
  }

  //Swapping the pointers between the two Matrices
  pointer_swap(Gptr,newMat);



}

//! Ising model evolution
/*!

  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]

  NOTE: Both matrices G and w are stored in row-major format.
*/
void ising_V1( int *G, double *w, int k, int n){

  //The second Matrix We use
  int * newG= (int *)malloc((size_t)sizeof(int)*n*n);

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){
    step_V1<<<1,n*n>>>(&G,&newG,w,n);
  }

  //Getting the right values to the initial Lattice matrix for odd number of spots
  //Freeing memory of the second matrix from the heap to avoid memory leaks.
  if((k%2)!=0){
    memcpy (newG, G, (size_t)sizeof(int)*n*n);
    free(G);
  }
  else
    free(newG);


}
