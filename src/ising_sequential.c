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
#include "essentials.h"

//Functions Declaration
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n);
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

  //The second Matrix We use
  int * secondG= (int *)malloc((size_t)sizeof(int)*n*n);

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){
    nextStateCalculation(G,secondG,w,n);

    //Swapping the pointers between the two Matrices
    pointer_swap(&G,&secondG);
  }

  //Getting the right values to the initial Lattice matrix for odd number of spots
  //Freeing memory of the second matrix from the heap to avoid memory leaks.
  if((k%2)!=0){
    memcpy (secondG, G, (size_t)sizeof(int)*n*n);
    free(G);
  }
  else
    free(secondG);





}

void nextStateCalculation(int *Gptr,int *newMat, double * w , int n){

  //scanning the Lattice
  for(int i=0; i<n; i++){
    for(int j=0;j<n;j++){
      getTheSpin(Gptr,newMat,w,n,i,j);
    }
  }

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
      idxR= (i  + n) % n ;
      idxC= (j  + n) % n ;

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
