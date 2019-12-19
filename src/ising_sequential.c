/*
*       Parallels and Distributed Systems Exercise 3
*       v0. Sequential version of Ising Model
*       Author:Michael Karatzas
*       AEM:9137
*/
#include <stdio.h>
#include <stdlib.h>
#include "ising.h"

//Functions Declaration
void step(int **Gptr,int **newMat, double * w , int n);
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int a , int b);
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

  //The second Matrix We use
  int * newG= (int *)malloc((size_t)sizeof(int)*n*n);

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){
    step(&G,&newG,w,n);
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

void step(int **Gptr,int **newMat, double * w , int n){

  //scanning the Lattice
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      getTheSpin(*Gptr,*newMat,w,n,i,j);
    }
  }

  //Swapping the pointers between the two Matrices
  pointer_swap(Gptr,newMat);



}

void getTheSpin(int * Lat,int * newLat, double * weights , int n, int a , int b){
  double total=0;
  int idxR,idxC;
  //Getting the total for a certain spin.
  for(int i=a-2;i<a+3;i++ ){
    for(int j=b-2;j<b+3;j++ ){
      if((i==a) && (j==b))
        continue;

      //using modulus arithmetic for handle the edges
      //Getting the modulus from the remainder in negative values of Cmodulus operator
      idxR= ((i%n)<0 )? (i%n +n) : i%n;
      idxC= ((j%n)<0 )? (j%n +n) : j%n;

      total+=Lat[ idxR*n + idxC] *weights[(2+i-a)*5 + (2+j-b)];
    }
  }

  //Checking the conditions
  //if (total ==0), with taking into account possible floating point errors
  if( (total<1e-6)  &&  (total>(-1e-6)) ){
    newLat[a*n+b]=Lat[a*n+b];
  }
  else if(total<0){
    newLat[a*n+b]=-1;
  }
  else if(total>0){
    newLat[a*n+b]=1;
  }

}

void pointer_swap(int **a , int **b){
  int * temp=*a;
  *a=*b;
  *b=temp;
}
