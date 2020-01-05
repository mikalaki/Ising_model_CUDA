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
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n, int *flag);
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,int colIndex, int * flag);

void ising( int *G, double *w, int k, int n){

  //The second Matrix We use
  int * secondG= (int *)malloc((size_t)sizeof(int)*n*n);
  if(!secondG){
    printf("Couldn't allocate memory!");
    exit(1);
  }

  //Flag for indicate if there was no changes in the lattice during a step,in order to terminate the evolving.
  int no_changes_flag;

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){
    //no_changes_flag=1, indicates nochange in the lattice, if there are changes , next function will update its value.
    no_changes_flag=1;

    nextStateCalculation(G,secondG,w,n,&no_changes_flag);

    //Swapping the pointers between the two Matrices
    pointer_swap(&G,&secondG);

    //If there are no changes in the lattice we stop evolving the model
    if(no_changes_flag){
      //getting the actual steps happened to handle the memcpy bellow
      k=i+1;
      break;
    }

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
//Function that scans the lattice
void nextStateCalculation(int *Gptr,int *newMat, double * w , int n,int * flag){

  //Scanning the Lattice
  for(int i=0; i<n; i++){
    for(int j=0;j<n;j++){
      getTheSpin(Gptr,newMat,w,n,i,j,flag);
    }
  }

}
//Function that computes the next spin of a certain point
void getTheSpin(int * Lat,int * newLat, double * weights , int n, int rowIndex,
  int colIndex, int *flag){
  // int rowIndex= index/n;
  // int colIndex= index%n;

  double total=0;
  int idxR,idxC;

  //Calculating the Total influence for a certain spot.
  for(int i=rowIndex-2;i<rowIndex+3;i++ ){
    for(int j=colIndex-2;j<colIndex+3;j++ ){
      if((i==rowIndex) && (j==colIndex))
        continue;

      //using modulus arithmetic for handle the boundaries' conditions
      //Getting the positive modulus
      idxR= (i  + n) % n ;
      idxC= (j  + n) % n ;

      //Total influence update
      total+=Lat[ idxR*n + idxC] *weights[(2+i-rowIndex)*5 + (2+j-colIndex)];
    }
  }

  //Checking the conditions in order to get the next state spin
  //if (total ==0), with taking into account possible floating point errors
  if( (total<1e-6)  &&  (total>(-1e-6)) ){
    newLat[rowIndex*n+colIndex]=Lat[rowIndex*n+colIndex];
  }
  else if(total<0){
    if(Lat[rowIndex*n+colIndex]!=-1)
      *flag=0;
    newLat[rowIndex*n+colIndex]=-1;
  }
  else if(total>0){
    if(Lat[rowIndex*n+colIndex]!=1)
      *flag=0;
    newLat[rowIndex*n+colIndex]=1;
  }

}
