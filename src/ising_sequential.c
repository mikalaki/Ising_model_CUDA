#include "stdio.h"
#include "stdlib.h"


void step(int * G, double * w , int n);
void evolve(int * Lat, double * weights , int n, int a , int b);

// static double weights[25]={ 0.004,  0.016,  0.026,  0.016,   0.004,
//                             0.016,  0.071,  0.117,  0.071,   0.016,
//                             0.026,  0.117,    0  ,  0.117,   0.026,
//                             0.016,  0.071,  0.117,  0.071,   0.016,
//                             0.004,  0.016,  0.026,  0.016,   0.004};
//! Ising model evolution
/*!

  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]

  NOTE: Both matrices G and w are stored in row-major format.
*/




void ising( int *G, double *w, int k, int n){
  // w= (double)malloc(25*sizeof(double));
  // w={ 0.004,  0.016,  0.026,  0.016,   0.004,
  //                             0.016,  0.071,  0.117,  0.071,   0.016,
  //                             0.026,  0.117,    0  ,  0.117,   0.026,
  //                             0.016,  0.071,  0.117,  0.071,   0.016,
  //                             0.004,  0.016,  0.026,  0.016,   0.004};

  //Evolving the model for k steps
  for(int i=0 ; i<k ;i++){
    step(G,w,n);
  }


}

void step(int * G, double * w , int n){

  // /*We will use step equals to 3 in the for loops which scan the G matrix(lattice)
  //   in order sequential controls aren't
  //   dependent on immediate previous-ones (we scan spins by 3 ex. if we check the spin
  //   in G[1][1] in one iteration, in the next iteration we will check the spin in
  //   G[4][4], which hasn't G[1][1] in its "weight window").In order to scan all the
  //   spins of the lattice we use offsets.
  //   */

  for(int c_offset =0 ; c_offset<3 ; c_offset++ ){
    for (int r_offset = 0; r_offset <3; r_offset++) {
      //scanning the Lattice
      for(int i=c_offset; i<n; i+=3){
        for(int j=r_offset; j<n; j+=3){
          evolve(G,w,n,i,j);
        }
      }
    }
  }


}

void evolve(int * Lat, double * weights , int n, int a , int b){
  double total=0;
  int idxR,idxC;
  //Getting the total for a certain spin.
  for(int i=a-2;i<a+3;i++ ){
    for(int j=b-2;j<b+3;j++ ){
      if((i==a) && (j==b))
        continue;
      idxR= ((i%n)<0 )? (i%n +n) : i%n;
      idxC= ((j%n)<0 )? (j%n +n) : j%n;

      total+=Lat[ idxR*n + idxC] *weights[(2+i-a)*5 + (2+j-b)];
    }
  }

  double dE=2*Lat[a*n+b]*total;
  if(dE <= 0){
    Lat[a*n+b]=Lat[a*n+b]*(-1);
  }
     // for i in range(n-1, n+2):
   //     for j in range(m-1, m+2):

   //         if i == n and j == m:
   //             continue
   //         total += field[i % N, j % M]
   // dE = 2 * field[n, m] * total
   // if dE <= 0:
   //     field[n, m] *= -1
   // elif np.exp(-dE * beta) > np.random.rand():
   //     field[n, m] *= -1
}
