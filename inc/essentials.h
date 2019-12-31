/*
*       Parallels and Distributed Systems Exercise 3
*       Essentials functions found in different implementations of the program.
*       Author:Michael Karatzas
*       AEM:9137
*/

void pointer_swap(int **a , int **b){
  int * temp=*a;
  *a=*b;
  *b=temp;
}
