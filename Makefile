#      Parallels and Distributed Systems Exercise 3
#			 Makefile
#      Author:Michael Karatzas
#      AEM:9137

SHELL := /bin/bash
CC = gcc-5
#NVCC = nvcc  -gencode=arch=compute_20,code=sm_21 -Wno-deprecated-gpu-targets -lcudart  -Xptxas=-v#use this flag as I own an old gpu with 2.1 compute capability
NVCC =nvcc    -Xptxas=-v
CFLAGS =  #-Wall -O3
INCLUDES = -I ./inc

clean:
	find ./ -name "*.a" -o -name "*.o" -o -executable -a -type f | xargs rm -f



lib: ising_sequential.o ising_V1.o ising_V2.o
	ar rcs lib/ising_sequential.a lib/ising_sequential.o
	ar rcs lib/ising_V1.a lib/ising_V1.o #lib/knnring_sequential.o
	ar rcs lib/ising_V2.a lib/ising_V2.o
# #	#### COMMENT OUT THE LINE BELLOW TO RUN SYNCHRONOYS MPI (ln 19 must be commented) ####
# #	rm -f lib/knnring_mpi.a; ar rcs lib/knnring_mpi.a lib/knnring_synchronous.o lib/knnring_sequential.o
#
# #	#### COMMENT OUT THE LINE BELLOW TO RUN ASYNCHRONOYS MPI (ln 16 should be commented) ####
# 	rm -f lib/knnring_mpi.a ;ar rcs lib/knnring_mpi.a lib/knnring_asynchronous.o lib/knnring_sequential.o

	rm lib/*.o

# knnring_synchronous.o:
# 	$(MPICC) $(CFLAGS) $(INCLUDES) -c src/knnring_synchronous.c -o lib/knnring_synchronous.o
#
# knnring_asynchronous.o:
# 	$(MPICC) $(CFLAGS) $(INCLUDES) -c src/knnring_asynchronous.c -o lib/knnring_asynchronous.o

ising_sequential.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c src/ising_sequential.c -o lib/ising_sequential.o

ising_V1.o:
#	$(NVCC)  $(INCLUDES) -c src/ising_V1.cu -o lib/ising_V1.o
	$(NVCC)  $(INCLUDES)  -c src/ising_V1.cu -o lib/ising_V1.o

ising_V2.o:
#	$(NVCC)  $(INCLUDES) -c src/ising_V1.cu -o lib/ising_V1.o
	$(NVCC)  $(INCLUDES)  -c src/ising_V2.cu -o lib/ising_V2.o
