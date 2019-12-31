#      Parallels and Distributed Systems Exercise 3
#			 Makefile
#      Author:Michael Karatzas
#      AEM:9137

SHELL := /bin/bash
CC = gcc#-5
#Set your system's Cuda installation's gcc version and the flags needed for your system.
#NVCC = nvcc  -gencode=arch=compute_20,code=sm_21 -Wno-deprecated-gpu-targets -lcudart  -Xptxas=-v#use this flags as I own an old gpu with 2.1 compute capability
NVCC =nvcc   # -Xptxas=-v
CFLAGS = #-Wall -O3
INCLUDES = -I ./inc

clean:
	find ./ -name "*.a" -o -name "*.o" -o -executable -a -type f | xargs rm -f



lib: ising_sequential.o ising_V1.o ising_V2.o ising_V3.o
	ar rcs lib/ising_sequential.a lib/ising_sequential.o
	ar rcs lib/ising_V1.a lib/ising_V1.o #lib/knnring_sequential.o
	ar rcs lib/ising_V2.a lib/ising_V2.o
	ar rcs lib/ising_V3.a lib/ising_V3.o

	rm lib/*.o


ising_sequential.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c src/ising_sequential.c -o lib/ising_sequential.o -pg

ising_V1.o:
#	$(NVCC)  $(INCLUDES) -c src/ising_V1.cu -o lib/ising_V1.o
	$(NVCC)  $(INCLUDES)  -c src/ising_V1.cu -o lib/ising_V1.o

ising_V2.o:
#	$(NVCC)  $(INCLUDES) -c src/ising_V1.cu -o lib/ising_V1.o
	$(NVCC)  $(INCLUDES)  -c src/ising_V2.cu -o lib/ising_V2.o

ising_V3.o:
#	$(NVCC)  $(INCLUDES) -c src/ising_V1.cu -o lib/ising_V1.o
	$(NVCC)  $(INCLUDES)  -c src/ising_V3.cu -o lib/ising_V3.o
