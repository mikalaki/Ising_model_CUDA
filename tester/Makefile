#      Parallels and Distributed Systems Exercise 3
#			 Makefile for tester
#      Author:Michael Karatzas
#      AEM:9137

SHELL := /bin/bash
#Set your system's Cuda installation's gcc version and the flags needed for your system.
#I use gcc-5 because is the gcc version my CUDA version(8.0) supports,in order all of the project's C code get compiled by the same compiler.
CC = gcc#-5
#NVCC = nvcc  -gencode=arch=compute_20,code=sm_21 -Wno-deprecated-gpu-targets -lcudart  -Xptxas=-v#use this flag as I own an old gpu with 2.1 compute capability
NVCC = nvcc


test_sequential:
	cd ../; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(CC) test.c ising_sequential.a -o $@
	./test_sequential
	rm test_sequential


test_v1:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) test_cuda.cu ising_V1.a -o $@ -lm
	#nvprof ./test_v1
	./test_v1
	rm test_v1

test_v2:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) test_cuda.cu ising_V2.a -o $@ -lm
	#nvprof ./test_v2
	./test_v2
	rm test_v2

test_v3:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) test_cuda.cu ising_V3.a -o $@ -lm
	#nvprof ./test_v3
	./test_v3
	rm test_v3
