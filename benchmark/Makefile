#      Parallels and Distributed Systems Exercise 3
#			 Makefile for benchmarks
#      Author:Michael Karatzas
#      AEM:9137

SHELL := /bin/bash
#Set your system's Cuda installation's gcc version and the flags needed for your system.
#I use gcc-5 because is the gcc version my CUDA version(8.0) supports,in order all of the project's C code get compiled by the same compiler.
CC = gcc #-5 -pg #
#NVCC = nvcc  -gencode=arch=compute_20,code=sm_21 -Wno-deprecated-gpu-targets -lcudart  -Xptxas=-v#use this flag as I own an old gpu with 2.1 compute capability
NVCC = nvcc   # -Xptxas=-v #use this flag as I own an old gpu with 2.1 compute capability


bench_sequential:
	cd ../; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(CC) bench.c ising_sequential.a -o $@ -lm
	./bench_sequential
	rm bench_sequential


bench_v1:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) bench.cu ising_V1.a -o $@ -lm
	nvprof ./bench_v1
	rm bench_v1

bench_v2:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) bench.cu ising_V2.a -o $@ -lm
	nvprof ./bench_v2
	rm bench_v2

bench_v3:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) bench.cu ising_V3.a -o $@ -lm
	nvprof ./bench_v3
	rm bench_v3
