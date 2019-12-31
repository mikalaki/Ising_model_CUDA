#      Parallels and Distributed Systems Exercise 3
#			 Makefile
#      Author:Michael Karatzas
#      AEM:9137

SHELL := /bin/bash
#Set your system's Cuda installation's gcc version and the flags needed for your system.
CC = gcc-5 -pg
#NVCC = nvcc  -gencode=arch=compute_20,code=sm_21 -Wno-deprecated-gpu-targets -lcudart  -Xptxas=-v#use this flag as I own an old gpu with 2.1 compute capability
NVCC = nvcc   # -Xptxas=-v #use this flag as I own an old gpu with 2.1 compute capability


bench_sequential:
	cd ../; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(CC) bench.c ising_sequential.a -o $@
	./bench_sequential | tee stdoutV0.txt
	rm bench_sequential


bench_v1:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) bench.cu ising_V1.a -o $@ -lm
	#nvprof ./bench_v1
	nvprof ./bench_v1 | tee stdoutV1.txt
	rm bench_v1

bench_v2:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) bench.cu ising_V2.a -o $@ -lm
	#nvprof ./bench_v2
	nvprof ./bench_v2 | tee stdoutV2.txt
	rm bench_v2

bench_v3:
	cd ..; make clean; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) bench.cu ising_V3.a -o $@ -lm
	#nvprof ./bench_v3
	nvprof ./bench_v3 | tee stdoutV3.txt
	rm bench_v3
