SHELL := /bin/bash
CC = gcc-7
NVCC = nvcc
# MPICC = mpicc
# MPIRUN = mpirun --oversubscribe -np 4

test_sequential:
	cd ../; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(CC) test.c ising_sequential.a -o $@
	./test_sequential
	rm test_sequential


test_mpi:
	cd ..; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) tester_mpi.c ising_V1.a -o $@ -lm 
	./test_mpi
	rm test_mpi
