SHELL := /bin/bash
CC = gcc-7
NVCC = nvcc


test_sequential:
	cd ../; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(CC) test.c ising_sequential.a -o $@
	./test_sequential
	rm test_sequential


test_v1:
	cd ..; make lib;
	cp ../lib/*.a ../inc/ising.h ./;
	$(NVCC) test_v1.c ising_V1.a -o $@ -lm
	./test_v1
	rm test_v1
