CC = gcc
NVCC = nvcc
CFLAGS = -O3 -Wall
NVCCFLAGS = -O3
CUDA_PATH = /usr/local/cuda
CUDA_INCLUDE = -I$(CUDA_PATH)/include
CUDA_LIB = -L$(CUDA_PATH)/lib64 -lcudart
MPICC = ~/opt/openmpi/bin/mpicc
MPICXX = ~/opt/openmpi/bin/mpicxx
MPI_INCLUDE = -I$(HOME)/opt/openmpi/include

# Updated NVCCFLAGS with MPI include path
NVCCFLAGS = -O3 $(MPI_INCLUDE)

# Main targets
all: flouri_cpu flouri_cuda benchmark flouri_openmpi_cpu flouri_openmpi_gpu

# CPU implementation
flouri_cpu: flouri_cpu.o
	$(CC) $(CFLAGS) -o $@ $^

flouri_cpu.o: flouri_cpu.c flouri_cpu.h
	$(CC) $(CFLAGS) -c $<

# Object file with CPU functions but without main
flouri_cpu_lib.o: flouri_cpu.c flouri_cpu.h
	$(CC) $(CFLAGS) -c -DSKIP_MAIN -o $@ $<

# CUDA implementation
flouri_cuda: flouri_cuda.o flouri_cpu.o
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(CUDA_LIB)

flouri_cuda.o: flouri_cuda.cu flouri_cpu.h
	$(NVCC) $(NVCCFLAGS) $(CUDA_INCLUDE) -c $<

# Benchmark
benchmark: benchmark.o flouri_cuda.o flouri_cpu_lib.o
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(CUDA_LIB)

benchmark.o: benchmark.cu flouri_cpu.h
	$(NVCC) $(NVCCFLAGS) $(CUDA_INCLUDE) -c $<

# OpenMPI CPU implementation
flouri_openmpi_cpu: flouri_openmpi_cpu.c flouri_cpu_lib.o
	$(MPICC) -o $@ $^ -Wall -O3

# OpenMPI GPU implementation
flouri_openmpi_gpu.o: flouri_openmpi_gpu.cu flouri_cpu.h
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

flouri_openmpi_gpu: flouri_openmpi_gpu.o flouri_cuda.o flouri_cpu_lib.o
	$(MPICXX) -o $@ $^ $(CUDA_LIB)

# Utility targets
clean:
	rm -f *.o flouri_cpu flouri_cuda benchmark flouri_cpu_lib.o flouri_openmpi_cpu flouri_openmpi_gpu

.PHONY: all clean 