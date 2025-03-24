
## Overview  
This code provides the implementation of LCP_k table by Flouri in O(n^2/p) time using OpenMPI in two settings: (1) CPU-only (2) CPU-GPU
## Steps to Compile and Run  

### 1. Compile the Code  

```bash
gcc -c flouri_cpu.c -o flouri_cpu.o

nvcc -c flouri_cuda.cu -o flouri_cuda.o

mpicc -o flouri_openmpi \
      flouri_openmpi.c \
      flouri_cpu.o \
      flouri_cuda.o \
      -I/usr/local/cuda/include \
      -L/usr/local/cuda/lib64 \
      -lcudart
```

### 2. Run the Code on the Cluster  

```bash
time mpirun -np 4 ./flouri_openmpi data/10.fasta 3 -c
```
-c refers to CPU-GPU computation

