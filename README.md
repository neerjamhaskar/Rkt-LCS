
# OpenMPI Installation and Usage Guide

## Overview  
This guide provides instructions for installing and using OpenMPI on the CAS server, where users may not have access to `sudo`. There are two options for installing OpenMPI:  

1. **Request Installation:** Ask Derek to install OpenMPI system-wide.  
2. **Local Installation:** Install OpenMPI in your local directory.  

## Steps to Compile and Run  

### 1. Compile the Code  
Use the following command to compile the code with `mpicc` provided in the dependencies:  
```bash
./dependencies/bin/mpicc -o flouri_openmpi.o flouri_openmpi.c
```  

### 2. Run the Code on the Cluster  
Execute the compiled program using `mpirun` with the desired number of processes:  
```bash
./dependencies/bin/mpirun -np 4 ./flouri_openmpi.o fake_reads.fasta 3
```  
### 2. Debug runtime memory issues

```bash
valgrind --leak-check=full --track-origins=yes ./flouri_openmpi.o fake_reads.fasta 3
```

## Notes  
- Ensure that the paths to the `mpicc` and `mpirun` binaries in the `dependencies/bin` directory are correct.  
- Replace `fake_reads.fasta` and `3` with your own input file and parameters, if needed.
- Code for k-errata Trie can be found in: https://bitbucket.org/srirampc/bhavani/src/master/
