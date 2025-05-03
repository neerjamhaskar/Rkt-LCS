# Longest Common Substring (LCS) Implementation

A high-performance implementation of the Approximate Longest Common Substring (ALCS) algorithm of multiple strings with support for k-mismatches. This project provides implementations for CPU (sequential), CUDA (GPU), and OpenMPI (parallel CPU and GPU).

## Features

- Efficient LCS computation with k-mismatches support
- Multiple implementation options:
  - Sequential CPU implementation
  - CUDA-accelerated GPU implementation
  - MPI-based distributed CPU implementation
  - Hybrid MPI+CUDA implementation
- Benchmarking tools for performance comparison
- FASTA file format support for biological sequence data
- Configurable parameters:
  - k: Maximum number of mismatches allowed
  - t: Minimum number of strings which have ALCS occurances
  - tau: Minimum LCS length threshold

## Requirements

- C/C++ compiler (GCC recommended)
- CUDA Toolkit (for GPU implementations)
- OpenMPI (for distributed implementations)
- Make

## Building

Use the provided Makefile to build all implementations:

```bash
make all
```

Or build specific implementations:

```bash
make flouri_cpu        # CPU-only implementation
make flouri_cuda       # CUDA implementation
make flouri_openmpi_cpu  # OpenMPI CPU implementation 
make flouri_openmpi_gpu  # OpenMPI + CUDA implementation
make benchmark         # Benchmarking tool
```

## Usage

### Sequential CPU Version

```bash
./flouri_cpu <fasta_file> <candidate_num> <k> <t> <tau>
```

### CUDA GPU Version

```bash
./flouri_cuda <fasta_file> <candidate_num> <k> <t> <tau>

```

### MPI CPU Version

```bash
mpirun -np <num_processes> ./flouri_openmpi_cpu <fasta_file> <k> <t> <tau>
time mpirun -np 32 ./flouri_openmpi_cpu data/10000.fasta 2 4 5
```

### MPI GPU Version

```bash
mpirun -np <num_processes> ./flouri_openmpi_gpu <fasta_file> <k> <t> <tau>
time mpirun -np 32 ./flouri_openmpi_gpu data/10000.fasta 2 4 5
```

### Benchmarking Tool

```bash
./benchmark <fasta_file> <k> <t> <tau> [verbose]
```

### Parameters:

- `fasta_file`: Input FASTA file containing DNA sequences
- `k`: Maximum number of mismatches allowed
- `t`: Minimum number of strings which have ALCS occurances
- `tau`: Minimum match length
- `verbose`: (Optional) Display detailed output

## Performance

The CUDA implementation offers significant speedup over the sequential CPU version, especially for large datasets. Benchmark results can be found in the `benchmark_results.csv` file.

## Algorithm Details

The implementation is based on the MaxLCP_k algorithm. The core components include:

1. **MaxLCP_k Computation**: Calculating the longest common prefix with up to k mismatches
2. **Multi-Sequence Search**: Finding common substrings across multiple sequences
3. **Rkt-LCS Algorithm**: Given finite integers $k, t, m \in \mathbb{N}$, $t \leq m$, and a set $\s{S} = \{s_1, s_2, \dots, s_m\}$ of strings, find the longest substring $u$ of any string in $S$ such that there exists a subset $S' \subseteq S$ of size $t$ of substrings $u'_i$ of length $|u|$ satisfying $d_H(u, u'_i) \leq k$.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).

## Citation

If you use this software in your research or refer to our algorithm, please cite:

```
@inproceedings{rktLCS2023,
  title     = {Restricted $kt$-Longest Common Substring of multiple strings: Experimental Evaluation},
  booktitle = {Prague Stringology Conference},
  author    = {},
  year      = {?},
  publisher = {Prague Stringology Conference}
}
```

## Contact

For questions or feedback, please contact the repository maintainer.

