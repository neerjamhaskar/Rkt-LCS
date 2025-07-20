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

## Algorithm Details


**Rkt-LCS Algorithm**: Given finite integers $k, t, m \in \mathbb{N}$, $t \leq m$, and a set $\s{S} = \{s_1, s_2, \dots, s_m\}$ of strings, find the longest substring $u$ of any string in $S$ such that there exists a subset $S' \subseteq S$ of size $t$ of substrings $u'_i$ of length $|u|$ satisfying $d_H(u, u'_i) \leq k$.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).

## Citation

If you use this software in your research or refer to our algorithm, please cite:

```
@inproceedings{rktLCS2023,
  title     = {Approximate Longest Common Substring of Multiple Strings: Experimental Evaluation},
  author    = {Hasibi, Hamed and Mhaskar, Neerja and Smyth, WF},
  booktitle = {Proceedings of the Prague Stringology Conference},
  year      = {2025}
}
```

## Contact

For questions or feedback, please contact the repository maintainer.

