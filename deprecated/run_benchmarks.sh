#!/bin/bash

# Array of string lengths
lengths=(100 1000 2000 3000 10000)

# Array of k values
k_values=(1 10 20 30)

# Output file
output_file="benchmark_results.csv"

# Write header to CSV file
echo "Length,k,CPU_Time,CUDA_Transfer_Time,CUDA_Kernel_Time,CUDA_Total_Time,Speedup" > $output_file

# Run benchmarks for each combination
for l in "${lengths[@]}"; do
    for k in "${k_values[@]}"; do
        echo "Running benchmark for length=$l, k=$k"
        # Run benchmark and capture output
        result=$(./flouri_benchmark $l $k 2>&1)
        
        # Extract times using grep and sed
        cpu_time=$(echo "$result" | grep "CPU Time:" | sed 's/.*: \([0-9.]*\) .*/\1/')
        transfer_time=$(echo "$result" | grep "Memory Transfer:" | sed 's/.*: \([0-9.]*\) .*/\1/')
        kernel_time=$(echo "$result" | grep "Kernel Execution:" | sed 's/.*: \([0-9.]*\) .*/\1/')
        total_cuda_time=$(echo "$result" | grep "Total CUDA Time:" | sed 's/.*: \([0-9.]*\) .*/\1/')
        speedup=$(echo "$result" | grep "Speedup:" | sed 's/.*: \([0-9.]*\).*/\1/')
        
        # Write to CSV
        echo "$l,$k,$cpu_time,$transfer_time,$kernel_time,$total_cuda_time,$speedup" >> $output_file
        
        # Debug output
        echo "Extracted times:"
        echo "  CPU: $cpu_time"
        echo "  Transfer: $transfer_time"
        echo "  Kernel: $kernel_time"
        echo "  Total: $total_cuda_time"
        echo "  Speedup: $speedup"
    done
done

echo "Benchmark results saved to $output_file" 