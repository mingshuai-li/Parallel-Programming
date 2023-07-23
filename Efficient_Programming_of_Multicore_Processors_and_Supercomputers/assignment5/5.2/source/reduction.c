/*
Group7: MPI Array Reduction
Build: $ mpicc -o reduction reduction.c
Run: $ sbatch job.scp
where the script job.scp comprises:
=========================================================
#!/bin/bash
#SBATCH ...
# other SLURM directives

module load slurm_setup
mpiexec -n <processesNumber> ./reduction <arrayLength> 
=========================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

bool isPowerOfTwo(int number) {
    while (number > 1 && number % 2 == 0) {
        number /= 2;
    }
    return number == 1;
}

int main(int argc, char** argv) {
    
    if (argc != 2) {
        printf("Please provide exactly one integral argument to specify the array length.\n");
        return 1; 
    }
    
    int n = atoi(argv[1]);
    int rank, size;
    int local_sum = 0;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (!isPowerOfTwo(size)) {
        printf("Number of MPI processes is not a power of 2. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        MPI_Finalize();
        return 1;
    }

    // Allocate memory for the local part of array A
    int local_size = n / size;
    int remaining = n % size;
    if (rank < remaining) {
        local_size += 1;
    }
    int* local_A = (int*)malloc(local_size * sizeof(int));
    
    // Initialize local part of array A
    if (rank < remaining){
        for (int i = 0; i < local_size; i++) {
            local_A[i] = rank * local_size + i;
        }
    }
    else{
        for (int i = 0; i < local_size; i++) {
            local_A[i] = remaining * (local_size + 1) + (rank - remaining) * local_size + i;
        }
    }

    // Local Reduction
    for (int i = 0; i < local_size; i++) {
        local_sum += local_A[i];
    }
    
    // Global Reduction
    for (int mid = size/2; mid > 0; mid /= 2) {
        if (rank < mid) {
            int received_sum;
            MPI_Recv(&received_sum, 1, MPI_INT, rank + mid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_sum += received_sum;
        }
        else {
            MPI_Send(&local_sum, 1, MPI_INT, rank - mid, 0, MPI_COMM_WORLD);
            break;
        }
    }

    // Process 0 prints the global sum
    if (rank == 0) {
        printf("Global sum: %d\n", local_sum);
    }
    
    free(local_A);
    
    MPI_Finalize();
    
    return 0;
}
