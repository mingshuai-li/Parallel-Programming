#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NUM_MESSAGES 1000

int main(int argc, char** argv) {
    int rank, size;
    MPI_Status status;
    double start_time, end_time, total_time;
    int message_sizes[25];
    int i, j, size_idx;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0) {
            printf("This program requires exactly 2 processes.\n");
        }
        MPI_Finalize();
        return 0;
    }

    for (i = 0; i < 25; i++) {
        message_sizes[i] = 1;
        for (j = 0; j < i; j++){
            message_sizes[i] *= 2;
        }
    }

    // Perform ping-pong communication
    for (size_idx = 0; size_idx < 25; size_idx++) {
        int message_size = message_sizes[size_idx];
        char* message = (char*)malloc(message_size);

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            // Process A sends message to Process B
            start_time = MPI_Wtime();
            for (i = 0; i < NUM_MESSAGES; i++) {
                MPI_Send(message, message_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(message, message_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
            }
            end_time = MPI_Wtime();
        } else if (rank == 1) {
            // Process B receives message from Process A and sends it back
            for (i = 0; i < NUM_MESSAGES; i++) {
                MPI_Recv(message, message_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
                MPI_Send(message, message_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }

        total_time = end_time - start_time;

        if (rank == 0) {
            double startup_time = total_time / (2 * NUM_MESSAGES);
            double bandwidth = (message_size * 1e-6) / startup_time;
            
            printf("Message Size: %d bytes\n", message_size);
            printf("Startup Time: %.6f seconds\n", startup_time);
            printf("Bandwidth: %.6f MB/s\n", bandwidth);
            printf("===========================================\n");
        }

        free(message);
    }

    MPI_Finalize();
    return 0;
}
