#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <immintrin.h>

void dgemm(float alpha, const float *a, const float *b, float beta, float *c) {
    int endvec = MATRIX_SIZE - MATRIX_SIZE % 8;
    for (int i = 0; i < MATRIX_SIZE; i++) { 

        for (int j = 0; j < MATRIX_SIZE; j++) {
            
            c[i * MATRIX_SIZE + j] *= beta;

            float partial_sum_array[8] = {0};

            __m256 partial_sum = _mm256_setzero_ps();
            float sum = 0;

            for (int k = 0; k < endvec; k += 8) {
                __m256 avec = _mm256_loadu_ps(a + i * MATRIX_SIZE + k);
                __m256 bvec = _mm256_loadu_ps(b + j * MATRIX_SIZE + k);
                __m256 prod = _mm256_mul_ps(avec, bvec);
                partial_sum = _mm256_add_ps(partial_sum, prod);
            }

            _mm256_storeu_ps(partial_sum_array, partial_sum);
            

            sum = partial_sum_array[0] + partial_sum_array[1] + partial_sum_array[2] + partial_sum_array[3] + partial_sum_array[4] + partial_sum[5] + partial_sum[6] + partial_sum[7];


            for (int k = endvec; k < MATRIX_SIZE; k++) {
                sum += a[i * MATRIX_SIZE + k] * b[j * MATRIX_SIZE + k];
            }

            c[i * MATRIX_SIZE + j] += alpha * sum;
        }
    }
}

int main(int, char **) {
    float alpha, beta;

    // mem allocations
    int mem_size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
    auto a = (float *) malloc(mem_size);
    auto b = (float *) malloc(mem_size);
    auto c = (float *) malloc(mem_size);

    // check if allocated
    /*if (nullptr == a || nullptr == b || nullptr == c) {
        printf("Memory allocation failed\n");
        if (nullptr != a) free(a);
        if (nullptr != b) free(b);
        if (nullptr != c) free(c);
        return 0;
    }*/

    generateProblemFromInput(alpha, a, b, beta, c);

    std::cerr << "Launching dgemm step." << std::endl;
    // matrix-multiplication
    dgemm(alpha, a, b, beta, c);

    outputSolution(c);

    free(a);
    free(b);
    free(c);
    return 0;
}
