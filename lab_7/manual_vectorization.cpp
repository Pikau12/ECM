#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>
#include <stdint.h>


typedef struct {
    double* base_ptr;
    double* A;
    double* AT;
    double* B;
    double* BA;
    double* BA_neg;
    double* I;
    double* R;
    double* S;
    double* new_S;
    double* R_power;
    double* A_inv;
    double* product;
    double* Identity;
    double* temp;
} MatrixSet;

void allocateMatrices(int N, MatrixSet* matrices) {

    matrices->A = (double*)calloc(N * N, sizeof(double));
    matrices->AT = (double*)calloc(N * N, sizeof(double));
    matrices->B = (double*)calloc(N * N, sizeof(double));
    matrices->BA = (double*)calloc(N * N, sizeof(double));
    matrices->BA_neg = (double*)calloc(N * N, sizeof(double));
    matrices->I = (double*)calloc(N * N, sizeof(double));
    matrices->R = (double*)calloc(N * N, sizeof(double));
    matrices->S = (double*)calloc(N * N, sizeof(double));
    matrices->new_S = (double*)calloc(N * N, sizeof(double));
    matrices->R_power = (double*)calloc(N * N, sizeof(double));
    matrices->A_inv = (double*)calloc(N * N, sizeof(double));
    matrices->product = (double*)calloc(N * N, sizeof(double));
    matrices->Identity = (double*)calloc(N * N, sizeof(double));
    matrices->temp = (double*)calloc(N * N, sizeof(double));
}

void freeMatrices(MatrixSet* matrices) {
    free(matrices->A);

    matrices->A = NULL;
    matrices->AT = NULL;
    matrices->B = NULL;
    matrices->BA = NULL;
    matrices->BA_neg = NULL;
    matrices->I = NULL;
    matrices->R = NULL;
    matrices->S = NULL;
    matrices->R_power = NULL;
    matrices->A_inv = NULL;
    matrices->product = NULL;
    matrices->Identity = NULL;
    matrices->temp = NULL;
    matrices->new_S = NULL;
}

void createIdentityMatrix(double* result, int N) { 
    for (int i = 0; i < N; ++i) {
        result[i * N + i] = 1.0;
    }
}

void transpose(double* A, double* B, int N) { 
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[j * N + i] = B[i * N + j];
        }
    }
}

void addMatrix(double* C, double* A, double* B, int N) {
    int size = N * N;
    int i = 0;
    for (; i <= size - 4; i += 4) {
        __m256d a_vec = _mm256_loadu_pd(A + i);
        __m256d b_vec = _mm256_loadu_pd(B + i);
        __m256d c_vec = _mm256_add_pd(a_vec, b_vec);
        _mm256_storeu_pd(C + i, c_vec);
    }
    for (; i < size; ++i) {
        C[i] = A[i] + B[i];
    }
}

void scalarMultiply(double scalar, double* B, double* A, int N) {
    int size = N * N;
    int i = 0;
    __m256d scalar_vec = _mm256_set1_pd(scalar);
    for (; i <= size - 4; i += 4) {
        __m256d a_vec = _mm256_loadu_pd(A + i);
        __m256d b_vec = _mm256_mul_pd(scalar_vec, a_vec);
        _mm256_storeu_pd(B + i, b_vec);
    }
    for (; i < size; ++i) {
        B[i] = scalar * A[i];
    }
}

void multiplyMatrix(double* result, double* A, double* B, int n) {
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            double temp = A[i * n + k];
            for (int j = 0; j < n; j += 4) {
                __m256d result_vec = _mm256_loadu_pd(&result[i * n + j]); 
                __m256d B_vec = _mm256_loadu_pd(&B[k * n + j]);         
                result_vec = _mm256_add_pd(result_vec, _mm256_mul_pd(_mm256_set1_pd(temp), B_vec)); 
                _mm256_storeu_pd(&result[i * n + j], result_vec); 
            }
        }
    }
}

void copyMatrix(double* A,  double* B, int N) {
    for (int i = 0; i < N * N; ++i) {
        A[i] = B[i];
    }
}

double norm1(double* A, int N) { // норма по стобцу 
    double maxSum = 0.0;
    for (int j = 0; j < N; ++j) {
        double columnSum = 0.0;
        for (int i = 0; i < N; ++i) {
            columnSum += fabs(A[i * N + j]);
        }
        if (columnSum > maxSum) {
            maxSum = columnSum;
        }
    }
    return maxSum;
}

double normInf(double* A, int N) { // норма по строке 
    double maxSum = 0.0;
    for (int i = 0; i < N; ++i) {
        double rowSum = 0.0;
        for (int j = 0; j < N; ++j) {
            rowSum += fabs(A[i * N + j]);
        }
        if (rowSum > maxSum) {
            maxSum = rowSum;
        }
    }
    return maxSum;
}

void printMatrix(double* A, int N, const char* description) {
    printf("%s:\n", description);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%10.6lf ", A[i * N + j]);
        }
        printf("\n");
    }
}

double* generateRandomMatrix(double* A, int N) {
    for (int i = 0; i < N * N; ++i) {
        A[i] = ((double)rand() / RAND_MAX) * 100.0;
    }
    return A;
}

void dataEntry(int* N, int* M) { 
    printf("Введите N: ");
    if (scanf("%d", N) != 1 || *N <= 0) {
        fprintf(stderr, "Некорректный ввод N\n");
        exit(EXIT_FAILURE); 
    }

    printf("Введите M: ");
    if (scanf("%d", M) != 1 || *M <= 0) {
        fprintf(stderr, "Некорректный ввод M\n");
        exit(EXIT_FAILURE); 
    }
}

void generateMatrix(MatrixSet matrices, int N, int M, bool generate) { 
    if (generate == 1) {
        generateRandomMatrix(matrices.A, N);
        if (N <= 20){
            printMatrix(matrices.A, N, "Матрица A");
        }
    } else {
        printf("Введите матрицу A (%d элементов):\n", N * N);
        for (int i = 0; i < N * N; ++i) {
            if (scanf("%lf", &(matrices.A)[i]) != 1) { 
                freeMatrices(&matrices);
                exit(EXIT_FAILURE);
            }
        }
    }
}

int main(int argc, char* argv[]) {
    int N, M;
    bool generate = false;
    MatrixSet matrices;

    if (argc > 1 && strcmp(argv[1], "1") == 0) {
        generate = true;
    }

    dataEntry(&N, &M);
    allocateMatrices(N, &matrices);
    generateMatrix(matrices, N, M, generate);

    struct timespec start, end; 
    if (clock_gettime(CLOCK_MONOTONIC, &start) == -1) { // получаем время началы работы
        freeMatrices(&matrices);
        return EXIT_FAILURE;
    }

    double normA1 = norm1(matrices.A, N);
    double normAInf = normInf(matrices.A, N);
    double s = normA1 * normAInf;

    //B = (1/s) * A^T
    transpose(matrices.AT, matrices.A, N);
    scalarMultiply(1.0 / s, matrices.B, matrices.AT, N);

    //R = I - BA
    multiplyMatrix(matrices.BA, matrices.B, matrices.A, N);
    scalarMultiply(-1.0, matrices.BA_neg, matrices.BA, N);
    createIdentityMatrix(matrices.I, N);
    addMatrix(matrices.R, matrices.I, matrices.BA_neg, N);

    //S = I + R + R^2 + ... + R^(M-1)
    createIdentityMatrix(matrices.S, N);
    createIdentityMatrix(matrices.R_power, N);

    for (int k = 1; k < M; ++k) {
        //R_power = R_power * R
        multiplyMatrix(matrices.temp, matrices.R_power, matrices.R, N);
        //S = S + R_power
        addMatrix(matrices.new_S, matrices.S, matrices.temp, N);
        copyMatrix(matrices.S, matrices.new_S, N); 
    }

    multiplyMatrix(matrices.A_inv, matrices.S, matrices.B, N);

    if (clock_gettime(CLOCK_MONOTONIC, &end) == -1) {
        freeMatrices(&matrices);
        return EXIT_FAILURE;
    }

    double elapsed_time = (end.tv_sec - start.tv_sec) * 1000.0;
    elapsed_time += (end.tv_nsec - start.tv_nsec) / 1.0e6;

    if (N <= 20){
        printMatrix(matrices.A_inv, N, "A_reverse");
    }

    printf("Время вычислений: %.3lf миллисекунд\n", elapsed_time);

// типо тесты
    multiplyMatrix(matrices.product, matrices.A_inv, matrices.A, N);
    if (N <= 20){
        printMatrix(matrices.product, N, "A_reverse * A");
    }
    createIdentityMatrix(matrices.Identity, N);

    int ind1 = 0;
    int ind2 = 0;
    double maxDifference = 0.0;
    for (int i = 0; i < N * N; ++i) {
        double diff = fabs(matrices.product[i] - matrices.Identity[i]);
        if (diff > maxDifference) {
            maxDifference = diff;
            ind1 = i / N;
            ind2 = i % N;
        }
    }

    printf("Max модуль разности (A_reverse * A) и I: %.6lf\n", maxDifference);
    printf("Ячейка: [%d, %d]\n", ind1, ind2);

    freeMatrices(&matrices);
    return 0;
}
