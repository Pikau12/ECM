#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include </usr/include/x86_64-linux-gnu/cblas.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>

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
    size_t count_matrices = 14;
    size_t total_size = count_matrices * N * N * sizeof(double);

    matrices->base_ptr = (double*)calloc(total_size, 1);

    matrices->A = matrices->base_ptr + 0 * N * N;
    matrices->AT = matrices->base_ptr + 1 * N * N;
    matrices->B = matrices->base_ptr + 2 * N * N;
    matrices->BA = matrices->base_ptr + 3 * N * N;
    matrices->I = matrices->base_ptr + 4 * N * N;
    matrices->R = matrices->base_ptr + 5 * N * N;
    matrices->S = matrices->base_ptr + 6 * N * N;
    matrices->R_power  = matrices->base_ptr + 7 * N * N;
    matrices->A_inv    = matrices->base_ptr + 8 * N * N;
    matrices->product  = matrices->base_ptr + 9 * N * N;
    matrices->Identity = matrices->base_ptr + 10 * N * N;
    matrices->BA_neg = matrices->base_ptr + 11 * N * N;
    matrices->temp = matrices->base_ptr + 12 * N * N;
    matrices->new_S = matrices->base_ptr + 13 * N * N;
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

void multiplyMatrix(double* C, double* A, double* B,  int N) {
    const double alpha = 1.0;
    const double beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, alpha,
                A, N,
                B, N,
                beta,
                C, N);
}

void addMatrix(double* C, double* A, double* B, int N) {
    cblas_dcopy(N * N, A, 1, C, 1);
    cblas_daxpy(N * N, 1.0, B, 1, C, 1);
}

void scalarMultiply(double scalar, double* B,  double* A, int N) {
    cblas_dcopy(N * N, A, 1, B, 1);
    cblas_dscal(N * N, scalar, B, 1);
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

void dataEntry(int* N, int* M) { //+
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

void generateMatrix(MatrixSet matrices, int N, int M, bool generate) { //+
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
if (N <= 20) {
    printMatrix(matrices.product, N, "A_reverse * A");
}
    createIdentityMatrix(matrices.Identity, N);

    int ind1 = 0;
    int ind2 = 0;
    double maxDifference = 0.0;
    for (int i = 0; i < N * N; ++i) {
        
        double diff = fabs(matrices.product[i] - matrices.Identity[i]); // Исправлено!
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
