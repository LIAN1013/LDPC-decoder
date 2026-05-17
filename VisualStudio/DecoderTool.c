#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

#include "DecoderTool.h"


int N = -1, K = -1, M = -1;
RowInfo* H_rows = NULL;
ColInfo* H_cols = NULL;

int* H_row_places_flat = NULL;
int* H_col_places_flat = NULL;

long long GLOBAL_TRUEBIT = 0, GLOBAL_FALSEBIT = 0;



// Read H Matrix
void Read_H_MatrixByAlist(void) {
    FILE* fp = fopen(H_MATRIX_FILE, "r");

    if (!fp) {
        printf("[SYSTEM] Open file failed: %s\n", H_MATRIX_FILE);
        return;
    }

    int max_col_weight = 0, max_row_weight = 0;
    fscanf(fp, "%d %d", &N, &M);
    fscanf(fp, "%d %d", &max_col_weight, &max_row_weight);

    H_rows = (RowInfo*)malloc(M * sizeof(RowInfo));
    H_cols = (ColInfo*)malloc(N * sizeof(ColInfo));

    // 分配連續的一維空間 (大小為 總行數*最大權重)
    H_col_places_flat = (int*)malloc(N * max_col_weight * sizeof(int));
    H_row_places_flat = (int*)malloc(M * max_row_weight * sizeof(int));

    // 讀取權重並計算偏移量
    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d", &H_cols[i].column_weight);
        H_cols[i].offset = i * max_col_weight;
    }
    for (int i = 0; i < M; i++) {
        fscanf(fp, "%d", &H_rows[i].row_weight);
        H_rows[i].offset = i * max_row_weight;
    }

    int temp = 0;
    // 讀取 Column 位置到一維陣列
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < H_cols[i].column_weight; j++) {
            fscanf(fp, "%d", &temp);
            H_col_places_flat[H_cols[i].offset + j] = temp - 1;
        }
    }
    // 讀取 Row 位置到一維陣列
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < H_rows[i].row_weight; j++) {
            fscanf(fp, "%d", &temp);
            H_row_places_flat[H_rows[i].offset + j] = temp - 1;
        }
    }
    fclose(fp);
}


// Generate Gaussian Test Data
double generate_gaussian(double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 <= 0.0) u1 = 1e-9;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    return z0 * sigma;
}

void Generate_TestData(const int count, int N, int M, double EbN0_dB) {
    FILE* f_out = fopen(TEST_DATA_FILE, "w");
    if (!f_out) {
        perror("[System] Failed to open Gaussian TestData file\n");
        return;
    }

    double rate = (double)(N - M) / N;
    double snr = pow(10.0, EbN0_dB / 10.0);
    double sigma = sqrt(1.0 / (2.0 * rate * snr));

    srand((unsigned int)time(NULL));
    int* codeword = (int*)calloc(N, sizeof(int));

    for (int c = 0; c < count; c++) {
        for (int i = 0; i < N; i++) fprintf(f_out, "%d%s", codeword[i], (i == N - 1) ? "" : " ");
        fprintf(f_out, "\n%lf\n", sigma);
        for (int i = 0; i < N; i++) {
            double s = (codeword[i] == 0) ? 1.0 : -1.0;
            double y = s + generate_gaussian(sigma);
            fprintf(f_out, "%lf%s", y, (i == N - 1) ? "" : " ");
        }
        fprintf(f_out, "\n");
    }
    free(codeword);
    fclose(f_out);
    printf("Successfully generated %d sets of data with Sigma: %lf\n", count, sigma);
}

// --- SPA 解碼核心副程式 ---
double phi(double x) {
    // Version 1 : only limit maxinum minimum
    double value = 0;
    if (value > 30.0)
        return 1e-13;
    else if (value < 1e-13)
        return 30.0;
    else
        return -log(tanh(x / 2.0));
    // version 2 : Split to X(EX:8) point

}



