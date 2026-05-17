#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

#include "DecoderTool.h"

/*
    Update Record
    2026/05/17 18:07 Add Path Prefix
    2026/05/17 19:36

*/










// --- 主程式 ---
int main(void) 
{
    Read_H_MatrixByAlist();

    const int num_test_cases = 100;
    const int decoder_iteration_max = 50;
    const double EbN0_dB = 5.0;

    int* TEST_countEachIteration = (int*)malloc(num_test_cases * sizeof(int));
    for (int i = 0; i < num_test_cases; i++) TEST_countEachIteration[i] = -1;

    int totalSucessTime = 0, totalFailTime = 0;
    GLOBAL_TRUEBIT = 0; GLOBAL_FALSEBIT = 0;

    Generate_TestData(num_test_cases, N, M, EbN0_dB);

    FILE* TEST_DATA = fopen(TEST_DATA_FILE, "r");
    if (!TEST_DATA) return 1;

    double* channel_output = (double*)malloc(N * sizeof(double));
    double* channel_sigma = (double*)malloc(sizeof(double));
    int* decode_iteration = (int*)malloc(sizeof(int));

    for (int testtime = 0; testtime < num_test_cases; testtime++) {
        double dummy;
        for (int j = 0; j < N; j++) fscanf(TEST_DATA, "%lf", &dummy);
        fscanf(TEST_DATA, "%lf", channel_sigma);
        for (int j = 0; j < N; j++) fscanf(TEST_DATA, "%lf", &channel_output[j]);

        bool isThisIterationSucess = Decode_SumProduct(channel_sigma, decode_iteration, decoder_iteration_max, channel_output);

        if (isThisIterationSucess) {
            totalSucessTime += 1;
        }
        else {
            totalFailTime += 1;
        }
        TEST_countEachIteration[testtime] = *decode_iteration;
    }

    printf("Total Test Time = %d\n", num_test_cases);
    printf("MAIN_SucessTime = %d   Failed Time = %d\n", totalSucessTime, totalFailTime);
    printf("MAIN_CorrectRate = %lf \n", (float)totalSucessTime / (float)num_test_cases);
    printf("MAIN_In Eb/N0 %f dB\n", EbN0_dB);

    double GLOBAL_ERR_RATE = (double)GLOBAL_FALSEBIT / (double)(GLOBAL_TRUEBIT + GLOBAL_FALSEBIT);
    printf("MAIN_TRUEBIT = %lld  FALSEBIT = %lld\nMAIN_ERROR_RATE = %g\n\n", GLOBAL_TRUEBIT, GLOBAL_FALSEBIT, GLOBAL_ERR_RATE);

    // 釋放記憶體
    free(H_row_places_flat); free(H_col_places_flat);
    free(H_rows); free(H_cols);
    free(channel_output); free(channel_sigma); free(decode_iteration); free(TEST_countEachIteration);
    fclose(TEST_DATA);

    return 0;
}
