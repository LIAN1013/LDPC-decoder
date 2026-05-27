#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

typedef struct
{
    int row_weight;
    int row_offset; // 紀錄這一列在 H_rows_data 中的起點位置
} RowInfo;

typedef struct
{
    int column_weight;
    int column_offset; // 紀錄這一欄在 H_cols_data 中的起點位置
} ColInfo;

typedef struct
{
    double minBetaValue;
    int minBetaPosition;
} minBeta;

int N = -1, K = -1, M = -1;
RowInfo* H_rows = NULL;
ColInfo* H_cols = NULL;

// 集中管理所有連邊資訊的一維大陣列
int* H_rows_data = NULL;
int* H_cols_data = NULL;

long long GLOBAL_TRUEBIT, GLOBAL_FALSEBIT;
const int for_test_only_true_time_max = 0;

void Read_H_MatrixByAlist(void)
{
    //char filename[100] = { "C:\\Users\\user\\Documents\\project\\Gallager_3_6.txt" };
    char filename[100] = { "C:\\Users\\USER\\Desktop\\Min sum decoder\\Gallager_3_6.txt" };
    FILE* fp = fopen(filename, "r");

    if (!fp)
    {
        printf("Open file failed\n");
        return;
    }

    int max_row_weight = 0;
    int max_col_weight = 0;

    fscanf(fp, "%d %d", &N, &M);
    fscanf(fp, "%d %d", &max_row_weight, &max_col_weight);

    H_rows = (RowInfo*)malloc(M * sizeof(RowInfo));
    H_cols = (ColInfo*)malloc(N * sizeof(ColInfo));

    // column
    for (int i = 0; i < N; i++)
    {
        fscanf(fp, "%d", &H_cols[i].column_weight);
    }

    H_cols[0].column_offset = 0;
    for (int i = 0; i < N - 1; i++)
    {
        H_cols[i + 1].column_offset = H_cols[i].column_offset + H_cols[i].column_weight;
    }
    int total_col_edges = H_cols[N - 1].column_offset + H_cols[N - 1].column_weight;

    // row
    for (int i = 0; i < M; i++)
    {
        fscanf(fp, "%d", &H_rows[i].row_weight);
    }
    H_rows[0].row_offset = 0;
    for (int i = 0; i < M - 1; i++)
    {
        H_rows[i + 1].row_offset = H_rows[i].row_offset + H_rows[i].row_weight;
    }
    int total_row_edges = H_rows[M - 1].row_offset + H_rows[M - 1].row_weight;

    H_cols_data = (int*)malloc(total_col_edges * sizeof(int));
    H_rows_data = (int*)malloc(total_row_edges * sizeof(int));
    
    // dim2
    int temp = 0;
    for (int i = 0; i < N; i++)
    {
        int start_index = H_cols[i].column_offset;
        for (int j = 0; j < H_cols[i].column_weight; j++)
        {
            fscanf(fp, "%d", &temp);
            H_cols_data[start_index + j] = temp - 1;
        }
    }
    for (int i = 0; i < M; i++)
    {
        int start_idx = H_rows[i].row_offset;
        for (int j = 0; j < H_rows[i].row_weight; j++)
        {
            fscanf(fp, "%d", &temp);
            H_rows_data[start_idx + j] = temp - 1;
        }
    }

    fclose(fp);
    printf("Read Finish\n");
    return;
}

double generate_gaussian(double sigma)
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 <= 0.0)
        u1 = 1e-9; 
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    return z0 * sigma;
}

void Generate_TestData(const int count, int N, int M, double EbN0_dB)
{
    //char output_name[] = "C:\\Users\\user\\Documents\\project\\test_data.txt";
    char output_name[] = "C:\\Users\\USER\\Desktop\\Min sum decoder\\test_data.txt";
    FILE* f_out = fopen(output_name, "w");
    if (!f_out)
    {
        perror("Failed to open output file");
        return;
    }

    double rate = (double)(N - M) / N;
    double snr = pow(10.0, EbN0_dB / 10.0);
    double sigma = sqrt(1.0 / (2.0 * rate * snr));

    srand((unsigned int)time(NULL));

    int* codeword = (int*)calloc(N, sizeof(int));

    for (int c = 0; c < count; c++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(f_out, "%d%s", codeword[i], (i == N - 1) ? "" : " ");
        }
        fprintf(f_out, "\n");

        fprintf(f_out, "%lf\n", sigma);

        for (int i = 0; i < N; i++)
        {
            double s = (codeword[i] == 0) ? 1.0 : -1.0;
            double noise = generate_gaussian(sigma);
            double y = s + noise;
            fprintf(f_out, "%lf%s", y, (i == N - 1) ? "" : " ");
        }
        fprintf(f_out, "\n");
    }
    free(codeword);
    fclose(f_out);
    printf("Successfully generated %d sets of data with Sigma: %lf\n", count, sigma);
}

bool isDecodeSucess(double* channel_sigma, int* decoder_iteration, int decoder_max_iteration, double* LLR_Current)
{
    if (N == -1 || M == -1) {
        printf("[SYSTEM] READ Alist Failed\n");
        return false;
    }

    double* LLR_FinalResult = (double*)malloc(N * sizeof(double));
    int* Ci = (int*)malloc(N * sizeof(int));
    double* ReturnTotal = (double*)malloc(N * sizeof(double));

    int total_row_edges = H_rows[M - 1].row_offset + H_rows[M - 1].row_weight;
    double* R_c2v = (double*)calloc(total_row_edges, sizeof(double));

    double* Q_v2c = (double*)malloc(total_row_edges * sizeof(double));

    double sigma_sq = (*channel_sigma) * (*channel_sigma);
    
    for (int i = 0; i < M; i++)
    {
        int start = H_rows[i].row_offset;
        int end = start + H_rows[i].row_weight;
        for (int k = start; k < end; k++)
        {
            int v_node = H_rows_data[k];
            Q_v2c[k] = 2.0 * LLR_Current[v_node] / sigma_sq;
        }
    }

    *decoder_iteration = 0;
    bool isAllZero = false;

    // --- 解碼主迴圈 ---
    while (*decoder_iteration < decoder_max_iteration)
    {
        for (int i = 0; i < N; i++) ReturnTotal[i] = 0.0;

        /// --- STEP 1: CHECK NODE UPDATE (Pure Min-Sum) ---
        for (int i = 0; i < M; i++)
        {
            int start = H_rows[i].row_offset;
            int end = start + H_rows[i].row_weight;

            int TotalNegativeCount = 0;
            minBeta first, second;
            // 初始化最小與次小值
            first.minBetaValue = second.minBetaValue = 1e100;

            // 尋找最小與次小絕對值
            for (int k = start; k < end; k++)
            {
                int v_node = H_rows_data[k];
                double v_message = Q_v2c[k]; 

                if (v_message < 0) TotalNegativeCount++;

                double tempValue = fabs(v_message);
                
                if (tempValue < first.minBetaValue) {
                    second = first;
                    first.minBetaValue = tempValue; 
                    first.minBetaPosition = v_node;
                }
                else if (tempValue < second.minBetaValue) {
                    second.minBetaValue = tempValue; 
                    second.minBetaPosition = v_node;
                }
            }

            // 更新 Check to Variable Node 訊息
            for (int k = start; k < end; k++)
            {
                int tempPosition = H_rows_data[k];
                double v_message = Q_v2c[k]; 

                int self_sign = (v_message < 0) ? -1 : 1;
                int alpha = (TotalNegativeCount % 2 == 0) ? 1 : -1;
                int exclusion_alpha = alpha * self_sign;

                // 純 Min-Sum 邏輯：如果你是最小值，那就用次小值；否則就用最小值
                double pureMinValues = (tempPosition == first.minBetaPosition) ? second.minBetaValue : first.minBetaValue;

                double new_R = exclusion_alpha * pureMinValues;
                R_c2v[k] = new_R;
                ReturnTotal[tempPosition] += new_R; 
            }
        }

        /// --- STEP 2: VARIABLE NODE UPDATE & HARD DECISION ---
        for (int i = 0; i < N; i++)
        {
            double LLR_intrinsic = 2.0 * LLR_Current[i] / sigma_sq;
            LLR_FinalResult[i] = LLR_intrinsic + ReturnTotal[i]; 

            Ci[i] = (LLR_FinalResult[i] >= 0.0) ? 0 : 1;
        }

        for (int i = 0; i < M; i++)
        {
            int start = H_rows[i].row_offset;
            int end = start + H_rows[i].row_weight;
            for (int k = start; k < end; k++)
            {
                int v_node = H_rows_data[k];
                Q_v2c[k] = LLR_FinalResult[v_node] - R_c2v[k]; 
            }
        }

        /// --- STEP 3: PARITY CHECK (H * c^T == 0) ---
        isAllZero = true;
        for (int i = 0; i < M; i++)
        {
            int decision = 0;
            int start = H_rows[i].row_offset;
            int end = start + H_rows[i].row_weight;
            for (int k = start; k < end; k++) {
                decision += Ci[H_rows_data[k]];
            }
            if (decision % 2 != 0) {
                isAllZero = false;
                break;
            }
        }

        if (isAllZero) break; 

        *decoder_iteration += 1;
    }

    /// --- 終點線記分板 ---
    long long this_block_true = 0;
    long long this_block_false = 0;

    for (int i = 0; i < N; i++)
    {
        if (Ci[i] == 0) {
            this_block_true++;  
        }
        else {
            this_block_false++; 
        }
    }

    GLOBAL_TRUEBIT += this_block_true;
    GLOBAL_FALSEBIT += this_block_false;

    free(LLR_FinalResult); free(Ci); free(ReturnTotal); free(R_c2v); free(Q_v2c);

    return isAllZero; 
}

int main(void)
{
    clock_t start_read, end_read, start_run, end_run;
    start_read = clock();
    Read_H_MatrixByAlist();
    const int num_test_cases = 3000;
    const int deocder_iteration_max = 16;
    const double EbN0_dB = 1.6;

    int* TEST_countEachIteration = (int*)malloc(num_test_cases * sizeof(int));

    for (int i = 0; i < num_test_cases; i++)
    {
        TEST_countEachIteration[i] = -1;
    }
    int totalSucessTime = 0, totalFailTime = 0;
    GLOBAL_TRUEBIT = 0, GLOBAL_FALSEBIT = 0;

    Generate_TestData(num_test_cases, N, M, EbN0_dB);

    //char filename[100] = { "C:\\Users\\user\\Documents\\project\\test_data.txt" };
    char filename[100] = { "C:\\Users\\USER\\Desktop\\Min sum decoder\\test_data.txt" };
    FILE* TEST_DATA = fopen(filename, "r");
    if (!TEST_DATA)
    {
        printf("[SYSTEM]MAIN_Open file failed\n");
        return 0; 
    }
    end_read = clock();
    printf("Read time: %lf\n", (double)(end_read - start_read) / CLOCKS_PER_SEC);

    double* channel_output = (double*)malloc(N * sizeof(double));
    double* channel_sigma = (double*)malloc(sizeof(double));
    int* decode_iteration = (int*)malloc(sizeof(int));

    start_run = clock();
    for (int testtime = 0; testtime < num_test_cases; testtime++)
    {
        int dummy;
        for (int j = 0; j < N; j++)
        {
            fscanf(TEST_DATA, "%d", &dummy);
        }
        fscanf(TEST_DATA, "%lf", channel_sigma);

        for (int j = 0; j < N; j++)
        {
            fscanf(TEST_DATA, "%lf", &channel_output[j]);
        }
        bool isThisIterationSucess = isDecodeSucess(channel_sigma, decode_iteration, deocder_iteration_max, channel_output);

        if (isThisIterationSucess == true)
        {
            totalSucessTime += 1;
        }
        else
        {
            totalFailTime += 1;
        }
    }

    end_run = clock();
    printf("Run time: %lf\n", (double)(end_run - start_run) / CLOCKS_PER_SEC);
    printf("Total Test Time = %d\n", num_test_cases);
    printf("deocder_iteration_max = %d\n", deocder_iteration_max);
    printf("MAIN_SucessTime = %d   Failed Time = %d\n", totalSucessTime, totalFailTime);
    printf("MAIN_CorrectRate = %lf \n", (float)totalSucessTime / (float)num_test_cases);

    printf("\n");
    printf("MAIN_In Eb/N0 %f dB\n", EbN0_dB);
    double GLOBAL_ERR_RATE = (double)GLOBAL_FALSEBIT / (double)(GLOBAL_TRUEBIT + GLOBAL_FALSEBIT);
    printf("MIAN_TRUEBIT = %lld  FALSEBIT = %lld\nMAIN_ERROR_RATE = %g\n\n", GLOBAL_TRUEBIT, GLOBAL_FALSEBIT, GLOBAL_ERR_RATE);

    printf("---------------------------------------------\n");
    printf("Bit Error Rate (BER)  : %e\n", GLOBAL_ERR_RATE);
    printf("Frame Error Rate (FER): %e\n", (double)totalFailTime / (double)num_test_cases);

    free(TEST_countEachIteration);
    free(channel_output);
    free(channel_sigma);
    free(decode_iteration);

    return 0;
}
