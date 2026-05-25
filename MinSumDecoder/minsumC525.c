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

// 【全新登場】集中管理所有連邊資訊的一維大陣列（告別零碎記憶體）
int* H_rows_data = NULL; // 集中存放所有 Check Node 連接的 V-node 編號
int* H_cols_data = NULL; // 集中存放所有 Variable Node 連接的 C-node 編號

long long GLOBAL_TRUEBIT, GLOBAL_FALSEBIT;
const int for_test_only_true_time_max = 0;

void Read_H_MatrixByAlist(void)
{
    char filename[100] = { "C:\\Users\\user\\Documents\\project\\Gallager_3_6.txt" };
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

// (DO NOT CHANGE!!!)ͦT
double generate_gaussian(double sigma)
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 <= 0.0)
        u1 = 1e-9; // קK log(0)
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    return z0 * sigma;
}

// (DO NOT CHANGE!!!)
void Generate_TestData(const int count, int N, int M, double EbN0_dB)
{
    char output_name[] = "C:\\Users\\user\\Documents\\project\\test_data.txt";
    FILE* f_out = fopen(output_name, "w");
    if (!f_out)
    {
        perror("Failed to open output file");
        return;
    }

    // p Sigma
    double rate = (double)(N - M) / N;
    double snr = pow(10.0, EbN0_dB / 10.0);
    double sigma = sqrt(1.0 / (2.0 * rate * snr));

    srand((unsigned int)time(NULL));

    int* codeword = (int*)calloc(N, sizeof(int));

    for (int c = 0; c < count; c++)
    {
        // --- 1. ͦ Codeword (o̼ȥΥsAY G x}iאּHX) ---
        // ܡGYnդPXrAiHͦ줸ATO Hc^T = 0
        for (int i = 0; i < N; i++)
        {
            fprintf(f_out, "%d%s", codeword[i], (i == N - 1) ? "" : " ");
        }
        fprintf(f_out, "\n");

        // --- 2. gJ Sigma ---
        fprintf(f_out, "%lf\n", sigma);

        // --- 3. ͦLqD᪺ y (BPSK: 0 -> +1.0, 1 -> -1.0) ---
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


double s_tilde(double correctionValue, double x, double y)
{
    if (fabs(x + y) < 2 && fabs(x - y) > 2 * fabs(x + y))
        return correctionValue;
    else if (fabs(x - y) < 2 && fabs(x + y) > 2 * fabs(x - y))
        return -correctionValue;
    else
        return 0;
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

    // 計算總邊數 (動態對齊你的 H_rows 結構)
    int total_row_edges = H_rows[M - 1].row_offset + H_rows[M - 1].row_weight;
    double* R_c2v = (double*)calloc(total_row_edges, sizeof(double));
    double* Q_v2c = (double*)malloc(total_row_edges * sizeof(double));

    // 初始化 LLR
    double sigma_sq = (*channel_sigma) * (*channel_sigma);

    // 💡 優化 1：初始化改用線性遞增指標 edge_idx，讓記憶體寫入連續化
    int edge_idx = 0;
    for (int i = 0; i < M; i++)
    {
        int weight = H_rows[i].row_weight;
        for (int j = 0; j < weight; j++)
        {
            int v_node = H_rows_data[edge_idx];
            Q_v2c[edge_idx] = 2.0 * LLR_Current[v_node] / sigma_sq;
            edge_idx++;
        }
    }

    *decoder_iteration = 0;
    double correctionValue = 0.1;
    bool isAllZero = true;

    // --- 解碼主迴圈 ---
    while (*decoder_iteration < decoder_max_iteration)
    {
        for (int i = 0; i < N; i++) ReturnTotal[i] = 0.0;

        /// --- STEP 1: CHECK NODE UPDATE (全面改用線性遞增指標) ---
        edge_idx = 0; // 每一輪疊代開始，重新從第 0 條邊出發
        for (int i = 0; i < M; i++)
        {
            int weight = H_rows[i].row_weight;

            int TotalNegativeCount = 0;
            minBeta first, second, third;
            first.minBetaValue = second.minBetaValue = third.minBetaValue = 1e100;

            // 第一次 Loop：直接依序讀取連續的記憶體區段
            for (int j = 0; j < weight; j++)
            {
                int k = edge_idx + j; // 計算這條邊在全域陣列中的絕對位置
                int v_node = H_rows_data[k];
                double v_message = Q_v2c[k];

                if (v_message < 0) TotalNegativeCount++;

                double tempValue = fabs(v_message);
                if (tempValue < first.minBetaValue) {
                    third = second; second = first;
                    first.minBetaValue = tempValue; first.minBetaPosition = v_node;
                }
                else if (tempValue < second.minBetaValue) {
                    third = second;
                    second.minBetaValue = tempValue; second.minBetaPosition = v_node;
                }
                else if (tempValue < third.minBetaValue) {
                    third.minBetaValue = tempValue;
                }
            }

            // 第二次 Loop：加入修正項並發送
            for (int j = 0; j < weight; j++)
            {
                int k = edge_idx + j;
                int tempPosition = H_rows_data[k];
                double v_message = Q_v2c[k];

                int self_sign = (v_message < 0) ? -1 : 1;
                int alpha = (TotalNegativeCount % 2 == 0) ? 1 : -1;
                int exclusion_alpha = alpha * self_sign;

                double pureMinValues = 0.0;
                double s_correction = 0.0;

                if (tempPosition == first.minBetaPosition) {
                    s_correction = s_tilde(correctionValue, second.minBetaValue, third.minBetaValue);
                    pureMinValues = second.minBetaValue;
                }
                else if (tempPosition == second.minBetaPosition) {
                    s_correction = s_tilde(correctionValue, first.minBetaValue, third.minBetaValue);
                    pureMinValues = first.minBetaValue;
                }
                else {
                    s_correction = s_tilde(correctionValue, first.minBetaValue, second.minBetaValue);
                    pureMinValues = first.minBetaValue;
                }

                double new_R = exclusion_alpha * (pureMinValues + s_correction);
                R_c2v[k] = new_R;
                ReturnTotal[tempPosition] += new_R;
            }

            edge_idx += weight; // 💡 關鍵：處理完一個 Check Node，直接跳過該列的權重長度
        }

        /// --- STEP 2: VARIABLE NODE UPDATE & HARD DECISION ---
        for (int i = 0; i < N; i++)
        {
            double LLR_intrinsic = 2.0 * LLR_Current[i] / sigma_sq;
            LLR_FinalResult[i] = LLR_intrinsic + ReturnTotal[i];
            Ci[i] = (LLR_FinalResult[i] >= 0.0) ? 0 : 1;
        }

        // 💡 優化 2：更新 Q_v2c 邊訊息時，同步改用線性累加指標 edge_idx
        edge_idx = 0;
        for (int i = 0; i < M; i++)
        {
            int weight = H_rows[i].row_weight;
            for (int j = 0; j < weight; j++)
            {
                int v_node = H_rows_data[edge_idx];
                Q_v2c[edge_idx] = LLR_FinalResult[v_node] - R_c2v[edge_idx];
                edge_idx++;
            }
        }

        /// --- STEP 3: PARITY CHECK (同步改用線性累加指標優化) ---
        //isAllZero = true;
        edge_idx = 0;
        for (int i = 0; i < M; i++)
        {
            int decision = 0;
            int weight = H_rows[i].row_weight;
            for (int j = 0; j < weight; j++) {
                decision += Ci[H_rows_data[edge_idx++]];
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
    printf("LDPC Decoder Test Start\n");
    clock_t start_read, end_read, start_run, end_run;
    // Read H Matrix
    start_read = clock();
    Read_H_MatrixByAlist();
    end_read = clock();
    const int num_test_cases = 100;
    const int deocder_iteration_max = 16;
    const double EbN0_dB = 1.6;

    int* TEST_countEachIteration = (int*)malloc(num_test_cases * sizeof(int));

    // initial
    for (int i = 0; i < num_test_cases; i++)
    {
        TEST_countEachIteration[i] = -1;
    }
    int totalSucessTime = 0, totalFailTime = 0;
    GLOBAL_TRUEBIT = 0, GLOBAL_FALSEBIT = 0;

    // Generate TestData
    Generate_TestData(num_test_cases, N, M, EbN0_dB);

    // open test data file
    char filename[100] = { "C:\\Users\\user\\Documents\\project\\test_data.txt" };
    FILE* TEST_DATA = fopen(filename, "r");
    if (!TEST_DATA)
    {
        printf("[SYSTEM]MAIN_Open file failed\n");
        return 0; // 💡 順手幫你將 main 的錯誤回傳改成標準 int 值 0
    }

    printf("Read time: %lf\n", (double)(end_read - start_read) / CLOCKS_PER_SEC);

    // allocate memory for y(after AWGN channel)
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
            // printf("%lf ",channel_output[j]);
        }
        fscanf(TEST_DATA, "%lf", channel_sigma);
        // printf("MAIN Channel_sigma = %lf \n",*channel_sigma);

        for (int j = 0; j < N; j++)
        {
            fscanf(TEST_DATA, "%lf", &channel_output[j]);
            // printf("%lf ",channel_output[j]);
        }
        bool isThisIterationSucess = isDecodeSucess(channel_sigma, decode_iteration, deocder_iteration_max, channel_output);

        if (isThisIterationSucess == true)
        {
            // printf("MAIN_Decode Sucess In %d Iterations\n",*decode_iteration);
            totalSucessTime += 1;
        }
        else
        {
            // printf("MAIN_Decode Failed In %d Iterations\n",*decode_iteration);
            totalFailTime += 1;
        }
        //ST_countEachIteration[testtime] = *decode_iteration;

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
    // BER 實際上就是你原本算出來的 GLOBAL_ERR_RATE
    printf("Bit Error Rate (BER)  : %e\n", GLOBAL_ERR_RATE);
    // FER = 解碼失敗的訊框數 (totalFailTime) / 總測試訊框數 (num_test_cases)
    printf("Frame Error Rate (FER): %e\n", (double)totalFailTime / (double)num_test_cases);


    // printf("MAIN_Run More %d TIMEs\n", for_test_only_true_time_max);

    // 釋放動態記憶體避免 Memory Leak
    free(TEST_countEachIteration);
    free(channel_output);
    free(channel_sigma);
    free(decode_iteration);

    return 0;
}
