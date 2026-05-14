#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

// --- 定義與全域變數 ---
typedef struct {
    int row_weight;
    int offset; // 改為紀錄在一維陣列中的起始索引
} RowInfo;

typedef struct {
    int column_weight;
    int offset; // 改為紀錄在一維陣列中的起始索引
} ColInfo;

int N = -1, K = -1, M = -1;
RowInfo *H_rows = NULL;
ColInfo *H_cols = NULL;

// --- 【新增】一維儲存大陣列 ---
int *H_row_places_flat = NULL;
int *H_col_places_flat = NULL;

long long GLOBAL_TRUEBIT = 0, GLOBAL_FALSEBIT = 0;
const int for_test_only_true_time_max = 500;

// --- 檔案讀取 (Alist) 改為一維存儲 ---
void Read_H_MatrixByAlist(void) {
    char filename[100] = {"C:\\Users\\user\\Documents\\project\\Gallager_3_6.txt"};
    FILE *fp = fopen(filename, "r");

    if (!fp) {
        printf("[SYSTEM] Open file failed: %s\n", filename);
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

// --- 產生測資與高斯雜訊 (維持原樣) ---
double generate_gaussian(double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 <= 0.0) u1 = 1e-9;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    return z0 * sigma;
}

void Generate_TestData(const int count, int N, int M, double EbN0_dB) {
    char output_name[] = "C:\\Users\\user\\Documents\\project\\test_data.txt";
    FILE *f_out = fopen(output_name, "w");
    if (!f_out) {
        perror("Failed to open output file");
        return;
    }

    double rate = (double)(N - M) / N;
    double snr = pow(10.0, EbN0_dB / 10.0);
    double sigma = sqrt(1.0 / (2.0 * rate * snr));

    srand((unsigned int)time(NULL));
    int *codeword = (int *)calloc(N, sizeof(int));

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
    if (x < 1e-12) return 30.0;
    double value = -log(tanh(x / 2.0));
    if (value > 30.0) value = 30.0;
    if (value < 1e-12) value = 1e-12;
    return value;
}

bool isDecodeSucess(double* channel_sigma, int* decoder_iteration, int decoder_max_iteration, double* channel_output) {
    double* LLR_Channel = (double*)malloc(N * sizeof(double));
    double* LLR_Current = (double*)malloc(N * sizeof(double));
    double* phibeta = (double*)malloc(N * sizeof(double));
    int* Ci = (int*)malloc(N * sizeof(int));

    for (int i = 0; i < N; i++) {
        LLR_Channel[i] = 2.0 * channel_output[i] / ((*channel_sigma) * (*channel_sigma));
        LLR_Current[i] = LLR_Channel[i];
    }

    *decoder_iteration = 0;
    static int for_test_only_true_time = 0;

    while (*decoder_iteration < decoder_max_iteration) {
        long long truebit = 0, falsebit = 0;
        double ReturnTotalExceptI[N];
        for (int i = 0; i < N; i++) ReturnTotalExceptI[i] = 0;

        for (int j = 0; j < N; j++) phibeta[j] = phi(fabs(LLR_Current[j]));

        // Check Node Update (使用一維陣列存取)
        for (int i = 0; i < M; i++) {
            int TotalNegativeCount = 0;
            double TotalPhi = 0;
            int offset = H_rows[i].offset;

            for (int j = 0; j < H_rows[i].row_weight; j++) {
                int col_idx = H_row_places_flat[offset + j]; // 改成一維索引
                if (LLR_Current[col_idx] < 0) TotalNegativeCount++;
                TotalPhi += phibeta[col_idx];
            }

            for (int j = 0; j < H_rows[i].row_weight; j++) {
                int col_idx = H_row_places_flat[offset + j];
                int sign = ((LLR_Current[col_idx] < 0 ? TotalNegativeCount - 1 : TotalNegativeCount) % 2 == 0) ? 1 : -1;
                double except_phi = TotalPhi - phibeta[col_idx];
                ReturnTotalExceptI[col_idx] += sign * phi(except_phi);
            }
        }

        bool isAllZero = true;
        for (int i = 0; i < N; i++) {
            double LLR_FinalResult = LLR_Channel[i] + ReturnTotalExceptI[i];
            Ci[i] = (LLR_FinalResult > 0) ? 0 : 1;
            if (Ci[i] == 0) truebit++; else falsebit++;
            LLR_Current[i] = LLR_FinalResult; 
        }

        // 校驗 (使用一維陣列存取)
        for (int i = 0; i < M; i++) {
            int decision = 0;
            int offset = H_rows[i].offset;
            for (int j = 0; j < H_rows[i].row_weight; j++) {
                decision += Ci[H_row_places_flat[offset + j]];
            }
            if (decision % 2 != 0) isAllZero = false;
        }

        if (isAllZero == true) {
            GLOBAL_TRUEBIT += truebit; GLOBAL_FALSEBIT += falsebit;
            for_test_only_true_time++;
            free(LLR_Channel); free(LLR_Current); free(phibeta); free(Ci);
            return true;
        } else {
            if (*decoder_iteration == decoder_max_iteration - 1) {
                GLOBAL_TRUEBIT += truebit; GLOBAL_FALSEBIT += falsebit;
                free(LLR_Channel); free(LLR_Current); free(phibeta); free(Ci);
                return false;
            }
        }
        *decoder_iteration += 1;
    }
    free(LLR_Channel); free(LLR_Current); free(phibeta); free(Ci);
    return false;
}

// --- 主程式 ---
int main(void) {
    Read_H_MatrixByAlist();

    const int num_test_cases = 100;
    const int decoder_iteration_max = 50; 
    const double EbN0_dB = 5.0;

    int *TEST_countEachIteration = (int*)malloc(num_test_cases * sizeof(int));
    for (int i = 0; i < num_test_cases; i++) TEST_countEachIteration[i] = -1;
    
    int totalSucessTime = 0, totalFailTime = 0;
    GLOBAL_TRUEBIT = 0; GLOBAL_FALSEBIT = 0;

    Generate_TestData(num_test_cases, N, M, EbN0_dB);

    char filename[100] = {"C:\\Users\\user\\Documents\\project\\test_data.txt"};
    FILE *TEST_DATA = fopen(filename, "r");
    if (!TEST_DATA) return 1;

    double *channel_output = (double*)malloc(N * sizeof(double));
    double *channel_sigma = (double*)malloc(sizeof(double));
    int *decode_iteration = (int*)malloc(sizeof(int));

    for (int testtime = 0; testtime < num_test_cases; testtime++) {
        double dummy;
        for (int j = 0; j < N; j++) fscanf(TEST_DATA, "%lf", &dummy);
        fscanf(TEST_DATA, "%lf", channel_sigma);
        for (int j = 0; j < N; j++) fscanf(TEST_DATA, "%lf", &channel_output[j]);

        bool isThisIterationSucess = isDecodeSucess(channel_sigma, decode_iteration, decoder_iteration_max, channel_output);

        if (isThisIterationSucess) {
            totalSucessTime += 1;
        } else {
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