#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

// 生成高斯雜訊 (DO NOT CHANGE!!!)
double generate_gaussian(double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 <= 0.0) u1 = 1e-9;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    return z0 * sigma;
}

// 副函式：讀取 TXT 檔案並驗算統計特性
void verify_gaussian_distribution(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("錯誤：無法開啟檔案 %s\n", filename);
        return;
    }

    char line[256];
    long long total_samples = 0;
    double bin_start, count;

    // 建立暫存結構來儲存讀取的數據
    // 因為我們知道 num_bins 大約是 4000 (或更多)，可以動態分配
    int capacity = 12800;
    double *bins = malloc(capacity * sizeof(double));
    double *counts = malloc(capacity * sizeof(double));
    int n = 0;

    // 跳過標頭，尋找數據起始點
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "總執行次數:")) {
            sscanf(line, "總執行次數: %lld", &total_samples);
        }
        if (strstr(line, "區間下限,出現次數")) {
            break; // 找到數據開頭
        }
    }

    // 讀取數據點
    while (fscanf(fp, "%lf,%lf", &bin_start, &count) == 2) {
        bins[n] = bin_start + 0.025; // 使用區間中點進行計算 (0.05/2)
        counts[n] = count;
        n++;
    }
    fclose(fp);

    // --- 開始統計計算 ---
    double mean = 0, var = 0, skew = 0, kurt = 0;
    double sum_counts = 0;

    // 1. 計算平均值 (1st Moment)
    for (int i = 0; i < n; i++) {
        mean += bins[i] * counts[i];
        sum_counts += counts[i];
    }
    mean /= sum_counts;

    // 2. 計算變異數 (2nd), 偏度 (3rd), 峰度 (4th)
    for (int i = 0; i < n; i++) {
        double diff = bins[i] - mean;
        var  += counts[i] * pow(diff, 2);
        skew += counts[i] * pow(diff, 3);
        kurt += counts[i] * pow(diff, 4);
    }
    var  /= sum_counts;
    skew /= (sum_counts * pow(var, 1.5));
    kurt /= (sum_counts * pow(var, 2));

    // --- 輸出驗證結果 ---
    printf("\n========== 檔案驗算報告 ==========\n");
    printf("檔案名稱: %s\n", filename);
    printf("樣本總數: %.0f\n", sum_counts);
    printf("----------------------------------\n");
    printf("1. 平均值 (Mean)   : %10.6f (理論值: 0.0)\n", mean);
    printf("2. 變異數 (Variance): %10.6f (理論值: 1.0)\n", var);
    printf("3. 偏度   (Skewness): %10.6f (理論值: 0.0)\n", skew);
    printf("4. 峰度   (Kurtosis): %10.6f (理論值: 3.0)\n", kurt);
    printf("----------------------------------\n");

    // 簡單邏輯判定
    if (fabs(mean) < 0.01 && fabs(var - 1.0) < 0.01 && fabs(skew) < 0.05) {
        printf("結果判定：符合高斯分佈 (PASS)\n");
    } else {
        printf("結果判定：不符合預期 (FAIL)\n");
    }
    printf("==================================\n");

    free(bins);
    free(counts);
}

int main(void) {
    // --- 參數設定 ---
    const double sigma = 1.0;
    const long long iterations = 1e7;
    const double min_limit = -10.0;
    const double max_limit = 10.0;
    const double step = 0.05;
    const int num_bins = (int)((max_limit - min_limit) / step);

    double min_val = DBL_MAX;
    double max_val = -DBL_MAX;

    unsigned int *histogram = (unsigned int *)calloc(num_bins, sizeof(unsigned int));
    if (histogram == NULL) {
        printf("記憶體配置失敗！\n");
        return 1;
    }

    srand((unsigned int)time(NULL));
    printf("計算中並寫入檔案...\n");

    clock_t start_time = clock();

    // 1. 進行模擬與統計
    for (long long i = 0; i < iterations; i++) {
        double val = generate_gaussian(sigma);
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;

        if (val >= min_limit && val < max_limit) {
            int index = (int)((val - min_limit) / step);
            histogram[index]++;
        }
    }

    // 2. 輸出至 TXT 檔案
    FILE *fp = fopen("gaussian_results.txt", "w");
    if (fp == NULL) {
        printf("無法建立檔案！\n");
        free(histogram);
        return 1;
    }

    // 在檔案開頭寫入摘要資訊
    fprintf(fp, "=== 高斯雜訊統計結果 ===\n");
    fprintf(fp, "總執行次數: %lld\n", iterations);
    fprintf(fp, "實際最小值: %f\n", min_val);
    fprintf(fp, "實際最大值: %f\n\n", max_val);
    fprintf(fp, "區間下限,出現次數\n");

    // 逐行寫入矩陣數據
    for (int i = 0; i < num_bins; i++) {
        double bin_start = min_limit + (i * step);
        fprintf(fp, "%.2f,%u\n", bin_start, histogram[i]);
    }

    fclose(fp); // 務必關閉檔案

    clock_t end_time = clock();
    double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("成功！結果已儲存至 gaussian_results.txt\n");
    printf("總耗時: %f 秒\n", cpu_time_used);

    free(histogram);


    verify_gaussian_distribution("gaussian_results.txt");

    return 0;
}
