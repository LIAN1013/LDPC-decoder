#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//Alist輸入格式(儲存H矩陣)
typedef struct {
    int N, M, max_col, max_row;
    int *col_weight, *row_weight;
    int **col_index, **row_index;
} Alist;

//通道初始值函式(1/1+exp(2*yi/sigma^2))=>qij(1)=Pi(1)的初始值
double initial(double yi, double sigma);

//計算從檢查節點到變數節點的訊息函式
double cal_r_for_fnode(Alist *matrix, int c_idx, int exclude_v, int M_size, double *q1_matrix);

int syndrome(Alist* matrix, int c[]);

//讀取Alist格式的H矩陣
Alist* readalist(const char *matrix) {
    FILE *fp = fopen(matrix, "r");

    if (!fp) return NULL;

    Alist *alist = (Alist *)malloc(sizeof(Alist));

    // 讀取N(column)、M(row)、max_col(最多1的行數)、max_row(最多1的列數)
    fscanf(fp, "%d %d", &alist->N, &alist->M);
    fscanf(fp, "%d %d", &alist->max_col, &alist->max_row);

    //記錄每一行和每一列1的權重
    alist->col_weight = (int *)malloc(alist->N * sizeof(int));
    alist->row_weight = (int *)malloc(alist->M * sizeof(int));

    for (int i = 0; i < alist->N; i++) {
        fscanf(fp, "%d", &alist->col_weight[i]);
    }

    for (int i = 0; i < alist->M; i++) {
        fscanf(fp, "%d", &alist->row_weight[i]);
    }

    //二維陣列儲存每一行1的位置(從0開始)
    alist->col_index = (int **)malloc(alist->N * sizeof(int *));

    for (int i = 0; i < alist->N; i++) {
        //alist->max_col => alist->col_weight[i]，因為每行的1的數量不一定，所以分配col_weight[i]的空間即可
        alist->col_index[i] = (int *)malloc(alist->col_weight[i] * sizeof(int));
        for (int j = 0; j < alist->col_weight[i]; j++) {
            fscanf(fp, "%d", &alist->col_index[i][j]);
        }
    }

    //二維陣列儲存每一列1的位置(從0開始)
    alist->row_index = (int **)malloc(alist->M * sizeof(int *));

    for (int i = 0; i < alist->M; i++) {
        //alist->max_row => alist->row_weight[i]，因為每列的1的數量不一定，所以分配row_weight[i]的空間即可
        alist->row_index[i] = (int *)malloc(alist->row_weight[i] * sizeof(int));
        for (int j = 0; j < alist->row_weight[i]; j++) {
            fscanf(fp, "%d", &alist->row_index[i][j]);
        }
    }

    fclose(fp);
    return alist;
}

int main() {

    clock_t start_read, end_read, start_run, end_run;
    //讀取計時開始
    start_read = clock();

    //讀取Alist格式的H矩陣
    Alist *Hmatrix = readalist("AlistFile01.txt");
    end_read = clock();
    if (!Hmatrix) {
        printf("failed to read file\n");
        return 1;
    }

    //定義sigma^2、接收的碼字長度(N)、接收的碼字陣列、計算用的陣列、矩陣大小、最大迭代次數(15)、解碼後的碼字陣列
    //移除50的額外空間，因為在這裡我們不需要額外的空間來存儲接收的碼字或計算用的陣列，直接使用N大小的陣列即可，這樣可以節省記憶體並避免不必要的複雜性。
    double sigma;
    int received_len = Hmatrix->N;
    double y[received_len], caly[received_len];
    int N_size = Hmatrix->N, M_size = Hmatrix->M;
    int max_iter = 15;
    int c[Hmatrix->N];

    printf("enter sigma^2:");
    scanf("%lf", &sigma);

    int count = 0;
    printf("enter received codeword:");
    while (count < received_len) {
        //移除 y[count] == -1 除非有特定需求，否則不應該在這裡使用 -1 作為結束輸入的標誌，因為 -1 可能是有效的接收值。建議改為使用 EOF 或其他非數字的輸入來結束輸入。
        if (scanf("%lf", &y[count]) != 1) break;

        //計算初始值並存入caly陣列
        caly[count] = initial(y[count], sigma);
        count++;
    }

    //讀取計時結束
    

    //迭代計時開始
    

    //分配記憶體給q1_matrix並初始化為0 (一維陣列模擬二維陣列，大小為N_size * M_size)
    double *q1_matrix = (double *)calloc(N_size * M_size, sizeof(double));

    start_run = clock();
    for (int i = 0; i < Hmatrix->N; i++) {
        for (int k = 0; k < Hmatrix->col_weight[i]; k++) {
            // col_index[i][k] 是從1開始的，所以要減1來對應到0-based的索引
            int temp = Hmatrix->col_index[i][k] - 1;
            if (temp >= 0) {
                // 存取方式：[i * M_size + temp] 二維轉一維
                q1_matrix[i * M_size + temp] = caly[i];
            }
        }
    }

    //分配記憶體給r_matrix，並初始化為0(calloc)
    double *r_matrix = (double *)calloc(M_size * N_size, sizeof(double));
    memset(r_matrix, 0, M_size * N_size * sizeof(double));
    int iter = 0;
    while (iter < max_iter) {
        //每次迭代開始前，將r_matrix重置為0 
        //memset的參數說明：第一個參數是要填充的記憶體區域的指標，第二個參數是要設定的值（0表示將所有位元設為0），第三個參數是要設定的位元數量（這裡是M_size * N_size * sizeof(double)，表示整個r_matrix的大小）。
        printf("%d iteration\n", iter);

        for (int j = 0; j < Hmatrix->M; j++) {
            for (int k = 0; k < Hmatrix->row_weight[j]; k++) {
                // row_index[j][k] 是從1開始的，所以要減1來對應到0-based的索引
                int i_vnode = Hmatrix->row_index[j][k] - 1;
                if (i_vnode >= 0) {
                    // 存取方式：[j * N_size + i_vnode] (rij(0)的值)，二維轉一維
                    r_matrix[j * N_size + i_vnode] = cal_r_for_fnode(Hmatrix, j, i_vnode, M_size, q1_matrix);
                }
            }
        }

        //找Qi來進行判決
        for (int i = 0; i < Hmatrix->N; i++) {
            double Q1_all = caly[i];//Pi(1)=>Pi
            double Q0_all = 1.0 - caly[i];//Pi(0)=>1-Pi

            for (int k = 0; k < Hmatrix->col_weight[i]; k++) {
                int row_j = Hmatrix->col_index[i][k] - 1;
                if (row_j >= 0) {
                    Q0_all *= r_matrix[row_j * N_size + i];//(1-Pi)*rji(0)的值
                    Q1_all *= (1.0 - r_matrix[row_j * N_size + i]);//Pi*rji(1)的值
                }
            }

            //hard decision，判定為1的機率大於0.5則判定為1，否則為0
            double Q1_final = Q1_all / (Q0_all + Q1_all);
            c[i] = (Q1_final > 0.5) ? 1 : 0;
            printf("位元 %d: P(1)=%f, 判定為 %d\n", i, Q1_final, c[i]);

            //排除掉對應的rji(0)和rji(1)來更新qij(1)的值，這裡的qij(1)是從變數節點i傳給檢查節點j的訊息，更新後的qij(1)不包含來自檢查節點j的訊息（Extrinsic information），這樣可以避免迴圈中的資訊循環問題。
            //可透過Q0_all和Q1_all來計算更新後的qij(1)值，因為Q0_all和Q1_all已經包含了所有rji(0)和rji(1)的連乘結果，所以我們可以直接使用這兩個值來計算更新後的qij(1)，而不需要再進行一次迴圈來排除掉對應的rji(0)和rji(1)。這樣可以簡化程式碼並提高效率。
            //可簡化
            for (int k = 0; k < Hmatrix->col_weight[i]; k++) {
                int target_row = Hmatrix->col_index[i][k] - 1;

                //可省略?
                if (target_row < 0) continue;

                double q1_prod = caly[i];
                double q0_prod = 1.0 - caly[i];

                for (int u = 0; u < Hmatrix->col_weight[i]; u++) {
                    int neighbor_row = Hmatrix->col_index[i][u] - 1;
                    if (neighbor_row >= 0 && neighbor_row != target_row) {
                        q0_prod *= r_matrix[neighbor_row * N_size + i];
                        q1_prod *= (1.0 - r_matrix[neighbor_row * N_size + i]);
                    }
                }
                q1_matrix[i * M_size + target_row] = q1_prod / (q0_prod + q1_prod);
            }
        }

        int judge = syndrome(Hmatrix, c);
        if (judge == 0) break;

        printf("\n");
        iter++;
    }
    end_run = clock();

    free(q1_matrix);
    free(r_matrix);

    printf("讀取與初始化時間: %f 秒\n", (double)(end_read - start_read) / CLOCKS_PER_SEC);
    printf("解碼運算執行時間: %f 秒\n", (double)(end_run - start_run) / CLOCKS_PER_SEC);
    return 0;
}

//通道初始值函式(1/1+exp(2*yi/sigma^2))=>qij(1)=Pi(1)的初始值
double initial(double yi, double sigma) {
    double b = 2.0 * yi / sigma;
    double a = 1.0 / (1.0 + exp(b));
    return a;
}

//計算從檢查節點傳給變數節點的訊息函式
double cal_r_for_fnode(Alist *matrix, int c_idx, int exclude_v, int M_size, double *q1_matrix) {
    double b, c = 1.0;
    int row_degree = matrix->row_weight[c_idx];

    for (int i = 0; i < row_degree; i++) {
        int v_neighbor = matrix->row_index[c_idx][i] - 1;

        //排除exclude_v，因為這是我們要計算的目標變數節點，不應該包含在計算中(Extrinsic information)
        //移除 v_neighbor >= 0 的檢查，因為在Alist格式中，v_neighbor應該永遠不會是負數，如果出現負數可能表示資料錯誤或索引問題，建議在讀取Alist時就確保資料的正確性，而不是在這裡進行檢查。
        if (v_neighbor != exclude_v) {
            // 一維存取：v_neighbor 為列索引，c_idx 為欄索引
            //計算1-2qi'j的連乘，這裡的q1_matrix[v_neighbor * M_size + c_idx]是從變數節點v_neighbor傳給檢查節點c_idx的訊息
            b = 1.0 - (2.0 * q1_matrix[v_neighbor * M_size + c_idx]);
            c *= b;
        }
    }
    //rji(0)=0.5+0.5*c，rji(1)=0.5-0.5*c，這裡我們只需要rji(0)的值來更新變數節點的訊息，所以返回0.5+0.5*c即可，c是所有1-2qi'j的連乘結果
    return (0.5 + 0.5 * c);
}

//計算syndrome的值，判斷是否解碼成功，syndrome為0表示解碼成功，為1表示解碼失敗
int syndrome(Alist* matrix, int c[]) {
    for (int j = 0; j < matrix->M; j++) {
        int row_sum = 0;
        for (int i = 0; i < matrix->row_weight[j]; i++) {
            int v_idx = matrix->row_index[j][i] - 1;
            if (v_idx >= 0) {
                row_sum += c[v_idx];
            }
        }
        if (row_sum % 2 != 0) return 1;
    }
    return 0;
}
