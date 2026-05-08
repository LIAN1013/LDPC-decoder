#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {
    int N, M, max_col, max_row;
    int *col_weight, *row_weight;
    int **col_index, **row_index;
} Alist;


double initial(double yi, double sigma);
double cal_r_for_fnode(Alist *matrix, int c_idx, int exclude_v, int M_size, double q1_matrix[][M_size]);
int syndrome(Alist* matrix, int c[]);

Alist* readalist(const char *matrix) {
    FILE *fp = fopen(matrix, "r");
    if (!fp) return NULL;

    Alist *alist = (Alist *)malloc(sizeof(Alist));
    fscanf(fp, "%d %d", &alist->N, &alist->M);
    fscanf(fp, "%d %d", &alist->max_col, &alist->max_row);

    alist->col_weight = (int *)malloc(alist->N * sizeof(int));
    alist->row_weight = (int *)malloc(alist->M * sizeof(int));

    for (int i = 0; i < alist->N; i++) {
        fscanf(fp, "%d", &alist->col_weight[i]);
    }

    for (int i = 0; i < alist->M; i++) {
        fscanf(fp, "%d", &alist->row_weight[i]);
    }

    alist->col_index = (int **)malloc(alist->N * sizeof(int *));
    for (int i = 0; i < alist->N; i++) {
        alist->col_index[i] = (int *)malloc(alist->max_col * sizeof(int));
        for (int j = 0; j < alist->max_col; j++) {
            fscanf(fp, "%d", &alist->col_index[i][j]);
        }
    }

    alist->row_index = (int **)malloc(alist->M * sizeof(int *));
    for (int i = 0; i < alist->M; i++) {
        alist->row_index[i] = (int *)malloc(alist->max_row * sizeof(int));
        for (int j = 0; j < alist->max_row; j++) {
            fscanf(fp, "%d", &alist->row_index[i][j]);
        }
    }

    fclose(fp);
    return alist;
}

int main() {
    Alist *Hmatrix = readalist("test.txt");
    if (!Hmatrix) return 1;

    double sigma;
    int received_len = Hmatrix->N;
    double y[received_len + 50], caly[received_len + 50];
    int N_size = Hmatrix->N + 50, M_size = Hmatrix->M + 50;
    int max_iter = 15;
    int c[Hmatrix->N];

    printf("enter sigma^2:");
    scanf("%lf", &sigma);

    int count = 0;
    printf("enter received codeword:");
    while (count < received_len) {
        if (scanf("%lf", &y[count]) != 1 || y[count] == -1) break;
        caly[count] = initial(y[count], sigma);
        count++;
    }

    //不能用 {0} 初始化，改用 memset
    double q1_matrix[N_size][M_size];
    memset(q1_matrix, 0, sizeof(q1_matrix));

    for (int i = 0; i < Hmatrix->N; i++) {
        for (int k = 0; k < Hmatrix->col_weight[i]; k++) {
            int temp = Hmatrix->col_index[i][k] - 1;
            if (temp >= 0) {
                q1_matrix[i][temp] = caly[i];
            }
        }
    }
    //printf("caly done\n");
    int iter = 0;
    while (iter < max_iter) {
        double r_matrix[M_size][N_size];
        memset(r_matrix, 0, sizeof(r_matrix));
        //printf("in while\n");
        for (int j = 0; j < Hmatrix->M; j++) {
            for (int k = 0; k < Hmatrix->row_weight[j]; k++) {
                int i_vnode = Hmatrix->row_index[j][k] - 1;
                if (i_vnode >= 0) {
                    r_matrix[j][i_vnode] = cal_r_for_fnode(Hmatrix, j, i_vnode, M_size, q1_matrix);
                }
               // printf("in for\n");
            }
            //printf("in 2for\n");
        }
    //printf("r done\n");
        for (int i = 0; i < Hmatrix->N; i++) {
            double Q1_all = caly[i];
            double Q0_all = 1.0 - caly[i];

            for (int k = 0; k < Hmatrix->col_weight[i]; k++) {
                int row_j = Hmatrix->col_index[i][k] - 1;
                if (row_j >= 0) {
                    Q0_all *= r_matrix[row_j][i];
                    Q1_all *= (1.0 - r_matrix[row_j][i]);
                }
            }
                //printf("Q done\n");
            double Q1_final = Q1_all / (Q0_all + Q1_all);
            c[i] = (Q1_final > 0.5) ? 1 : 0;
            printf("位元 %d: P(1)=%f, 判定為 %d\n", i, Q1_final, c[i]);

            for (int k = 0; k < Hmatrix->col_weight[i]; k++) {
                int target_row = Hmatrix->col_index[i][k] - 1;
                if (target_row < 0) continue;

                double q1_prod = caly[i];
                double q0_prod = 1.0 - caly[i];

                for (int u = 0; u < Hmatrix->col_weight[i]; u++) {
                    int neighbor_row = Hmatrix->col_index[i][u] - 1;
                    if (neighbor_row >= 0 && neighbor_row != target_row) {
                        q0_prod *= r_matrix[neighbor_row][i];
                        q1_prod *= (1.0 - r_matrix[neighbor_row][i]);
                    }
                }
                q1_matrix[i][target_row] = q1_prod / (q0_prod + q1_prod);
            }
        }

        int judge = syndrome(Hmatrix, c);
        if (judge == 0) break;

        printf("\n");
        iter++;
    }

    printf("End\n");
    return 0;
}

double initial(double yi, double sigma) {
    double b = 2.0 * yi / sigma;
    double a = 1.0 / (1.0 + exp(b));
    return a;
}

double cal_r_for_fnode(Alist *matrix, int c_idx, int exclude_v, int M_size, double q1_matrix[][M_size]) {
    double b, c = 1.0;
    int row_degree = matrix->row_weight[c_idx];

    for (int i = 0; i < row_degree; i++) {
        int v_neighbor = matrix->row_index[c_idx][i] - 1;
        if (v_neighbor >= 0 && v_neighbor != exclude_v) {
            b = 1.0 - (2.0 * q1_matrix[v_neighbor][c_idx]);
            c *= b;
        }
    }
    return (0.5 + 0.5 * c);
}

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
