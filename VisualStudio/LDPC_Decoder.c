#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

#include "DecoderTool.h"

// initial LLR=2yi/(sigma^2) data
void initial(double* channel_sigma, double* channel_output, double* LLR_ChannelCurrent, double* LLR_FinalResult)
{
    for (int i = 0; i < N; i++)
    {
        LLR_ChannelCurrent[i] = 2.0 * channel_output[i] / ((*channel_sigma) * (*channel_sigma));
        LLR_FinalResult[i] = LLR_ChannelCurrent[i];
    }
}

// --- log domain decoder main function ---
bool Decode_SumProduct(double* channel_sigma, int* decoder_iteration, int decoder_max_iteration, double* channel_output) {
    double* LLR_Channel = (double*)malloc(N * sizeof(double));
    double* LLR_Current = (double*)malloc(N * sizeof(double));
    double* phibeta = (double*)malloc(N * sizeof(double));
    int* Ci = (int*)malloc(N * sizeof(int));

    
    initial(channel_sigma, channel_output, LLR_Channel, LLR_Current);

    *decoder_iteration = 0;

    while (*decoder_iteration < decoder_max_iteration) 
    {
        long long truebit = 0, falsebit = 0;
        double* ReturnTotalExceptI = (double*)malloc(N * sizeof(double));
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
            free(LLR_Channel); free(LLR_Current); free(phibeta); free(Ci);
            return true;
        }
        else {
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


// --- minsum-c decoder main function ---



