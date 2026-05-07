#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

/**

    V 一切以H矩陣優先，先讀取H矩陣
    X 根據H矩陣生成G矩陣
    X 根據G矩陣生成多組合法codeword
    V 根據合法codeword生成多筆經過AWGN通道的測試資訊y(通道雜訊變異數不同)
    V 根據y解碼，檢查解碼是否成功
    X 計算時間效率

**/


/**
    INFORMATION:
    1. L(qij)=L(ci)=2yi/sigma^2 in BI-AWANG Channel
    2. Row 列(水平) Column 行(垂直)
    3. Example use (8,4) Min-Sum Decoder example in RYAN_LDPC
        (ppt example in page 36-37)
    4. PRINTF add -><- to prevent unexpected blank or \n
    5.
**/



// store H matrix
typedef struct{
    int row_weight;
    int *row_place;
}RowInfo;
typedef struct{
    int column_weight;
    int *column_place;
}ColInfo;

typedef struct{
    double minBetaValue;
    int minBetaPosition;
}minBeta;

int N=-1,K=-1,M=-1;
RowInfo *H_rows=NULL;
ColInfo *H_cols=NULL;



// Alist Structure
/**

    Alist structure

    N M (N=Vnode, M=Cnode, K=N-M)
    max_col_weight max_row_weight
    Col1_RowWeight Col2_RowWeight... ...Col20_RowWeight
    Row1_ColWeight Row2_ColWeight... ...Row16_ColWeight
    Col1_Row1 Col1_Row2... ...Col1_RowWeight
    Col2_Row1 Col2_Row2... ...Col2_RowWeight
    ...
    Col20_Row1 Col20_Row2... ...Col20_RowWeight
    Row1_Col1 Row1_Col2... ...Row1_ColWeight
    Row2_Col1 Row2_Col2... ...Row2_ColWeight
    ...
    Row16_Col1 Row16_Col2... ...Row16_ColWeight


 **/

// OK
void Read_H_MatrixByAlist(void)
{
    char filename[100]={"C:\\Users\\09043\\OneDrive\\Python\\LDPC\\HammingCodeEX.txt"};
    FILE *fp=fopen(filename,"r");

    if(!fp){
        printf("Open file failed\n");
        return;
    }

    int max_row_weight=0;
    int max_col_weight=0;

    fscanf(fp,"%d %d",&N,&M);
    fscanf(fp,"%d %d",&max_row_weight,&max_col_weight);

    //printf("->N=%d<- ->M=%d<-\n",N,M);
    //printf("->%d<- ->%d<-\n",max_row_weight,max_col_weight);

    H_rows=(RowInfo*)malloc(M*sizeof(RowInfo));
    H_cols=(ColInfo*)malloc(N*sizeof(ColInfo));

    // dim1
    for(int i=0;i<N;i++){
        fscanf(fp,"%d",&H_cols[i].column_weight);
        //printf("->%d<- ",H_cols[i].column_weight);
    }
    printf("\n");
    for(int i=0;i<M;i++){
        fscanf(fp,"%d",&H_rows[i].row_weight);
        //printf("->%d<- ",H_rows[i].row_weight);
    }
    //printf("\n");

    // dim2
    int temp=0;
    for(int i=0;i<N;i++){
        H_cols[i].column_place=(int*)malloc(H_cols[i].column_weight*sizeof(int));
        for(int j=0;j<H_cols[i].column_weight;j++){
            fscanf(fp,"%d",&temp);
            H_cols[i].column_place[j]=temp-1;
            //printf("->%d<- ",H_cols[i].column_place[j]);
        }
        //printf("\n");
    }

    for(int i=0;i<M;i++){
        H_rows[i].row_place=(int*)malloc(H_rows[i].row_weight*sizeof(int));
        for(int j=0;j<H_rows[i].row_weight;j++){
            fscanf(fp,"%d",&temp);
            H_rows[i].row_place[j]=temp-1;
            //printf("->%d<- ",H_rows[i].row_place[j]);
        }
        //printf("\n");
    }
    fclose(fp);
    //printf("Read Finish\n");
    return;
}


// 輔助函式：將 Dense 矩陣寫入 Alist 檔案
void Write_Dense_To_Alist(const char* filename, int **matrix, int rows, int cols) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Failed to create file: %s\n", filename);
        return;
    }

    int *col_weights = (int *)calloc(cols, sizeof(int));
    int *row_weights = (int *)calloc(rows, sizeof(int));
    int max_col_w = 0, max_row_w = 0;

    // 1. 統計權重
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] == 1) {
                row_weights[i]++;
                col_weights[j]++;
            }
        }
    }
    for (int j = 0; j < cols; j++) if (col_weights[j] > max_col_w) max_col_w = col_weights[j];
    for (int i = 0; i < rows; i++) if (row_weights[i] > max_row_w) max_row_w = row_weights[i];

    // 2. 寫入 N M 與最大權重
    fprintf(fp, "%d %d\n%d %d\n", cols, rows, max_col_w, max_row_w);

    // 3. 寫入各欄與各列權重
    for (int j = 0; j < cols; j++) fprintf(fp, "%d ", col_weights[j]);
    fprintf(fp, "\n");
    for (int i = 0; i < rows; i++) fprintf(fp, "%d ", row_weights[i]);
    fprintf(fp, "\n");

    // 4. 寫入 Column-wise 非零位置 (1-based)
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            if (matrix[i][j] == 1) fprintf(fp, "%d ", i + 1);
        }
        fprintf(fp, "\n");
    }

    // 5. 寫入 Row-wise 非零位置 (1-based)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] == 1) fprintf(fp, "%d ", j + 1);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    free(col_weights);
    free(row_weights);
}

// 核心副函式：優化 H 矩陣並寫入檔案
void Optimize_H_To_Alist() {
    // 1. 讀取原始 H
    Read_H_MatrixByAlist();
    if (N <= 0 || M <= 0) return;

    int K = N - M;
    int i, j, r, c;

    // 2. 建立 Dense 矩陣進行高斯消去
    int **h_dense = (int **)malloc(M * sizeof(int *));
    for (i = 0; i < M; i++) {
        h_dense[i] = (int *)calloc(N, sizeof(int));
        for (j = 0; j < H_rows[i].row_weight; j++) {
            h_dense[i][H_rows[i].row_place[j]] = 1;
        }
    }

    // 3. 高斯-約旦消去法：將右側 M x M 區域化為單位矩陣 I
    for (i = 0; i < M; i++) {
        int pivot_col = K + i; // 目標是化簡最後 M 個 Columns
        int pivot_row = i;

        // A. 尋找該 Column 中的 Pivot (1)
        while (pivot_row < M && h_dense[pivot_row][pivot_col] == 0) pivot_row++;

        if (pivot_row < M) {
            // 交換 Rows
            int *temp = h_dense[i];
            h_dense[i] = h_dense[pivot_row];
            h_dense[pivot_row] = temp;

            // B. 消去該 Column 中其他所有 Row 的 1 (GF(2) XOR)
            for (r = 0; r < M; r++) {
                if (r != i && h_dense[r][pivot_col] == 1) {
                    for (c = 0; c < N; c++) {
                        h_dense[r][c] ^= h_dense[i][c];
                    }
                }
            }
        } else {
            printf("Warning: Right part of H is not Full Rank at Col %d. Gaussian failed.\n", pivot_col);
        }
    }

    // 4. 寫入新的 Alist 檔案
    char output_path[] = "C:\\Users\\09043\\OneDrive\\Python\\LDPC\\H_Optimized_Alist.txt";
    Write_Dense_To_Alist(output_path, h_dense, M, N);

    printf("Optimized H [P|I] has been written to: %s\n", output_path);

    // 5. 清理記憶體
    for (i = 0; i < M; i++) free(h_dense[i]);
    free(h_dense);
}




// (DO NOT CHANGE!!!)生成高斯雜訊
double generate_gaussian(double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 <= 0.0) u1 = 1e-9; // 避免 log(0)
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    return z0 * sigma;
}
// (DO NOT CHANGE!!!)
void Generate_TestData(const int count, int N, int M, double EbN0_dB) {
    char output_name[] = "C:\\Users\\09043\\OneDrive\\Python\\LDPC\\test_data.txt";
    FILE *f_out = fopen(output_name, "w");
    if (!f_out) {
        perror("Failed to open output file");
        return;
    }

    // 計算 Sigma
    double rate = (double)(N - M) / N;
    double snr = pow(10.0, EbN0_dB / 10.0);
    double sigma = sqrt(1.0 / (2.0 * rate * snr));

    srand((unsigned int)time(NULL));

    int *codeword = (int *)calloc(N, sizeof(int));
    for (int c = 0; c < count; c++) {
        // --- 1. 生成 Codeword (這裡暫用全零，若有 G 矩陣可改為隨機碼) ---
        // 提示：若要測試不同碼字，可隨機生成位元，但須確保滿足 Hc^T = 0
        for (int i = 0; i < N; i++) {
            fprintf(f_out, "%d%s", codeword[i], (i == N - 1) ? "" : " ");
        }
        fprintf(f_out, "\n");

        // --- 2. 寫入 Sigma ---
        fprintf(f_out, "%lf\n", sigma);

        // --- 3. 生成過通道後的 y (BPSK: 0 -> +1.0, 1 -> -1.0) ---
        for (int i = 0; i < N; i++) {
            double s = (codeword[i] == 0) ? 1.0 : -1.0;
            double noise = generate_gaussian(sigma);
            double y = s + noise;
            fprintf(f_out, "%lf%s", y, (i == N - 1) ? "" : " ");
        }
        fprintf(f_out, "\n");


        // 每一組中間可以加個空行方便肉眼閱讀（可選）
        // fprintf(f_out, "\n");
    }
    free(codeword);

    fclose(f_out);
    printf("Successfully generated %d sets of data with Sigma: %lf\n", count, sigma);
}






// 解碼步驟
bool isDecodeSucess(double* channel_sigma,int* decoder_iteration,int decoder_max_iteration, double* channel_output)
{
    /** read file **/
    //Read_H_MatrixByAlist();

    /** decoder **/
    if(N==-1 || M==-1)printf("[SYSTEM] READ Alist Failed\n");

    double LLR_Current[N],LLR_FinalResult[N];
    int Ci[N];

    //printf("ITERATION %d LLR_Current: \n",*decoder_iteration);
    for(int i=0;i<N;i++)
    {
        LLR_Current[i]=2.0*channel_output[i]/( (*channel_sigma) * (*channel_sigma) );
        LLR_FinalResult[i]=LLR_Current[i];
        //printf("%llf ",LLR_Current[i]);
    }
    //printf("\n");

    *decoder_iteration=0;
    while(*decoder_iteration<decoder_max_iteration)
    {
        /// VARIABLE
        double ReturnTotal[N],ReturnTotalExcpetI[N];
        /// INITIAL
        for(int i=0;i<N;i++){ReturnTotal[i]=0;ReturnTotalExcpetI[i]=0;}

        for(int i=0;i<M;i++)
        {
            /// VARIABLE
            // alpha
            int TotalNegativeCount=0, isNegativeI[H_rows[i].row_weight];
            // beta
            minBeta first,second;

            /// INITIAL
            for(int j=0;j<H_rows[i].row_weight;j++){isNegativeI[j]=1;}
            first.minBetaValue=1e100;
            second.minBetaValue=1e100;

            // 第一遍遍歷每一row，先統計每一row中的負號位置和最小值與次小值數值大小與位置
            //printf("ITERATION %d-%d-%d FIRST SETP\n",testtime,*decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)// EX 2,3,4,5 in Hamming code first line
            {
                // alpha
                if(LLR_Current[H_rows[i].row_place[j]]<0){
                    TotalNegativeCount+=1;
                    isNegativeI[j]=-1;
                }
                // beta
                double tempValue=fabs(LLR_Current[H_rows[i].row_place[j]]);
                if(tempValue<first.minBetaValue && tempValue<second.minBetaValue)
                {
                    second=first;

                    first.minBetaValue=tempValue;
                    first.minBetaPosition=H_rows[i].row_place[j];
                    //printf("1.Value < First < Second -> Replace First By Value\n");
                    //printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
                else if(tempValue>=first.minBetaValue && tempValue<second.minBetaValue)
                {
                    second.minBetaValue=tempValue;
                    second.minBetaPosition=H_rows[i].row_place[j];
                    //printf("2.First < Value < Second -> Replace Second By Value\n");
                    //printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
                else{
                    //printf("3.First < second < Value\n");
                    //printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
            }
            //printf("\n");

            // 第二遍遍歷計算出用於最終結果的回傳值(以總負號和總phi數值計算)和繼續迭代的LLR(以總負號扣除該位置負號數量判斷正負，以總phi值扣除該位置phi數值判斷數值大小)
            //printf("ITERATION %d-%d-%d SECOND SETP\n",testtime,*decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)
            {
                ///
                // alpha
                int alpha=TotalNegativeCount%2==0?1:-1;
                // beta
                int tempPosition=H_rows[i].row_place[j];


                if(first.minBetaPosition!=tempPosition){
                    ReturnTotalExcpetI[tempPosition]+=alpha*isNegativeI[j]*first.minBetaValue;
                    ReturnTotal[tempPosition]+=alpha*first.minBetaValue;
                }
                else{
                    ReturnTotalExcpetI[tempPosition]+=alpha*isNegativeI[j]*second.minBetaValue;
                    ReturnTotal[tempPosition]+=alpha*first.minBetaValue;
                }
                //printf("At %d Point -> alpha=%d  isNegative=%d  firstVal=%lf  secondVal=%lf\n",tempPosition,alpha,isNegativeI[j],first.minBetaValue,second.minBetaValue);
                //printf("At %d Point -> ReturnTotal=%lf  ReturnTotalExceptI=%lf \n",tempPosition,ReturnTotal[tempPosition],ReturnTotalExcpetI[tempPosition]);
            }
            //for(int j=0;j<N;j++){printf("%d-> ReturnTotal->%llf  ReturnTotalExceptI->%llf\n",j,ReturnTotal[j],ReturnTotalExcpetI[j]);}
        }

        /// step 3 check decoder finish or not
        //printf("STEP3 in %d-%d iteration\n",testtime,*decoder_iteration);
        for(int i=0;i<N;i++)
        {
            // 不可交換
            LLR_FinalResult[i]=LLR_Current[i]+ReturnTotal[i];
            LLR_Current[i]+=ReturnTotalExcpetI[i];
            //printf("LLR After->%llf\n",LLR_Current[i]);
            Ci[i]=LLR_FinalResult[i]>0?0:1;
            //printf("RESULT:%d  FINALRESULT:LLR(Qi)=%llf  LLR(qij)=%llf Ci[%d]=%d\n",i,LLR_FinalResult[i],LLR_Current[i],i,Ci[i]);
        }

        /// FINAL DECISION
        bool isAllZero=true;

        // check H*c^T==0
        for(int i=0;i<M;i++){
            bool decision=0;
            for(int j=0;j<H_rows[i].row_weight;j++){
                decision+=Ci[H_rows[i].row_place[j]];
            }
            isAllZero=decision!=0?false:true;
        }
        if(isAllZero==true){
            return true;
        }
        else{
            //printf("DECODE FAILED in %d iteration\n",*decoder_iteration);
            if(*decoder_iteration==decoder_max_iteration-1)return false;
        }
        *decoder_iteration+=1;
    }
}

int main(void)
{
    // Read H Matrix
    Read_H_MatrixByAlist();
    //Optimize_H_To_Alist();

    // u can change the constant below >_<
    const int num_test_cases = 1000000;
    const int deocder_iteration_max = 500;
    const double EbN0_dB = -10.0;

    int *TEST_countEachIteration=(int*)malloc(num_test_cases * sizeof(int));

    // initial
    for(int i=0;i<num_test_cases;i++){TEST_countEachIteration[i]=-1;}
    int totalSucessTime=0,totalFailTime=0;

    // Generate TestData
    Generate_TestData(num_test_cases, N, M, EbN0_dB);

    // open test data file
    char filename[100]={"C:\\Users\\09043\\OneDrive\\Python\\LDPC\\test_data.txt"};
    FILE *TEST_DATA=fopen(filename,"r");
    if(!TEST_DATA){
        printf("[SYSTEM]MAIN_Open file failed\n");
        return false;
    }

    // allocate memory for y(after AWGN channel)
    double *channel_output=(double*)malloc(N*sizeof(double));
    double *channel_sigma=(double*)malloc(sizeof(double));
    int *decode_iteration=(int*)malloc(sizeof(int));

    //
    for(int testtime=0;testtime<num_test_cases;testtime++)
    {
        //
        for(int j=0;j<N;j++){
            fscanf(TEST_DATA,"%lf",&channel_output[j]);
            //printf("%lf ",channel_output[j]);
        }
        fscanf(TEST_DATA,"%lf",channel_sigma);
        //printf("MAIN Channel_sigma = %lf \n",*channel_sigma);

        for(int j=0;j<N;j++){
            fscanf(TEST_DATA,"%lf",&channel_output[j]);
            //printf("%lf ",channel_output[j]);
        }
        bool isThisIterationSucess=isDecodeSucess(channel_sigma,decode_iteration,deocder_iteration_max,channel_output);


        if(isThisIterationSucess==true){
            //printf("MAIN_Decode Sucess In %d Iterations\n",*decode_iteration);
            totalSucessTime+=1;
        }
        else{
            //printf("MAIN_Decode Failed In %d Iterations\n",*decode_iteration);
            totalFailTime+=1;
        }
        TEST_countEachIteration[testtime]=*decode_iteration;

        // print each decoder result
        // printf("MAIN_ ");
        // for(int j=0;j<N;j++){printf(" %lf ",channel_output[j]);}
        // printf("\n");
        if(testtime%10000==0)printf("FINISH %8d TOTAL %d\n",testtime,num_test_cases);
    }

    /// print final result
    //printf("\n\n\nMAIN_IF TIME==-1 -> Decode Failed!\n");
    //for(int g=0;g<num_test_cases;g++){printf("MAIN_In %d Time Try %d Iteration\n",g,TEST_countEachIteration[g]);}
    printf("MAIN_SucessTime = %d   Failed Time = %d\n",totalSucessTime,totalFailTime);
    printf("MAIN_CorrectRate = %lf \n",(float)totalSucessTime/(float)num_test_cases);

    return 0;
}
