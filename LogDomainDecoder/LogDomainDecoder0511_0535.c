#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

/**
    INFORMATION:
    1. L(qij)=L(ci)=2yi/sigma^2 in BI-AWGN Channel
    2. Row 列(水平) Column 行(垂直)
    3. Example use (8,4) Min-Sum Decoder example in RYAN_LDPC
        (ppt example in page 36-37)
    4. PRINTF add -><- to prevent unexpected blank or \n
    5.
**/


// store H matrix
typedef struct{
    int row_weight; //該列包含幾個1
    int *row_place; //記錄這些1在哪些行(column索引)
}RowInfo;

typedef struct{
    int column_weight; //該行包含幾個1
    int *column_place; //記錄這些1在哪些列(row索引)
}ColInfo;

//用於minSum decoder的beta值紀錄結構，紀錄最小值和次小值的數值和位置
/*typedef struct{
    double minBetaValue;
    int minBetaPosition;
}minBeta;*/

int N=-1,K=-1,M=-1;    // H(8,4) hamming code
                       //N=Vnode=8(總位元數/行數), M=Cnode=4(校驗方程式數量/列數)
                       //K=N-M=4(資料位元數)
RowInfo *H_rows=NULL;
ColInfo *H_cols=NULL;

// 目前接收端暫時寫死，可改成副函式開檔讀取
//double channel_output[20]={0.2,0.2,-0.9,0.6,0.5,-1.1,-0.4,-1.2};//proability example in RYAN_LDPC 60-61
double channel_output[20]={-1.5,0.8,-0.9,0.7,0.5,-1.1,-0.4,-1.2};//min sum decoder example in RYAN_LDPC 74-75
//double channel_output[20]={0.6,-2,2.87,1.45,-3.4,-2.6,-0.7,3.2};


// test {}
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
    char filename[100]={"AlistFile01.txt"};//相對路徑，放在專案資料夾內
    FILE *fp=fopen(filename,"r");

    if(!fp){
        printf("Open file failed\n");
        return;
    }

    int max_col_weight=0;
    int max_row_weight=0;

    fscanf(fp,"%d %d",&N,&M);
    fscanf(fp,"%d %d",&max_col_weight,&max_row_weight);//col <-> row

    // 測試是否正確讀取
    printf("->N=%d<- ->M=%d<-\n",N,M);
    printf("->%d<- ->%d<-\n",max_col_weight,max_row_weight);//col <-> row

    H_cols=(ColInfo*)malloc(N*sizeof(ColInfo));
    H_rows=(RowInfo*)malloc(M*sizeof(RowInfo));

    // dim1(讀取每一行包含幾個1的資訊)
    for(int i=0;i<N;i++){
        fscanf(fp,"%d",&H_cols[i].column_weight);
        printf("->%d<- ",H_cols[i].column_weight);
    }
    printf("\n");
    // dim1(讀取每一列包含幾個1的資訊)
    for(int i=0;i<M;i++){
        fscanf(fp,"%d",&H_rows[i].row_weight);
        printf("->%d<- ",H_rows[i].row_weight);
    }
    printf("\n");

    // dim2(讀取每一行包含1的位置資訊(row索引))
    int temp=0;
    for(int i=0;i<N;i++){
        H_cols[i].column_place=(int*)malloc(H_cols[i].column_weight*sizeof(int));
        for(int j=0;j<H_cols[i].column_weight;j++){
            fscanf(fp,"%d",&temp);
            // Alist中位置是從1開始的，但C中的矩陣從0開始存，所以要-1轉成從0開始的索引
            H_cols[i].column_place[j]=temp-1;
            printf("->%d<- ",H_cols[i].column_place[j]);
        }
        printf("\n");
    }
    // dim2(讀取每一列包含1的位置資訊(column索引))
    for(int i=0;i<M;i++){
        H_rows[i].row_place=(int*)malloc(H_rows[i].row_weight*sizeof(int));
        for(int j=0;j<H_rows[i].row_weight;j++){
            fscanf(fp,"%d",&temp);
            // Alist中位置是從1開始的，但C中的矩陣從0開始存，所以要-1轉成從0開始的索引
            H_rows[i].row_place[j]=temp-1;
            printf("->%d<- ",H_rows[i].row_place[j]);
        }
        printf("\n");
    }
    fclose(fp);
    printf("Read Finish\n");
    return;
}
// phi(x)=-ln(tanh(x/2))
double phi(double x)
{
    /*//1.方法一
    double value=x;//進入前已經取絕對值了，所以這裡直接用x就好
    value=-log(tanh(value/2));

    //防止數值過小導致的計算不穩定
    if(x<1e-12)
        value=1e-12;
    //*/

    // 直接計算，x 在外部應保證為正數
    double value = -log(tanh(x / 2.0));

    // 防護 1: 避免數值發散到無限大 (當 x 非常接近 0 時)
    if (value > 30.0) {
        value = 30.0; // 設定一個合理的上限值 (Clipping)
    }

    // 防護 2: 避免數值完全為 0 (當 x 非常大時)
    if (value < 1e-12) {
        value = 1e-12;
    }

    return value;
}

int main(void)
{
    // sigma 設為0.5是為了讓LLR的數值不要太大，方便觀察
    const double channel_sigma=sqrt(0.5);
    int decoder_iteration=0;
    const int decoder_max_iteration=30;

    // --- 【新增】計時用變數 ---
    clock_t start_time, end_time;
    double read_time_used, decode_time_used;
    // -------------------------

    /** read file **/
    start_time = clock(); // 開始記錄讀取時間
    Read_H_MatrixByAlist();
    end_time = clock();   // 結束記錄讀取時間

    // 計算讀取時間 (轉換為秒並帶有小數點)
    read_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    // 計算讀取時間 (轉換為秒並帶有小數點)
    printf("[計時] 讀取 Alist 檔案耗時: %llf 秒\n", read_time_used);

    /** decoder **/ //讀取失敗
    if(N==-1 || M==-1)printf("READ Alist Failed\n");

    double LLR_Current[N],LLR_FinalResult[N];
    int Ci[N];

    printf("ITERATION %d LLR_Current: \n",decoder_iteration);

    start_time = clock(); // 開始記錄解碼時間
    for(int i=0;i<N;i++)
    {
        // sigma -> sigma^2
        LLR_Current[i]=2.0*channel_output[i]/(channel_sigma*channel_sigma);
        LLR_FinalResult[i]=LLR_Current[i];
        printf("%llf ",LLR_Current[i]);
    }printf("\n");


    while(decoder_iteration<decoder_max_iteration)
    {
        /// VARIABLE
        // step 1(Check Node update) &2(Variable Node update)
        double ReturnTotal[N],ReturnTotalExceptI[N],phibeta[N];

        // INITIAL(calloc?最後可釋放空間)
        for(int i=0;i<N;i++){ReturnTotal[i]=0,ReturnTotalExceptI[i]=0;}

        // Count phi beta i'j before summation
        // 由於phi數在一次迭代中不會更動，因此再迭代前先提前計算
        printf("ITERATION %d PHIBETA Pre Calculate\n",decoder_iteration);
        for(int j=0;j<N;j++){
            phibeta[j]=phi(fabs(LLR_Current[j]));
            printf("phi%d->%llf  ",j,phibeta[j]);
        }printf("\n");



        for(int i=0;i<M;i++)//EX 3 iterations for hamming code
        {
            /// VARIABLE
            // alpha // 統計負號數量用於判斷回傳值正負，isNegativeI紀錄每個位置是否為負號(1 or -1)
            int TotalNegativeCount=0, isNegativeI[H_rows[i].row_weight];
            // beta
            double TotalPhi=0;

            /// INITIAL
            for(int j=0;j<H_rows[i].row_weight;j++){isNegativeI[j]=1;}

            // 第一遍遍歷每一row，先統計每一row中的負號位置和計算總負號數量與總phi數值
            printf("ITERATION %d-%d FIRST SETP\n",decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)// EX 2,3,4,5 in Hamming code first line
            {
                // 紀錄第i列第j行的值是否為負號，並統計總負號數量
                if(LLR_Current[H_rows[i].row_place[j]]<0){
                    TotalNegativeCount+=1;
                    isNegativeI[j]=-1;
                }
                // 總phi數值計算(所有資訊)
                TotalPhi+=phibeta[H_rows[i].row_place[j]];
                printf("%dTotalphi=%llf  ",j,TotalPhi);
            }printf("\n");

            // 第二遍遍歷計算出用於最終結果的回傳值(以總負號和總phi數值計算)和繼續迭代的LLR alpha(以總負號數量-該位置正負號)->判斷正負，以(總phi值-該位置phi數值)->判斷數值大小
            printf("ITERATION %d-%d SECOND SETP\n",decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)
            {
                // 總負號數量
                int alpha=TotalNegativeCount%2==0?1:-1;

                // 以(總phi值-該位置phi數值)->判斷數值大小
                double totalPhiExceptI=TotalPhi-phibeta[H_rows[i].row_place[j]];

                double returnValue=alpha*phi(TotalPhi);

                // 紀錄每個位置的回傳值和繼續迭代的LLR值 L(rij)
                double returnValueExceptI=alpha*isNegativeI[j]*phi(totalPhiExceptI);
                printf("%d TotalPhiExceptI->%llf  ",j,totalPhiExceptI);

                ReturnTotal[H_rows[i].row_place[j]]+=returnValue;
                ReturnTotalExceptI[H_rows[i].row_place[j]]+=returnValueExceptI;
                //新增紀錄row和column位置(均由矩陣序列+1)的回傳值，方便debug
                printf("Row %d column %d ReturnValue->%llf Row %d column %d ReturnValueExcept%d->%llf\n",i+1,H_rows[i].row_place[j]+1,returnValue,i+1,H_rows[i].row_place[j]+1,returnValueExceptI);
            }
            for(int j=0;j<N;j++)
            {
                //j -> column(已校正)
                printf("column=%d TOTAL:%llf  TOTALEXCEPT:%llf\n",j+1,ReturnTotal[j],ReturnTotalExceptI[j]);
            }
            printf("\n");
        }

        // step 3 check decoder finish or not
        printf("STEP3 in %d iteration\n",decoder_iteration);
        for(int i=0;i<N;i++)
        {
            printf("LLR Before->%llf   ",LLR_Current[i]);

            // 不可交換
            // L(Qi)=L(ci)+sum(L(rji)) 作為硬性判定用的最終結果
            LLR_FinalResult[i]=LLR_Current[i]+ReturnTotal[i];

            // L(qij)=L(ci)+sum(L(rj'i)) 作為繼續迭代用的LLR值
            LLR_Current[i]+=ReturnTotalExceptI[i];
            printf("LLR After->%llf\n",LLR_Current[i]);
            Ci[i]=LLR_FinalResult[i]>0?0:1;
            printf("RESULT:%d  FINALRESULT:LLR(Qi)=%llf  LLR(qij)=%llf Ci[%d]=%d\n",i,LLR_FinalResult[i],LLR_Current[i],i,Ci[i]);
        }
        printf("\n\n");
        // decision
        int decision[M];
        bool isAllZero=true;

        // Hc^T(校驗檢查)來判斷是否成功解碼，若成功(x*H^T=0)則停止迭代，失敗(x*H^T=1)則繼續迭代
        for(int i=0;i<M;i++){
            // initial
            decision[i]=0;
            for(int j=0;j<H_rows[i].row_weight;j++){
                decision[i]+=Ci[H_rows[i].row_place[j]];
            }
            if(decision[i]%2==1)isAllZero=false;
            //新增紀錄每一列的decision值，方便debug
            printf("row %d:%d ",i+1,decision[i]);
        }
        if(isAllZero==true){
            printf("DECODE SUCESS\n");
            break;
        }
        else{
            printf("DECODE FAILED in %d iteration\n",decoder_iteration);
        }
        decoder_iteration+=1;
    }
    end_time = clock();   // 結束記錄解碼時間

    // 計算並印出解碼時間
    decode_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("\n[計時] 解碼運算耗時: %llf 秒\n", decode_time_used);
    return 0;
}
