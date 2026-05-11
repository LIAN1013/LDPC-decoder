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
    int row_weight;//該列包含幾個1
    int *row_place;//記錄這些1在哪些行(column索引)
}RowInfo;

typedef struct{
    int column_weight;//該行包含幾個1
    int *column_place;//記錄這些1在哪些列(row索引)
}ColInfo;

//用於minSum decoder的beta值紀錄結構，紀錄最小值和次小值的數值和位置
typedef struct{
    double minBetaValue;
    int minBetaPosition;
}minBeta;

int N=-1,K=-1,M=-1;    // H(8,4) hamming code
                       //N=Vnode=8(總位元數/行數), M=Cnode=4(校驗方程式數量/列數)
                       //K=N-M=4(資料位元數)
RowInfo *H_rows=NULL;
ColInfo *H_cols=NULL;

// 目前接收端暫時寫死，可改成副函式開檔讀取
double channel_output[20]={-1.5,0.8,-0.9,0.7,0.5,-1.1,-0.4,-1.2};
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
    char filename[100]={"C:\\Users\\09043\\OneDrive\\Python\\LDPC\\HammingCodeEX.txt"};
    //char filename[100]={"HammingCodeEX.txt"};//相對路徑，放在專案資料夾內
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


    H_rows=(RowInfo*)malloc(M*sizeof(RowInfo));
    H_cols=(ColInfo*)malloc(N*sizeof(ColInfo));

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

int main(void)
{
    // 時間設置
    clock_t StartReadTime,EndReadTime,StartDecodeTime,EndDecodeTime;

    // sigma 設為0.5是為了讓LLR的數值不要太大，方便觀察，實際上可以根據需要調整
    const double channel_sigma=sqrt(0.5);
    const int decoder_max_iteration=1;
    int decoder_iteration=0;

    /** read file **/
    StartReadTime=clock();
    Read_H_MatrixByAlist();
    EndReadTime=clock();

    /** decoder **/ //讀取失敗
    if(N==-1 || M==-1)printf("READ Alist Failed\n");

    double LLR_Current[N],LLR_FinalResult[N];
    int Ci[N];

    StartDecodeTime=clock();
    printf("ITERATION %d LLR_Current: \n",decoder_iteration);
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
        double ReturnTotal[N],ReturnTotalExceptI[N];

        /// INITIAL(calloc?最後可釋放空間)
        for(int i=0;i<N;i++){ReturnTotal[i]=0;ReturnTotalExceptI[i]=0;}

        for(int i=0;i<M;i++)//EX 3 iterations for hamming code
        {
            /// VARIABLE
            // alpha// 統計負號數量用於判斷回傳值正負，isNegativeI紀錄每個位置是否為負號(1 or -1)
            int TotalNegativeCount=0, isNegativeI[H_rows[i].row_weight];
            // min beta
            minBeta first,second;

            /// INITIAL
            for(int j=0;j<H_rows[i].row_weight;j++){isNegativeI[j]=1;}

            // 初始化最小值和次小值為一個很大的數，確保在比較時能被替換掉
            first.minBetaValue=1e100;
            second.minBetaValue=1e100;

            // 第一遍遍歷每一row，先統計每一row中的負號位置和最小值與次小值數值大小與位置
            printf("ITERATION %d-%d FIRST SETP\n",decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)// EX 2,3,4,5 in Hamming code first line
            {
                // alpha->紀錄第i列第j行的值是否為負號，並統計總負號數量
                if(LLR_Current[H_rows[i].row_place[j]]<0){
                    TotalNegativeCount+=1;
                    isNegativeI[j]=-1;
                }
                // beta
                double tempValue=fabs(LLR_Current[H_rows[i].row_place[j]]);
                // 更新最小值和次小值，如果當前值比最小值還小，就把最小值更新為當前值，並把原來的最小值更新為次小值；
                if(tempValue<first.minBetaValue && tempValue<second.minBetaValue)
                {
                    second=first;

                    first.minBetaValue=tempValue;
                    first.minBetaPosition=H_rows[i].row_place[j];
                    printf("1.Value < First < Second -> Replace First By Value\n");
                    printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
                // 如果當前值比最小值大但比次小值還小，就把次小值更新為當前值，位置也更新為當前值的位置；
                else if(tempValue>=first.minBetaValue && tempValue<second.minBetaValue)
                {
                    second.minBetaValue=tempValue;
                    second.minBetaPosition=H_rows[i].row_place[j];
                    printf("2.First < Value < Second -> Replace Second By Value\n");
                    printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
                // 如果當前值比最小值和次小值都大，則不更新最小值和次小值。(這裡的else是可省略的，因為如果不滿足前兩個條件，就不需要做任何操作)
                else{
                    printf("3.First < second < Value\n");
                    printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
            }printf("\n");

            // 第二遍遍歷計算出用於最終結果的回傳值(以總負號和總phi數值計算)和繼續迭代的LLR(以總負號扣除該位置負號數量判斷正負，以總phi值扣除該位置phi數值判斷數值大小)
            printf("ITERATION %d-%d SECOND SETP\n",decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)
            {
                // alpha(總負號數量)
                int alpha=TotalNegativeCount%2==0?1:-1;
                // beta
                int tempPosition=H_rows[i].row_place[j];
                // 如果當前位置不是最小值的位置，則LLR回傳值的數值為最小值
                if(first.minBetaPosition!=tempPosition)
                {
                    // 繼續迭代的LLR值 L(rij)
                    ReturnTotalExceptI[tempPosition]+=alpha*isNegativeI[j]*first.minBetaValue;
                    // 最終結果的回傳值
                    ReturnTotal[tempPosition]+=alpha*first.minBetaValue;
                }
                // 如果當前位置是最小值的位置，則LLR回傳值的數值為次小值
                else
                {
                    // 繼續迭代的LLR值 L(rij)
                    ReturnTotalExceptI[tempPosition]+=alpha*isNegativeI[j]*second.minBetaValue;
                    // 最終結果的回傳值
                    ReturnTotal[tempPosition]+=alpha*first.minBetaValue;
                }
                printf("At %d Point -> alpha=%d  isNegative=%d  firstVal=%lf  secondVal=%lf\n",tempPosition,alpha,isNegativeI[j],first.minBetaValue,second.minBetaValue);
                printf("At %d Point -> ReturnTotal=%lf  ReturnTotalExceptI=%lf \n",tempPosition,ReturnTotal[tempPosition],ReturnTotalExceptI[tempPosition]);
            }
            for(int j=0;j<N;j++){printf("%d-> ReturnTotal->%llf  ReturnTotalExceptI->%llf\n",j,ReturnTotal[j],ReturnTotalExceptI[j]);}
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
            printf("row %d: %d ",i,decision[i]);
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


    EndDecodeTime=clock();

    printf("Read H Matrix Time = %g second\n",(double)(EndReadTime-StartReadTime)/CLOCKS_PER_SEC);
    printf("Decode Time = %g second\n",(double)(EndDecodeTime-StartDecodeTime)/CLOCKS_PER_SEC);

    return 0;
}
