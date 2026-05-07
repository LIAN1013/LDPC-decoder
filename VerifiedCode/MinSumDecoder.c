#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

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

int N=-1,K=-1,M=-1;    // H(8,4) hamming code
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

// 這邊還沒做到
// FAILED to save now!!!!!
void Save_H_MatrixByAlist(void)
{
    char filename[]={"C:\\C\\TESTZONE\\testfile.txt"};
    FILE *fp=fopen(filename,"w");

    if(!fp){
        printf("Open file failed\n");
        return;
    }

    int max_row_weight=0;
    int max_col_weight=0;

    fprintf(fp,"%d %d\n",N,M);
    fprintf(fp,"%d %d\n",max_row_weight,max_col_weight);


    // dim1
    for(int i=0;i<N;i++){
        fprintf(fp,"%d",H_cols[i].column_weight);
    }
    fprintf(fp,"\n");
    for(int i=0;i<M;i++){
        fprintf(fp,"%d",H_rows[i].row_weight);
    }
    fprintf(fp,"\n");

    // dim2
    for(int i=0;i<N;i++){
        for(int j=0;j<H_cols[i].column_weight;j++){
            fprintf(fp,"%d",H_cols[i].column_place[j]);
        }
        fprintf(fp,"\n");
    }

    for(int i=0;i<M;i++){
        for(int j=0;j<H_rows[i].row_weight;j++){
            fprintf(fp,"%d",H_rows[i].row_place[j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    printf("WRITE Finish\n");
    return;
}

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

    printf("->N=%d<- ->M=%d<-\n",N,M);
    printf("->%d<- ->%d<-\n",max_row_weight,max_col_weight);

    H_rows=(RowInfo*)malloc(M*sizeof(RowInfo));
    H_cols=(ColInfo*)malloc(N*sizeof(ColInfo));

    // dim1
    for(int i=0;i<N;i++){
        fscanf(fp,"%d",&H_cols[i].column_weight);
        printf("->%d<- ",H_cols[i].column_weight);
    }
    printf("\n");
    for(int i=0;i<M;i++){
        fscanf(fp,"%d",&H_rows[i].row_weight);
        printf("->%d<- ",H_rows[i].row_weight);
    }
    printf("\n");

    // dim2
    int temp=0;
    for(int i=0;i<N;i++){
        H_cols[i].column_place=(int*)malloc(H_cols[i].column_weight*sizeof(int));
        for(int j=0;j<H_cols[i].column_weight;j++){
            fscanf(fp,"%d",&temp);
            H_cols[i].column_place[j]=temp-1;
            printf("->%d<- ",H_cols[i].column_place[j]);
        }
        printf("\n");
    }

    for(int i=0;i<M;i++){
        H_rows[i].row_place=(int*)malloc(H_rows[i].row_weight*sizeof(int));
        for(int j=0;j<H_rows[i].row_weight;j++){
            fscanf(fp,"%d",&temp);
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
    double channel_sigma=0.5;
    int decoder_iteration=0;
    const int decoder_max_iteration=1;

    /** read file **/
    Read_H_MatrixByAlist();

    /** decoder **/
    if(N==-1 || M==-1)printf("READ Alist Failed\n");

    double LLR_Current[N],LLR_FinalResult[N];
    int Ci[N];

    printf("ITERATION %d LLR_Current: \n",decoder_iteration);
    for(int i=0;i<N;i++)
    {
        LLR_Current[i]=2.0*channel_output[i]/channel_sigma;
        LLR_FinalResult[i]=LLR_Current[i];
        printf("%llf ",LLR_Current[i]);
    }printf("\n");


    while(decoder_iteration<decoder_max_iteration)
    {
        /// VARIABLE
        // step 1&2
        double ReturnTotal[N],ReturnTotalExcpetI[N];
        /// INITIAL
        for(int i=0;i<N;i++){ReturnTotal[i]=0;ReturnTotalExcpetI[i]=0;}

        for(int i=0;i<M;i++)//EX 3 iterations for hamming code
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
            printf("ITERATION %d-%d FIRST SETP\n",decoder_iteration,i);
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
                    printf("1.Value < First < Second -> Replace First By Value\n");
                    printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
                else if(tempValue>=first.minBetaValue && tempValue<second.minBetaValue)
                {
                    second.minBetaValue=tempValue;
                    second.minBetaPosition=H_rows[i].row_place[j];
                    printf("2.First < Value < Second -> Replace Second By Value\n");
                    printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
                else{
                    printf("3.First < second < Value\n");
                    printf("FIRST: %lf  SECOND: %lf  Value: %lf\n",first.minBetaValue,second.minBetaValue,tempValue);
                }
            }printf("\n");

            // 第二遍遍歷計算出用於最終結果的回傳值(以總負號和總phi數值計算)和繼續迭代的LLR(以總負號扣除該位置負號數量判斷正負，以總phi值扣除該位置phi數值判斷數值大小)
            printf("ITERATION %d-%d SECOND SETP\n",decoder_iteration,i);
            for(int j=0;j<H_rows[i].row_weight;j++)
            {
                // alpha
                int alpha=TotalNegativeCount%2==0?1:-1;
                // beta
                int tempPosition=H_rows[i].row_place[j];
                if(first.minBetaPosition!=tempPosition)
                {
                    ReturnTotalExcpetI[tempPosition]+=alpha*isNegativeI[j]*first.minBetaValue;
                    ReturnTotal[tempPosition]+=alpha*first.minBetaValue;
                }
                else
                {
                    ReturnTotalExcpetI[tempPosition]+=alpha*isNegativeI[j]*second.minBetaValue;
                    ReturnTotal[tempPosition]+=alpha*first.minBetaValue;
                }
                printf("At %d Point -> alpha=%d  isNegative=%d  firstVal=%lf  secondVal=%lf\n",tempPosition,alpha,isNegativeI[j],first.minBetaValue,second.minBetaValue);
                printf("At %d Point -> ReturnTotal=%lf  ReturnTotalExceptI=%lf \n",tempPosition,ReturnTotal[tempPosition],ReturnTotalExcpetI[tempPosition]);
            }
            for(int j=0;j<N;j++){printf("%d-> ReturnTotal->%llf  ReturnTotalExceptI->%llf\n",j,ReturnTotal[j],ReturnTotalExcpetI[j]);}
        }

        // step 3 check decoder finish or not
        printf("STEP3 in %d iteration\n",decoder_iteration);
        for(int i=0;i<N;i++)
        {
            printf("LLR Before->%llf   ",LLR_Current[i]);
            // 不可交換
            LLR_FinalResult[i]=LLR_Current[i]+ReturnTotal[i];
            LLR_Current[i]+=ReturnTotalExcpetI[i];
            printf("LLR After->%llf\n",LLR_Current[i]);
            Ci[i]=LLR_FinalResult[i]>0?0:1;
            printf("RESULT:%d  FINALRESULT:LLR(Qi)=%llf  LLR(qij)=%llf Ci[%d]=%d\n",i,LLR_FinalResult[i],LLR_Current[i],i,Ci[i]);
        }
        printf("\n\n");
        // decision
        int decision[M];
        bool isAllZero=true;

        // Hc^T
        for(int i=0;i<M;i++){
            // initial
            decision[i]=0;
            for(int j=0;j<H_rows[i].row_weight;j++){
                decision[i]+=Ci[H_rows[i].row_place[j]];
            }
            if(decision[i]%2==1)isAllZero=false;
            printf("%d ",decision[i]);
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

    return 0;
}
