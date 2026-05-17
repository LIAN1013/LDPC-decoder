#pragma once
#ifndef DECODERTOOL_H
#define DECODERTOOL_H


#ifndef PATH_PREFIX
#define PATH_PREFIX "C:\\Users\\09043\\OneDrive\\Python\\LDPC\\"
#endif

#define TEST_DATA_FILE PATH_PREFIX "test_data.txt"
#define H_MATRIX_FILE  PATH_PREFIX "GallagerH.txt"




// Defination
typedef struct {
    int row_weight;
    int offset; // 改為紀錄在一維陣列中的起始索引
} RowInfo;

typedef struct {
    int column_weight;
    int offset; // 改為紀錄在一維陣列中的起始索引
} ColInfo;

// 用 extern 讓其他.c檔案看見
extern int N, K, M;
extern RowInfo* H_rows;
extern ColInfo* H_cols;

extern int* H_row_places_flat;
extern int* H_col_places_flat;

extern long long GLOBAL_TRUEBIT, GLOBAL_FALSEBIT;



#endif // ! DECODERTOOL_H
