HammingCodeEX是H矩陣的ALIST檔案  

LogDominSPADecoder是SPA原始檔案，MinSum也是，先保留不要更動，如果要修正上傳新的檔案不要直接覆蓋(我電腦沒有存，直接開始嘗試優化了  

MinSum_AutoCreateData包含了minsum和自動創建合法codeword(目前侷限在全0合法code，已知會有正確率過高問題)，可以調整雜訊SNR->Eb 表示合法碼字通過AWGN的待解碼碼字，可以調整運行次數看錯誤比例大概是多少(目前是設定1M次，會生成一個大約100MB的測試檔案，要調整檔案路徑)，目前測試在我的專題PPT內H(8,4)矩陣在codeword=全零下設置SNR>=4左右幾乎有98%可以成功解碼

