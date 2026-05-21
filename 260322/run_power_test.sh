#!/bin/bash
# =====================================================================
# run_power_test.sh (全自動CPU cycle收集腳本)
# =====================================================================

# 定義要測試的 Dataset 和 Method
DATASETS=("email" "fetal" "gene" "mush" "loan" "gesture" "arcene" "gas")
METHODS=("seq" "naive" "mcp" "lsh" "radix" "fptree" "union")

LOG_FILE="power_results_log.txt"

# 清空舊的 Log 檔
echo "==========================================================" > $LOG_FILE
echo " Power Measurement Results " >> $LOG_FILE
echo " Date: $(date)" >> $LOG_FILE
echo "==========================================================" >> $LOG_FILE

# 雙重迴圈：跑遍所有 Dataset 與所有 Method
for ds in "${DATASETS[@]}"; do
    echo "=========================================================="
    echo "▶️ Testing Dataset: $ds"
    echo "=========================================================="
    
    for method in "${METHODS[@]}"; do
        echo "   -> Running $method ..."
        
        # 在 Log 中加上標題
        echo "" >> $LOG_FILE
        echo "[Dataset: $ds | Method: $method]" >> $LOG_FILE
        
        # 執行 perf 並把標準錯誤輸出 (stderr) 附加到 Log 檔中
        # 注意：perf 的輸出是在 stderr，所以用 2>> 重導向
        sudo perf stat -e cycles:u Rscript measure_packing.R $ds $method 2>> $LOG_FILE
        
    done
    echo "Done with $ds."
done

echo ""
echo "🎉 All tests completed! Please check '$LOG_FILE' for the power consumption results."
