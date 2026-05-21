# 260412

In this week, we changed the simulate function, using `simulate_hardware_inference` to replace the old LCA-based simulate method `simulate_cam_array_mismatch` to evaluate the performance of the testing dataset in ACAM tiles. 


## 🏗️ 核心結構差異 (Structural Differences)

* **舊版 (`simulate_cam_array_mismatch`) - 邏輯拓樸模擬：**
    * **機制：** 依賴 Random Forest 的決策樹拓樸結構，透過尋找測試資料落點與目標路徑的「最低共同祖先 (LCA, Lowest Common Ancestor)」來判斷 Mismatch 節點。
* **新版 (`simulate_hardware_inference`) - 實體硬體模擬：**
    * **機制：** 徹底捨棄 LCA 邏輯。直接讀取依照演算法排佈好的「實體硬體映射表 (`hw_state_table`)」。
    * **優勢：** 透過 R 語言底層的矩陣向量化運算 (Vectorization, 如 `sweep` 與 `max.col`)，嚴格依照硬體的物理順序（由左至右的 Tile 排列）對測試資料進行掃描。這能完美捕捉硬體真實的「提早關閉 (Early Deactivation)」行為，運算速度更快且符合真實物理現象。

## 📦 新增的 Return 內容與對應作用 (New Returns & Roles)

新版本的 `main_simulation` 流程中，除了原本的 `row_tiles` 與 `physical_tiles` 外，我們新增並豐富了以下關鍵的輸出結構，為後續的功耗計算鋪路：

### 1. `res_static$hardware_state_tables` (硬體實體映射表)
* **內容：** 將全域的 ACAM 狀態映射到真實硬體位置上。紀錄了在特定 Packing 方法下，每一個 Array (`array_id`)、每一個實體列 (`row_in_array`) 具體存放了哪個特徵，以及該特徵對應的 ACAM 物理寫入區間 (`min_state`, `max_state`)。
* **作用：** **用於計算「編程功耗 (Programming Power)」。** 透過這張表，可以直接套用單一 Cell 的寫入耗能公式，精準算出整顆晶片載入模型時的總寫入成本。

### 2. `res_dynamic$mismatch_matrix` (全域 Mismatch 追蹤矩陣)
* **內容：** 一個覆蓋所有樣本 (Samples) 與所有路徑 (Paths) 的巨型追蹤矩陣。裡面的數值代表該筆資料在「第幾個橫向 Tile」發生了物理上的第一次 Mismatch；若完美 Match 則顯示為 `NA`。
* **作用：** **用於建立「逐週期動態功耗模型 (Cycle-by-cycle Dynamic Power Model)」。** 它能精確指出每條 Path 的 WL/BL 在哪個 Tile 之後會處於絕對的 Deactive 狀態，從而計算極度精確的漏電流 (Leakage Power) 與動態功耗，完全排除了舊版 LCA 造成的功耗高估。

### 3. `res_dynamic$stats` (校正後的 WL/BL 狀態統計)
* **內容：** 包含 `cell_oo`, `cell_ox`, `cell_xo`, `cell_xx` 的數值。
* **作用：** **用於評估「Mask 污染程度與最佳化效果」。** 這些數據現在是基於物理 Mismatch 矩陣動態結算出來的，能 100% 真實反映 `radix_union` 等不同演算法在關閉不必要電路上的實際省電效益。
