# =====================================================================
# measure_packing.R (獨立測量腳本：每次只跑單一 Dataset 與 單一 Method)
# =====================================================================

# 載入所有打包函數
source("function.R")

# 1. 接收終端機傳來的參數: [Dataset名稱] [Method名稱]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript measure_packing.R <dataset_name> <method>\nExample: Rscript measure_packing.R email radix")
}

dataset_name <- args[1]
method       <- args[2]

# 2. 組合檔名並讀取前處理好的資料
file_name <- paste0("paths_ranked_", dataset_name, ".rds")
if (!file.exists(file_name)) {
  stop(sprintf("File '%s' not found! Please run prep_data.R first.\n", file_name))
}

loaded_data <- readRDS(file_name)
paths_ranked <- loaded_data$paths
n_features <- loaded_data$n_features
n_paths <- length(paths_ranked)
tile_rows <- 128L

# 3. 根據參數執行對應的演算法
cat(sprintf("\n[Task] Dataset: %s (N=%d) | Method: %s\n", dataset_name, n_paths, method))


# --- 開始執行打包 ---
if (method == "seq") {
  res <- pack_all_tiles_sequential(paths_ranked, tile_rows)
} else if (method == "naive") {
  res <- pack_all_tiles_naive(paths_ranked, tile_rows)
} else if (method == "mcp") {
  res <- pack_all_tiles_mcp(paths_ranked, tile_rows)
} else if (method == "radix") {
  res <- pack_all_tiles_radix_union(paths_ranked, tile_rows, window_size = 512L)
} else if (method == "lsh") {
  res <- pack_all_tiles_lsh(paths_ranked, tile_rows)
} else if (method == "fptree") {
  res <- pack_all_tiles_fptree_union(paths_ranked, tile_rows, window_size = 512L)
} else if (method == "union") {
  res <- pack_all_tiles_union(paths_ranked, tile_rows)
} else {
  stop("Unknown method. Please specify seq, naive, mcp, radix, lsh, fptree, or union.\n")
}

cat("Packing finished successfully.\n")
