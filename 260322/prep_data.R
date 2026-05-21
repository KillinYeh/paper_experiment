# =====================================================================
# prep_data.R (復現版本：讀取模型檔 -> 萃取 -> 排序 -> 存檔)
# =====================================================================
source("function.R")

# 定義要處理的 Dataset 名稱
dataset_names <- c("email", "fetal", "gene", "mush", "loan", "gesture", "arcene", "gas")

cat("==================================================\n")
cat("Starting Reproducible Path-Extraction Process...\n")
cat("==================================================\n\n")

for (ds_name in dataset_names) {
  cat(sprintf("Processing Dataset: [%s] ...\n", ds_name))
  
  # 1. 讀取預先訓練好的模型檔案
  model_file <- sprintf("fit_%s.rds", ds_name)
  if (!file.exists(model_file)) {
    cat(sprintf("   ❌ Skip: Model file '%s' not found.\n", model_file))
    next
  }
  
  cat(sprintf("   -> Loading model from '%s'...\n", model_file))
  current_fit <- readRDS(model_file)
  
  # 2. 獲取參數並萃取路徑
  n_features <- current_fit$num.independent.variables
  cat("   -> Extracting paths from forest...\n")
  extracted <- extract_forest_info_optimized(current_fit, n_features)
  paths_list <- extracted$paths
  
  # 3. 特徵頻率排序與 Rank 映射
  cat("   -> Reordering features by frequency...\n")
  ord <- feature_frequency_order(paths_list, n_features, one_based = TRUE)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = TRUE)
  attr(paths_ranked, "n_features") <- n_features
  
  # 4. 存檔提供給 measure_packing.R 使用
  var_name <- paste0("paths_ranked_", ds_name)
  save_data <- list(
    paths = paths_ranked,
    n_features = n_features
  )
  
  file_name <- paste0(var_name, ".rds")
  saveRDS(save_data, file = file_name)
  cat(sprintf("   ✅ Processed and saved to '%s'\n\n", file_name))
}

cat("==================================================\n")
cat("Reproducible Prep Complete!\n")
cat("==================================================\n")
