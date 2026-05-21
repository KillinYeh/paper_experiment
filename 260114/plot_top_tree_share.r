calculate_tile_purity <- function(res_static) {
  # 檢查輸入
  if (is.null(res_static) || is.null(res_static$row_tiles)) {
    stop("Input error: res_static must contain 'row_tiles'.")
  }
  
  df <- res_static$row_tiles
  
  # 我們要對每一列 (每一個 Tile) 進行計算
  # 回傳一個詳細的 Data Frame
  
  purity_list <- lapply(seq_len(nrow(df)), function(i) {
    
    # 1. 取得該 Tile 內所有 Path 所屬的 Tree ID
    # tree_ids_in_tile 是一個 list-column，所以要用 [[i]] 取出向量
    t_ids <- unlist(df$tree_ids_in_tile[[i]])
    
    # 防呆：如果 Tile 是空的 (雖然不應該發生)
    if (length(t_ids) == 0) {
      return(data.frame(
        method = df$method[i],
        tile_id = df$tile_id[i],
        dominant_tree = NA,
        max_count = 0,
        total_count = 0,
        ratio = 0
      ))
    }
    
    # 2. 統計頻率 (Frequency Count)
    # table() 會計算每個 Tree ID 出現幾次
    counts <- table(t_ids)
    
    # 3. 找出最大值 (Most Frequent)
    max_count <- max(counts)
    
    # 找出是哪一棵樹 (如果有平手，取第一個)
    dom_tree <- names(counts)[which.max(counts)]
    
    # 4. 計算總數 (其實就是 length(t_ids) 或 rows_in_tile)
    total_count <- sum(counts)
    
    # 5. 計算比例
    ratio <- max_count / total_count
    
    data.frame(
      method = df$method[i],
      tile_id = df$tile_id[i],
      dominant_tree = as.integer(dom_tree),
      max_count = as.integer(max_count),
      total_count = as.integer(total_count),
      ratio = ratio
    )
  })
  
  # 合併成大表格
  purity_df <- do.call(rbind, purity_list)
  
  # --- 輸出統計報告 ---
  cat("========================================\n")
  cat("       Tile Purity Analysis Report      \n")
  cat("   (Ratio of Dominant Tree in Tile)     \n")
  cat("========================================\n")
  
  methods <- unique(purity_df$method)
  
  for (m in methods) {
    sub_df <- purity_df[purity_df$method == m, ]
    avg_ratio <- mean(sub_df$ratio)
    median_ratio <- median(sub_df$ratio)
    
    cat(sprintf("Method: %-10s\n", m))
    cat(sprintf("  - Mean Purity   : %.4f\n", avg_ratio))
    cat(sprintf("  - Median Purity : %.4f\n", median_ratio))
    cat(sprintf("  - Min / Max     : %.4f / %.4f\n", min(sub_df$ratio), max(sub_df$ratio)))
    cat("----------------------------------------\n")
  }
  
  return(purity_df)
}



# how to using ?
purity_stats <- calculate_tile_purity(res_static)

# 如果你想畫圖 (Boxplot 比較不同方法的純度分布)
library(ggplot2)
ggplot(purity_stats, aes(x = method, y = ratio, fill = method)) +
  geom_boxplot() +
  labs(title = "Tile Purity Distribution",
       y = "Dominant Tree Ratio (Purity)",
       x = "Method") +
  theme_minimal()
