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


progress_bar <- function(current, total, len = 50) {
  percent <- current / total
  filled <- round(percent * len)
  bar <- paste(rep("=", filled), collapse = "")
  empty <- paste(rep(" ", len - filled), collapse = "")
  cat(sprintf("\r[%s] %.1f%%", paste0(bar, empty), percent * 100))
  if (current == total) cat("\n")
}

feature_frequency_order <- function(paths_list, n_features, one_based = TRUE) {
  freq <- integer(n_features)
  for (p in paths_list) {
    f <- if (one_based) p else (p + 1L)
    if (length(f)) freq[f] <- freq[f] + 1L
  }
  order(freq, decreasing = TRUE)
}

remap_paths_to_rank <- function(paths_list, order_idx, one_based = TRUE) {
  F <- length(order_idx)
  inv <- integer(F)
  inv[order_idx] <- seq_len(F) 
  
  out <- lapply(paths_list, function(v) {
    f <- if (one_based) v else (v + 1L)
    if (!length(f)) return(integer(0L))
    inv[f]
  })
  attr(out, "n_features") <- F
  out
}

show_resources <- function(res_static) {
  # 檢查輸入是否合法
  if (is.null(res_static) || is.null(res_static$physical_tiles) || is.null(res_static$paths_ranked)) {
    stop("Input error: The input must be the result list (res_static) containing 'physical_tiles' and 'paths_ranked'.")
  }

  # 1. 計算總 Path 數
  # paths_ranked 是一個 list，每個 element 代表一條從 Root 到 Leaf 的路徑
  total_paths <- length(res_static$paths_ranked)

  # 2. 計算總 Tile 數
  # physical_tiles 是一個 data.frame
  # 它的每一列 (row) 都代表一個實體的硬體 Tile (包含 row grouping 和 column splitting 的結果)
  total_physical_tiles <- nrow(res_static$physical_tiles)

  # --- 額外資訊 (Optional) ---
  # 邏輯陣列數 (Logical Arrays): 也就是縱向切了幾刀 (Row Grouping)
  # 這代表我們用了幾個 Array ID
  n_logical_arrays <- length(unique(res_static$physical_tiles$array_id))

  # 橫向切分最大值 (Column Splits): 也就是特徵太多時，橫向切了幾刀
  n_col_splits <- max(res_static$physical_tiles$tile_in_array)

  # --- 輸出報告 ---
  cat("========================================\n")
  cat(sprintf("1. Total Paths (Rules)     : %d\n", total_paths))
  cat(sprintf("2. Total Physical Tiles    : %d\n", total_physical_tiles))
  cat("----------------------------------------\n")
  cat(sprintf("   - Logical Row Groups    : %d (Arrays)\n", n_logical_arrays))
  cat(sprintf("   - Max Column Splits     : %d (Tiles wide)\n", n_col_splits))
  cat("========================================\n")

  # 回傳一個 data frame 方便後續畫圖或紀錄
  return(data.frame(
    total_paths = total_paths,
    total_tiles = total_physical_tiles,
    logical_arrays = n_logical_arrays,
    max_col_splits = n_col_splits
  ))
}

############################################################
## 1. 靜態 WL/BL 統計
############################################################

wlbl_stats_for_window <- function(paths_ranked, rows_idx, start, end) {
  n_rows <- length(rows_idx) 
  win_len <- end - start + 1L
  if (win_len <= 0L || n_rows == 0L) {
    return(data.frame(wl_open=0L, wl_closed=0L, bl_open=0L, bl_closed=0L, 
                      cell_oo=0L, cell_ox=0L, cell_xo=0L, cell_xx=0L))
  }
  
  wl_open_vec <- logical(n_rows)
  bl_used <- rep(FALSE, win_len)
  
  for (i in seq_len(n_rows)) {
    r <- rows_idx[i]	#rows_idx[i] = the ith path in a tile
    v <- paths_ranked[[r]]
    if (!length(v)) next
    idx <- v[v >= start & v <= end]
    if (length(idx)) {
      wl_open_vec[i] <- TRUE
      bl_used[idx - start + 1L] <- TRUE
    }
  }
  
  wl_open   <- sum(wl_open_vec)
  wl_closed <- n_rows - wl_open
  bl_open   <- sum(bl_used)
  bl_closed <- win_len - bl_open
  
  data.frame(
    wl_open = wl_open, wl_closed = wl_closed,
    bl_open = bl_open, bl_closed = bl_closed,
    cell_oo = wl_open * bl_open,   cell_ox = wl_open * bl_closed,
    cell_xo = wl_closed * bl_open, cell_xx = wl_closed * bl_closed
  )
}

build_physical_tiles <- function(paths_ranked, row_tiles, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L || nrow(row_tiles) == 0L) return(NULL)
  
  res <- list()
  idx <- 0L
  
  for (i in seq_len(nrow(row_tiles))) {
    rows_idx <- row_tiles$rows_idx[[i]]
    array_id <- row_tiles$tile_id[i] # 在這裡array指的是完整的path , tile_in_array指的是在一個array中以tile_cols為單位切割出來的真實的Tile
    method   <- row_tiles$method[i]
    
    tile_in_array <- 0L
    for (start in seq.int(1L, n_features, by = tile_cols)) {
      end <- min(start + tile_cols - 1L, n_features)
      tile_in_array <- tile_in_array + 1L
      win_stats <- wlbl_stats_for_window(paths_ranked, rows_idx, start, end)
      
      idx <- idx + 1L
      res[[idx]] <- cbind(
        data.frame(method=method, array_id=array_id, tile_in_array=tile_in_array,
                   start_col=start, end_col=end, win_len=(end - start + 1L)),
        win_stats
      )
    }
  }
  do.call(rbind, res)
}

############################################################
## 2. Row Grouping (radix sort / LSH)
############################################################
pack_all_tiles_tries <- function(paths_ranked, tile_rows = 128L) {
  n_paths <- length(paths_ranked)
  n_features <- attr(paths_ranked, "n_features")
  
  # 記憶體與速度優化：只取最高頻的前 N 個特徵來做前綴排序
  # 由於 Tile 寬度就是 128，更後面的特徵對「前綴」的影響力微乎其微
  max_prefix_len <- min(n_features, 128L) 
  
  cat("    -> Building Bit Matrix for Radix Sort...\n")
  mat <- matrix(0L, nrow = n_paths, ncol = max_prefix_len)
  
  for (i in seq_len(n_paths)) {
    feats <- paths_ranked[[i]]
    # 只標記在前 max_prefix_len 範圍內的特徵
    valid_feats <- feats[feats <= max_prefix_len]
    if (length(valid_feats) > 0) {
      mat[i, valid_feats] <- 1L
    }
  }
  
  df_mat <- as.data.frame(mat)
  
  # 執行全域字典序排序 (Radix Sort)
  cat("    -> Executing fast Radix Sort...\n")
  # decreasing = TRUE 代表優先將有使用該特徵(1)的排在一起
  order_idx <- do.call(order, c(df_mat, list(decreasing = TRUE, method = "radix")))
  
  # 依序切分成 Tiles
  group_idx <- (seq_len(n_paths) - 1L) %/% tile_rows + 1L
  chunks <- split(order_idx, group_idx)
  
  do.call(rbind, lapply(names(chunks), function(nm) {
    rows <- chunks[[nm]]
    df <- data.frame(tile_id = as.integer(nm), rows_in_tile = length(rows))
    df$rows_idx <- I(list(rows))
    df
  }))
}

############################################################
## 2.2 Row Grouping (LSH / MinHash)
############################################################
pack_all_tiles_lsh <- function(paths_ranked, tile_rows = 128L, num_hash = 15L) {
  n_paths <- length(paths_ranked)
  n_features <- attr(paths_ranked, "n_features")
  # 自動適應邏輯：特徵越多 (越稀疏)，條件要越寬鬆 (num_hash 越小)
  if (n_features > 1000) {
    num_hash <- 1L      # Gene, Email
  } else if (n_features > 50) {
    num_hash <- 2L      # Gas
  } else {
    num_hash <- 4L      # Gesture, Fetal, Loan, Mush
  }
  
  cat("    -> Generating MinHash Signatures...\n")
  set.seed(42) # 固定 seed 以利實驗重現
  
  # 產生 num_hash 組隨機排列的特徵 ID (Hash Functions)
  hash_funcs <- replicate(num_hash, sample.int(n_features))
  
  signatures <- character(n_paths)
  for (i in seq_len(n_paths)) {
    feats <- paths_ranked[[i]]
    if (length(feats) == 0) {
      signatures[i] <- "empty"
    } else {
      # MinHash 核心邏輯：找出當前路徑的特徵中，在每個隨機排列裡最靠前的 Index
      sig <- apply(hash_funcs[feats, , drop = FALSE], 2, min)
      signatures[i] <- paste(sig, collapse = "_")
    }
  }
  
  cat("    -> Grouping by LSH Buckets...\n")
  # 根據 Hash Bucket (Signature) 進行全域排序，讓相似的路徑排在一起
  order_idx <- order(signatures)
  
  # 依序切分成 Tiles
  group_idx <- (seq_len(n_paths) - 1L) %/% tile_rows + 1L
  chunks <- split(order_idx, group_idx)
  
  do.call(rbind, lapply(names(chunks), function(nm) {
    rows <- chunks[[nm]]
    df <- data.frame(tile_id = as.integer(nm), rows_in_tile = length(rows))
    df$rows_idx <- I(list(rows))
    df
  }))
}

############################################################
## 2.3 Row Grouping (Radix Sort + Union Sliding Window)
############################################################

pack_all_tiles_radix_union <- function(paths_ranked, tile_rows = 128L, window_size = 512L) {
  n_paths <- length(paths_ranked)
  n_features <- attr(paths_ranked, "n_features")
  
  # --- Step 1: Radix Sort (全域前綴排序) ---
  # 為了速度與記憶體，只取最高頻的前 128 個特徵來做前綴排序
  max_prefix_len <- min(n_features, 128L) 
  
  cat("    -> Building Bit Matrix for Radix Sort...\n")
  mat <- matrix(0L, nrow = n_paths, ncol = max_prefix_len)
  
  for (i in seq_len(n_paths)) {
    feats <- paths_ranked[[i]]
    valid_feats <- feats[feats <= max_prefix_len]
    if (length(valid_feats) > 0) mat[i, valid_feats] <- 1L
  }
  
  df_mat <- as.data.frame(mat)
  
  cat("    -> Executing fast Radix Sort...\n")
  # 進行全域字典序排序，得到排列好的 index
  order_idx <- do.call(order, c(df_mat, list(decreasing = TRUE, method = "radix")))
  
  # --- Step 2: Sliding Window + Union Best-Fit ---
  cat(sprintf("    -> Packing tiles with Sliding Window (size=%d)...\n", window_size))
  
  unassigned_idx <- order_idx
  res_tiles <- list()
  tile_counter <- 1L
  
  # 使用 progress bar 追蹤進度
  total_tiles_estimate <- ceiling(n_paths / tile_rows)
  
  while (length(unassigned_idx) > 0) {
    # 1. 圈出當前的滑動視窗 (最多 window_size 個)
    current_window_size <- min(length(unassigned_idx), window_size)
    window_pool <- unassigned_idx[1:current_window_size]
    
    if (length(window_pool) <= tile_rows) {
      # 如果剩下的路徑不夠裝滿一個 Tile，就全部打包
      selected_for_tile <- window_pool
      unassigned_idx <- unassigned_idx[-(1:length(selected_for_tile))]
    } else {
      # 2. 挑選 Base Path (視窗內長度最短的 Path)
      path_lengths <- lengths(paths_ranked[window_pool])
      base_local_idx <- which.min(path_lengths)
      base_idx <- window_pool[base_local_idx]
      
      selected_for_tile <- integer(tile_rows)
      selected_for_tile[1] <- base_idx
      current_mask <- paths_ranked[[base_idx]]
      
      # 將 Base 移出候選池
      window_pool <- window_pool[-base_local_idx]
      
      # 3. 貪婪尋找剩下的 127 條路徑 (Best-Fit)
      for (k in 2:tile_rows) {
        # 計算視窗內每個 candidate 加入後，會新增多少個 1 (Cost)
        # sum(!(candidate %in% mask)) 就是增加的特徵數量
        costs <- vapply(window_pool, function(idx) {
          sum(!(paths_ranked[[idx]] %in% current_mask))
        }, numeric(1))
        
        # 找出 Cost 最小的 candidate (若平手會自動選視窗中最前面的，也就是前綴最像的)
        best_local_idx <- which.min(costs)
        best_idx <- window_pool[best_local_idx]
        
        # 加入 Tile 並更新 Mask
        selected_for_tile[k] <- best_idx
        current_mask <- unique(c(current_mask, paths_ranked[[best_idx]]))
        
        # 將選中的人移出當前候選池
        window_pool <- window_pool[-best_local_idx]
      }
      
      # 4. 將這 128 個被選中的人，從全域的未分配名單中移除
      # 注意：R 的 setdiff 會自動保留 unassigned_idx 原本的排序，完美維持 Radix Sort 的效果
      unassigned_idx <- setdiff(unassigned_idx, selected_for_tile)
    }
    
    # 紀錄這個 Tile 的資訊
    res_tiles[[tile_counter]] <- data.frame(
      tile_id = tile_counter, 
      rows_in_tile = length(selected_for_tile)
    )
    res_tiles[[tile_counter]]$rows_idx <- I(list(selected_for_tile))
    
    if (tile_counter %% 10 == 0) {
       progress_bar(min(tile_counter, total_tiles_estimate), total_tiles_estimate)
    }
    tile_counter <- tile_counter + 1L
  }
  progress_bar(total_tiles_estimate, total_tiles_estimate)
  
  do.call(rbind, res_tiles)
}





############################################################
## 3. 核心比較邏輯 
############################################################

compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE, path_tree_id, leaf_ids_list) {
  n_features <- as.integer(n_features)
  
  # 1. Feature Reordering (全域優化，必須保留)
  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features
  
  cat("Packing tiles (tries (radix sort) / LSH)...\n")
  progress_bar(2, 5)
  cat("[Method 2] Tries (Radix Sort)...\n")
  tries_df <- pack_all_tiles_tries(paths_ranked)
  tries_df$method <- "tries_radix"

  cat("[Method 3] LSH (MinHash)...\n")
  lsh_df <- pack_all_tiles_lsh(paths_ranked, num_hash = 15L)
  lsh_df$method <- "lsh_minhash"
  
  cat("\n[Method 2] Radix + Union (Window=512)...\n")
  radix_union_df <- pack_all_tiles_radix_union(paths_ranked, tile_rows = 128L, window_size = 512L)
  radix_union_df$method <- "radix_union"



  
  # 附加 Meta Info
  append_meta <- function(df) {
    df$tree_ids_in_tile <- lapply(df$rows_idx, function(idx) path_tree_id[idx])
    df$leaf_ids_in_tile <- lapply(df$rows_idx, function(idx) leaf_ids_list[idx])
    df
  }
  
  
  all_row_tiles <- rbind(
    append_meta(tries_df),
    append_meta(lsh_df),
    append_meta(radix_union_df)
  )

  # 在這裡之後，每一個df包含以下成員：
  # tile_id = array id = 第幾個array
  # rows_in_tile = 一個array有幾條path
  # rows_idx = 這 $rows_in_tile 條path完整的路徑
  # tree_ids_in_tile = 每一個path來自第幾棵tree
  # leaf_ids_in_tile = 每一個path的leaf node ID
 
  
  cat("build physical tile\n")
  progress_bar(3, 5)
  physical_tiles <- build_physical_tiles(paths_ranked, all_row_tiles, tile_cols = 128L)
  
  list(
    columns_order = ord,          # feature after sorted by feature frequency
    paths_ranked = paths_ranked,  # path reorder by column_order
    n_features = n_features,
    path_tree_id = path_tree_id,
    row_tiles = all_row_tiles,
    physical_tiles = physical_tiles
  )
}

############################################################
## 4. 提取 Forest 拓樸與 ACAM 區間 (瘦身版：移除 LCA)
############################################################
extract_forest_info_optimized <- function(fit, n_features) {
  all_paths_list   <- list()
  all_path_tree_id <- integer(0L)
  all_leaf_ids     <- integer(0L)
  leaf_info_list   <- list()
  
  # 儲存全域特徵的切割點
  all_split_vals_per_feat <- vector("list", n_features)
  
  cat("Extracting paths and ACAM cell intervals...\n")
  
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    
    paths_t <- list()
    leaf_ids_t <- integer(0L)
    path_intervals_t <- list()
    
    non_terms <- info[!info$terminal, ]
    if(nrow(non_terms) > 0) {
      for(i in 1:nrow(non_terms)) {
        f_id <- non_terms$splitvarID[i] + 1L 
        s_val <- non_terms$splitval[i]
        all_split_vals_per_feat[[f_id]] <- c(all_split_vals_per_feat[[f_id]], s_val)
      }
    }
    
    left_child <- info$leftChild + 1
    right_child <- info$rightChild + 1
    split_var <- info$splitvarID
    is_term <- info$terminal
    node_ids <- info$nodeID
    
    dfs <- function(idx, current_feats, lower_b, upper_b) {
      curr_node_id <- node_ids[idx]
      
      if (is_term[idx]) {
        # 修正：不需 + 1L，確保 current_feats 維持 1-based
        p_final <- if(length(current_feats)) unique(sort(current_feats)) else integer(0L)
        paths_t[[length(paths_t) + 1L]] <<- p_final
        leaf_ids_t <<- c(leaf_ids_t, curr_node_id)
        
        used_feats <- p_final
        if (length(used_feats) > 0) {
          intervals <- data.frame(
            feature_id = used_feats,
            min_val = lower_b[used_feats],
            max_val = upper_b[used_feats]
          )
          path_intervals_t[[as.character(curr_node_id)]] <<- intervals
        } else {
          path_intervals_t[[as.character(curr_node_id)]] <<- data.frame()
        }
      } else {
        feat <- split_var[idx] + 1L
        s_val <- info$splitval[idx]
        
        new_upper <- upper_b
        new_upper[feat] <- min(new_upper[feat], s_val)
        dfs(left_child[idx], c(current_feats, feat), lower_b, new_upper)
        
        new_lower <- lower_b
        new_lower[feat] <- max(new_lower[feat], s_val)
        dfs(right_child[idx], c(current_feats, feat), new_lower, upper_b)
      }
    }
    
    if (nrow(info) > 0) {
      dfs(1, integer(0L), rep(-Inf, n_features), rep(Inf, n_features))
    }
    
    all_paths_list <- c(all_paths_list, paths_t)
    all_path_tree_id <- c(all_path_tree_id, rep.int(t, length(paths_t)))
    all_leaf_ids <- c(all_leaf_ids, leaf_ids_t)
    
    leaf_info_list[[t]] <- list(
      leaf_nodes = leaf_ids_t,
      path_intervals = path_intervals_t
    )
  }
  
  feature_split_ranges <- lapply(all_split_vals_per_feat, function(x) {
    if(is.null(x)) return(numeric(0))
    sort(unique(x))
  })
  
  progress_bar(fit$num.trees, fit$num.trees)
  
  list(
    paths = all_paths_list,
    tree_ids = all_path_tree_id,
    leaf_ids = all_leaf_ids,
    leaf_info_list = leaf_info_list,
    feature_split_ranges = feature_split_ranges
  )
}

############################################################
## 4.2：建立全域 ACAM State 映射表 (用於計算 Program 功耗)
############################################################
generate_acam_state_table <- function(extracted_info) {
  cat("Generating Global ACAM State Table...\n")
  state_list <- list()
  idx <- 1L
  
  for (t in seq_along(extracted_info$leaf_info_list)) {
    info <- extracted_info$leaf_info_list[[t]]
    intervals_list <- info$path_intervals
    
    for (leaf_id_str in names(intervals_list)) {
      df <- intervals_list[[leaf_id_str]]
      if (nrow(df) > 0) {
        
        # 計算一個path中所有 internal node 落在哪個 State (利用 findInterval 與全域 feature_split_ranges 比較)
        # min_val == -Inf 代表無下限 (State 0)
        min_states <- mapply(function(val, f_id) {
          if (val == -Inf) 0L else findInterval(val, extracted_info$feature_split_ranges[[f_id]])
        }, df$min_val, df$feature_id)
        
        # max_val == Inf 代表無上限 (最大 State)
        max_states <- mapply(function(val, f_id) {
          if (val == Inf) length(extracted_info$feature_split_ranges[[f_id]]) else findInterval(val, extracted_info$feature_split_ranges[[f_id]])
        }, df$max_val, df$feature_id)
        
        df$min_state <- as.integer(min_states)
        df$max_state <- as.integer(max_states)
        df$tree_id <- as.integer(t)
        df$leaf_id <- as.integer(leaf_id_str)
        
        state_list[[idx]] <- df
        idx <- idx + 1L
      }
    }
  }
  
  final_df <- do.call(rbind, state_list)
  # 重排欄位順序使其清晰易讀
  final_df <- final_df[, c("tree_id", "leaf_id", "feature_id", "min_val", "max_val", "min_state", "max_state")]
  
  return(final_df)
}


precompute_bl_info <- function(physical_tiles) {
  physical_tiles$key <- paste(physical_tiles$method, physical_tiles$array_id, physical_tiles$tile_in_array, sep="_")
  bl_vec <- physical_tiles[, c("bl_open", "bl_closed", "win_len")]
  row.names(bl_vec) <- physical_tiles$key
  bl_vec
}

############################################################
## 4. 提取 Forest 拓樸與 ACAM 區間 (瘦身版：移除 LCA)
############################################################
extract_forest_info_optimized <- function(fit, n_features) {
  all_paths_list   <- list()
  all_path_tree_id <- integer(0L)
  all_leaf_ids     <- integer(0L)
  leaf_info_list   <- list()
  
  # 儲存全域特徵的切割點
  all_split_vals_per_feat <- vector("list", n_features)
  
  cat("Extracting paths and ACAM cell intervals...\n")
  
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    
    paths_t <- list()
    leaf_ids_t <- integer(0L)
    path_intervals_t <- list()
    
    non_terms <- info[!info$terminal, ]
    if(nrow(non_terms) > 0) {
      for(i in 1:nrow(non_terms)) {
        f_id <- non_terms$splitvarID[i] + 1L 
        s_val <- non_terms$splitval[i]
        all_split_vals_per_feat[[f_id]] <- c(all_split_vals_per_feat[[f_id]], s_val)
      }
    }
    
    left_child <- info$leftChild + 1
    right_child <- info$rightChild + 1
    split_var <- info$splitvarID
    is_term <- info$terminal
    node_ids <- info$nodeID
    
    dfs <- function(idx, current_feats, lower_b, upper_b) {
      curr_node_id <- node_ids[idx]
      
      if (is_term[idx]) {
        # 修正：不需 + 1L，確保 current_feats 維持 1-based
        p_final <- if(length(current_feats)) unique(sort(current_feats)) else integer(0L)
        paths_t[[length(paths_t) + 1L]] <<- p_final
        leaf_ids_t <<- c(leaf_ids_t, curr_node_id)
        
        used_feats <- p_final
        if (length(used_feats) > 0) {
          intervals <- data.frame(
            feature_id = used_feats,
            min_val = lower_b[used_feats],
            max_val = upper_b[used_feats]
          )
          path_intervals_t[[as.character(curr_node_id)]] <<- intervals
        } else {
          path_intervals_t[[as.character(curr_node_id)]] <<- data.frame()
        }
      } else {
        feat <- split_var[idx] + 1L
        s_val <- info$splitval[idx]
        
        new_upper <- upper_b
        new_upper[feat] <- min(new_upper[feat], s_val)
        dfs(left_child[idx], c(current_feats, feat), lower_b, new_upper)
        
        new_lower <- lower_b
        new_lower[feat] <- max(new_lower[feat], s_val)
        dfs(right_child[idx], c(current_feats, feat), new_lower, upper_b)
      }
    }
    
    if (nrow(info) > 0) {
      dfs(1, integer(0L), rep(-Inf, n_features), rep(Inf, n_features))
    }
    
    all_paths_list <- c(all_paths_list, paths_t)
    all_path_tree_id <- c(all_path_tree_id, rep.int(t, length(paths_t)))
    all_leaf_ids <- c(all_leaf_ids, leaf_ids_t)
    
    leaf_info_list[[t]] <- list(
      leaf_nodes = leaf_ids_t,
      path_intervals = path_intervals_t
    )
  }
  
  feature_split_ranges <- lapply(all_split_vals_per_feat, function(x) {
    if(is.null(x)) return(numeric(0))
    sort(unique(x))
  })
  
  progress_bar(fit$num.trees, fit$num.trees)
  
  list(
    paths = all_paths_list,
    tree_ids = all_path_tree_id,
    leaf_ids = all_leaf_ids,
    leaf_info_list = leaf_info_list,
    feature_split_ranges = feature_split_ranges
  )
}

############################################################
## 4.1：將 Global ACAM State 映射到實際的硬體排列 (依據 Packing Method)
############################################################
generate_hardware_mapped_state_table <- function(row_tiles, global_state_table) {
  cat("Mapping ACAM states to physical hardware layouts (by method)...\n")
  
  # 1. 展開 row_tiles：把每一種 method 的 array 佈局攤平
  # 建立一個基礎的「硬體位置對應表 (Hardware Position Map)」
  hw_map_list <- lapply(seq_len(nrow(row_tiles)), function(i) {
    method <- row_tiles$method[i]
    array_id <- row_tiles$tile_id[i]
    
    tree_ids <- unlist(row_tiles$tree_ids_in_tile[[i]])
    leaf_ids <- unlist(row_tiles$leaf_ids_in_tile[[i]])
    n_rows <- length(tree_ids)
    
    if (n_rows == 0) return(NULL)
    
    data.frame(
      method = method,
      array_id = array_id,
      row_in_array = seq_len(n_rows), # 1 ~ 128 (Tile 內的實體 Row Index)
      tree_id = as.integer(tree_ids),
      leaf_id = as.integer(leaf_ids),
      stringsAsFactors = FALSE
    )
  })
  
  hw_map_df <- do.call(rbind, hw_map_list)
  
  # 2. 進行 Left Join (合併)：
  # 根據 tree_id 和 leaf_id，把 global_state_table 裡面的 feature 和 state 貼到硬體位置上
  hw_state_table <- merge(hw_map_df, global_state_table, by = c("tree_id", "leaf_id"), all.x = TRUE)
  
  # 3. 重新排序：讓表格看起來就像真實硬體中的排列順序
  # 順序：方法 -> 陣列 ID -> 陣列內的 Row -> Feature ID
  hw_state_table <- hw_state_table[order(hw_state_table$method, hw_state_table$array_id, hw_state_table$row_in_array, hw_state_table$feature_id), ]
  
  # 4. 重排欄位，讓閱讀更直覺
  final_cols <- c("method", "array_id", "row_in_array", "tree_id", "leaf_id", "feature_id", "min_val", "max_val", "min_state", "max_state")
  hw_state_table <- hw_state_table[, final_cols]
  
  # 去除可能因為某些 Path 沒有用到任何 Feature 而產生的 NA row
  hw_state_table <- hw_state_table[!is.na(hw_state_table$feature_id), ]

  rownames(hw_state_table) <- NULL # 按照所有方法產生的 path 順序重新做編號
  
  return(hw_state_table)
}

############################################################
## 5. 整合與執行
############################################################

simulate_hardware_inference <- function(data_test, res_static, tile_cols = 128L) {
  
  hw_state_table <- res_static$hardware_state_tables
  row_tiles <- res_static$row_tiles
  
  # 建立 Feature ID 轉換至實體排列 Rank 的 Lookup Table
  ord <- res_static$columns_order
  feat_rank_lookup <- integer(res_static$n_features)
  feat_rank_lookup[ord] <- seq_along(ord)
  
  # 預算所有實體 Tile 的 BL 開關資訊
  bl_lookup <- precompute_bl_info(res_static$physical_tiles)
  
  mat_test <- as.matrix(data_test)
  n_samples <- nrow(mat_test)
  methods <- unique(hw_state_table$method)
  out <- list()
  
  for (method in methods) {
    cat("Simulating hardware inference for method:", method, "\n")
    
    sub_hw <- hw_state_table[hw_state_table$method == method, ]
    total_paths_method <- sum(row_tiles$rows_in_tile[row_tiles$method == method])
    
    total_oo <- 0; total_ox <- 0; total_xo <- 0; total_xx <- 0
    mat_mis_tile_all <- matrix(NA_integer_, nrow=n_samples, ncol=total_paths_method)
    col_offset <- 0L
    
    # 按照 Array_ID 將硬體排佈表切塊
    array_list <- split(sub_hw, sub_hw$array_id, drop = TRUE)
    
    for (a_id_str in names(array_list)) {
      a_id <- as.integer(a_id_str)
      arr_chunk <- array_list[[a_id_str]]
      
      # 取出該 Array 在 physical_tiles 中的 meta info 以計算迴圈長度
      pt_sub <- res_static$physical_tiles[res_static$physical_tiles$method == method & 
                                            res_static$physical_tiles$array_id == a_id, ]
      n_tiles_wide <- nrow(pt_sub)
      n_rows_in_array <- max(arr_chunk$row_in_array)
      
      mat_mis_tile <- matrix(NA_integer_, nrow=n_samples, ncol=n_rows_in_array)
      
      # 按照 Row 繼續切塊
      row_list <- split(arr_chunk, arr_chunk$row_in_array, drop = TRUE)
      
      for (r_id_str in names(row_list)) {
        r_id <- as.integer(r_id_str)
        df_path <- row_list[[r_id_str]]
        
        # 【物理正確性關鍵】：必須確保 df_path 裡面的 feature 是依照實體陣列從左到右 (Rank) 排序的
        df_path <- df_path[order(feat_rank_lookup[df_path$feature_id]), ]
        
        feats <- df_path$feature_id
        test_vals <- mat_test[, feats, drop = FALSE]
        
        # 向量化比對區間
        pass_mat <- sweep(test_vals, 2, df_path$min_val, ">") & 
                    sweep(test_vals, 2, df_path$max_val, "<=")
        
        is_match <- rowSums(!pass_mat) == 0
        
        # 找出物理上最左邊 (第一個) 產生 Mismatch 的 Feature
        first_miss_idx <- max.col(!pass_mat, ties.method = "first")
        miss_f_id <- feats[first_miss_idx]
        
        # 計算該 Mismatch 發生在第幾個橫向 Tile
        tiles_idx <- (feat_rank_lookup[miss_f_id] - 1L) %/% tile_cols + 1L
        
        # 將 Mismatch 的 Tile_ID 寫入矩陣中 (Match 者保持 NA)
        mat_mis_tile[!is_match, r_id] <- tiles_idx[!is_match]
      }
      
      # 更新至全域追蹤矩陣
      col_range <- (col_offset + 1L):(col_offset + n_rows_in_array)
      mat_mis_tile_all[, col_range] <- mat_mis_tile
      col_offset <- col_offset + n_rows_in_array
      
      # =========================================================
      # 動態功耗 WL/BL 計算邏輯 (完美無縫接軌)
      # =========================================================
      for (k in seq_len(n_tiles_wide)) {
        key <- paste(method, a_id, k, sep="_")
        if (!key %in% rownames(bl_lookup)) next
        
        bl_info <- bl_lookup[key, ]
        b_open  <- bl_info[["bl_open"]]
        b_closed <- bl_info[["bl_closed"]]
        
        # 如果是 NA (Match) 或 Mismatch 發生在第 k 個 Tile 或之後，代表在當前 k Tile 仍須保持 Active
        wl_is_active_mat <- is.na(mat_mis_tile) | (mat_mis_tile >= k)
        wl_open_counts <- rowSums(wl_is_active_mat)
        wl_closed_counts <- n_rows_in_array - wl_open_counts
        
        sum_wl_open <- sum(wl_open_counts)
        sum_wl_closed <- sum(wl_closed_counts)
        
        total_oo <- total_oo + sum_wl_open * b_open
        total_ox <- total_ox + sum_wl_open * b_closed
        total_xo <- total_xo + sum_wl_closed * b_open
        total_xx <- total_xx + sum_wl_closed * b_closed
      }
    }
    
    out[[method]] <- list(
      stats = c(cell_oo=total_oo, cell_ox=total_ox, cell_xo=total_xo, cell_xx=total_xx),
      mismatch_matrix = mat_mis_tile_all # 用於未來深度 Debug 查閱
    )
  }
  return(out)
}

main_simulation <- function(fit,  data_test) {
  library(ranger)
  
  feat_names <- fit$forest$independent.variable.names
  n_features <- length(feat_names)
  
  # 1. 提取 (瘦身版 DFS)
  extracted <- extract_forest_info_optimized(fit, n_features)
  
  # 2. 靜態 Mapping 與物理切割
  res_static <- compare_proposed_vs_naive(
    paths_list = extracted$paths,
    n_features = n_features,
    one_based = TRUE,
    path_tree_id = extracted$tree_ids,
    leaf_ids_list = extracted$leaf_ids
  )
  show_resources(res_static)
  
  # 3. 建立 Global State Table 與 Hardware State Table
  res_static$acam_state_table <- generate_acam_state_table(extracted)
  res_static$hardware_state_tables <- generate_hardware_mapped_state_table(
    row_tiles = res_static$row_tiles,
    global_state_table = res_static$acam_state_table
  )
  
  purity_stats <- calculate_tile_purity(res_static)
  
  # 4. 極速實體動態模擬 (取代原本 LCA 的 simulate_cam_array_mismatch)
  cat("Running physical dynamic simulation...\n")
  res_dynamic <- simulate_hardware_inference(
    data_test = data_test,
    res_static = res_static,
    tile_cols = 128L
  )
  
  cat("Finish the simulation!\n")
  progress_bar(5, 5)
  list(static = res_static, dynamic = res_dynamic)
}

#---------------------------- experiment -------------------------
sim_email <- main_simulation(email_forest100,email_test_x)
  
sim_fetal <- main_simulation(fetal_forest100,fetal_x)
  
sim_gene  <- main_simulation(gene_forest100,gene_test_x)
  
sim_mush  <- main_simulation(mush_forest100,mush_test_x)
  
sim_loan  <- main_simulation(loan_forest100,loan_test_x)
  
sim_gesture <- main_simulation(gesture_forest100,gesture_test_x)
   
sim_arcene <- main_simulation(arcene_forest100,arcene_test_x)
  
sim_gas <- main_simulation(gas_forest100,gas_test_x)

save.image(file = "backup_testing_dataset_performance_0412.RData")



#sim_cov <- main_simulation(cov_forest_100,cov_test_x)
