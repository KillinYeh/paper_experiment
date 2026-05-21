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
    array_id <- row_tiles$tile_id[i]
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
  #cat("[Method 2] Tries (Radix Sort)...\n")
  #tries_df <- pack_all_tiles_tries(paths_ranked)
  #tries_df$method <- "tries_radix"

  #cat("[Method 3] LSH (MinHash)...\n")
  #lsh_df <- pack_all_tiles_lsh(paths_ranked, num_hash = 15L)
  #lsh_df$method <- "lsh_minhash"
  
  cat("\n[Method 2] Radix + Union (Window=512)...\n")
  radix_union_df <- pack_all_tiles_radix_union(paths_ranked, tile_rows = 128L, window_size = 512L)
  radix_union_df$method <- "radix_union"



  
  # 附加 Meta Info
  append_meta <- function(df) {
    df$tree_ids_in_tile <- lapply(df$rows_idx, function(idx) path_tree_id[idx])
    df$leaf_ids_in_tile <- lapply(df$rows_idx, function(idx) leaf_ids_list[idx])
    df
  }
  
  # 這裡只會有 sequential 的結果
  all_row_tiles <- rbind(
  #  append_meta(tries_df),
  #  append_meta(lsh_df),
     append_meta(radix_union_df)
  )
  
  cat("build physical tile\n")
  progress_bar(3, 5)
  physical_tiles <- build_physical_tiles(paths_ranked, all_row_tiles, tile_cols = 128L)
  
  list(
    columns_order = ord,          #feature after sorted by feature frequency
    paths_ranked = paths_ranked,  #path reorder by column_order
    n_features = n_features,
    path_tree_id = path_tree_id,
    row_tiles = all_row_tiles,
    physical_tiles = physical_tiles
  )
}

############################################################
## 4. LCA Map 與 Mismatch 模擬 (優化 DFS)
############################################################


build_tree_lca_maps_optimized <- function(fit, leaf_info_list) {
  lca_maps <- list()

  for (t in seq_along(leaf_info_list)) {
    info <- leaf_info_list[[t]]
    leaves <- info$leaf_nodes
    leaf_paths <- info$node_paths

    n_leaves <- length(leaves)
    mat <- matrix(NA_integer_, nrow=n_leaves, ncol=n_leaves)
    rownames(mat) <- colnames(mat) <- as.character(leaves)

    for (i in 1:n_leaves) {
      lid_i <- as.character(leaves[i])
      path_i <- leaf_paths[[lid_i]]
      mat[i, i] <- -1L

      if (i < n_leaves) {
        for (j in (i+1):n_leaves) {
          lid_j <- as.character(leaves[j])
          path_j <- leaf_paths[[lid_j]]

          len <- min(length(path_i), length(path_j))
          same <- (path_i[1:len] == path_j[1:len])
          last_common_idx <- if(all(same)) len else which.min(same) - 1L

          if (last_common_idx > 0) {
            # 【修改】：直接將 LCA 的 node_id 存入矩陣中
            lca_node <- path_i[last_common_idx]
            mat[i, j] <- lca_node
            mat[j, i] <- lca_node
          }
        }
      }
    }
    lca_maps[[t]] <- mat
  }
  lca_maps
}

precompute_bl_info <- function(physical_tiles) {
  physical_tiles$key <- paste(physical_tiles$method, physical_tiles$array_id, physical_tiles$tile_in_array, sep="_")
  bl_vec <- physical_tiles[, c("bl_open", "bl_closed", "win_len")]
  row.names(bl_vec) <- physical_tiles$key
  bl_vec
}


simulate_cam_array_mismatch <- function(fit, data_test, res, leaf_info_list, feature_split_ranges, tile_cols = 128L) {
  row_tiles    <- res$row_tiles
  paths_ranked <- res$paths_ranked
  ord <- res$columns_order
  feat_rank_lookup <- integer(res$n_features)
  feat_rank_lookup[ord] <- seq_along(ord)

  cat("Predicting terminal nodes for all samples...\n")
  pred_leaf <- predict(fit, data = data_test, type = "terminalNodes")$predictions
  n_samples <- nrow(pred_leaf)
  
  # 【新增】：將 data_test 轉為 matrix 以利後續極速提取座標值
  mat_data_test <- as.matrix(data_test)

  cat("Building LCA Maps...\n")
  progress_bar(4, 5)
  lca_maps <- build_tree_lca_maps_optimized(fit, leaf_info_list)
  bl_lookup <- precompute_bl_info(res$physical_tiles)
  methods <- unique(row_tiles$method)
  out <- list()

  for (method in methods) {
    cat("Simulating method:", method, "\n")
    sub <- row_tiles[row_tiles$method == method, ]

    total_oo <- 0; total_ox <- 0; total_xo <- 0; total_xx <- 0

    for (i in seq_len(nrow(sub))) {
      rows_idx <- sub$rows_idx[[i]]
      tree_ids <- unlist(sub$tree_ids_in_tile[[i]])
      leaf_ids <- unlist(sub$leaf_ids_in_tile[[i]])
      array_id <- sub$tile_id[i]
      n_rows_in_array <- length(rows_idx)

      pt_sub <- res$physical_tiles[res$physical_tiles$method == method &
                                     res$physical_tiles$array_id == array_id, ]
      n_tiles_wide <- nrow(pt_sub)
      mat_mis_tile <- matrix(Inf, nrow=n_samples, ncol=n_rows_in_array)

      for (r in seq_len(n_rows_in_array)) {
        t_id <- tree_ids[r]
        l_tgt <- as.character(leaf_ids[r])
        l_true_vec <- as.character(pred_leaf[, t_id])

        diff_indices <- which(l_true_vec != l_tgt)
        if (length(diff_indices) > 0) {
          mismatch_true_leaves <- l_true_vec[diff_indices]
          map <- lca_maps[[t_id]]
          
          # 這裡現在取得的是 lca_node_id
          lca_nodes <- map[cbind(mismatch_true_leaves, l_tgt)] 
          lca_nodes_str <- as.character(lca_nodes)
          
          # 查表取得 feat_id
          mis_feats <- vapply(lca_nodes_str, function(n) leaf_info_list[[t_id]]$node_split_map[[n]]$feat_id, numeric(1))
          
          # =========================================================
          # 【關鍵新增：提取 Threshold 與 測試資料落點 State】
          # 如果您需要將這些資訊存下來給未來的 Quantile Simulation 用，
          # 您可以在這裡建立一個 log data.frame 把它存起來。
          # =========================================================
          # 1. 取得這顆節點真實的 Split Threshold
          mis_thresh <- vapply(lca_nodes_str, function(n) leaf_info_list[[t_id]]$node_split_map[[n]]$splitval, numeric(1))
          
          # 2. 取得這筆測試資料在該 Mismatch 特徵上的「真實數值」
          actual_test_vals <- mat_data_test[cbind(diff_indices, mis_feats)]
          
          # 3. 算出這個值落在該 Feature 總體 Range 的第幾個區間 (State/Bin)
          # 使用 findInterval 瞬間找出它屬於哪個 State Index
          test_state_bins <- mapply(function(val, f_id) {
            findInterval(val, feature_split_ranges[[f_id]])
          }, actual_test_vals, mis_feats)
          
          # 您可以將 mis_thresh, actual_test_vals, test_state_bins 匯出
          # 例如存進全域的 logger list 中...
          # =========================================================

          ranks <- feat_rank_lookup[mis_feats]
          tiles_idx <- (ranks - 1L) %/% tile_cols + 1L
          mat_mis_tile[diff_indices, r] <- tiles_idx
        }
      }

      # ... (下方 BL/WL 邏輯不變，保持原狀) ...
      for (k in seq_len(n_tiles_wide)) {
        key <- paste(method, array_id, k, sep="_")
        bl_info <- bl_lookup[key, ]
        b_open  <- bl_info[["bl_open"]]
        b_closed <- bl_info[["bl_closed"]]

        wl_is_active_mat <- (mat_mis_tile >= k)
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
    out[[method]] <- c(cell_oo=total_oo, cell_ox=total_ox, cell_xo=total_xo, cell_xx=total_xx)
  }
  out
}

############################################################
## 5. 整合與執行
############################################################

extract_forest_info_optimized <- function(fit, n_features) {
  all_paths_list   <- list()
  all_path_tree_id <- integer(0L)
  all_leaf_ids     <- integer(0L)
  leaf_info_list   <- list()

  # 用來儲存整個森林中，每個特徵所有出現過的 split_val
  all_split_vals_per_feat <- vector("list", n_features)

  cat("Extracting paths, topology, and ACAM cell intervals...\n")

  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)

    paths_t <- list()
    leaf_ids_t <- integer(0L)
    node_paths_t <- list()
    node_split_map <- list()
    path_intervals_t <- list() # 【新增】：儲存每條 Path 的真實區間

    non_terms <- info[!info$terminal, ]
    if(nrow(non_terms) > 0) {
      for(i in 1:nrow(non_terms)) {
        nid_str <- as.character(non_terms$nodeID[i])
        f_id <- non_terms$splitvarID[i] + 1L # 轉為 1-based
        s_val <- non_terms$splitval[i]

        node_split_map[[nid_str]] <- list(feat_id = f_id, splitval = s_val)
        all_split_vals_per_feat[[f_id]] <- c(all_split_vals_per_feat[[f_id]], s_val)
      }
    }

    left_child <- info$leftChild + 1
    right_child <- info$rightChild + 1
    split_var <- info$splitvarID
    is_term <- info$terminal
    node_ids <- info$nodeID

    # 【新增】：DFS 攜帶 lower_b 與 upper_b 陣列
    dfs <- function(idx, current_feats, current_nodes, lower_b, upper_b) {
      curr_node_id <- node_ids[idx]
      new_nodes <- c(current_nodes, curr_node_id)

      if (is_term[idx]) {
        p_final <- if(length(current_feats)) unique(sort(current_feats + 1L)) else integer(0L)
        paths_t[[length(paths_t) + 1L]] <<- p_final
        leaf_ids_t <<- c(leaf_ids_t, curr_node_id)
        node_paths_t[[as.character(curr_node_id)]] <<- new_nodes

        # =========================================================
        # 【核心新增】：紀錄抵達此 Leaf 時，各 Feature 最終收斂的區間
        # =========================================================
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

        # Ranger 邏輯：左子樹是 x <= splitval
        new_upper <- upper_b
        new_upper[feat] <- min(new_upper[feat], s_val)
        dfs(left_child[idx], c(current_feats, feat), new_nodes, lower_b, new_upper)

        # Ranger 邏輯：右子樹是 x > splitval
        new_lower <- lower_b
        new_lower[feat] <- max(new_lower[feat], s_val)
        dfs(right_child[idx], c(current_feats, feat), new_nodes, new_lower, upper_b)
      }
    }

    # 初始呼叫 DFS，所有 Feature 的上下限皆為 -Inf 到 Inf
    if (nrow(info) > 0) {
      dfs(1, integer(0L), integer(0L), rep(-Inf, n_features), rep(Inf, n_features))
    }

    all_paths_list <- c(all_paths_list, paths_t)
    all_path_tree_id <- c(all_path_tree_id, rep.int(t, length(paths_t)))
    all_leaf_ids <- c(all_leaf_ids, leaf_ids_t)

    leaf_info_list[[t]] <- list(
      leaf_nodes = leaf_ids_t,
      node_paths = node_paths_t,
      node_split_map = node_split_map,
      path_intervals = path_intervals_t # 將區間表包裝匯出
    )
  }

  # 全局 Feature States 排序去重
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


inspect_single_path_acam_states <- function(data_test, extracted_info, tree_id, leaf_id) {
  # 1. 取得該路徑的 ACAM 區間定義與全局 State Range
  leaf_info <- extracted_info$leaf_info_list[[tree_id]]
  path_intervals <- leaf_info$path_intervals[[as.character(leaf_id)]]
  feature_split_ranges <- extracted_info$feature_split_ranges
  
  if (nrow(path_intervals) == 0) {
    stop("This path has no feature splits (Root is Leaf).")
  }
  
  n_samples <- nrow(data_test)
  mat_test <- as.matrix(data_test)
  
  # 準備存放結果的 list
  results <- vector("list", n_samples)
  
  # 2. 針對每一筆 Testing Data 進行模擬比對
  for (i in seq_len(n_samples)) {
    is_match <- TRUE
    missed_f_id <- NA
    actual_val <- NA
    actual_state <- NA
    exp_min <- NA
    exp_max <- NA
    exp_min_state <- NA
    exp_max_state <- NA
    
    # 依序檢查這條 Path 上的每個 Feature (ACAM Cell)
    for (k in seq_len(nrow(path_intervals))) {
      f_id <- path_intervals$feature_id[k]
      min_v <- path_intervals$min_val[k]
      max_v <- path_intervals$max_val[k]
      
      val <- mat_test[i, f_id]
      
      # 判斷是否落入 ACAM Cell 允許的區間：(min_v, max_v]
      cond_pass <- (val > min_v) && (val <= max_v)
      
      if (!cond_pass) {
        is_match <- FALSE
        missed_f_id <- f_id
        actual_val <- val
        # 計算測試資料真實落點 State
        actual_state <- findInterval(val, feature_split_ranges[[f_id]])
        
        exp_min <- min_v
        exp_max <- max_v
        
        # 計算 ACAM 期望的 State 區間
        exp_min_state <- if(min_v == -Inf) 0 else findInterval(min_v, feature_split_ranges[[f_id]])
        exp_max_state <- if(max_v == Inf) length(feature_split_ranges[[f_id]]) else findInterval(max_v, feature_split_ranges[[f_id]])
        
        break # ACAM 只要有一個 Cell 不亮 (Miss)，後面的就不管了
      }
    }
    
    # 紀錄該筆 Sample 的比對結果
    results[[i]] <- data.frame(
      Sample_ID = i,
      Status = ifelse(is_match, "Matched ✅", "Missed ❌"),
      Miss_Feature = missed_f_id,
      ACAM_Range_Val = ifelse(is_match, "All Passed", sprintf("(%.3f, %.3f]", exp_min, exp_max)),
      ACAM_Range_State = ifelse(is_match, "All Passed", sprintf("[%d, %d]", exp_min_state, exp_max_state)),
      Test_Actual_Val = actual_val,
      Test_Actual_State = actual_state
    )
  }
  
  # 合併成一個完整的 DataFrame
  final_report <- do.call(rbind, results)
  return(final_report)
}


main_simulation <- function(fit,  data_test) {
  library(ranger)
  
  feat_names <- fit$forest$independent.variable.names
  n_features <- length(feat_names)
  
  # 1. 提取 (一次 DFS)
  extracted <- extract_forest_info_optimized(fit, n_features)
  
  # 2. 靜態 Mapping (Sequential Only)
  # 保持函數名稱為 compare_proposed_vs_naive
  res_static <- compare_proposed_vs_naive(
    paths_list = extracted$paths,
    n_features = n_features,
    one_based = TRUE,
    path_tree_id = extracted$tree_ids,
    leaf_ids_list = extracted$leaf_ids
  )
  show_resources(res_static)
  
  purity_stats <- calculate_tile_purity(res_static)

  # 如果你想畫圖 (Boxplot 比較不同方法的純度分布)
  library(ggplot2)
  ggplot(purity_stats, aes(x = method, y = ratio, fill = method)) +
    geom_boxplot() +
      labs(title = "Tile Purity Distribution",
         y = "Dominant Tree Ratio (Purity)",
         x = "Method") +
    theme_minimal()

  #3. 動態模擬
  cat("Running dynamic simulation...\n")
  res_dynamic <- simulate_cam_array_mismatch(
    fit, data_test, res_static, 
    leaf_info_list = extracted$leaf_info_list,
    tile_cols = 128L
  )
  cat("finish the simulation\n")
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

#save.image(file = "backup_tries_lsh_size128_20260301.RData")
#save.image(file = "backup_lsh_with_adaptive_hash_size_20260301.RData")
save.image(file = "backup_union_radix_20260301.RData")



sim_cov <- main_simulation(cov_forest_100,cov_test_x)
