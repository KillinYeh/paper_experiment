library(ranger)
library(ggplot2)

############################################################
## 1. 基礎工具與統計函式 (Utilities & Stats)
############################################################

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

precompute_bl_info <- function(physical_tiles) {
  physical_tiles$key <- paste(physical_tiles$method, physical_tiles$array_id, physical_tiles$tile_in_array, sep="_")
  bl_vec <- physical_tiles[, c("bl_open", "bl_closed", "win_len")]
  row.names(bl_vec) <- physical_tiles$key
  bl_vec
}

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
    r <- rows_idx[i]
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

calculate_tile_purity <- function(res_static) {
  if (is.null(res_static) || is.null(res_static$row_tiles)) {
    stop("Input error: res_static must contain 'row_tiles'.")
  }
  df <- res_static$row_tiles
  
  purity_list <- lapply(seq_len(nrow(df)), function(i) {
    t_ids <- unlist(df$tree_ids_in_tile[[i]])
    if (length(t_ids) == 0) {
      return(data.frame(
        method = df$method[i], tile_id = df$tile_id[i], dominant_tree = NA,
        max_count = 0, total_count = 0, ratio = 0
      ))
    }
    
    counts <- table(t_ids)
    max_count <- max(counts)
    dom_tree <- names(counts)[which.max(counts)]
    total_count <- sum(counts)
    ratio <- max_count / total_count
    
    data.frame(
      method = df$method[i], tile_id = df$tile_id[i], dominant_tree = as.integer(dom_tree),
      max_count = as.integer(max_count), total_count = as.integer(total_count), ratio = ratio
    )
  })
  
  purity_df <- do.call(rbind, purity_list)
  
  cat("\n========================================\n")
  cat("        Tile Purity Analysis Report      \n")
  cat("   (Ratio of Dominant Tree in Tile)     \n")
  cat("========================================\n")
  
  methods <- unique(purity_df$method)
  for (m in methods) {
    sub_df <- purity_df[purity_df$method == m, ]
    cat(sprintf("Method: %-10s | Mean Purity: %.4f | Median Purity: %.4f | Min/Max: %.4f/%.4f\n", 
                m, mean(sub_df$ratio), median(sub_df$ratio), min(sub_df$ratio), max(sub_df$ratio)))
  }
  cat("----------------------------------------\n")
  
  return(purity_df)
}

show_resources <- function(res_static) {
  if (is.null(res_static) || is.null(res_static$physical_tiles) || is.null(res_static$paths_ranked)) {
    stop("Input error: The input must be the result list (res_static).")
  }
  
  total_paths <- length(res_static$paths_ranked)
  total_physical_tiles <- nrow(res_static$physical_tiles)
  n_logical_arrays <- length(unique(res_static$physical_tiles$array_id))
  n_col_splits <- max(res_static$physical_tiles$tile_in_array)
  
  cat("\n========================================\n")
  cat(sprintf("1. Total Paths (Rules)     : %d\n", total_paths))
  cat(sprintf("2. Total Physical Tiles    : %d\n", total_physical_tiles))
  cat("----------------------------------------\n")
  cat(sprintf("   - Logical Row Groups    : %d (Arrays)\n", n_logical_arrays))
  cat(sprintf("   - Max Column Splits     : %d (Tiles wide)\n", n_col_splits))
  cat("========================================\n")
  
  return(invisible(data.frame(
    total_paths = total_paths, total_tiles = total_physical_tiles,
    logical_arrays = n_logical_arrays, max_col_splits = n_col_splits
  )))
}

############################################################
## 2. 路徑分裝演算法 (Packing Algorithms)
############################################################


# 1. Naive / X-Time Baseline (根據 Path 長度進行排序後直接切分)
pack_all_tiles_naive <- function(paths_ranked, tile_rows = 128L) {
  n_paths <- length(paths_ranked)

  cat("    -> Sorting paths by length (X-Time baseline)...\n")
  # 取得每條 Path 的長度 (特徵數量) 並進行排序
  path_lengths <- lengths(paths_ranked)
  order_idx <- order(path_lengths)

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

# 2. Global Union (全域暴力搜尋，時間複雜度 O(N^2))
pack_all_tiles_union <- function(paths_ranked, tile_rows = 128L) {
  # 1. 預處理：將所有 path 轉為整數向量，方便計算
  idx_full <- lapply(paths_ranked, function(v) if (!length(v)) integer(0L) else unique(as.integer(v)))

  # 2. 優先順序：短路徑優先 (Length Ascending)
  # 理由：短路徑通常是長路徑的子集 (Subset)，先定下來當基底，長路徑進來時增量可能為 0
  lens <- vapply(idx_full, length, integer(1))

  # 建立一個全局的處理順序 (尚未被使用的 indices)
  # 我們不直接改變 paths_ranked 的順序，而是用一個 indices 列表來控制
  sorted_indices <- order(lens, decreasing = FALSE)

  n_paths <- length(paths_ranked)
  used_flags <- rep(FALSE, n_paths) # 用來標記原始 index 是否已使用

  tiles <- list()
  tile_id <- 0L

  # 為了加速全域掃描，我們只在還沒用過的 indices 中搜尋

  while(any(!used_flags)) {
    tile_id <- tile_id + 1L
    rows <- integer(0L)

    # --- Step 1: 選種子 (Seed) ---
    # 從排序好的清單中，挑選第一個還沒被使用的 (即最短的)
    seed_ptr <- which(!used_flags[sorted_indices])[1]

    if (is.na(seed_ptr)) break # 全部做完

    seed_idx <- sorted_indices[seed_ptr]

    # 初始化 Tile
    rows <- c(rows, seed_idx)
    used_flags[seed_idx] <- TRUE

    # 目前 Tile 開啟的 BL (Target Mask)
    current_mask <- idx_full[[seed_idx]]

    # --- Step 2: 填滿 Tile (Greedy Union Increment) ---
    while(length(rows) < tile_rows && any(!used_flags)) {

      # 找出所有還沒被使用的候選人
      cand_indices <- which(!used_flags)
      #if( length(cand_indices) > 5000) {cand_indices <- cand_indices[1:5000]}

      # 全域掃描：計算每個候選人對 current_mask 的"增量"
      # Increment = length(union(cand, mask)) - length(mask)
      #           = length(setdiff(cand, mask))
      # 即：候選人有多少特徵是目前 Mask 沒有的

      increments <- vapply(cand_indices, function(i) {
        p <- idx_full[[i]]
        # 計算 p 中有多少個元素不在 current_mask 中
        sum(! (p %in% current_mask) )
      }, integer(1))

      # 策略：找增量最小的
      min_inc <- min(increments)

      # 如果有多個最小增量，挑選其中最短的 (Secondary Sort Key: Length)
      # 這有助於保持 Tile 內的特徵數不要暴增太快
      candidates_with_min_inc <- which(increments == min_inc)

      # 在這些最佳候選人中，找原始長度最短的
      # 我們需要回頭查 idx_full 的長度
      best_cand_local_idx <- candidates_with_min_inc[which.min(lens[cand_indices[candidates_with_min_inc]])]

      best_cand_global_idx <- cand_indices[best_cand_local_idx]

      # 加入 Tile
      rows <- c(rows, best_cand_global_idx)
      used_flags[best_cand_global_idx] <- TRUE

      # 更新 Mask (將新特徵加入)
      if (min_inc > 0) {
        new_feats <- idx_full[[best_cand_global_idx]]
        current_mask <- unique(c(current_mask, new_feats))
      }
    }

    # 存檔
    tiles[[tile_id]] <- data.frame(tile_id = tile_id, rows_in_tile = length(rows))
    tiles[[tile_id]]$rows_idx <- I(list(rows))

  }
  do.call(rbind, tiles)
}


# 3. MCP (Most Common Prefix / 純 Radix 前綴排序，無 Union 視窗)
pack_all_tiles_mcp <- function(paths_ranked, tile_rows = 128L) {
  n_paths <- length(paths_ranked)
  n_features <- attr(paths_ranked, "n_features")
  max_prefix_len <- min(n_features, 128L)

  mat <- matrix(0L, nrow = n_paths, ncol = max_prefix_len)
  for (i in seq_len(n_paths)) {
    valid_feats <- paths_ranked[[i]][paths_ranked[[i]] <= max_prefix_len]
    if (length(valid_feats) > 0) mat[i, valid_feats] <- 1L
  }

  df_mat <- as.data.frame(mat)
  order_idx <- do.call(order, c(df_mat, list(decreasing = TRUE, method = "radix")))

  group_idx <- (seq_len(n_paths) - 1L) %/% tile_rows + 1L
  chunks <- split(order_idx, group_idx)

  do.call(rbind, lapply(names(chunks), function(nm) {
    rows <- chunks[[nm]]
    df <- data.frame(tile_id = as.integer(nm), rows_in_tile = length(rows))
    df$rows_idx <- I(list(rows))
    df
  }))
}

# 4. Seq(sequential)
pack_all_tiles_sequential <- function(paths_ranked, tile_rows = 128L) {
  n_paths <- length(paths_ranked)
  # 直接按 1..N 切分，不進行排序
  order_rows <- seq_len(n_paths)

  group_idx <- (order_rows - 1L) %/% tile_rows + 1L
  chunks <- split(order_rows, group_idx)

  do.call(rbind, lapply(names(chunks), function(nm) {
    rows <- chunks[[nm]]
    df <- data.frame(tile_id = as.integer(nm), rows_in_tile = length(rows))
    df$rows_idx <- I(list(rows))
    df
  }))
}

# 5. Radix + Union
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

  # --- Step 5: Sliding Window + Union Best-Fit ---
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

# 6. LSH(adaptive)
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

# 7. fptree + union
pack_all_tiles_fptree_union <- function(paths_ranked, tile_rows = 128L, window_size = 512L) {
  n_paths <- length(paths_ranked)

  # 根據要求：只使用前 128 個 Feature 建樹，後面的直接 dropout
  max_prefix_len <- 128L

  cat("    -> Constructing Real FP-Tree (Max Feature = 128)...\n")

  # --- 步驟 1: 建立真實的 FP-Tree ---
  # 使用 environment 來模擬指標 (pointers)，建立真實的樹狀節點
  root <- new.env(parent = emptyenv())
  root$children <- new.env(parent = emptyenv())
  root$path_ids <- integer(0)

  for (i in seq_len(n_paths)) {
    feats <- paths_ranked[[i]]
    # Dropout 後面的特徵
    valid_feats <- feats[feats <= max_prefix_len]

    curr_node <- root

    if (length(valid_feats) > 0) {
      for (f in valid_feats) {
        fname <- as.character(f)
        # 如果子節點不存在，則創建新的真實節點
        if (!exists(fname, envir = curr_node$children, inherits = FALSE)) {
          new_child <- new.env(parent = emptyenv())
          new_child$children <- new.env(parent = emptyenv())
          new_child$path_ids <- integer(0)
          assign(fname, new_child, envir = curr_node$children)
        }
        # 指標往下移動 (Traversal)
        curr_node <- get(fname, envir = curr_node$children, inherits = FALSE)
      }
    }
    # 抵達該路徑在 FP-Tree 的終點，將 Path ID 存入節點
    curr_node$path_ids <- c(curr_node$path_ids, i)
  }

  cat("    -> Traversing FP-Tree (DFS)...\n")

  # --- 步驟 2: 對真實 FP-Tree 進行深度優先搜尋 (DFS) ---
  order_idx <- integer(n_paths)
  ptr <- 1L

  # 定義 DFS 遞迴函數
  dfs <- function(node) {
    # 如果這個節點有存放路徑，將其依序取出
    if (length(node$path_ids) > 0) {
      len <- length(node$path_ids)
      order_idx[ptr:(ptr + len - 1L)] <<- node$path_ids
      ptr <<- ptr + len
    }

    # 取得所有子節點的名字 (即 Feature ID)
    cnames <- ls(node$children, all.names = TRUE)
    if (length(cnames) > 0) {
      # 將名字轉回整數進行排序，確保高頻特徵 (ID 較小) 優先走訪
      cnames_int <- as.integer(cnames)
      cnames <- as.character(sort(cnames_int))

      for (cn in cnames) {
        child_node <- get(cn, envir = node$children, inherits = FALSE)
        dfs(child_node)
      }
    }
  }

  # 從 Root 啟動 DFS
  dfs(root)
  order_idx <- order_idx[1:(ptr - 1L)]

  # --- 步驟 3: Sliding Window + Union Best-Fit ---
  cat(sprintf("    -> Packing tiles with Sliding Window (size=%d)...\n", window_size))

  unassigned_idx <- order_idx
  res_tiles <- list()
  tile_counter <- 1L
  total_tiles_estimate <- ceiling(n_paths / tile_rows)

  while (length(unassigned_idx) > 0) {
    current_window_size <- min(length(unassigned_idx), window_size)
    window_pool <- unassigned_idx[1:current_window_size]

    if (length(window_pool) <= tile_rows) {
      selected_for_tile <- window_pool
      unassigned_idx <- unassigned_idx[-(1:length(selected_for_tile))]
    } else {
      # 挑選 Base Path (長度最短)
      path_lengths <- lengths(paths_ranked[window_pool])
      base_local_idx <- which.min(path_lengths)
      base_idx <- window_pool[base_local_idx]

      selected_for_tile <- integer(tile_rows)
      selected_for_tile[1] <- base_idx
      current_mask <- paths_ranked[[base_idx]]
      window_pool <- window_pool[-base_local_idx]

      # Union 精細挑選
      for (k in 2:tile_rows) {
        costs <- vapply(window_pool, function(idx) {
          sum(!(paths_ranked[[idx]] %in% current_mask))
        }, numeric(1))

        best_local_idx <- which.min(costs)
        best_idx <- window_pool[best_local_idx]

        selected_for_tile[k] <- best_idx
        current_mask <- unique(c(current_mask, paths_ranked[[best_idx]]))
        window_pool <- window_pool[-best_local_idx]
      }
      unassigned_idx <- setdiff(unassigned_idx, selected_for_tile)
    }

    res_tiles[[tile_counter]] <- data.frame(
      tile_id = tile_counter,
      rows_in_tile = length(selected_for_tile)
    )
    res_tiles[[tile_counter]]$rows_idx <- I(list(selected_for_tile))

    if (tile_counter %% 10 == 0) progress_bar(min(tile_counter, total_tiles_estimate), total_tiles_estimate)
    tile_counter <- tile_counter + 1L
  }
  progress_bar(total_tiles_estimate, total_tiles_estimate)

  do.call(rbind, res_tiles)
}


############################################################
## 3. 靜態硬體佈局邏輯 (Static Layout)
############################################################

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
## 3.1 核心比較邏輯 (執行所有 7 種 Packing 策略)
############################################################
compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE, path_tree_id, leaf_ids_list) {
  n_features <- as.integer(n_features)
  
  # 1. Feature Reordering (全域頻率排序，這是所有進階方法的基礎)
  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features
  
  cat("Packing tiles using multiple strategies...\n")
  
  # --- Method 1: Sequential (不排序，照原始樹順序切) ---
  cat("  - Running [Method 1] Sequential...\n")
  seq_df <- pack_all_tiles_sequential(paths_ranked, tile_rows = 128L)
  seq_df$method <- "sequential"
  
  # --- Method 2: Naive / X-Time (依照長度排序) ---
  cat("  - Running [Method 2] Naive / X-Time Baseline...\n")
  naive_df <- pack_all_tiles_naive(paths_ranked, tile_rows = 128L)
  naive_df$method <- "naive"
  
  # --- Method 3: LSH / MinHash ---
  cat("  - Running [Method 3] LSH (MinHash)...\n")
  lsh_df <- pack_all_tiles_lsh(paths_ranked, tile_rows = 128L, num_hash = 15L)
  lsh_df$method <- "lsh_minhash"
  
  # --- Method 4: MCP (Most Common Prefix / 純 Radix) ---
  cat("  - Running [Method 4] MCP (Radix Sort Only)...\n")
  mcp_df <- pack_all_tiles_mcp(paths_ranked, tile_rows = 128L)
  mcp_df$method <- "mcp_radix"
  
  # --- Method 5: Radix + Union ---
  cat("\n  - Running [Method 5] Radix + Union (Window=512)...\n")
  radix_union_df <- pack_all_tiles_radix_union(paths_ranked, tile_rows = 128L, window_size = 512L)
  radix_union_df$method <- "radix_union"
  
  # --- Method 6: FP-Tree + Union ---
  cat("\n  - Running [Method 6] FP-Tree + Union (Window=512)...\n")
  fptree_union_df <- pack_all_tiles_fptree_union(paths_ranked, tile_rows = 128L, window_size = 512L)
  fptree_union_df$method <- "fptree_union"

  # --- Method 7: Global Union (暴力搜尋) ---
  cat("\n  - Running [Method 7] Global Union (O(N^2) Warning!)...\n")
  global_union_df <- pack_all_tiles_union(paths_ranked, tile_rows = 128L)
  global_union_df$method <- "global_union"

  # 附加 Meta Info (將對應的 Tree ID 與 Leaf ID 貼回)
  append_meta <- function(df) {
    df$tree_ids_in_tile <- lapply(df$rows_idx, function(idx) path_tree_id[idx])
    df$leaf_ids_in_tile <- lapply(df$rows_idx, function(idx) leaf_ids_list[idx])
    df
  }
  
  # 將 7 種方法的結果綁定成一個巨大的 Data Frame
  all_row_tiles <- rbind(
    append_meta(seq_df),
    append_meta(naive_df),
    append_meta(lsh_df),
    append_meta(mcp_df),
    append_meta(radix_union_df),
    append_meta(fptree_union_df),
    append_meta(global_union_df)
  )
  
  cat("\n  - Building Physical Tiles for all methods...\n")
  physical_tiles <- build_physical_tiles(paths_ranked, all_row_tiles, tile_cols = 128L)
  
  list(
    columns_order = ord,          
    paths_ranked = paths_ranked,  
    n_features = n_features,
    path_tree_id = path_tree_id,
    row_tiles = all_row_tiles,
    physical_tiles = physical_tiles
  )
}
############################################################
## 4. 森林拓樸與狀態映射 (Forest Info & State Mapping)
############################################################

extract_forest_info_optimized <- function(fit, n_features) {
  all_paths_list   <- list()
  all_path_tree_id <- integer(0L)
  all_leaf_ids     <- integer(0L)
  leaf_info_list   <- list()
  
  all_split_vals_per_feat <- vector("list", n_features)
  
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
        p_final <- if(length(current_feats)) unique(sort(current_feats)) else integer(0L)
        paths_t[[length(paths_t) + 1L]] <<- p_final
        leaf_ids_t <<- c(leaf_ids_t, curr_node_id)
        
        used_feats <- p_final
        if (length(used_feats) > 0) {
          intervals <- data.frame(
            feature_id = used_feats, min_val = lower_b[used_feats], max_val = upper_b[used_feats]
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
    
    leaf_info_list[[t]] <- list(leaf_nodes = leaf_ids_t, path_intervals = path_intervals_t)
    
    if (t %% 10 == 0) progress_bar(t, fit$num.trees)
  }
  progress_bar(fit$num.trees, fit$num.trees)
  
  feature_split_ranges <- lapply(all_split_vals_per_feat, function(x) {
    if(is.null(x)) return(numeric(0))
    sort(unique(x))
  })
  
  list(
    paths = all_paths_list, tree_ids = all_path_tree_id, leaf_ids = all_leaf_ids,
    leaf_info_list = leaf_info_list, feature_split_ranges = feature_split_ranges
  )
}

generate_acam_state_table <- function(extracted_info) {
  state_list <- list()
  idx <- 1L
  
  for (t in seq_along(extracted_info$leaf_info_list)) {
    info <- extracted_info$leaf_info_list[[t]]
    intervals_list <- info$path_intervals
    
    for (leaf_id_str in names(intervals_list)) {
      df <- intervals_list[[leaf_id_str]]
      if (nrow(df) > 0) {
        min_states <- mapply(function(val, f_id) {
          if (val == -Inf) 0L else findInterval(val, extracted_info$feature_split_ranges[[f_id]])
        }, df$min_val, df$feature_id)
        
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
  final_df <- final_df[, c("tree_id", "leaf_id", "feature_id", "min_val", "max_val", "min_state", "max_state")]
  return(final_df)
}

generate_hardware_mapped_state_table <- function(row_tiles, global_state_table) {
  hw_map_list <- lapply(seq_len(nrow(row_tiles)), function(i) {
    method <- row_tiles$method[i]
    array_id <- row_tiles$tile_id[i]
    tree_ids <- unlist(row_tiles$tree_ids_in_tile[[i]])
    leaf_ids <- unlist(row_tiles$leaf_ids_in_tile[[i]])
    n_rows <- length(tree_ids)
    
    if (n_rows == 0) return(NULL)
    
    data.frame(
      method = method, array_id = array_id, row_in_array = seq_len(n_rows), 
      tree_id = as.integer(tree_ids), leaf_id = as.integer(leaf_ids), stringsAsFactors = FALSE
    )
  })
  
  hw_map_df <- do.call(rbind, hw_map_list)
  hw_state_table <- merge(hw_map_df, global_state_table, by = c("tree_id", "leaf_id"), all.x = TRUE)
  hw_state_table <- hw_state_table[order(hw_state_table$method, hw_state_table$array_id, hw_state_table$row_in_array, hw_state_table$feature_id), ]
  hw_state_table <- hw_state_table[, c("method", "array_id", "row_in_array", "tree_id", "leaf_id", "feature_id", "min_val", "max_val", "min_state", "max_state")]
  hw_state_table <- hw_state_table[!is.na(hw_state_table$feature_id), ]
  
  rownames(hw_state_table) <- NULL 
  return(hw_state_table)
}

############################################################
## 5. 極速實體動態模擬 (Dynamic Hardware Simulation)
############################################################

simulate_hardware_inference <- function(data_test, res_static, tile_cols = 128L) {
  hw_state_table <- res_static$hardware_state_tables
  row_tiles <- res_static$row_tiles
  
  ord <- res_static$columns_order
  feat_rank_lookup <- integer(res_static$n_features)
  feat_rank_lookup[ord] <- seq_along(ord)
  
  bl_lookup <- precompute_bl_info(res_static$physical_tiles)
  mat_test <- as.matrix(data_test)
  n_samples <- nrow(mat_test)
  methods <- unique(hw_state_table$method)
  out <- list()
  
  for (method in methods) {
    cat(sprintf("  - Simulating [%s]...\n", method))
    
    sub_hw <- hw_state_table[hw_state_table$method == method, ]
    total_paths_method <- sum(row_tiles$rows_in_tile[row_tiles$method == method])
    
    total_oo <- 0; total_ox <- 0; total_xo <- 0; total_xx <- 0
    mat_mis_tile_all <- matrix(NA_integer_, nrow=n_samples, ncol=total_paths_method)
    col_offset <- 0L
    
    array_list <- split(sub_hw, sub_hw$array_id, drop = TRUE)
    
    for (a_id_str in names(array_list)) {
      a_id <- as.integer(a_id_str)
      arr_chunk <- array_list[[a_id_str]]
      
      pt_sub <- res_static$physical_tiles[res_static$physical_tiles$method == method & 
                                            res_static$physical_tiles$array_id == a_id, ]
      n_tiles_wide <- nrow(pt_sub)
      n_rows_in_array <- max(arr_chunk$row_in_array)
      
      mat_mis_tile <- matrix(NA_integer_, nrow=n_samples, ncol=n_rows_in_array)
      row_list <- split(arr_chunk, arr_chunk$row_in_array, drop = TRUE)
      
      for (r_id_str in names(row_list)) {
        r_id <- as.integer(r_id_str)
        df_path <- row_list[[r_id_str]]
        
        df_path <- df_path[order(feat_rank_lookup[df_path$feature_id]), ]
        feats <- df_path$feature_id
        test_vals <- mat_test[, feats, drop = FALSE]
        
        pass_mat <- sweep(test_vals, 2, df_path$min_val, ">") & 
                    sweep(test_vals, 2, df_path$max_val, "<=")
        
        is_match <- rowSums(!pass_mat) == 0
        first_miss_idx <- max.col(!pass_mat, ties.method = "first")
        miss_f_id <- feats[first_miss_idx]
        
        tiles_idx <- (feat_rank_lookup[miss_f_id] - 1L) %/% tile_cols + 1L
        mat_mis_tile[!is_match, r_id] <- tiles_idx[!is_match]
      }
      
      col_range <- (col_offset + 1L):(col_offset + n_rows_in_array)
      mat_mis_tile_all[, col_range] <- mat_mis_tile
      col_offset <- col_offset + n_rows_in_array
      
      for (k in seq_len(n_tiles_wide)) {
        key <- paste(method, a_id, k, sep="_")
        if (!key %in% rownames(bl_lookup)) next
        
        bl_info <- bl_lookup[key, ]
        b_open  <- bl_info[["bl_open"]]
        b_closed <- bl_info[["bl_closed"]]
        
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
      mismatch_matrix = mat_mis_tile_all
    )
  }
  return(out)
}

############################################################
## 6. 主程式整合與執行 (Main Execution)
############################################################

main_simulation <- function(fit,  data_test) {
  feat_names <- fit$forest$independent.variable.names
  n_features <- length(feat_names)
  
  cat("\n[Step 1/5] Extracting Forest Topology & Intervals...\n")
  extracted <- extract_forest_info_optimized(fit, n_features)
  
  cat("\n[Step 2/5] Static Mapping & Hardware Packing...\n")
  res_static <- compare_proposed_vs_naive(
    paths_list = extracted$paths,
    n_features = n_features,
    one_based = TRUE,
    path_tree_id = extracted$tree_ids,
    leaf_ids_list = extracted$leaf_ids
  )
  show_resources(res_static)
  
  cat("\n[Step 3/5] Generating ACAM State Mapping Tables...\n")
  res_static$acam_state_table <- generate_acam_state_table(extracted)
  res_static$hardware_state_tables <- generate_hardware_mapped_state_table(
    row_tiles = res_static$row_tiles,
    global_state_table = res_static$acam_state_table
  )
  
  cat("\n[Step 4/5] Calculating Tile Purity...\n")
  purity_stats <- calculate_tile_purity(res_static)
  
#  # 畫出 Purity 盒鬚圖並顯示
#  p <- ggplot(purity_stats, aes(x = method, y = ratio, fill = method)) +
#    geom_boxplot() +
#    labs(title = "Tile Purity Distribution", y = "Dominant Tree Ratio (Purity)", x = "Method") +
#    theme_minimal()
#  print(p)
  
  gc(verbose = FALSE)

  cat("\n[Step 5/5] Running Physical Dynamic Simulation...\n")
  res_dynamic <- simulate_hardware_inference(
    data_test = data_test,
    res_static = res_static,
    tile_cols = 128L
  )
  
  cat("\n>>> Simulation Completed Successfully! <<<\n")
  list(static = res_static, dynamic = res_dynamic)
}
#---------------------------- Experiment -------------------------

# 建立一個輔助函數來執行、存檔並清空記憶體
run_and_save <- function(model, test_data, save_name) {
  cat(sprintf("\n========== Starting %s ==========\n", save_name))
  
  # 執行模擬
  result <- main_simulation(model, test_data)
  
  # 存成獨立的 RDS 檔案 (比 RData 更乾淨、讀取更彈性)
  saveRDS(result, file = paste0(save_name, "_0412.rds"))
  
  # 徹底刪除巨大結果，並呼叫 GC
  rm(result)
  gc(verbose = FALSE)
  cat(sprintf("========== Saved and Cleared %s ==========\n", save_name))
}

# 依序執行，保證記憶體永遠只承載一個 Dataset 的重量
run_and_save(email_forest100, email_test_x, "sim_email")
run_and_save(fetal_forest100, fetal_x,      "sim_fetal")
run_and_save(gene_forest100,  gene_test_x,  "sim_gene")
run_and_save(loan_forest100,  loan_test_x,  "sim_loan")
run_and_save(gesture_forest100, gesture_test_x, "sim_gesture")

# 警告：以下幾個因為 Rule 數量極大，建議睡前跑
run_and_save(mush_forest100,  mush_test_x,  "sim_mush")
run_and_save(gas_forest100,   gas_test_x,   "sim_gas")


