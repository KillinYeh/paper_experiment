# basic function
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

extract_forest_info_optimized <- function(fit, n_features) {
  all_paths_list   <- list()
  all_path_tree_id <- integer(0L)
  all_leaf_ids     <- integer(0L)
  leaf_info_list <- list()

  cat("Extracting paths and topology (Single Pass)...\n")

  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)

    paths_t <- list()
    leaf_ids_t <- integer(0L)
    node_paths_t <- list()
    node_split_map <- as.list(rep(NA, max(info$nodeID) + 1))
    names(node_split_map) <- 0:max(info$nodeID)

    non_terms <- info[!info$terminal, ]
    if(nrow(non_terms) > 0) {
      for(i in 1:nrow(non_terms)) {
        node_split_map[[as.character(non_terms$nodeID[i])]] <- non_terms$splitvarID[i]
      }
    }

    left_child <- info$leftChild + 1
    right_child <- info$rightChild + 1
    split_var <- info$splitvarID
    is_term <- info$terminal
    node_ids <- info$nodeID

    dfs <- function(idx, current_feats, current_nodes) {
      curr_node_id <- node_ids[idx]
      new_nodes <- c(current_nodes, curr_node_id)

      if (is_term[idx]) {
        p_final <- if(length(current_feats)) unique(sort(current_feats + 1L)) else integer(0L)
        paths_t[[length(paths_t) + 1L]] <<- p_final
        leaf_ids_t <<- c(leaf_ids_t, curr_node_id)
        node_paths_t[[as.character(curr_node_id)]] <<- new_nodes  # using leafID to index the full path
      } else {
        feat <- split_var[idx]
        dfs(left_child[idx], c(current_feats, feat), new_nodes)
        dfs(right_child[idx], c(current_feats, feat), new_nodes)
      }
    }

    if (nrow(info) > 0) dfs(1, integer(0L), integer(0L))

    all_paths_list <- c(all_paths_list, paths_t)
    all_path_tree_id <- c(all_path_tree_id, rep.int(t, length(paths_t)))
    all_leaf_ids <- c(all_leaf_ids, leaf_ids_t)

    leaf_info_list[[t]] <- list(
      leaf_nodes = leaf_ids_t,
      node_paths = node_paths_t,
      node_split_map = node_split_map
    )

  }
  progress_bar(fit$num.trees, fit$num.trees)

  list(
    paths = all_paths_list,
    tree_ids = all_path_tree_id,
    leaf_ids = all_leaf_ids,
    leaf_info_list = leaf_info_list
  )
}



############################################################
## 為了實驗補齊的輔助分組函數 (Naive/X-Time, Global Union, MCP)
############################################################

# 2. Naive / X-Time Baseline (根據 Path 長度進行排序後直接切分)
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

# 3. Global Union (全域暴力搜尋，時間複雜度 O(N^2))
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


# 4. MCP (Most Common Prefix / 純 Radix 前綴排序，無 Union 視窗)
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

# 5. Seq(sequential)
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

# 6. Radix + Union
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

# 7. LSH(adaptive)
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

# 8. fptree + union
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
## 核心時間測量實驗 (Benchmark Execution Time)
############################################################

run_packing_time_experiment <- function(paths_list, n_features, tile_rows = 128L) {
  n_features <- as.integer(n_features)
  n_paths <- length(paths_list)

  cat(sprintf("\n==================================================\n"))
  cat(sprintf("Starting Time Benchmark for %d paths...\n", n_paths))
  if (n_paths > 10000) {
    cat("WARNING: 'Union' (Global Best-Fit) is O(N^2) and might take hours!\n")
    cat("It will be automatically skipped if N > 10000 to prevent system hang.\n")
  }
  cat(sprintf("==================================================\n"))

  # --- 1. 前置處理 (不計入打包時間) ---
  cat("Preparing Frequency Reordering...\n")
  ord <- feature_frequency_order(paths_list, n_features, one_based = TRUE)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = TRUE)
  attr(paths_ranked, "n_features") <- n_features

  # 建立一個 Data Frame 來記錄時間
  time_results <- data.frame(
    Method = character(),
    Elapsed_Time_sec = numeric(),
    stringsAsFactors = FALSE
  )

  # 測量輔助函數
  measure_time <- function(method_name, func_call) {
    cat(sprintf("\nRunning %s ...\n", method_name))
    # system.time 會回傳 user, system, elapsed time，我們取真實流逝時間
    t_res <- system.time({ res <- func_call })
    cat(sprintf("-> Finished in %.3f seconds.\n", t_res["elapsed"]))
    return(t_res["elapsed"])
  }

  # --- 2. 開始測量各個方法的「打包時間」 ---

  # (1) Sequential (頻率排序後直切)
  t_seq <- measure_time("Sequential", pack_all_tiles_sequential(paths_ranked, tile_rows))
  time_results <- rbind(time_results, data.frame(Method="Seq", Elapsed_Time_sec=t_seq))

  # (2) Naive / X-Time (依 Path 長度排序後直切)
  t_naive <- measure_time("Naive (X-Time)", pack_all_tiles_naive(paths_ranked, tile_rows))
  time_results <- rbind(time_results, data.frame(Method="Naive (X-Time)", Elapsed_Time_sec=t_naive))

  # (3) MCP (Radix 前綴排序，無 Union)
  t_mcp <- measure_time("MCP", pack_all_tiles_mcp(paths_ranked, tile_rows))
  time_results <- rbind(time_results, data.frame(Method="MCP", Elapsed_Time_sec=t_mcp))

  # (4) Radix + Union (Window=512)
  t_radix_union <- measure_time("Radix + Union", pack_all_tiles_radix_union(paths_ranked, tile_rows, window_size = 512L))
  time_results <- rbind(time_results, data.frame(Method="Radix + Union", Elapsed_Time_sec=t_radix_union))

  # (5) LSH (Adaptive MinHash)
  # 假設您已經將前面的 adaptive LSH 實作加入了程式碼中
  t_lsh <- measure_time("LSH (Adaptive)", pack_all_tiles_lsh(paths_ranked, tile_rows))
  time_results <- rbind(time_results, data.frame(Method="LSH (Adaptive)", Elapsed_Time_sec=t_lsh))

  # (6) FP-Tree + Union (Window=512)
  t_fptree <- measure_time("FP-Tree + Union", pack_all_tiles_fptree_union(paths_ranked, tile_rows, window_size = 512L))
  time_results <- rbind(time_results, data.frame(Method="FP-Tree + Union", Elapsed_Time_sec=t_fptree))

  # (7) Global Union (O(N^2) 暴力搜尋)
  t_union <- measure_time("Global Union", pack_all_tiles_union(paths_ranked, tile_rows))
  time_results <- rbind(time_results, data.frame(Method="Global Union", Elapsed_Time_sec=t_union))

  # --- 3. 輸出結果與繪圖 ---
  cat("\n==================================================\n")
  cat("Execution Time Summary:\n")
  print(time_results)
  cat("==================================================\n")

  #library(ggplot2)
  # 為了讓圖表好看，我們將 Method 根據耗時由短到長排序
  #p <- ggplot(time_results, aes(x = reorder(Method, Elapsed_Time_sec), y = Elapsed_Time_sec, fill = Method)) +
  #  geom_bar(stat = "identity", color="black", alpha=0.8) +
  #  coord_flip() +
  #  labs(title = "Execution Time Comparison for Tile Packing",
  #       subtitle = sprintf("Dataset Size: %d Paths", n_paths),
  #       x = "Packing Method",
  #       y = "Elapsed Time (Seconds)") +
  #  theme_minimal() +
  #  theme(legend.position = "none")
  #print(p)

  return(time_results)
}


main_simulation <- function(fit)
{
	feat_names <- fit$forest$independent.variable.names
	n_features <- length(feat_names)

	cat("Step 1: Extracting paths from the Random Forest model...\n")
# 2. 呼叫之前的萃取函數，把樹的結構轉成路徑清單 (paths_list)
	extracted_info <- extract_forest_info_optimized(fit, n_features)
	paths_list <- extracted_info$paths

	cat("Step 2: Starting the Packing Time Benchmark...\n")
# 3. 把萃取出來的路徑，丟進您剛剛寫好的 Benchmark 函數中
	time_results <- run_packing_time_experiment(
  	paths_list = paths_list,
  	n_features = n_features,
  	tile_rows = 128L
	)

# 查看最終的時間比較表
	print(time_results)
}


#----------- experiment -----------------
time_email <- main_simulation(email_forest100)

time_fetal <- main_simulation(fetal_forest100)

time_gene  <- main_simulation(gene_forest100)

time_mush  <- main_simulation(mush_forest100)

time_loan  <- main_simulation(loan_forest100)

time_gesture <- main_simulation(gesture_forest100)

time_arcene <- main_simulation(arcene_forest100)

time_gas <- main_simulation(gas_forest100)

#time_cov <- main_simulation(cov_forest_100,cov_test_x)
save.image(file = "backup_packing_time_comparison.RData")


