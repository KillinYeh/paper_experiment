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

wlbl_stats_for_window <- function(paths_ranked, rows_idx, start, end, feature_bit_lengths_ranked) {
  n_rows <- length(rows_idx)
  win_len <- sum(feature_bit_lengths_ranked[start:end])

  if (win_len <= 0L || n_rows == 0L) {
    return(data.frame(
      wl_open=0L, wl_closed=0L, bl_open=0L, bl_closed=0L,
      cell_oo=0L, cell_ox=0L, cell_xo=0L, cell_xx=0L,
      win_len=win_len,
      bl_open_cols = I(list(integer(0L))),
      bl_closed_cols = I(list(integer(0L))),
      wl_open_rows = I(list(integer(0L))),
      wl_closed_rows = I(list(integer(0L)))
    ))
  }

  wl_open_vec <- logical(n_rows)
  bl_used <- rep(FALSE, end - start + 1L)

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

  cols <- start:end
  bit_lengths_in_window <- feature_bit_lengths_ranked[cols]

  bl_open   <- sum(bit_lengths_in_window[bl_used])
  bl_closed <- win_len - bl_open

  data.frame(
    wl_open = wl_open, wl_closed = wl_closed,
    bl_open = bl_open, bl_closed = bl_closed,
    cell_oo = wl_open * bl_open,   cell_ox = wl_open * bl_closed,
    cell_xo = wl_closed * bl_open, cell_xx = wl_closed * bl_closed,
    win_len = win_len,
    bl_open_cols = I(list(cols[bl_used])),
    bl_closed_cols = I(list(cols[!bl_used])),
    wl_open_rows = I(list(which(wl_open_vec))),
    wl_closed_rows = I(list(which(!wl_open_vec)))
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
      max_count = as.integer(max_count), total_count = as.integer(total_count),ratio = ratio
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
m, mean(sub_df$ratio), median(sub_df$ratio), min(sub_df$ratio),max(sub_df$ratio)))
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
      #= length(setdiff(cand, mask))
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
      # 注意：R 的 setdiff 會自動保留 unassigned_idx 原本的排序，完美維持 RadixSort 的效果
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

build_physical_tiles <- function(paths_ranked, row_tiles, feature_bit_lengths_ranked, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L || nrow(row_tiles) == 0L) return(NULL)

  cut_points <- list()
  curr_start <- 1L
  curr_sum <- 0L

  for (f in seq_len(n_features)) {
    bit_len <- feature_bit_lengths_ranked[f]

    if (curr_sum + bit_len > tile_cols && curr_sum > 0L) {
      cut_points[[length(cut_points) + 1L]] <- c(curr_start, f - 1L)
      curr_start <- f
      curr_sum <- bit_len
    } else {
      curr_sum <- curr_sum + bit_len
    }
  }
  cut_points[[length(cut_points) + 1L]] <- c(curr_start, n_features)

  res <- list()
  idx <- 0L

  for (i in seq_len(nrow(row_tiles))) {
    rows_idx <- row_tiles$rows_idx[[i]]
    array_id <- row_tiles$tile_id[i]
    method   <- row_tiles$method[i]

    tile_in_array <- 0L
    for (cp in cut_points) {
      start <- cp[1]
      end <- cp[2]
      tile_in_array <- tile_in_array + 1L
      win_stats <- wlbl_stats_for_window(
        paths_ranked = paths_ranked,
        rows_idx = rows_idx,
        start = start,
        end = end,
        feature_bit_lengths_ranked = feature_bit_lengths_ranked
      )

      idx <- idx + 1L
      res[[idx]] <- cbind(
        data.frame(method=method, array_id=array_id, tile_in_array=tile_in_array,
start_col=start, end_col=end),
        win_stats
      )
    }
  }
  do.call(rbind, res)
}

############################################################
## 3.1 核心比較邏輯 (執行所有 7 種 Packing 策略)
############################################################
compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE,path_tree_id, leaf_ids_list, feature_bit_lengths) {
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
  fptree_union_df <- pack_all_tiles_fptree_union(paths_ranked, tile_rows = 128L,window_size = 512L)
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
  feature_bit_lengths_ranked <- feature_bit_lengths[ord]

  cat("\n  - Building Physical Tiles for all methods...\n")
  physical_tiles <- build_physical_tiles(paths_ranked, all_row_tiles, feature_bit_lengths_ranked, tile_cols= 128L)

  list(
    ord = ord,
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

generate_acam_state_table <- function(extracted_info, ord = NULL) {
  state_list <- list()
  idx <- 1L

  original_to_ranked <- NULL
  if (!is.null(ord)) {
    original_to_ranked <- integer(length(ord))
    original_to_ranked[ord] <- seq_along(ord)
  }

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

        if (!is.null(original_to_ranked)) {
          df$original_feature_id <- as.integer(df$feature_id)
          df$feature_id <- as.integer(original_to_ranked[df$original_feature_id])
        }

        state_list[[idx]] <- df
        idx <- idx + 1L
      }
    }
  }

  final_df <- do.call(rbind, state_list)
  if ("original_feature_id" %in% names(final_df)) {
    final_df <- final_df[, c("tree_id", "leaf_id", "feature_id", "original_feature_id", "min_val", "max_val", "min_state", "max_state")]
  } else {
    final_df <- final_df[, c("tree_id", "leaf_id", "feature_id", "min_val", "max_val", "min_state", "max_state")]
  }
  return(final_df)
}

generate_hardware_mapped_state_table <- function(row_tiles, global_state_table){
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
  keep_cols <- c("method", "array_id", "row_in_array", "tree_id", "leaf_id", "feature_id")
  if ("original_feature_id" %in% names(hw_state_table)) {
    keep_cols <- c(keep_cols, "original_feature_id")
  }
  keep_cols <- c(keep_cols, "min_val", "max_val", "min_state", "max_state")
  hw_state_table <- hw_state_table[, keep_cols]
  hw_state_table <- hw_state_table[!is.na(hw_state_table$feature_id), ]

  rownames(hw_state_table) <- NULL
  return(hw_state_table)
}


############################################################
## 5. 功耗模型與靜態能量分析 (Power Model & Static Energy)
############################################################

# 向量化版本的 bit length 計算
get_bit_len_vec <- function(states, radix_base = 8L) {
  states <- as.numeric(states)
  out <- rep(1L, length(states))
  valid <- !is.na(states) & is.finite(states) & states >= radix_base
  out[valid] <- as.integer(ceiling(log(states[valid] + 1, base = radix_base)))
  out[is.infinite(states)] <- NA_integer_
  out
}

get_feature_bit_lengths <- function(extracted_info, n_features, radix_base = 8L) {
  max_states <- vapply(seq_len(n_features), function(f_id) {
    length(extracted_info$feature_split_ranges[[f_id]])
  }, integer(1))

  bit_lengths <- get_bit_len_vec(max_states, radix_base)

  data.frame(
    feature_id = seq_len(n_features),
    max_global_state = max_states,
    bit_length = bit_lengths
  )
}

split_state_to_digits <- function(state, bit_length, radix_base = 8L) {
  if (is.na(state) || bit_length <= 0L) return(integer(0L))

  digits <- integer(bit_length)
  value <- as.integer(state)
  for (i in seq_len(bit_length)) {
    digits[i] <- value %% radix_base
    value <- value %/% radix_base
  }
  digits
}

calculate_masked_column_cells <- function(paths_ranked, rows_idx, start, end, feature_bit_lengths_ranked) {
  if (length(rows_idx) == 0L || end < start) {
    return(list(masked_cells = 0, masked_columns = 0, masked_bit_columns = 0))
  }

  cols <- start:end
  used <- rep(FALSE, length(cols))

  for (r in rows_idx) {
    feats <- paths_ranked[[r]]
    if (!length(feats)) next
    idx <- feats[feats >= start & feats <= end]
    if (length(idx)) used[idx - start + 1L] <- TRUE
  }

  masked_cols <- cols[!used]
  if (!length(masked_cols)) {
    return(list(masked_cells = 0, masked_columns = 0, masked_bit_columns = 0))
  }

  masked_bit_columns <- sum(feature_bit_lengths_ranked[masked_cols])
  list(
    masked_cells = length(rows_idx) * masked_bit_columns,
    masked_columns = length(masked_cols),
    masked_bit_columns = masked_bit_columns
  )
}

calculate_active_cell_expansion <- function(hw_state_table, radix_base = 8L) {
  if (is.null(hw_state_table) || nrow(hw_state_table) == 0L) {
    return(list(active_summary = data.frame(), expanded_active_cells = 0))
  }

  # 使用每一筆 rule 的 max_state 決定該 cell 需要幾個 base-8 digit。
  # 例如 range = 50~900 時，即使 M1 的 min_state 只寫 50，
  # 仍需要能表示 max_state = 900 的 bit length。
  hw_state_table$bit_length <- get_bit_len_vec(hw_state_table$max_state, radix_base)

  active_summary <- aggregate(bit_length ~ method, data = hw_state_table, FUN = sum)
  names(active_summary)[2] <- "expanded_active_cells"

  list(
    active_summary = active_summary,
    expanded_active_cells = sum(active_summary$expanded_active_cells)
  )
}

calculate_aCAM_write_energy <- function(prev_state, curr_state,
                                        write_down_pj = 7.5,
                                        write_up_pj = 4.0) {
  ifelse(
    curr_state == prev_state,
    0,
    ifelse(curr_state < prev_state, write_down_pj, write_up_pj)
  )
}

make_default_state_matrix <- function(n_rows, n_cols, radix_base = 8L) {
  list(
    # m1 / m2 儲存每一個 physical aCAM cell 中兩個 RRAM 的 digit state。
    # 預設值就是 physical Don't Care: M1 = HRS/state 0, M2 = LRS/state 7。
    m1 = matrix(0L, nrow = n_rows, ncol = n_cols),
    m2 = matrix(radix_base - 1L, nrow = n_rows, ncol = n_cols),

    # active_mask 表示這個 physical cell 是否屬於某個 path 使用到的 feature。
    # 若某一側 bound 不存在，build_tile_state_matrix() 會直接把該側 RRAM
    # 寫成對應的 Don't Care state：M1 = state 0 / HRS，M2 = state 7 / LRS。
    # 因此 write energy 只需要根據 active_mask 同時計算 M1 與 M2 的 transition，
    # 不需要再另外維護 partial Don't Care mask。
    active_mask = matrix(FALSE, nrow = n_rows, ncol = n_cols)
  )
}

align_prev_state_matrix <- function(prev_state, n_rows, n_cols, radix_base = 8L) {
  aligned <- make_default_state_matrix(n_rows, n_cols, radix_base)

  if (is.null(prev_state)) return(aligned)

  common_rows <- min(n_rows, nrow(prev_state$m1))
  common_cols <- min(n_cols, ncol(prev_state$m1))

  if (common_rows > 0L && common_cols > 0L) {
    aligned$m1[seq_len(common_rows), seq_len(common_cols)] <-
      prev_state$m1[seq_len(common_rows), seq_len(common_cols)]
    aligned$m2[seq_len(common_rows), seq_len(common_cols)] <-
      prev_state$m2[seq_len(common_rows), seq_len(common_cols)]
  }

  aligned
}

build_tile_state_matrix <- function(pt, rt, hw_state_table_tile,
                                    feature_bit_lengths_ranked, radix_base = 8L) {
  rows_idx <- rt$rows_idx[[1]]
  n_rows <- length(rows_idx)
  n_cols <- as.integer(pt$win_len)

  tile_state <- make_default_state_matrix(n_rows, n_cols, radix_base)

  if (n_rows == 0L || n_cols == 0L || nrow(hw_state_table_tile) == 0L) {
    return(tile_state)
  }

  ranked_cols <- as.integer(pt$start_col):as.integer(pt$end_col)
  bit_lengths <- feature_bit_lengths_ranked[ranked_cols]

  # 假設 bit length 為 c(2, 1, 3)，則 col_starts 為 c(1, 3, 4)，
  # 代表這個 feature 在 expanded physical tile 中的起始 physical column。
  col_starts <- cumsum(c(1L, head(bit_lengths, -1L)))

  for (i in seq_len(nrow(hw_state_table_tile))) {
    # 這裡的 feature_id 已經在 generate_acam_state_table() 中被轉成 ranked feature ID。
    ord <- as.integer(hw_state_table_tile$feature_id[i])

    if (is.na(ord) || ord < pt$start_col || ord > pt$end_col) next

    local_feature_idx <- ord - pt$start_col + 1L
    physical_col_start <- col_starts[local_feature_idx]

    # 依照目前討論，這筆 rule 的實體 digit 長度主要由 max_state 決定。
    # 例如 range = 50~900 時，即使 M1/min_state 只需要較少 digit，
    # 這個 aCAM macro-cell 仍需要能表示 max_state = 900 的長度。
    #
    # 注意：max_val = Inf 代表 upper-bound Don't Care。正常情況下，
    # generate_acam_state_table() 已經會把 Inf 轉成該 feature 的最大 finite state，
    # 所以 max_state 不應該再是 Inf。若舊資料中仍出現 Inf/NA，這裡退回使用
    # 該 ranked feature 的 global bit length，避免 log(Inf) 造成 bit_len 錯誤。
    max_state_i <- hw_state_table_tile$max_state[i]
    bit_len <- get_bit_len_vec(max_state_i, radix_base)
    if (is.na(bit_len) || bit_len <= 0L) {
      bit_len <- bit_lengths[local_feature_idx]
    }

    if (bit_len <= 0L) next

    row_id <- as.integer(hw_state_table_tile$row_in_array[i])
    if (row_id < 1L || row_id > n_rows) next

    # 進一步拆開左右兩側 RRAM 的 active 狀態。
    # lower_has_rule = FALSE 代表 lower bound 是 -Inf，M1 是 lower-bound Don't Care，
    # 應固定為 State 0 / HRS。
    # upper_has_rule = FALSE 代表 upper bound 是 Inf，M2 是 upper-bound Don't Care，
    # 應固定為 State 7 / LRS。
    #
    # 若 min_state 的 digit 長度短於 max_state，M1 高位 digit 不是「不存在」，
    # 而是補成 State 0 / HRS。這樣才符合逐 digit 的 base-8 range encoding。
    lower_has_rule <- !is.infinite(hw_state_table_tile$min_val[i])
    upper_has_rule <- !is.infinite(hw_state_table_tile$max_val[i])

    m1_digits <- if (lower_has_rule) {
      split_state_to_digits(hw_state_table_tile$min_state[i], bit_len, radix_base)
    } else {
      rep(0L, bit_len)
    }

    m2_digits <- if (upper_has_rule) {
      split_state_to_digits(hw_state_table_tile$max_state[i], bit_len, radix_base)
    } else {
      rep(radix_base - 1L, bit_len)
    }

    physical_col_end <- physical_col_start + bit_len - 1L
    physical_cols <- physical_col_start:physical_col_end

    tile_state$m1[row_id, physical_cols] <- m1_digits
    tile_state$m2[row_id, physical_cols] <- m2_digits

    # cell_active_mask: 這些 physical columns 屬於這個 row/path 使用到的 feature。
    # 即使其中一側 bound 不存在，該側 RRAM 也已經在上面被寫成 Don't Care state，
    # 所以這個 active cell 的 M1 / M2 都會在 write energy 中一起計算。
    tile_state$active_mask[row_id, physical_cols] <- TRUE
  }

  tile_state
}

calculate_tile_write_energy <- function(prev_state, curr_state,
                                        write_down_pj = 7.5,
                                        write_up_pj = 4.0,
                                        radix_base = 8L) {
  n_rows <- nrow(curr_state$m1)
  n_cols <- ncol(curr_state$m1)
  prev_aligned <- align_prev_state_matrix(prev_state, n_rows, n_cols, radix_base)

  m1_write_pj <- sum(calculate_aCAM_write_energy(
    prev_state = prev_aligned$m1,
    curr_state = curr_state$m1,
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj
  ))

  m2_write_pj <- sum(calculate_aCAM_write_energy(
    prev_state = prev_aligned$m2,
    curr_state = curr_state$m2,
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj
  ))

  list(
    m1_write_pj = m1_write_pj,
    m2_write_pj = m2_write_pj,
    total_write_pj = m1_write_pj + m2_write_pj
  )
}

calculate_tile_write_energy_for_rows <- function(prev_state, curr_state, rows_to_write,
                                                 write_down_pj = 7.5,
                                                 write_up_pj = 4.0,
                                                 radix_base = 8L) {
  n_rows <- nrow(curr_state$m1)
  n_cols <- ncol(curr_state$m1)
  prev_aligned <- align_prev_state_matrix(prev_state, n_rows, n_cols, radix_base)

  rows_to_write <- rows_to_write[rows_to_write >= 1L & rows_to_write <= n_rows]

  if (!length(rows_to_write) || n_cols == 0L) {
    return(list(m1_write_pj = 0, m2_write_pj = 0, total_write_pj = 0))
  }

  m1_write_pj <- sum(calculate_aCAM_write_energy(
    prev_state = prev_aligned$m1[rows_to_write, , drop = FALSE],
    curr_state = curr_state$m1[rows_to_write, , drop = FALSE],
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj
  ))

  m2_write_pj <- sum(calculate_aCAM_write_energy(
    prev_state = prev_aligned$m2[rows_to_write, , drop = FALSE],
    curr_state = curr_state$m2[rows_to_write, , drop = FALSE],
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj
  ))

  list(
    m1_write_pj = m1_write_pj,
    m2_write_pj = m2_write_pj,
    total_write_pj = m1_write_pj + m2_write_pj
  )
}

calculate_tile_active_write_energy_for_rows <- function(prev_state, curr_state, rows_to_write,
                                                        write_down_pj = 7.5,
                                                        write_up_pj = 4.0,
                                                        radix_base = 8L) {
  n_rows <- nrow(curr_state$m1)
  n_cols <- ncol(curr_state$m1)
  prev_aligned <- align_prev_state_matrix(prev_state, n_rows, n_cols, radix_base)

  rows_to_write <- rows_to_write[rows_to_write >= 1L & rows_to_write <= n_rows]

  if (!length(rows_to_write) || n_cols == 0L) {
    return(list(m1_write_pj = 0, m2_write_pj = 0, total_write_pj = 0, active_rule_cells = 0))
  }
  # active_mask 代表 Tile 裡面實際有使用到的 row , column
  active_mask <- curr_state$active_mask
  if (is.null(active_mask)) {
    active_mask <- matrix(TRUE, nrow = n_rows, ncol = n_cols)
  }

  row_mask <- matrix(FALSE, nrow = n_rows, ncol = n_cols)
  row_mask[rows_to_write, ] <- TRUE

  # active_cell_mask 是真正有 feature rule 的 physical cells。
  # 對每個 active cell，M1 與 M2 都是實體存在且需要被維持在正確狀態的 RRAM。
  # 若 lower bound 不存在，M1 已在 curr_state$m1 中設為 state 0 / HRS。
  # 若 upper bound 不存在，M2 已在 curr_state$m2 中設為 state 7 / LRS。（透過 rep(radix_base - 1L, bit_len) ）
  # 因此這裡直接用同一個 active_cell_mask 計算 M1 與 M2 的 write transition，
  # 不再另外拆 partial Don't Care RRAM。
  active_cell_mask <- active_mask & row_mask
  active_rule_cells <- sum(active_cell_mask)

  if (active_rule_cells == 0L) {
    return(list(m1_write_pj = 0, m2_write_pj = 0, total_write_pj = 0, active_rule_cells = 0))
  }

  m1_write_pj <- sum(calculate_aCAM_write_energy(
    prev_state = prev_aligned$m1[active_cell_mask],
    curr_state = curr_state$m1[active_cell_mask],
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj
  ))

  m2_write_pj <- sum(calculate_aCAM_write_energy(
    prev_state = prev_aligned$m2[active_cell_mask],
    curr_state = curr_state$m2[active_cell_mask],
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj
  ))

  list(
    m1_write_pj = m1_write_pj,
    m2_write_pj = m2_write_pj,
    total_write_pj = m1_write_pj + m2_write_pj,
    active_rule_cells = active_rule_cells
  )
}

update_residual_state_for_rows <- function(prev_state, curr_state, rows_to_write,
                                           radix_base = 8L) {
  n_rows <- nrow(curr_state$m1)
  n_cols <- ncol(curr_state$m1)
  residual_state <- align_prev_state_matrix(prev_state, n_rows, n_cols, radix_base)

  rows_to_write <- rows_to_write[rows_to_write >= 1L & rows_to_write <= n_rows]

  if (length(rows_to_write) && n_cols > 0L) {
    residual_state$m1[rows_to_write, ] <- curr_state$m1[rows_to_write, , drop = FALSE]
    residual_state$m2[rows_to_write, ] <- curr_state$m2[rows_to_write, , drop = FALSE]
  }

  residual_state
}

calculate_write_energy_stats <- function(res_static, feature_bit_lengths_ranked,
                                         dont_care_write_cell_pj = 23,
                                         write_down_pj = 7.5,
                                         write_up_pj = 4.0,
                                         radix_base = 8L) {
  if (is.null(res_static$physical_tiles) || is.null(res_static$row_tiles) ||
      is.null(res_static$hardware_state_tables)) {
    stop("Input error: res_static must contain physical_tiles, row_tiles, and hardware_state_tables.")
  }

  physical_tiles <- res_static$physical_tiles
  physical_tiles$tile_order <- seq_len(nrow(physical_tiles))

  write_records <- list()
  rec_idx <- 0L

  for (m in unique(physical_tiles$method)) {
    sub_tiles <- physical_tiles[physical_tiles$method == m, ]
    sub_tiles <- sub_tiles[order(sub_tiles$array_id, sub_tiles$tile_in_array), ]
    prev_state <- NULL

    for (i in seq_len(nrow(sub_tiles))) {
      pt <- sub_tiles[i, ]
      rt <- res_static$row_tiles[
        res_static$row_tiles$method == pt$method &
          res_static$row_tiles$tile_id == pt$array_id,
      ]

      if (nrow(rt) == 0L) next

      hw_state_table_tile <- res_static$hardware_state_tables[
        res_static$hardware_state_tables$method == pt$method &
          res_static$hardware_state_tables$array_id == pt$array_id,
      ]

      curr_state <- build_tile_state_matrix(
        pt = pt,
        rt = rt,
        hw_state_table_tile = hw_state_table_tile,
        feature_bit_lengths_ranked = feature_bit_lengths_ranked,
        radix_base = radix_base
      )

      rows_in_tile <- length(rt$rows_idx[[1]])
      wl_open_rows <- if ("wl_open_rows" %in% names(sub_tiles)) {
        pt$wl_open_rows[[1]]
      } else {
        seq_len(rows_in_tile)
      }

      tile_write <- calculate_tile_active_write_energy_for_rows(
        prev_state = prev_state,
        curr_state = curr_state,
        rows_to_write = wl_open_rows,
        write_down_pj = write_down_pj,
        write_up_pj = write_up_pj,
        radix_base = radix_base
      )

      total_write_cells <- as.integer(pt$cell_oo + pt$cell_ox)
      masked_write_cells <- as.integer(pt$cell_ox)
      active_rule_cells <- as.integer(tile_write$active_rule_cells)
      scattered_write_cells <- max(0L, as.integer(pt$cell_oo) - active_rule_cells)

      scattered_write_pj <- scattered_write_cells * dont_care_write_cell_pj
      masked_write_pj <- masked_write_cells * dont_care_write_cell_pj
      fake_match_write_saved_pj <- masked_write_pj

      # tile_write$total_write_pj 指的是所有 cell_oo 且非處於 don't care 的那些 cell 貢獻的寫入能量
      # scattered_write_pj 指的是 cell_oo 且處於 don't care 的那些 cell 貢獻的寫入能量
      baseline_write_pj <- tile_write$total_write_pj + scattered_write_pj + masked_write_pj
      rewrite_write_pj <- baseline_write_pj
      fake_match_write_pj <- tile_write$total_write_pj + scattered_write_pj

      rec_idx <- rec_idx + 1L
      write_records[[rec_idx]] <- data.frame(
        method = pt$method,
        array_id = pt$array_id,
        tile_in_array = pt$tile_in_array,
        rows_in_tile = rows_in_tile,
        win_len = as.integer(pt$win_len),
        total_write_cells = total_write_cells,
        active_rule_cells_write = active_rule_cells,
        scattered_write_cells = scattered_write_cells,
        masked_write_cells = masked_write_cells,
        baseline_write_pj = baseline_write_pj,
        rewrite_write_pj = rewrite_write_pj,
        fake_match_write_saved_pj = fake_match_write_saved_pj,
        fake_match_write_pj = fake_match_write_pj,
        m1_write_pj = tile_write$m1_write_pj,
        m2_write_pj = tile_write$m2_write_pj,
        scattered_write_pj = scattered_write_pj,
        masked_write_pj = masked_write_pj
      )

      prev_state <- update_residual_state_for_rows(
        prev_state = prev_state,
        curr_state = curr_state,
        rows_to_write = wl_open_rows,
        radix_base = radix_base
      )
    }
  }

  write_tile_table <- do.call(rbind, write_records)

  write_summary <- aggregate(
    cbind(total_write_cells, active_rule_cells_write, scattered_write_cells,
          masked_write_cells, baseline_write_pj, rewrite_write_pj,
          fake_match_write_saved_pj, fake_match_write_pj,
          m1_write_pj, m2_write_pj, scattered_write_pj, masked_write_pj) ~ method,
    data = write_tile_table,
    FUN = sum
  )

  write_summary$baseline_write_avg_pj_per_cell <- ifelse(
    write_summary$total_write_cells == 0,
    NA_real_,
    write_summary$baseline_write_pj / write_summary$total_write_cells
  )
  write_summary$rewrite_write_avg_pj_per_cell <- ifelse(
    write_summary$total_write_cells == 0,
    NA_real_,
    write_summary$rewrite_write_pj / write_summary$total_write_cells
  )
  write_summary$fake_match_write_avg_pj_per_cell <- ifelse(
    write_summary$total_write_cells == 0,
    NA_real_,
    write_summary$fake_match_write_pj / write_summary$total_write_cells
  )

  # 對於每一個 cell，引入 Fake Match 平均寫入能量是原本的幾倍
  write_summary$fake_over_rewrite_write_ratio <- ifelse(
    write_summary$rewrite_write_avg_pj_per_cell == 0,
    NA_real_,
    write_summary$fake_match_write_avg_pj_per_cell / write_summary$rewrite_write_avg_pj_per_cell
  )
  # 引入 Fake Match 省下原本總寫入能量的比例
  write_summary$fake_write_reduction_ratio <- ifelse(
    write_summary$baseline_write_pj == 0,
    NA_real_,
    write_summary$fake_match_write_saved_pj / write_summary$baseline_write_pj
  )

  list(
    write_tile_table = write_tile_table,
    write_summary = write_summary[order(write_summary$method), ]
  )
}

calculate_column_mask_energy <- function(res_static, extracted_info, tile_cols = 128L,
                                         base_cell_fj = 13.33,
                                         dont_care_cell_fj = 0.25,
                                         fake_match_cell_fj = 10.88,
                                         radix_base = 8L,
                                         dont_care_write_cell_pj = 23,
                                         write_down_pj = 7.5,
                                         write_up_pj = 4.0) {
  if (is.null(res_static$physical_tiles) || is.null(res_static$row_tiles)) {
    stop("Input error: res_static must contain physical_tiles and row_tiles.")
  }

  n_features <- res_static$n_features
  feature_bit_table <- get_feature_bit_lengths(extracted_info, n_features, radix_base)
  feature_bit_lengths_ranked <- feature_bit_table$bit_length[res_static$ord]

  masked_tile_table <- res_static$physical_tiles
  masked_tile_table$masked_columns <- masked_tile_table$bl_closed
  masked_tile_table$masked_cells <- masked_tile_table$cell_ox
  masked_tile_table$rewrite_dont_care_fj <- masked_tile_table$masked_cells * dont_care_cell_fj
  masked_tile_table$fake_match_fj <- masked_tile_table$masked_cells * fake_match_cell_fj

  wlbl_summary <- aggregate(
    cbind(cell_oo, cell_ox, cell_xo, cell_xx, masked_columns, masked_cells,
          rewrite_dont_care_fj, fake_match_fj) ~ method,
    data = masked_tile_table,
    FUN = sum
  )

  active_expansion <- calculate_active_cell_expansion(
    hw_state_table = res_static$hardware_state_tables,
    radix_base = radix_base
  )

  method_paths <- aggregate(rows_in_tile ~ method, data = res_static$row_tiles, FUN = sum)
  method_tiles <- aggregate(tile_id ~ method, data = res_static$row_tiles, FUN = length)
  names(method_tiles)[2] <- "logical_tiles"

  method_summary <- merge(wlbl_summary, active_expansion$active_summary, by = "method", all = TRUE)
  method_summary <- merge(method_summary, method_paths, by = "method", all.x = TRUE)
  method_summary <- merge(method_summary, method_tiles, by = "method", all.x = TRUE)

  numeric_cols <- setdiff(names(method_summary), "method")
  for (nm in numeric_cols) method_summary[[nm]][is.na(method_summary[[nm]])] <- 0

  # 只有 WL open 的 cells 會參與本次 read/evaluate 平均功耗：
  # cell_oo: WL open + BL open
  # cell_ox: WL open + BL closed
  method_summary$total_cells_with_bit_slicing <-
    method_summary$cell_oo + method_summary$cell_ox

  # 真正可以被整條 column masking 的 read/evaluate cells 只有 cell_ox
  method_summary$masked_column_cells <- method_summary$cell_ox

  # 不是整條 column 都是 Don't Care 的零星 Don't Care cells，只作為觀察指標
  method_summary$scattered_dont_care_cells <- pmax(
    0,
    method_summary$cell_oo - method_summary$expanded_active_cells
  )

  # Baseline read/evaluate energy：
  # ACAMRF 的 base_cell_fj 視為平均每個 evaluated cell 的功耗
  method_summary$baseline_total_fj <-
    method_summary$total_cells_with_bit_slicing * base_cell_fj

  # Rewrite policy：將 cell_ox 從 baseline energy 替換成 physical Don't Care energy
  method_summary$rewrite_total_fj <-
    method_summary$baseline_total_fj -
    method_summary$masked_column_cells * (base_cell_fj - dont_care_cell_fj)

  # Fake Match policy：將 cell_ox 從 baseline energy 替換成 Fake Match energy
  method_summary$fake_match_total_fj <-
    method_summary$baseline_total_fj -
    method_summary$masked_column_cells * (base_cell_fj - fake_match_cell_fj)

  method_summary$baseline_avg_fj_per_cell <- ifelse(
    method_summary$total_cells_with_bit_slicing == 0,
    NA_real_,
    method_summary$baseline_total_fj / method_summary$total_cells_with_bit_slicing
  )

  method_summary$rewrite_avg_fj_per_cell <- ifelse(
    method_summary$total_cells_with_bit_slicing == 0,
    NA_real_,
    method_summary$rewrite_total_fj / method_summary$total_cells_with_bit_slicing
  )

  method_summary$fake_match_avg_fj_per_cell <- ifelse(
    method_summary$total_cells_with_bit_slicing == 0,
    NA_real_,
    method_summary$fake_match_total_fj / method_summary$total_cells_with_bit_slicing
  )

  method_summary$fake_over_rewrite_avg_ratio <- ifelse(
    method_summary$rewrite_avg_fj_per_cell == 0,
    NA_real_,
    method_summary$fake_match_avg_fj_per_cell / method_summary$rewrite_avg_fj_per_cell
  )
# ------------------ For write ----------------------
  write_energy <- calculate_write_energy_stats(
    res_static = res_static,
    feature_bit_lengths_ranked = feature_bit_lengths_ranked,
    dont_care_write_cell_pj = dont_care_write_cell_pj,
    write_down_pj = write_down_pj,
    write_up_pj = write_up_pj,
    radix_base = radix_base
  )

  method_summary <- merge(method_summary, write_energy$write_summary, by = "method", all.x = TRUE)

  numeric_cols <- setdiff(names(method_summary), "method")
  for (nm in numeric_cols) method_summary[[nm]][is.na(method_summary[[nm]])] <- 0

  list(
    method_summary = method_summary[order(method_summary$method), ],
    masked_tile_table = masked_tile_table,
    write_tile_table = write_energy$write_tile_table,
    feature_bit_table = feature_bit_table,
    parameters = list(
      base_cell_fj = base_cell_fj,
      dont_care_cell_fj = dont_care_cell_fj,
      fake_match_cell_fj = fake_match_cell_fj,
      dont_care_write_cell_pj = dont_care_write_cell_pj,
      write_down_pj = write_down_pj,
      write_up_pj = write_up_pj,
      radix_base = radix_base,
      tile_cols = tile_cols
    )
  )
}

print_energy_report <- function(energy_stats) {
  cat("\n========================================\n")
  cat("       Column Mask Energy Report        \n")
  cat("========================================\n")
  report_cols <- c(
    "method", "logical_tiles",
    "cell_oo", "cell_ox", "cell_xo", "cell_xx",
    "expanded_active_cells", "scattered_dont_care_cells", "masked_cells",
    "total_cells_with_bit_slicing",
    "baseline_total_fj", "rewrite_total_fj", "fake_match_total_fj",
    "baseline_avg_fj_per_cell", "rewrite_avg_fj_per_cell",
    "fake_match_avg_fj_per_cell", "fake_over_rewrite_avg_ratio",
    "total_write_cells", "active_rule_cells_write",
    "scattered_write_cells", "masked_write_cells",
    "baseline_write_pj", "rewrite_write_pj", "fake_match_write_pj",
    "baseline_write_avg_pj_per_cell", "rewrite_write_avg_pj_per_cell",
    "fake_match_write_avg_pj_per_cell", "fake_over_rewrite_write_ratio",
    "fake_write_reduction_ratio"
  )
  report_cols <- report_cols[report_cols %in% names(energy_stats$method_summary)]
  print(energy_stats$method_summary[, report_cols], row.names = FALSE)
  cat("----------------------------------------\n")
  invisible(energy_stats$method_summary)
}

############################################################
## 6. 主程式整合與執行 (Main Execution)
############################################################

main_simulation <- function(fit) {
  feat_names <- fit$forest$independent.variable.names
  n_features <- length(feat_names)

  cat("\n[Step 1/5] Extracting Forest Topology & Intervals...\n")
  extracted <- extract_forest_info_optimized(fit, n_features)

  feature_bit_table <- get_feature_bit_lengths(extracted, n_features, radix_base = 8L)

  cat("\n[Step 2/5] Static Mapping & Hardware Packing...\n")
  res_static <- compare_proposed_vs_naive(
    paths_list = extracted$paths,
    n_features = n_features,
    one_based = TRUE,
    path_tree_id = extracted$tree_ids,
    leaf_ids_list = extracted$leaf_ids,
    feature_bit_lengths = feature_bit_table$bit_length
  )
  show_resources(res_static)

  cat("\n[Step 3/5] Generating ACAM State Mapping Tables...\n")
  res_static$acam_state_table <- generate_acam_state_table(extracted, ord = res_static$ord)
  res_static$hardware_state_tables <- generate_hardware_mapped_state_table(
    row_tiles = res_static$row_tiles,
    global_state_table = res_static$acam_state_table
  )

  cat("\n[Step 4/5] Calculating Tile Purity...\n")
  purity_stats <- calculate_tile_purity(res_static)

  gc(verbose = FALSE)

  cat("\n[Step 5/5] Calculating Column Mask / Fake Match Energy...\n")
  energy_stats <- calculate_column_mask_energy(
    res_static = res_static,
    extracted_info = extracted,
    tile_cols = 128L,
    base_cell_fj = 13.33,
    dont_care_cell_fj = 0.25,
    fake_match_cell_fj = 10.88,
    radix_base = 8L
  )
  print_energy_report(energy_stats)

  cat("\n>>> Simulation Completed Successfully! <<<\n")
  list(static = res_static, purity = purity_stats, energy = energy_stats)
}

#---------------------------- Experiment -------------------------

run_and_save <- function(model, save_name) {
  cat(sprintf("\n========== Starting %s ==========\n", save_name))

  result <- main_simulation(model)

  saveRDS(result, file = paste0(save_name, "_read_write_avg_energy_0518.rds"))

  rm(result)
  gc(verbose = FALSE)
  cat(sprintf("========== Saved and Cleared %s ==========\n", save_name))
}

run_and_save(energy_y1_forest100, "sim_energy")
run_and_save(breast_forest100, "sim_breast")
run_and_save(heart_forest100,  "sim_heart")
run_and_save(parkinsons_motor_forest100,  "sim_parkinsons")

