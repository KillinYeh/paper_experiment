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
## 2.2 Row Grouping (Real FP-Tree Construction + Union Window)
############################################################

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
## 3. 核心比較邏輯 
############################################################

compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE, path_tree_id, leaf_ids_list) {
  n_features <- as.integer(n_features)
  
  # 1. Feature Reordering (全域優化，必須保留)
  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features
  
  cat("Packing tiles (FP-Tree)...\n")
  progress_bar(2, 5)
  cat("\n[Method 1] Real FP-Tree + Union (Window=512, Dropout > 128)...\n")
  fptree_union_df <- pack_all_tiles_fptree_union(paths_ranked, tile_rows = 128L, window_size = 512L)
  fptree_union_df$method <- "fptree_union"


  # 附加 Meta Info
  append_meta <- function(df) {
    df$tree_ids_in_tile <- lapply(df$rows_idx, function(idx) path_tree_id[idx])
    df$leaf_ids_in_tile <- lapply(df$rows_idx, function(idx) leaf_ids_list[idx])
    df
  }
  
  # 這裡只會有 sequential 的結果
  all_row_tiles <- rbind(
	append_meta(fptree_union_df)	
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

# 優化：直接利用 extract 階段產生的 leaf_info，不再重複 DFS
build_tree_lca_maps_optimized <- function(fit, leaf_info_list) {
  lca_maps <- list()
  
  for (t in seq_along(leaf_info_list)) {
    info <- leaf_info_list[[t]]
    leaves <- info$leaf_nodes
    leaf_paths <- info$node_paths
    node_split_map <- info$node_split_map
    
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
            lca_node <- path_i[last_common_idx]
            feat_id <- node_split_map[[as.character(lca_node)]] + 1L 
            mat[i, j] <- feat_id
            mat[j, i] <- feat_id
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



simulate_cam_array_mismatch <- function(fit, data_test, res, leaf_info_list, tile_cols = 128L) {
  row_tiles    <- res$row_tiles
  paths_ranked <- res$paths_ranked
  
  ord <- res$columns_order
  feat_rank_lookup <- integer(res$n_features)
  feat_rank_lookup[ord] <- seq_along(ord)
  
  cat("Predicting terminal nodes for all samples...\n")
  pred_leaf <- predict(fit, data = data_test, type = "terminalNodes")$predictions
  n_samples <- nrow(pred_leaf)
  
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
          mis_feats <- map[cbind(mismatch_true_leaves, l_tgt)]
          
          ranks <- feat_rank_lookup[mis_feats]
          tiles_idx <- (ranks - 1L) %/% tile_cols + 1L
          mat_mis_tile[diff_indices, r] <- tiles_idx
        }
      }
      
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

sim_cov <- main_simulation(cov_forest_100,cov_test_x)
save.image(file = "backup_fptree_radix_20260310.RData")
