############################################################
## 0. 基礎工具
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
## 2. Row Grouping (Only Sequential)
############################################################

# 這是唯一保留的分組邏輯：依序切分
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

############################################################
## 3. 核心比較邏輯 (修改為只跑 Sequential)
############################################################

compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE, path_tree_id, leaf_ids_list) {
  n_features <- as.integer(n_features)
  
  # 1. Feature Reordering (全域優化，必須保留)
  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features
  
  cat("Packing tiles (Sequential only)...\n")
  progress_bar(2, 5)
  
  # 2. Grouping: 僅執行 Sequential
  seq_df <- pack_all_tiles_sequential(paths_ranked)
  seq_df$method <- "sequential"
  
  # 附加 Meta Info
  append_meta <- function(df) {
    df$tree_ids_in_tile <- lapply(df$rows_idx, function(idx) path_tree_id[idx])
    df$leaf_ids_in_tile <- lapply(df$rows_idx, function(idx) leaf_ids_list[idx])
    df
  }
  
  # 這裡只會有 sequential 的結果
  all_row_tiles <- rbind(
    append_meta(seq_df)
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
        node_paths_t[[as.character(curr_node_id)]] <<- new_nodes  #record leaf -> full path
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
  
  # 3. 動態模擬
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
