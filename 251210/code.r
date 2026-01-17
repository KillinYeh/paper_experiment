############################################################
## 0. 基礎工具 (維持不變)
############################################################

# 頻率排序
feature_frequency_order <- function(paths_list, n_features, one_based = TRUE) {
  freq <- integer(n_features)
  for (p in paths_list) {
    f <- if (one_based) p else (p + 1L)
    if (length(f)) freq[f] <- freq[f] + 1L
  }
  order(freq, decreasing = TRUE)
}

# Remap Feature Index -> Global Rank
remap_paths_to_rank <- function(paths_list, order_idx, one_based = TRUE) {
  F <- length(order_idx)
  inv <- integer(F)
  inv[order_idx] <- seq_len(F) # feature ID -> rank

  out <- lapply(paths_list, function(v) {
    f <- if (one_based) v else (v + 1L)
    if (!length(f)) return(integer(0L))
    inv[f]
  })
  attr(out, "n_features") <- F
  out
}

hamming_sets <- function(a, b) {
  if ((!length(a)) && (!length(b))) return(0L)
  inter_len <- length(intersect(a, b))
  length(a) + length(b) - 2L * inter_len
}

############################################################
## 1. 靜態 WL/BL 統計 (維持不變)
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

closed_cols_across_windows_for_tile <- function(paths_ranked, rows_idx, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L) return(c(closed = 0L, total = 0L))
  
  total_closed <- 0L
  total_cols   <- 0L
  
  for (start in seq.int(1L, n_features, by = tile_cols)) {
    end <- min(start + tile_cols - 1L, n_features)
    win_len <- end - start + 1L
    total_cols <- total_cols + win_len
    
    used <- rep(FALSE, win_len)
    for (r in rows_idx) {
      v <- paths_ranked[[r]]
      if (!length(v)) next
      idx <- v[v >= start & v <= end]
      if (length(idx)) used[idx - start + 1L] <- TRUE
    }
    total_closed <- total_closed + sum(!used)
  }
  c(closed = total_closed, total = total_cols)
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
## 2. Row Grouping (Proposed / Naive)
############################################################

pack_all_tiles_proposed <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
  # 簡單實作：以全域長度當種子，用 Hamming 找最近鄰
  idx_full <- lapply(paths_ranked, function(v) if (!length(v)) integer(0L) else unique(as.integer(v)))
  len_full <- vapply(idx_full, length, integer(1))
  
  used  <- rep(FALSE, length(paths_ranked))
  tiles <- list()
  tile_id <- 0L
  
  while (any(!used)) {
    tile_id <- tile_id + 1L
    cand <- which(!used)
    seed <- cand[which.min(len_full[cand])]
    
    rows <- seed
    used[seed] <- TRUE
    
    if (length(cand) > 1L) {
      a <- idx_full[[seed]]
      others <- setdiff(cand, seed)
      if (length(others)) {
        dists <- vapply(others, function(j) {
          bj <- idx_full[[j]]
          hamming_sets(a, bj)
        }, integer(1))
        
        take <- head(others[order(dists, decreasing = FALSE)], max(0L, tile_rows - 1L))
        if (length(take)) {
          rows <- c(rows, take)
          used[take] <- TRUE
        }
      }
    }
    tiles[[tile_id]] <- data.frame(tile_id = tile_id, rows_in_tile = length(rows))
    tiles[[tile_id]]$rows_idx <- I(list(rows))
  }
  do.call(rbind, tiles)
}

pack_all_tiles_naive <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
  lens <- vapply(paths_ranked, length, integer(1))
  order_rows <- order(lens, decreasing = FALSE)
  
  tiles <- list()
  tile_id <- 0L
  for (i in seq(1, length(order_rows), by = tile_rows)) {
    tile_id <- tile_id + 1L
    rows <- order_rows[i : min(i + tile_rows - 1L, length(order_rows))]
    tiles[[tile_id]] <- data.frame(tile_id = tile_id, rows_in_tile = length(rows))
    tiles[[tile_id]]$rows_idx <- I(list(rows))
 
  }
  do.call(rbind, tiles)
}

############################################################
## 3. 核心比較邏輯 (靜態部分)
############################################################

compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE, path_tree_id, leaf_ids_list) {
  # leaf_ids_list: 與 paths_list 對應的 leaf node ID，用於後續動態模擬
  
  n_features <- as.integer(n_features)
  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features
  
  cat("using proposed/naive to pack the path\n")
  progress_bar(2,5)
  # Grouping
  prop_df  <- pack_all_tiles_proposed(paths_ranked)
  prop_df$method <- "proposed"
  naive_df <- pack_all_tiles_naive(paths_ranked)
  naive_df$method <- "naive"
  
  # 紀錄每個 array 包含的 tree 和 leaf 資訊 (為了動態模擬查表)
  append_meta <- function(df) {
    df$tree_ids_in_tile <- lapply(df$rows_idx, function(idx) path_tree_id[idx])
    df$leaf_ids_in_tile <- lapply(df$rows_idx, function(idx) leaf_ids_list[idx])
    df
  }
  prop_df  <- append_meta(prop_df)
  naive_df <- append_meta(naive_df)
  
  # 合併
  all_row_tiles <- rbind(prop_df, naive_df)
  cat("build physical tile\n")
  progress_bar(3,5)
  # 預先計算物理 tile 靜態結構 (BL)
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
## 4. 關鍵修正：透過 LCA Map 計算 Mismatch
############################################################

# 4.1 建立 Tree 結構的 LCA Map
# 回傳一個 list，key 為 tree ID，value 為 matrix [leaf_id_1, leaf_id_2] -> split_feature
build_tree_lca_maps <- function(fit, n_features) {
  lca_maps <- list()
  
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    leaves <- info$nodeID[info$terminal]
    
    # 用 path 記錄每個 leaf 的祖先 split feature
    # 為了快，我們先建立 node -> parent & split feature 的 lookup
    # 但 ranger treeInfo 由上而下比較好找
    
    # 建立 node -> path_features (包含 splitvarID)
    # 我們需要知道兩個 leaf 的分歧點是哪個 feature
    # 方法：DFS 遍歷，紀錄每個 leaf 的完整路徑 (nodes list)
    
    leaf_paths <- list() # key = leaf_id, val = vector of nodeIDs from root
    node_split <- integer(max(info$nodeID) + 1L) # nodeID -> splitVarID
    names(node_split) <- 0:max(info$nodeID)
    
    # 填 node_split
    non_term <- info[!info$terminal, ]
    node_split[as.character(non_term$nodeID)] <- non_term$splitvarID
    
    dfs <- function(nid, current_path) {
      row_idx <- which(info$nodeID == nid)
      if (info$terminal[row_idx]) {
        leaf_paths[[as.character(nid)]] <<- current_path
      } else {
        left <- info$leftChild[row_idx]
        right <- info$rightChild[row_idx]
        new_path <- c(current_path, nid)
        dfs(left, new_path)
        dfs(right, new_path)
      }
    }
    dfs(0, integer(0L)) # root is 0
    
    # 建立 Matrix
    # row/col 名稱是 leafID
    n_leaves <- length(leaves)
    mat <- matrix(NA_integer_, nrow=n_leaves, ncol=n_leaves)
    rownames(mat) <- colnames(mat) <- as.character(leaves)
    
    # 兩兩比較
    # 這是 O(L^2 * Depth)，對於 ranger 預設 leaf 數量通常可接受
    for (i in 1:n_leaves) {
      lid_i <- as.character(leaves[i])
      path_i <- leaf_paths[[lid_i]]
      
      mat[i, i] <- -1L # 自己跟自己沒有 mismatch
      
      if (i < n_leaves) {
        for (j in (i+1):n_leaves) {
          lid_j <- as.character(leaves[j])
          path_j <- leaf_paths[[lid_j]]
          
          # 找分歧點：路徑上最後一個相同的 node
          # 路徑是 root -> ... -> parent_of_leaf
          len <- min(length(path_i), length(path_j))
          same <- (path_i[1:len] == path_j[1:len])
          last_common_idx <- which.min(same) - 1L 
          # 如果全部 same (一長一短)，則 last_common_idx = len
          if (all(same)) last_common_idx <- len
          if (length(same) == 0) last_common_idx <- 0 # 其實不會發生因為都有 root
          
          # 分歧點 nodeID
          lca_node <- path_i[last_common_idx] 
          
          # 該 node 的 split feature 即為 mismatch feature
          # ranger splitvarID 是 0-based
          feat_id <- node_split[as.character(lca_node)] + 1L 
          
	  #mat[i,j]:path_i與path_j在feat_id的地方發生分歧
          mat[i, j] <- feat_id
          mat[j, i] <- feat_id
        }
      }
    }
    lca_maps[[t]] <- mat
  }
  lca_maps
}

# 4.2 預先計算 BL 結構 (加速模擬)
precompute_bl_info <- function(physical_tiles) {
  # 轉成 list of list 結構方便存取: bl_lookup[[method]][[array_id]][[tile_in_array]]
  # 這裡簡化：直接對 physical_tiles 建立 key
  physical_tiles$key <- paste(physical_tiles$method, physical_tiles$array_id, physical_tiles$tile_in_array, sep="_")
  
  # 只需要存 bl_open, bl_closed, win_len
  bl_vec <- physical_tiles[, c("bl_open", "bl_closed", "win_len")]
  row.names(bl_vec) <- physical_tiles$key
  bl_vec
}

############################################################
## 4.3 核心模擬函式 (修正版：支援 Mixed Trees in Tile)
############################################################

simulate_cam_array_mismatch <- function(fit, data_test, res, tile_cols = 128L) {
  row_tiles    <- res$row_tiles
  paths_ranked <- res$paths_ranked
  
  # Feature Rank Lookup (Feature ID -> Global Rank)
  # 用來將 LCA 找到的 Feature ID 轉換成 CAM 上的 Column Index
  ord <- res$columns_order 
  feat_rank_lookup <- integer(res$n_features)
  feat_rank_lookup[ord] <- seq_along(ord)
  
  # 1. 預測所有 Sample 在所有 Tree 的 True Leaf
  # pred_leaf: [n_samples x n_trees] 矩陣
  # 這代表了 Ground Truth
  cat("Predicting terminal nodes for all samples...\n")
  pred_leaf <- predict(fit, data = data_test, type = "terminalNodes")$predictions
  n_samples <- nrow(pred_leaf)
  
  # 2. 建立所有 Tree 的 LCA Map
  # lca_maps[[tree_id]] 是一個矩陣，存儲該樹任意兩個 leaf 之間的分歧 feature
  cat("Building LCA Maps for all trees...\n")
  progress_bar(4,5)
  lca_maps <- build_tree_lca_maps(fit, res$n_features)
  
  # 3. 準備 BL 結構 lookup (靜態結構)
  bl_lookup <- precompute_bl_info(res$physical_tiles)
  
  methods <- unique(row_tiles$method)
  out <- list()
  
  for (method in methods) {
    cat("Simulating method:", method, "\n")
    sub <- row_tiles[row_tiles$method == method, ]
    
    # 累加器 (Cell 開啟次數)
    total_oo <- 0; total_ox <- 0; total_xo <- 0; total_xx <- 0
    
    # === 遍歷每一個 Array (Row grouping) ===
    # 每個 Array 可能包含來自不同 Tree 的 Path
    for (i in seq_len(nrow(sub))) {
      
      # 取得此 Array 的基本資訊
      rows_idx <- sub$rows_idx[[i]]            # 這些 Row 對應原本 paths_list 的 index
      tree_ids <- unlist(sub$tree_ids_in_tile[[i]]) # 每一條 Row 來自哪棵樹 (Vector)
      leaf_ids <- unlist(sub$leaf_ids_in_tile[[i]]) # 每一條 Row 代表哪個 Target Leaf (Vector)
      array_id <- sub$tile_id[i]
      
      n_rows_in_array <- length(rows_idx)
      
      # 找出這個 Array 橫向切成了幾顆 Tile
      pt_sub <- res$physical_tiles[res$physical_tiles$method == method & 
                                     res$physical_tiles$array_id == array_id, ]
      n_tiles_wide <- nrow(pt_sub)
      
      # === 步驟 A: 計算 Mismatch Tile Index (針對所有 Sample, 所有 Row) ===
      # mat_mis_tile: [Sample, Row] -> 該 Row 在第幾顆 Tile 發生 Mismatch (之後可關閉)
      # 預設為 Inf (表示完全 Match 或 Mismatch 發生在最後之後)
      mat_mis_tile <- matrix(Inf, nrow=n_samples, ncol=n_rows_in_array)
      
      # 為了效率，我們這裡逐「Row」處理 (因為 Row 數通常較少，約 128)
      # 這樣可以利用 R 的向量化運算處理所有 Sample
      for (r in seq_len(n_rows_in_array)) {
        
        # 1. 鎖定當前這條 Row 的身份
        t_id <- tree_ids[r]          # 這條 Row 屬於第幾棵樹
        l_tgt <- as.character(leaf_ids[r]) # 這條 Row 的目標 Leaf ID (轉字串查表用)
        
        # 2. 找出所有 Sample 在「這棵樹」的 True Leaf
        # 這裡就是解決您疑問的關鍵：我們只看與這條 Row 相關的那棵樹的預測結果
        l_true_vec <- as.character(pred_leaf[, t_id])
        
        # 3. 找出 Mismatch 的 Sample (True Leaf != Target Leaf)
        # 如果 True Leaf == Target Leaf，代表完全命中，Mismatch Tile = Inf (預設值)，不用改
        diff_indices <- which(l_true_vec != l_tgt)
        
        if (length(diff_indices) > 0) {
          # 4. 取出這些 Mismatch Sample 的 True Leaf ID
          mismatch_true_leaves <- l_true_vec[diff_indices]
          
          # 5. 查這棵樹的 LCA Map
          # map[row, col] -> split feature ID
          # 我們一次查多個 (Vectorized Lookup)
          # 注意：lca_maps[[t_id]] 的 row/col names 是 leaf ID 字串
          map <- lca_maps[[t_id]]
          
          # 透過 matrix indexing 快速查找: map[cbind(row_names, col_names)]
          mis_feats <- map[cbind(mismatch_true_leaves, l_tgt)]
          
          # 6. Feature ID -> Global Rank -> Tile Index
          # 這是所有 Tree 共用的 Global Feature Mapping(origin feature_id -> reorder feature_id)
          ranks <- feat_rank_lookup[mis_feats]
          
          # 計算在哪一顆 Tile 發生 Mismatch (1-based index)
          # 例如: Rank 1~128 -> Tile 1, Rank 129~256 -> Tile 2
          tiles_idx <- (ranks - 1L) %/% tile_cols + 1L
          
          # 7. 填回矩陣
	  #mat_mis_tile[mis_sample,path_id]:所有在path_id為r的samples他的mismatch發生在哪一個tile
          mat_mis_tile[diff_indices, r] <- tiles_idx
        }
      }
      
      # === 步驟 B: 根據 Mismatch Tile 統計 Cell 狀態 ===
      # 遍歷橫向的每一顆 Tile (k = 1..N)
      for (k in seq_len(n_tiles_wide)) 
      {
        # 取得這顆 Tile 的靜態 BL 結構 (bl_open, bl_closed)
        key <- paste(method, array_id, k, sep="_")
        bl_info <- bl_lookup[key, ] 
        b_open  <- bl_info[["bl_open"]]
        b_closed <- bl_info[["bl_closed"]]
        
        # 判斷 WL 是否開啟 (WL Open Status)
        # 對於某個 Sample 的某條 Row：
        # 如果 Mismatch 發生在 Tile k 之後 (mis_tile > k)，WL 必須開著。
        # 如果 Mismatch 發生在 Tile k 本身 (mis_tile == k)，WL 也必須開著 (才能偵測到 Mismatch)。
        # 如果 Mismatch 發生在 Tile k 之前 (mis_tile < k)，WL 可以關閉 (Selective Precharge)。
        
        # [Sample, Row] 的"邏輯"矩陣 , return T/F
        wl_is_active_mat <- (mat_mis_tile >= k) 
        
        # 統計:
        # wl_open_counts[s] = 第 s 個 sample 在這顆 Tile 開啟了幾條 Row , rowSums = sum over a row.
        wl_open_counts <- rowSums(wl_is_active_mat) 
        
        # wl_closed_counts[s] = 第 s 個 sample 在這顆 Tile 關閉了幾條 Row
        wl_closed_counts <- n_rows_in_array - wl_open_counts
        
        # 累加到總數 (Sum over all samples)
        sum_wl_open <- sum(wl_open_counts)
        sum_wl_closed <- sum(wl_closed_counts)
        
        # 計算 Cell 狀態 (四象限)
        # Cell OO: WL開 & BL開
        total_oo <- total_oo + sum_wl_open * b_open
        # Cell OX: WL開 & BL關 (Match 不需要耗電，但 WL 還是充了電)
        total_ox <- total_ox + sum_wl_open * b_closed
        # Cell XO: WL關 & BL開 (BL 充了電，但 Row 沒開)
        total_xo <- total_xo + sum_wl_closed * b_open
        # Cell XX: WL關 & BL關 (最省電狀態)
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

extract_paths_and_leaves_ranger <- function(fit, n_features) {
  paths_list   <- list()
  path_tree_id <- integer(0L)
  leaf_ids_list <- integer(0L)
  
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    tmp  <- extract_paths_with_leaf(info) # Helper from user code or above
    
    paths_t <- lapply(tmp$paths, function(p1) {
      if (is.null(p1) || !length(p1)) return(integer(0L))
      p1 <- as.integer(p1) + 1L # ranger is 0-based
      p1 <- p1[p1 >= 1L & p1 <= n_features]
      if (length(p1)) unique(p1) else integer(0L)
    })
    
    paths_list   <- c(paths_list, paths_t)
    path_tree_id <- c(path_tree_id, rep.int(t, length(paths_t)))
    leaf_ids_list <- c(leaf_ids_list, tmp$leaf_ids)
  }
  list(paths=paths_list, tree_ids=path_tree_id, leaf_ids=leaf_ids_list)
}

# 輔助：從 treeInfo 遞迴抓路徑 (User 原始程式碼的變體)
extract_paths_with_leaf <- function(tree_info) {
  paths    <- list()
  leaf_ids <- integer(0L)
  
  # Ranger treeInfo format check
  if (!"splitvarID" %in% names(tree_info)) return(list(paths=list(), leaf_ids=integer(0)))
  
  # Build simple adjacency for speed
  left_child <- tree_info$leftChild + 1 # 1-based index for R
  right_child <- tree_info$rightChild + 1
  split_var <- tree_info$splitvarID
  is_term <- tree_info$terminal
  node_ids <- tree_info$nodeID
  
  dfs <- function(idx, current_feats) {
    if (is_term[idx]) {
      paths[[length(paths) + 1L]] <<- current_feats
      leaf_ids <<- c(leaf_ids, node_ids[idx])
    } else {
      # Current node split feature
      feat <- split_var[idx]
      # Go left
      dfs(left_child[idx], c(current_feats, feat))
      # Go right
      dfs(right_child[idx], c(current_feats, feat))
    }
  }
  
  if (nrow(tree_info) > 0) dfs(1, integer(0L))
  
  list(paths = paths, leaf_ids = leaf_ids)
}

main_simulation <- function(fit,  data_test) {
  library(ranger)
  
  # 1. 解析 Features
  feat_names <- fit$forest$independent.variable.names
  n_features <- length(feat_names)
  
  # 2. 提取路徑
  cat("Extracting paths...\n")
  progress_bar(1,5)
  extracted <- extract_paths_and_leaves_ranger(fit, n_features)
  
  # 3. 靜態比較與 Tile 建構
  cat("Running static mapping (Proposed vs Naive)...\n")
  res_static <- compare_proposed_vs_naive(
    paths_list = extracted$paths,
    n_features = n_features,
    one_based = TRUE,
    path_tree_id = extracted$tree_ids,
    leaf_ids_list = extracted$leaf_ids
  )
  
  # 4. 動態功耗模擬
  cat("Running dynamic simulation on test data...\n")
  res_dynamic <- simulate_cam_array_mismatch(fit, data_test, res_static, tile_cols = 128L)
  cat("finish the simulation\n")
  progress_bar(5,5)
  list(static = res_static, dynamic = res_dynamic)
}



#-------------------------experiment----------------------
sim_email <- main_simulation(email_forest100,email_test_x)

sim_fetal <- main_simulation(fetal_forest100,fetal_x)

sim_gene  <- main_simulation(gene_forest100,gene_test_x)

sim_mush  <- main_simulation(mush_forest100,mush_test_x)

sim_loan  <- main_simulation(loan_forest100,loan_test_x)

sim_gesture <- main_simulation(gesture_forest100,gesture_test_x)

sim_arcene <- main_simulation(arcene_forest100,arcene_test_x)

sim_gas <- main_simulation(gas_forest100,gas_test_x)

