## ============ Helpers：頻率排序與映射（沿用） ============
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
  inv <- integer(F); inv[order_idx] <- seq_len(F)
  out <- lapply(paths_list, function(v) {
    f <- if (one_based) v else (v + 1L)
    inv[f]
  })
  attr(out, "n_features") <- F
  out
}


## ============ WL/BL 四象限統計（單一視窗） ============
# rows_idx: 這個 tile 的 row 索引（固定，不重排）
# 回傳一列 data.frame：wl_open, wl_closed, bl_open, bl_closed, cell_oo, cell_ox, cell_xo, cell_xx
wlbl_stats_for_window <- function(paths_ranked, rows_idx, start, end) {
  n_rows <- length(rows_idx)
  win_len <- end - start + 1L
  if (win_len <= 0L || n_rows == 0L) {
    return(data.frame(
      wl_open=0L, wl_closed=0L, bl_open=0L, bl_closed=0L,
      cell_oo=0L, cell_ox=0L, cell_xo=0L, cell_xx=0L
    ))
  }
  # 每個 row 是否在視窗內用到至少一欄（WL 是否開）
  wl_open_vec <- logical(n_rows)
  # 視窗內每個 column 是否被至少一個 row 用到（BL 是否開）
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
  wl_open <- sum(wl_open_vec)
  wl_closed <- n_rows - wl_open
  bl_open <- sum(bl_used)
  bl_closed <- win_len - bl_open

  cell_oo <- wl_open  * bl_open
  cell_ox <- wl_open  * bl_closed
  cell_xo <- wl_closed* bl_open
  cell_xx <- wl_closed* bl_closed

  data.frame(
    wl_open=wl_open, wl_closed=wl_closed,
    bl_open=bl_open, bl_closed=bl_closed,
    cell_oo=cell_oo, cell_ox=cell_ox, cell_xo=cell_xo, cell_xx=cell_xx
  )
}

## ============ Proposed：最短＋（全 feature）Hamming 距離近鄰打包 ============
# 這裡 Hamming 距離用「整個 feature 空間」計算鄰居，不再只看前 cols_eff；
# 但 tile 的第一窗統計仍依 cols_eff = min(tile_cols, n_features)。
pack_all_tiles_proposed <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
    n_features <- attr(paths_ranked, "n_features")
    cols_eff <- min(as.integer(tile_cols), as.integer(n_features))
    
    # 視窗切片：只用來算第一個視窗的 active/closed
    idx_win <- lapply(paths_ranked, function(v) {
        if (!length(v)) integer(0L) else as.integer(v[v <= cols_eff])
    })
    # 全 feature 集合：用於 Hamming & 全域長度
    idx_full <- lapply(paths_ranked, function(v) {
        if (!length(v)) integer(0L) else unique(as.integer(v))
    })
    
    # ★ 全域 path 長度：一條 path 用到多少個 feature（不看只限前 128）
    len_full <- vapply(idx_full, length, integer(1))
    
    used  <- rep(FALSE, length(paths_ranked))
    tiles <- list()
    tile_id <- 0L
    
    while (any(!used)) {
        tile_id <- tile_id + 1L
        cand <- which(!used)
        
        # ★ 種子：全域最短 path（len_full），而不是視窗內最短
        seed <- cand[which.min(len_full[cand])]
        
        rows <- seed
        used[seed] <- TRUE
        
        # 2) 用「全 feature」Hamming 找最近的 127 條 path
        if (length(cand) > 1L) {
            a <- idx_full[[seed]]
            others <- setdiff(cand, seed)
            if (length(others)) {
                dists <- vapply(others, function(j) {
                    bj <- idx_full[[j]]
                    if (!length(a) && !length(bj)) 0L
                    else length(a) + length(bj) - 2L * length(intersect(a, bj))
                }, integer(1))
                
                take <- head(others[order(dists, decreasing = FALSE)], max(0L, tile_rows - 1L))
                if (length(take)) {
                    rows <- c(rows, take)
                    used[take] <- TRUE
                }
            }
        }
        
        # 3) 第一個視窗（1..cols_eff）的 active/closed
        cover <- rep(FALSE, cols_eff)
        for (r in rows) {
            if (length(idx_win[[r]])) cover[idx_win[[r]]] <- TRUE
        }
        active_cols <- sum(cover)
        closed_cols <- cols_eff - active_cols
        
        tiles[[tile_id]] <- data.frame(
            tile_id      = tile_id,
            rows_in_tile = length(rows),
            active_cols  = active_cols,
            closed_cols  = closed_cols
        )
        tiles[[tile_id]]$rows_idx <- I(list(rows))
    }
    
    do.call(rbind, tiles)
}
 
## ============ Naive：維持原本長度排序切塊，補上四象限統計 ============
pack_all_tiles_naive <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
  lens <- sapply(paths_ranked, length)
  order_rows <- order(lens, decreasing = FALSE)

  n_features <- attr(paths_ranked, "n_features")
  cols_eff <- min(as.integer(tile_cols), as.integer(n_features))

  tiles <- list()
  tile_id <- 0L
  for (i in seq(1, length(order_rows), by = tile_rows)) {
    tile_id <- tile_id + 1L
    rows <- order_rows[i : min(i + tile_rows - 1L, length(order_rows))]

    # 第一個視窗 active/closed
    cover <- rep(FALSE, cols_eff)
    for (r in rows) {
      v <- paths_ranked[[r]]
      if (length(v)) {
        idx <- v[v <= cols_eff]
        if (length(idx)) cover[idx] <- TRUE
      }
    }
    active_cols <- sum(cover)
    closed_cols <- cols_eff - active_cols


    tiles[[tile_id]] <- data.frame(
      tile_id      = tile_id,
      rows_in_tile = length(rows),
      active_cols  = active_cols,
      closed_cols  = closed_cols
    )
    tiles[[tile_id]]$rows_idx <- I(list(rows))
  }
  do.call(rbind, tiles)
}

## ============ 跨所有 column 視窗的總可關列 & 四象限總和（不重排 row） ============
# 仍沿用你之前的「跨視窗」統計，只是補上四象限的逐窗累加
closed_cols_across_windows_for_tile <- function(paths_ranked, rows_idx, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L) return(0L)
  total_closed <- 0L
  for (start in seq.int(1L, n_features, by = tile_cols)) {
    end <- min(start + tile_cols - 1L, n_features)
    win_len <- end - start + 1L
    used <- rep(FALSE, win_len)
    for (r in rows_idx) {
      v <- paths_ranked[[r]]
      if (!length(v)) next
      idx <- v[v >= start & v <= end]
      if (length(idx)) used[idx - start + 1L] <- TRUE
    }
    total_closed <- total_closed + sum(!used)
  }
  total_closed
}

# 四象限的跨視窗總和
wlbl_across_windows_for_tile <- function(paths_ranked, rows_idx, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L) {
    return(c(cell_oo=0L, cell_ox=0L, cell_xo=0L, cell_xx=0L))
  }
  acc <- c(cell_oo=0L, cell_ox=0L, cell_xo=0L, cell_xx=0L)
  for (start in seq.int(1L, n_features, by = tile_cols)) {
    end <- min(start + tile_cols - 1L, n_features)
    wlbl <- wlbl_stats_for_window(paths_ranked, rows_idx, start, end)
    acc <- acc + c(cell_oo=wlbl$cell_oo, cell_ox=wlbl$cell_ox,
                   cell_xo=wlbl$cell_xo, cell_xx=wlbl$cell_xx)
  }
  acc
}
compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE) {
    n_features <- as.integer(n_features)
    if (is.na(n_features) || n_features <= 0L) {
        stop("n_features must be a positive integer (> 0).")
    }
    
    ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
    progress_bar(1,5)
    
    paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
    attr(paths_ranked, "n_features") <- n_features
    progress_bar(2,5)
    
    # Proposed（哈明距離版，現在只決定每個 tile 的 rows）
    prop_df  <- pack_all_tiles_proposed(paths_ranked, tile_rows = 128L, tile_cols = 128L)
    prop_df$method <- "proposed"
    
    # Naive（按 path length 切塊）
    naive_df <- pack_all_tiles_naive(paths_ranked, tile_rows = 128L, tile_cols = 128L)
    naive_df$method <- "naive"
    
    progress_bar(3,5)
    
    
    # ★ 跨所有視窗：總可關列（真正你要比較的量）
    prop_df$closed_cols_all_windows  <- vapply(
        prop_df$rows_idx,
        function(x) closed_cols_across_windows_for_tile(paths_ranked, x[[1]], tile_cols = 128L),
        integer(1)
    )
    naive_df$closed_cols_all_windows <- vapply(
        naive_df$rows_idx,
        function(x) closed_cols_across_windows_for_tile(paths_ranked, x[[1]], tile_cols = 128L),
        integer(1)
    )
    
    progress_bar(4,5)
    
    # 合併所有 tile 資訊
    all <- rbind(prop_df, naive_df)
    
    # ★ 每種 method 的「全局總可關列」：你最關心的指標
    total_sum <- aggregate(
        closed_cols_all_windows ~ method,
        data = all,
        FUN = base::sum
    )
    colnames(total_sum)[2] <- "total_closed_cols_all_windows"
    
    # 順手給你每 method 的 tile-level closed_cols_all_windows 統計（平均、median、p90）
    summary_global <- aggregate(
        closed_cols_all_windows ~ method,
        data = all,
        FUN = function(x) c(
            mean   = base::mean(x),
            median = stats::median(x),
            p90    = stats::quantile(x, 0.9, names = FALSE)
        )
    )
    
    progress_bar(5,5)
    return(list(
        columns_order = ord,
        tiles   = all,                          # 每個 tile 的 rows_idx + closed_cols_all_windows
        total_closed_cols_all_windows = total_sum,
        summary_global = summary_global         # 原本叫 summary_first_window，現在變成全局關列的 summary
    ))
}

## ============ 你的 glue：把 fit+data 變成 paths_list 後丟進比較 ============
proposed_vs_naive_main <- function(fit, data, one_based = TRUE) {
  library(ranger)
  feat_names <- fit$forest$independent.variable.names
  n_features <- if (!is.null(feat_names) && length(feat_names) > 0L) length(feat_names) else ncol(data)
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features <= 0L) stop("n_features must be a positive integer (> 0).")

  paths_list <- list()
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    paths_t <- extract_paths(info)  # 0-based splitvarID
    paths_t <- lapply(paths_t, function(p0) {
      if (is.null(p0) || length(p0) == 0L) return(integer(0))
      p1 <- as.integer(p0) + 1L
      p1 <- p1[p1 >= 1L & p1 <= n_features]
      if (length(p1)) unique(p1) else integer(0)
    })
    paths_list <- c(paths_list, paths_t)
  }

  compare_proposed_vs_naive(paths_list, n_features, one_based = TRUE)
}

