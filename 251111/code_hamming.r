## ============ Helpers：頻率排序與映射（沿用） ============
# res: proposed_vs_naive_main() 的結果
# method: "proposed" 或 "naive"
# tile_id: 想看的 tile 編號；如果是 NULL，就對該 method 的所有 tile 做 summary
tile_tree_usage <- function(res, method = c("proposed", "naive"), tile_id = NULL) {
  method    <- match.arg(method)
  row_tiles <- res$row_tiles
  
  # 只取指定 method 的 tiles
  sub <- row_tiles[row_tiles$method == method, ]
  if (nrow(sub) == 0L) {
    stop("No tiles found for method = ", method)
  }
  
  # 檢查有沒有 tree_ids_in_tile 這個欄位
  if (is.null(sub$tree_ids_in_tile)) {
    stop("row_tiles 中沒有 'tree_ids_in_tile' 欄位，請確認 compare_proposed_vs_naive 有填入 path_tree_id。")
  }
  
  # 若有指定 tile_id：回傳這個 tile 的 tree 分布（table）
  if (!is.null(tile_id)) {
    idx <- which(sub$tile_id == tile_id)
    if (!length(idx)) {
      stop("No tile with tile_id = ", tile_id, " for method = ", method)
    }
    tree_ids <- sub$tree_ids_in_tile[[idx]]
    return(sort(table(tree_ids), decreasing = TRUE))
  }
  
  # 沒指定 tile_id：對這個 method 的所有 tile 做 summary，回傳 list
  out <- lapply(seq_len(nrow(sub)), function(i) {
    tree_ids <- sub$tree_ids_in_tile[[i]]
    list(
      method       = method,
      tile_id      = sub$tile_id[i],
      rows_in_tile = sub$rows_in_tile[i],
      tree_counts  = sort(table(tree_ids), decreasing = TRUE)
    )
  })
  names(out) <- paste0(method, "_tile", sub$tile_id)
  out
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
  inv <- integer(F); inv[order_idx] <- seq_len(F)
  out <- lapply(paths_list, function(v) {
    f <- if (one_based) v else (v + 1L)
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

## ============ WL/BL 四象限統計（單一視窗, 用來建 128x128 實體 tile） ============

# rows_idx: 這個 row-tile 的 row 索引（固定，不重排）
# 視窗範圍: [start, end] 的 feature index
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
  wl_open   <- sum(wl_open_vec)
  wl_closed <- n_rows - wl_open
  bl_open   <- sum(bl_used)
  bl_closed <- win_len - bl_open

  cell_oo <- wl_open  * bl_open
  cell_ox <- wl_open  * bl_closed
  cell_xo <- wl_closed* bl_open
  cell_xx <- wl_closed* bl_closed

  data.frame(
    wl_open  = wl_open,  wl_closed = wl_closed,
    bl_open  = bl_open,  bl_closed = bl_closed,
    cell_oo  = cell_oo,  cell_ox   = cell_ox,
    cell_xo  = cell_xo,  cell_xx   = cell_xx
  )
}

## ============ Proposed：最短＋（全 feature）Hamming 距離近鄰打包（只決定 row grouping） ============

pack_all_tiles_proposed <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")

  # 全 feature 集合：用於 Hamming & 全域長度
  idx_full <- lapply(paths_ranked, function(v) {
    if (!length(v)) integer(0L) else unique(as.integer(v))
  })

  # 全域 path 長度：一條 path 用到多少個 feature
  len_full <- vapply(idx_full, length, integer(1))

  used  <- rep(FALSE, length(paths_ranked))
  tiles <- list()
  tile_id <- 0L

  while (any(!used)) {
    tile_id <- tile_id + 1L
    cand <- which(!used)

    # 種子：全域最短 path（len_full）
    seed <- cand[which.min(len_full[cand])]

    rows <- seed
    used[seed] <- TRUE

    # 用「全 feature」Hamming 找最近的 127 條 path
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

    tiles[[tile_id]] <- data.frame(
      tile_id      = tile_id,
      rows_in_tile = length(rows)
    )
    tiles[[tile_id]]$rows_idx <- I(list(rows))
  }

  do.call(rbind, tiles)
}

## ============ Naive：維持原本長度排序切塊（只決定 row grouping） ============

pack_all_tiles_naive <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
  lens <- vapply(paths_ranked, length, integer(1))
  order_rows <- order(lens, decreasing = FALSE)

  tiles <- list()
  tile_id <- 0L
  for (i in seq(1, length(order_rows), by = tile_rows)) {
    tile_id <- tile_id + 1L
    rows <- order_rows[i : min(i + tile_rows - 1L, length(order_rows))]

    tiles[[tile_id]] <- data.frame(
      tile_id      = tile_id,
      rows_in_tile = length(rows)
    )
    tiles[[tile_id]]$rows_idx <- I(list(rows))
  }
  do.call(rbind, tiles)
}

## ============ 跨所有 column 視窗的「row-tile」總可關列（128 x n_features） ============

# 回傳 c(closed, total)，其中 total = Σ 每個 window 的 win_len
closed_cols_across_windows_for_tile <- function(paths_ranked, rows_idx, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L) {
    return(c(closed = 0L, total = 0L))
  }

  total_closed <- 0L
  total_cols   <- 0L

  for (start in seq.int(1L, n_features, by = tile_cols)) {
    end     <- min(start + tile_cols - 1L, n_features)
    win_len <- end - start + 1L

    total_cols <- total_cols + win_len

    used <- rep(FALSE, win_len)
    for (r in rows_idx) {
      v <- paths_ranked[[r]]
      if (!length(v)) next
      idx <- v[v >= start & v <= end]
      if (length(idx)) {
        used[idx - start + 1L] <- TRUE
      }
    }
    total_closed <- total_closed + sum(!used)
  }

  c(closed = total_closed, total = total_cols)
}

## ============ 建立「實體 128x128 tile」統計（row-tile × window） ============

build_physical_tiles <- function(paths_ranked, tiles_df, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L || nrow(tiles_df) == 0L) return(NULL)

  res <- list()
  idx <- 0L

  for (i in seq_len(nrow(tiles_df))) {
    rows_idx     <- tiles_df$rows_idx[[i]]
    logic_tile_id <- tiles_df$tile_id[i]
    method        <- tiles_df$method[i]

    window_id <- 0L
    for (start in seq.int(1L, n_features, by = tile_cols)) {
      end <- min(start + tile_cols - 1L, n_features)
      win_len <- end - start + 1L
      window_id <- window_id + 1L

      win_stats <- wlbl_stats_for_window(paths_ranked, rows_idx, start, end)

      idx <- idx + 1L
      res[[idx]] <- data.frame(
        method        = method,
        logic_tile_id = logic_tile_id,   # 哪一包 row-tile
        window_id     = window_id,       # 第幾個 column window
        start_col     = start,
        end_col       = end,
        win_len       = win_len,         # 這顆實體 tile 的 column 數（最後一個 window 可能 < 128）
        wl_open       = win_stats$wl_open,
        wl_closed     = win_stats$wl_closed,
        bl_open       = win_stats$bl_open,
        bl_closed     = win_stats$bl_closed,  # ★ 這顆實體 tile 關掉幾根 BL
        cell_oo       = win_stats$cell_oo,
        cell_ox       = win_stats$cell_ox,
        cell_xo       = win_stats$cell_xo,
        cell_xx       = win_stats$cell_xx
      )
    }
  }

  do.call(rbind, res)
}

## ============ 主比較函式：Proposed vs Naive ============

compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE,path_tree_id) 
{
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features <= 0L) {
    stop("n_features must be a positive integer (> 0).")
  }

  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  progress_bar(1, 5)

  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features
  progress_bar(2, 5)

  # 1) 決定 row grouping
  prop_df  <- pack_all_tiles_proposed(paths_ranked, tile_rows = 128L, tile_cols = 128L)
  prop_df$method <- "proposed"

  naive_df <- pack_all_tiles_naive(paths_ranked, tile_rows = 128L, tile_cols = 128L)
  naive_df$method <- "naive"
  prop_df$tree_ids_in_tile  <- lapply(prop_df$rows_idx,  function(idx) path_tree_id[idx])
  naive_df$tree_ids_in_tile <- lapply(naive_df$rows_idx, function(idx) path_tree_id[idx])


  progress_bar(3, 5)

  # 2) row-tile 視角：跨所有視窗的總可關列（128 x n_features）
  prop_closed_total <- t(vapply(
    prop_df$rows_idx,
    function(x) closed_cols_across_windows_for_tile(paths_ranked, x, tile_cols = 128L),
    integer(2L)
  ))
  colnames(prop_closed_total) <- c("closed_cols_all_windows", "total_cols_all_windows")
  prop_df <- cbind(prop_df, prop_closed_total)

  naive_closed_total <- t(vapply(
    naive_df$rows_idx,
    function(x) closed_cols_across_windows_for_tile(paths_ranked, x, tile_cols = 128L),
    integer(2L)
  ))
  colnames(naive_closed_total) <- c("closed_cols_all_windows", "total_cols_all_windows")
  naive_df <- cbind(naive_df, naive_closed_total)

  prop_df$active_cols_all_windows  <- prop_df$total_cols_all_windows  - prop_df$closed_cols_all_windows
  naive_df$active_cols_all_windows <- naive_df$total_cols_all_windows - naive_df$closed_cols_all_windows

  # 3) Hamming distance（row-tile 裡第一條 vs 最後一條 path）
  prop_df$ham_dis_in_tile <- vapply(
    seq_len(nrow(prop_df)),
    function(i) {
      rows <- prop_df$rows_idx[[i]]
      n_r  <- length(rows)
      if (n_r < 2L) {
        0L
      } else {
        a <- paths_ranked[[rows[1L]]]
        b <- paths_ranked[[rows[n_r]]]
        hamming_sets(a, b)
      }
    },
    integer(1)
  )

  naive_df$ham_dis_in_tile <- vapply(
    seq_len(nrow(naive_df)),
    function(i) {
      rows <- naive_df$rows_idx[[i]]
      n_r  <- length(rows)
      if (n_r < 2L) {
        0L
      } else {
        a <- paths_ranked[[rows[1L]]]
        b <- paths_ranked[[rows[n_r]]]
        hamming_sets(a, b)
      }
    },
    integer(1)
  )

  progress_bar(4, 5)

  # 合併 row-tile
  all_row_tiles <- rbind(prop_df, naive_df)

  # Hamming summary（以 row-tile 為單位）
  hamming_summary <- aggregate(
    ham_dis_in_tile ~ method,
    data = all_row_tiles,
    FUN = function(x) c(
      min    = base::min(x),
      median = stats::median(x),
      max    = base::max(x)
    )
  )

  # row-tile global closed stats（128 x n_features 視角）
  rowtile_closed_sum <- aggregate(
    closed_cols_all_windows ~ method,
    data = all_row_tiles,
    FUN = base::sum
  )
  rowtile_total_sum <- aggregate(
    total_cols_all_windows ~ method,
    data = all_row_tiles,
    FUN = base::sum
  )
  rowtile_global <- merge(rowtile_closed_sum, rowtile_total_sum, by = "method")
  colnames(rowtile_global)[colnames(rowtile_global) == "closed_cols_all_windows"] <- "total_closed_cols_all_windows"
  colnames(rowtile_global)[colnames(rowtile_global) == "total_cols_all_windows"]  <- "total_cols_all_windows"
  rowtile_global$closed_ratio_rowtiles <- with(
    rowtile_global,
    total_closed_cols_all_windows / total_cols_all_windows
  )

  progress_bar(5, 5)

  # 4) physical tiles（實體 128x128 視角）
  phys_prop  <- build_physical_tiles(paths_ranked, prop_df,  tile_cols = 128L)
  phys_naive <- build_physical_tiles(paths_ranked, naive_df, tile_cols = 128L)
  physical_tiles <- rbind(phys_prop, phys_naive)

  # 每 method：每顆 128x128 tile 的 BL 關閉分布
  summary_physical <- aggregate(
    bl_closed ~ method,
    data = physical_tiles,
    FUN = function(x) c(
      min    = base::min(x),
      mean   = base::mean(x),
      median = stats::median(x),
      max    = base::max(x)
    )
  )

  # 每 method：全部 128x128 tile 的總 BL / 關閉 BL & 關閉比例
  physical_global <- aggregate(
    cbind(bl_closed, win_len) ~ method,
    data = physical_tiles,
    FUN  = base::sum
  )
  colnames(physical_global)[colnames(physical_global) == "bl_closed"] <- "total_closed_bl_physical_tiles"
  colnames(physical_global)[colnames(physical_global) == "win_len"]   <- "total_bl_physical_tiles"
  physical_global$closed_ratio_physical <- with(
    physical_global,
    total_closed_bl_physical_tiles / total_bl_physical_tiles
  )

  # 每 method：cell 四象限總和（以實體 128x128 tile 為單位）
  wlbl_total <- aggregate(
    cbind(cell_oo, cell_ox, cell_xo, cell_xx) ~ method,
    data = physical_tiles,
    FUN  = base::sum
  )

  list(
    columns_order      = ord,
    row_tiles          = all_row_tiles,   # 每個 row-tile 的資訊（rows_idx, closed_cols_all_windows, ...）
    rowtile_global     = rowtile_global,  # 128 x n_features 視角的總關列/比例
    physical_tiles     = physical_tiles,  # 每一顆實體 128x128 tile 的 WL/BL/cell 統計
    summary_physical   = summary_physical,# 以實體 128x128 tile 為單位的 BL 關閉 min/mean/median/max
    physical_global    = physical_global, # 全部實體 tile 的總 BL / 關閉 BL / 比例
    wlbl_total         = wlbl_total,      # cell 四象限總和（物理 tile 累加）
    hamming_summary    = hamming_summary  # row-tile 內 Hamming 距離分布
  )
}

## ============ 你的 glue：把 fit+data 變成 paths_list 後丟進比較 ============

proposed_vs_naive_main <- function(fit, data, one_based = TRUE) {
  library(ranger)
  feat_names <- fit$forest$independent.variable.names
  n_features <- if (!is.null(feat_names) && length(feat_names) > 0L) length(feat_names) else ncol(data)
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features <= 0L) {
    stop("n_features must be a positive integer (> 0).")
  }

  paths_list <- list()
  path_tree_id <- integer(0L)

  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    paths_t <- extract_paths(info)  # 0-based splitvarID
    paths_t <- lapply(paths_t, function(p0) {
      if (is.null(p0) || length(p0) == 0L) return(integer(0L))
      p1 <- as.integer(p0) + 1L
      p1 <- p1[p1 >= 1L & p1 <= n_features]
      if (length(p1)) unique(p1) else integer(0L)
    })
    paths_list <- c(paths_list, paths_t)
    path_tree_id <- c(path_tree_id , rep.int(t, length(paths_t)))
  }

  compare_proposed_vs_naive(paths_list, n_features, one_based = one_based, path_tree_id = path_tree_id)
}



#-------------- experiment ----------------------
pvn_email <- proposed_vs_naive_main(email_forest100,email.x)
print(pvn_email$wlbl_total)
print(pvn_email$hamming_summary)
email_distri <- tile_tree_usage(pvn_email, method = "proposed")

pvn_fetal <- proposed_vs_naive_main(fetal_forest100,fetal_x)
print(pvn_fetal$wlbl_total)
print(pvn_fetal$hamming_summary)
fetal_distri <- tile_tree_usage(pvn_fetal, method = "proposed")

pvn_gene  <- proposed_vs_naive_main(gene_forest100,gene_train_x)
print(pvn_gene$wlbl_total)
print(pvn_gene$hamming_summary)
gene_distri <- tile_tree_usage(pvn_gene, method = "proposed")

pvn_mush  <- proposed_vs_naive_main(mush_forest100,mush.x)
print(pvn_mush$wlbl_total)
print(pvn_mush$hamming_summary)
mush_distri <- tile_tree_usage(pvn_mush, method = "proposed")

pvn_loan  <- proposed_vs_naive_main(loan_forest100,loan_train_.x)
print(pvn_loan$wlbl_total)
print(pvn_loan$hamming_summary)
loan_distri <- tile_tree_usage(pvn_loan, method = "proposed")

pvn_gesture <- proposed_vs_naive_main(gesture_forest100,gesture_train_x)
print(pvn_gesture$wlbl_total)
print(pvn_gesture$hamming_summary)
gesture_distri <- tile_tree_usage(pvn_gesture, method = "proposed")

#----------------------11/30 add 2 more dataset to distinguish "Many feature" ----------


pvn_arcene <- proposed_vs_naive_main(arcene_forest100,arcene_train_x)
print(pvn_arcene$wlbl_total)
print(pvn_arcene$hamming_summary)
arcene_distri <- tile_tree_usage(pvn_arcene, method = "proposed")


pvn_gas <- proposed_vs_naive_main(gas_forest100,gas_train_x)
print(pvn_gas$wlbl_total)
print(pvn_gas$hamming_summary)
gas_distri <- tile_tree_usage(pvn_gas, method = "proposed")
