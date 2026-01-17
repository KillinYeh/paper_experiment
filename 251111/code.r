## ============ Helpers: 資料結構與排序（同前） ============
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

## ============ 視窗內增量（用 cover 長度當視窗大小）（同前） ============

inc_in_front128 <- function(path_rank, cover128) {
  if (!length(path_rank) || !length(cover128)) return(0L)
  win <- length(cover128)
  idx <- path_rank[path_rank <= win]
  if (!length(idx)) return(0L)
  sum(!cover128[idx])
}

## ============ Proposed：保留演算法，但 tiles 加上 rows_idx ============

pack_one_tile_proposed <- function(paths_ranked, used, tile_rows = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  cols_eff <- min(128L, as.integer(n_features))

  cand <- which(!used)
  if (!length(cand)) return(list(rows = integer(0), active_cols = 0L))

  cover128 <- rep(FALSE, cols_eff)

  used_in_win <- sapply(cand, function(i) {
    v <- paths_ranked[[i]]
    if (!length(v)) return(0L)
    sum(v <= cols_eff)
  })
  seed <- cand[which.min(used_in_win)]

  selected <- seed
  used[seed] <- TRUE
  if (length(paths_ranked[[seed]])) {
    idx <- paths_ranked[[seed]]
    idx <- idx[idx <= cols_eff]
    if (length(idx)) cover128[idx] <- TRUE
  }

  eat_zero_inc <- function() {
    if (length(selected) >= tile_rows) return(invisible(0L))
    pool <- which(!used)
    if (!length(pool)) return(invisible(0L))
    incs <- sapply(pool, function(i) inc_in_front128(paths_ranked[[i]], cover128))
    zeros <- pool[incs == 0L]
    if (!length(zeros)) return(invisible(0L))
    room <- tile_rows - length(selected)
    zeros <- head(zeros, room)
    used[zeros] <<- TRUE
    selected <<- c(selected, zeros)
    invisible(length(zeros))
  }
  eat_zero_inc()

  while (length(selected) < tile_rows && any(!used)) {
    pool <- which(!used)
    incs <- sapply(pool, function(i) inc_in_front128(paths_ranked[[i]], cover128))
    min_inc <- min(incs)
    cand_min <- pool[incs == min_inc]

    if (length(cand_min) > 1L) {
      used_cnt <- sapply(cand_min, function(i) {
        v <- paths_ranked[[i]]
        if (!length(v)) return(0L)
        sum(v <= cols_eff)
      })
      pick <- cand_min[order(used_cnt, decreasing = FALSE)][1]
    } else {
      pick <- cand_min[1]
    }

    used[pick] <- TRUE
    selected <- c(selected, pick)

    v <- paths_ranked[[pick]]
    if (length(v)) {
      idx <- v[v <= cols_eff]
      if (length(idx)) cover128[idx] <- TRUE
    }

    eat_zero_inc()
    if (!any(!used)) break
  }

  active_cols <- sum(cover128)
  list(rows = selected, active_cols = active_cols)
}

pack_all_tiles_proposed <- function(paths_ranked, tile_rows = 128L, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  cols_eff <- min(as.integer(tile_cols), as.integer(n_features))

  used <- rep(FALSE, length(paths_ranked))
  tiles <- list()
  tile_id <- 0L

  while (any(!used)) {
    tile_id <- tile_id + 1L
    res <- pack_one_tile_proposed(paths_ranked, used, tile_rows = tile_rows)
    rows <- res$rows
    if (!length(rows)) break
    used[rows] <- TRUE

    active_cols <- res$active_cols
    closed_cols <- cols_eff - active_cols

    tiles[[tile_id]] <- data.frame(
      tile_id      = tile_id,
      rows_in_tile = length(rows),
      active_cols  = active_cols,
      closed_cols  = closed_cols
    )
    # 重要：存 row 索引，後面跨視窗計算要用（list-column）
    tiles[[tile_id]]$rows_idx <- I(list(rows))
  }
  do.call(rbind, tiles)
}

## ============ Naive：一樣加上 rows_idx ============

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

    cover128 <- rep(FALSE, cols_eff)
    for (r in rows) {
      idx <- paths_ranked[[r]]
      if (length(idx)) {
        idx <- idx[idx <= cols_eff]
        if (length(idx)) cover128[idx] <- TRUE
      }
    }

    active_cols <- sum(cover128)
    closed_cols <- cols_eff - active_cols

    tiles[[tile_id]] <- data.frame(
      tile_id = tile_id,
      rows_in_tile = length(rows),
      active_cols = active_cols,
      closed_cols = closed_cols
    )
    tiles[[tile_id]]$rows_idx <- I(list(rows))
  }
  do.call(rbind, tiles)
}

## ============ 新增：計算「跨所有 column 視窗」的可關欄數（不重排 row） ============

# 對一個 tile（其 rows_idx 已固定），在所有 column 視窗（1..n_features 以 128 等距分窗）上
# 計算可關欄數總和。條件：某視窗的某欄位，若該視窗內「所有這個 tile 的 rows」都沒用到，才視為可關。
closed_cols_across_windows_for_tile <- function(paths_ranked, rows_idx, tile_cols = 128L) {
  n_features <- attr(paths_ranked, "n_features")
  if (n_features <= 0L) return(0L)

  total_closed <- 0L
  # 每個視窗的左右界
  for (start in seq.int(1L, n_features, by = tile_cols)) {
    end <- min(start + tile_cols - 1L, n_features)
    win_len <- end - start + 1L
    # 這個視窗內，先假設全部可關（全 FALSE = 未使用）
    used <- rep(FALSE, win_len)
    # 把 tile 的 rows 在此視窗用到的 rank 打上 TRUE
    for (r in rows_idx) {
      v <- paths_ranked[[r]]
      if (!length(v)) next
      idx <- v[v >= start & v <= end]
      if (length(idx)) used[idx - start + 1L] <- TRUE
    }
    # 沒被使用的列才可關
    total_closed <- total_closed + sum(!used)
  }
  total_closed
}

## ============ 一鍵比較：加入跨視窗總可關欄數（不重排 row） ============

compare_proposed_vs_naive <- function(paths_list, n_features, one_based = TRUE) {
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features <= 0L) {
    stop("n_features must be a positive integer (> 0).")
  }

  ord <- feature_frequency_order(paths_list, n_features, one_based = one_based)
  paths_ranked <- remap_paths_to_rank(paths_list, ord, one_based = one_based)
  attr(paths_ranked, "n_features") <- n_features

  # 1) 先完成第一輪（第一個 column 視窗）的 tiles 結果
  prop_df  <- pack_all_tiles_proposed(paths_ranked, tile_rows = 128L, tile_cols = 128L)
  prop_df$method <- "proposed"

  naive_df <- pack_all_tiles_naive(paths_ranked, tile_rows = 128L, tile_cols = 128L)
  naive_df$method <- "naive"

  # 2) 在不重排 rows 的前提下，沿 column 視窗掃描到 n_features，計算每 tile 可關欄數總和
  prop_df$closed_cols_all_windows  <- vapply(prop_df$rows_idx,  function(x)
    closed_cols_across_windows_for_tile(paths_ranked, x[[1]], tile_cols = 128L), integer(1))

  naive_df$closed_cols_all_windows <- vapply(naive_df$rows_idx, function(x)
    closed_cols_across_windows_for_tile(paths_ranked, x[[1]], tile_cols = 128L), integer(1))

  # 3) 彙整輸出
  all <- rbind(prop_df, naive_df)

  # 只針對「跨所有視窗」做 proposed vs naive 的總可關欄數比較
  total_sum <- aggregate(closed_cols_all_windows ~ method, data = all, sum)
  colnames(total_sum)[2] <- "total_closed_cols_all_windows"

  # 也保留第一個視窗（原本）的摘要，方便對照
  summary_first_window <- aggregate(
    cbind(active_cols, closed_cols) ~ method,
    data = all,
    FUN = function(x) c(mean = mean(x), median = median(x), p90 = quantile(x, 0.9))
  )

  list(
    columns_order = ord,
    tiles   = all,                          # 含每 tile 的 rows_idx 與各窗總可關欄數
    total_closed_cols_all_windows = total_sum, # 兩法整體省下的欄數（跨所有視窗）
    summary_first_window = summary_first_window
  )
}

proposed_vs_naive_main <- function(fit, data, one_based = TRUE) {
  library(ranger)

  # 1) 取得模型實際使用的 feature 數
  feat_names <- fit$forest$independent.variable.names
  if (!is.null(feat_names) && length(feat_names) > 0L) {
    n_features <- length(feat_names)
  } else {
    n_features <- ncol(data)
  }
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features <= 0L) {
    stop("n_features must be a positive integer (> 0).")
  }

  # 2) 逐棵樹萃取所有 root->leaf paths（0-based），並轉為 1-based
  paths_list <- list()
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    paths_t <- extract_paths(info)  # 每條路徑是一個整數向量（0-based feature id）

    # 0-based -> 1-based；過濾到 1..n_features，並移除重複
    paths_t <- lapply(paths_t, function(p0) {
      if (is.null(p0) || length(p0) == 0L) return(integer(0))
      p1 <- as.integer(p0) + 1L
      p1 <- p1[p1 >= 1L & p1 <= n_features]
      if (length(p1)) unique(p1) else integer(0)
    })

    # 合併到總清單（list 用 c() 合併）
    paths_list <- c(paths_list, paths_t)
  }

  # 3) 丟進前面完成的比較主流程（會自動：以頻率排序欄、proposed/naive 打包、
  #    第一個 128 視窗評估，以及跨所有 column 視窗的總可關欄數）
  res <- compare_proposed_vs_naive(paths_list, n_features, one_based = TRUE)

  return(res)
}


#----------------------- experiment --------------------------------


pvn_email <- proposed_vs_naive_main(email_forest100,email.x)
pvn_email$total_closed_cols_all_windows
print(pvn_email$wlbl_total)
print(pvn_email$hamming_summary)

pvn_fetal <- proposed_vs_naive_main(fetal_forest100,fetal_x)
pvn_fetal$total_closed_cols_all_windows
print(pvn_fetal$wlbl_total)
print(pvn_fetal$hamming_summary)

pvn_gene  <- proposed_vs_naive_main(gene_forest100,gene_train_x)
pvn_gene$total_closed_cols_all_windows
print(pvn_gene$wlbl_total)
print(pvn_gene$hamming_summary)

pvn_mush  <- proposed_vs_naive_main(mush_forest100,mush.x)
pvn_mush$total_closed_cols_all_windows
print(pvn_mush$wlbl_total)
print(pvn_mush$hamming_summary)

pvn_loan  <- proposed_vs_naive_main(loan_forest100,loan_train_.x)
pvn_loan$total_closed_cols_all_windows
print(pvn_loan$wlbl_total)
print(pvn_loan$hamming_summary)

pvn_gesture <- proposed_vs_naive_main(gesture_forest100,gesture_train_x)
pvn_gesture$total_closed_cols_all_windows
print(pvn_gesture$wlbl_total)
print(pvn_gesture$hamming_summary)
