analyze_tile_tree_usage <- function(distri, outdir = "~/paper/ranger/R/ranger_proj_test/experiment/251204/chart") {

  # =============== 準備輸出資料夾 ===============
  outdir <- path.expand(outdir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # 抓 dataset 名稱（例如 gesture_distri）
  dataset_name <- deparse(substitute(distri))

  # =============== 基本資訊 ===============
  n_tiles <- length(distri)

  # 抓全部 tree id
  all_tree_ids <- unique(as.integer(unlist(
    lapply(distri, function(x) names(x$tree_counts))
  )))
  n_trees <- max(all_tree_ids)

  # =============== 建立 tile × tree 使用矩陣 ===============
  tile_tree_mat <- matrix(0L, nrow = n_tiles, ncol = n_trees)

  for (i in seq_len(n_tiles)) {
    tc <- distri[[i]]$tree_counts
    ids  <- as.integer(names(tc))
    vals <- as.integer(tc)
    tile_tree_mat[i, ids] <- vals
  }

  colnames(tile_tree_mat) <- paste0("tree_", seq_len(n_trees))
  rownames(tile_tree_mat) <- paste0("tile_", seq_len(n_tiles))

  # =============== 計算 Top-share & Entropy ===============

  tile_totals <- rowSums(tile_tree_mat)

  top_share <- apply(tile_tree_mat, 1, function(x) {
    if (sum(x) == 0) return(NA_real_)
    max(x) / sum(x)
  })

  entropy <- apply(tile_tree_mat, 1, function(x) {
    s <- sum(x)
    if (s == 0) return(NA_real_)
    p <- x[x > 0] / s
    -sum(p * log2(p))
  })

  tile_summary <- data.frame(
    tile_id   = seq_len(n_tiles),
    top_share = top_share,
    entropy   = entropy
  )

  # =============== 計算累積使用量 ===============
  cum_usage <- apply(tile_tree_mat, 2, cumsum)

  # =============== (1) Top-share plot ===============
  png(file.path(outdir, sprintf("%s_top_treeshare.png", dataset_name)),
      width = 1400, height = 900)
  plot(tile_summary$tile_id, tile_summary$top_share,
       type = "l", col = "blue", lwd = 1.8,
       xlab = "Tile Index",
       ylab = "Top Tree Share",
       main = sprintf("%s Top Tree Share per Tile", dataset_name))
  dev.off()

  # =============== (2) Entropy plot ===============
  png(file.path(outdir, sprintf("%s_entropy.png", dataset_name)),
      width = 1400, height = 900)
  plot(tile_summary$tile_id, tile_summary$entropy,
       type = "l", col = "darkgreen", lwd = 1.8,
       xlab = "Tile Index",
       ylab = "Entropy",
       main = sprintf("%s Tree Usage Entropy per Tile", dataset_name))
  dev.off()

  # =============== (3) Cumulative usage plot ===============
  trees_to_plot <- n_trees
  png(file.path(outdir, sprintf("%s_cumulative.png", dataset_name)),
      width = 1600, height = 1000)
  matplot(
    x = seq_len(n_tiles),
    y = cum_usage[, 1:trees_to_plot, drop = FALSE],
    type = "l", lty = 1, lwd = 1.8,
    xlab = "Tile Index",
    ylab = "Cumulative Usage",
    main = sprintf("%s Cumulative Tree Usage (first %d trees)",
                   dataset_name, trees_to_plot)
  )
  legend("topleft",
         legend = paste0("tree_", 1:trees_to_plot),
         col = 1:trees_to_plot,
         lty = 1, lwd = 1.8, cex = 0.7)
  dev.off()

  # =============== 回傳結果 ===============
  return(list(
    tile_tree_mat = tile_tree_mat,
    tile_summary  = tile_summary,
    cum_usage     = cum_usage,
    outdir = outdir
  ))
}

#-------------------experiment--------------------
email_analyze <-  analyze_tile_tree_usage(email_distri)
fetal_analyze <-  analyze_tile_tree_usage(fetal_distri)
gene_analyze <-  analyze_tile_tree_usage(gene_distri)
mush_analyze <-  analyze_tile_tree_usage(mush_distri)
loan_analyze <-  analyze_tile_tree_usage(loan_distri)
gesture_analyze <-  analyze_tile_tree_usage(gesture_distri)
arcene_analyze <-  analyze_tile_tree_usage(arcene_distri)
gas_analyze <-  analyze_tile_tree_usage(gas_distri)
