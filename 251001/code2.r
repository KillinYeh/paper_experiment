library(ranger)

# DFS 找所有 root→leaf path
extract_paths <- function(info) {
  left  <- setNames(info$leftChild,  info$nodeID)
  right <- setNames(info$rightChild, info$nodeID)
  term  <- setNames(info$terminal,   info$nodeID)
  split_id <- setNames(info$splitvarID, info$nodeID)  # feature index (0-based)

  paths <- list()
  dfs <- function(node, used_feats) {
    if (isTRUE(term[as.character(node)])) {
      paths[[length(paths) + 1]] <<- used_feats
    } else {
      fid <- split_id[as.character(node)]
      new_used <- unique(c(used_feats, fid))
      if (!is.na(left[as.character(node)])) {
        dfs(left[as.character(node)], new_used)
      }
      if (!is.na(right[as.character(node)])) {
        dfs(right[as.character(node)], new_used)
      }
    }
  }
  dfs(0L, integer(0))
  paths
}

# 計算 tile utilization (row=256, col=min(2^floor(log2(Nfeat)),256))
tile_utilization_one_tree <- function(info, Nfeat, row_tile_size = 256) {
  paths <- extract_paths(info)
  n_paths <- length(paths)

  # 動態決定 column tile size
  col_tile_size <- min(2^(floor(log2(Nfeat))), 256)

  n_row_tiles <- ceiling(n_paths / row_tile_size)
  n_col_tiles <- ceiling(Nfeat / col_tile_size)

  results <- data.frame()
  for (rt in 0:(n_row_tiles-1)) {
    row_start <- rt * row_tile_size + 1
    row_end   <- min((rt+1) * row_tile_size, n_paths)
    row_paths <- paths[row_start:row_end]

    for (ct in 0:(n_col_tiles-1)) {
      col_start <- ct * col_tile_size
      col_end   <- min((ct+1) * col_tile_size - 1, Nfeat-1)

      used <- 0
      for (p in row_paths) {
        used <- used + sum(p >= col_start & p <= col_end)
      }
      util <- used / (row_tile_size * col_tile_size)
      results <- rbind(results, data.frame(
        row_tile = rt,
        col_tile = ct,
        used_features = used,
        utilization = util,
        row_tile_size = row_tile_size,
        col_tile_size = col_tile_size
      ))
    }
  }
  results
}

# Forest 層級
tile_utilization_forest <- function(fit, data, row_tile_size = 256) {
  Nfeat <- ncol(data)
  res <- list()
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    df <- tile_utilization_one_tree(info, Nfeat, row_tile_size)
    df$tree_id <- t
    res[[t]] <- df
  }
  do.call(rbind, res)
}
#-------------------experiment-------------------

Tile_utilization_email = tile_utilization_forest(email_forest500 , email.x)
mean(Tile_utilization_email$utilization)

Tile_utilization_fetal = tile_utilization_forest(fetal_forest , fetal_x)
mean(Tile_utilization_fetal$utilization)

Tile_utilization_gene = tile_utilization_forest(gene_forest500 , gene_train_x)
mean(Tile_utilization_gene$utilization)

Tile_utilization_mush = tile_utilization_forest(mush_forest500 , mush.x)
mean(Tile_utilization_mush$utilization)

Tile_utilization_loan = tile_utilization_forest(loan_forest500 , loan_train_x)
mean(Tile_utilization_loan$utilization)

Tile_utilization_gesture = tile_utilization_forest(gesture_forest500 , gesture_train_x)
mean(Tile_utilization_gesture$utilization)


#------------------  result ----------------------
import pandas as pd
import matplotlib.pyplot as plt

# 整理你的數據
data = {
    "dataset": ["email", "fetal", "gene", "mush", "loan", "gesture"],
    "mean_utilization": [
        0.003601231,
        0.1048333,
        0.00001993928,
        0.7899027,
        0.1781936,
        0.3275146,
    ]
}

df = pd.DataFrame(data)

# 畫柱狀圖
plt.figure(figsize=(8,5))
bars = plt.bar(df["dataset"], df["mean_utilization"], color="skyblue", edgecolor="black")

plt.title("Average CAM tile Utilization with adaptive column size")
plt.ylabel("Mean Utilization (used features / total features)")
plt.xlabel("Dataset")

# 在柱子上標出數值
for bar, val in zip(bars, df["mean_utilization"]):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
             f"{val:.3f}", ha="center", va="bottom", fontsize=9)

plt.ylim(0, 1.1)  # 因為最大值接近 1
plt.tight_layout()
plt.show()

