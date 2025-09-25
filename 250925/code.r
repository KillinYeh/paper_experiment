library(ranger)

#------------------function------------------
# 計算一棵樹的所有 root→leaf path，並統計每條 path 的利用率
path_utilization_one_tree <- function(info, Nfeat) {
  # nodeID -> children
  left  <- setNames(info$leftChild,  info$nodeID)
  right <- setNames(info$rightChild, info$nodeID)
  term  <- setNames(info$terminal,   info$nodeID)
  split_name <- setNames(as.character(info$splitvarName), info$nodeID)

  paths <- list()

  # DFS 搜尋所有 root→leaf path
  dfs <- function(node, used_feats) {
    if (isTRUE(term[as.character(node)])) {
      # 到葉，紀錄這條 path 用到的 feature
      paths[[length(paths) + 1]] <<- used_feats
    } else {
      fname <- split_name[as.character(node)]
      new_used <- unique(c(used_feats, fname))
      if (!is.na(left[as.character(node)])) {
        dfs(left[as.character(node)], new_used)
      }
      if (!is.na(right[as.character(node)])) {
        dfs(right[as.character(node)], new_used)
      }
    }
  }
  dfs(0L, character(0))

  # 計算每條 path 的利用率
  k <- sapply(paths, length)
  utilization <- k / Nfeat

  data.frame(
    path_id = seq_along(k),
    used_features = k,
    utilization = utilization
  )
}

# 對整個 forest 計算
path_utilization_forest <- function(fit, data) {
  Nfeat <- ncol(data) - 1  # 假設 data 最後一欄是 target
  res <- list()
  for (t in seq_len(fit$num.trees)) {
    info <- ranger::treeInfo(fit, t)
    df <- path_utilization_one_tree(info, Nfeat)
    df$tree_id <- t
    res[[t]] <- df
  }
  do.call(rbind, res)
}

#---------------------experiment---------------------
utilization_email = path_utilization_forest(email_forest500 , email.x)
mean(utilization_email$utilization)

utilization_fetal = path_utilization_forest(fetal_forest , fetal_x)
mean(utilization_email$utilization)

utilization_gene = path_utilization_forest(gene_forest500 , gene_train_x)
mean(utilization_email$utilization)

utilization_mush = path_utilization_forest(mush_forest500 , mush.x)
mean(utilization_email$utilization)

utilization_loan = path_utilization_forest(loan_forest500 , loan_train_x)
mean(utilization_email$utilization)

utilization_gesture = path_utilization_forest(gesture_forest500 , gesture_train_x)
mean(utilization_email$utilization)



