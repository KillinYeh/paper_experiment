library(ranger)

# 由 treeInfo 算每個 node 的深度（root=0）
.node_depths_from_treeinfo <- function(info) {
  left  <- setNames(info$leftChild,  info$nodeID)
  right <- setNames(info$rightChild, info$nodeID)
  term  <- setNames(info$terminal,   info$nodeID)

  depth <- setNames(rep(NA_integer_, nrow(info)), info$nodeID)
  root <- 0L
  depth[as.character(root)] <- 0L
  q <- root

  while (length(q)) {
    nid <- q[1]; q <- q[-1]
    d <- depth[as.character(nid)]
    if (!isTRUE(term[as.character(nid)])) {
      lc <- left [as.character(nid)]
      rc <- right[as.character(nid)]
      if (!is.na(lc)) { depth[as.character(lc)] <- d + 1L; q <- c(q, lc) }
      if (!is.na(rc)) { depth[as.character(rc)] <- d + 1L; q <- c(q, rc) }
    }
  }
  as.integer(depth)
}

# 計算：在 caps=2^H 的設定下，單棵樹最壞需要幾次 load
max_loads_one_tree <- function(info, caps = 256L) {
  H <- as.integer(round(log2(caps)))
  if (2^H != caps) stop("caps 必須是 2 的冪（例如 2/4/8/.../256）")

  depths <- .node_depths_from_treeinfo(info)
  is_internal <- !info$terminal
  # 只在深度為 0, H, 2H, ... 的「內部節點」各觸發一次 load
  sum(is_internal & (depths %% H == 0L))
}

# 對整個 forest（ranger 模型）估計最壞 loads
max_loads_forest <- function(fit, caps = 256L) {
  stopifnot(inherits(fit, "ranger"), !is.null(fit$forest))
  per_tree <- sapply(seq_len(fit$num.trees), function(t) {
    info <- ranger::treeInfo(fit, t)
    max_loads_one_tree(info, caps = caps)
  })
  list(caps = caps, per_tree_loads = per_tree, total_loads = sum(per_tree))
}


#----------------experiment------------------

horizontal_email <- max_loads_forest(email_forest500 , caps=256)
print(horizontal_email$per_tree_loads)
print(mean(horizontal_email$per_tree_loads))
horizontal_fetal <- max_loads_forest(fetal_forest , caps=256)
print(horizontal_fetal$per_tree_loads)
print(mean(horizontal_fetal$per_tree_loads))
horizontal_gene <- max_loads_forest(gene_forest500 , caps=256)
print(horizontal_gene$per_tree_loads)
print(mean(horizontal_gene$per_tree_loads))
horizontal_mush <- max_loads_forest(mush_forest500 , caps=256)
print(horizontal_mush$per_tree_loads)
print(mean(horizontal_mush$per_tree_loads))
horizontal_loan <- max_loads_forest(loan_forest500 , caps=256)
print(horizontal_loan$per_tree_loads)
print(mean(horizontal_loan$per_tree_loads))
horizontal_gesture <- max_loads_forest(gesture_forest500 , caps=256)
print(horizontal_gesture$per_tree_loads)
print(mean(horizontal_gesture$per_tree_loads))
#horizontal_covtype <- max_loads_forest(covtype_forest500 , caps=256)
#print(horizontal_email$per_tree_loads)

