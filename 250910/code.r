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



#----------------plot chart------------------
import pandas as pd
import matplotlib.pyplot as plt

# 整理數據
data = {
    "dataset": ["email", "fetal", "gene", "mush", "loan", "gesture"],
    "avg_leaf_per_tree": [339.718, 115.19, 20.214, 1223.796, 112.474, 1649.476],
    "total_leafs": [169859, 57595, 10107, 611898, 56237, 824738],
    "horizontal_mean_loads": [43, 16, 2, 152, 15, 203],
    "vertical_mean_loads": [2, 1, 1, 5, 1, 7],
}

df = pd.DataFrame(data)

# 顯示表格
print(df)

# 畫長條圖比較 horizontal vs vertical
fig, ax = plt.subplots(figsize=(10,6))
width = 0.35
x = range(len(df))

ax.bar([i - width/2 for i in x], df["horizontal_mean_loads"], width, label="Horizontal")
ax.bar([i + width/2 for i in x], df["vertical_mean_loads"], width, label="Vertical")

ax.set_xticks(x)
ax.set_xticklabels(df["dataset"])
ax.set_ylabel("Mean load times")
ax.set_title("Horizontal vs Vertical CAM load times per dataset")
ax.legend()

plt.tight_layout()
plt.show()

# 也可以計算 horizontal/vertical 的比值
df["horiz_to_vert_ratio"] = df["horizontal_mean_loads"] / df["vertical_mean_loads"]
print("\nHorizontal / Vertical ratio:")
print(df[["dataset", "horiz_to_vert_ratio"]])

