# training dataset ,label distribution

# API -----------------------------------------------------------------
count_class <- function()
{
    counts <- table(email.data[[3002]])
    print("count total member of different class in email dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in email dataset")
    print(round(counts,3))
    
    counts <- table(fetal.data$fetal_health)
    print("count total member of different class in fetal dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in fetal dataset")
    print(round(counts,3))
    
    counts <- table(gene_y)
    print("count total member of different class in gene dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in gene dataset")
    print(round(counts,3))
    
    counts <- table(mush.data$class)
    print("count total member of different class in mush dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in mush dataset")
    print(round(counts,3))
    
    counts <- table(loan_data$loan_status)
    print("count total member of different class in loan dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in loan dataset")
    print(round(counts,3))
    
    # dataset provide from X-TIME 
    counts <- table(gesture_data[[33]])
    print("count total member of different class in gesture dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in gesture dataset")
    print(round(counts,3))
    
    counts <- table(cov_data$V55)
    print("count total member of different class in cov_type dataset")
    print(counts)
    counts <- prop.table(counts)
    print("count distribution of different class in cov_type dataset")
    print(round(counts,3))
}


# fit model , label distribution

library(ranger)

# 統計整個森林所有「葉節點標籤」的分佈
leaf_label_counts <- function(fit) {
  stopifnot(inherits(fit, "ranger"))
  labs <- unlist(lapply(seq_len(fit$num.trees), function(t) {
    info <- ranger::treeInfo(fit, t)
    as.character(info$prediction[info$terminal])  # 只取葉子
  }))
  counts <- sort(table(labs), decreasing = TRUE)
  props  <- prop.table(counts)
  list(counts = counts, proportion = props, leaf_labels_all = labs)
}

#（可選）同時回傳每棵樹各自的葉標籤向量
leaf_labels_per_tree <- function(fit) {
  lapply(seq_len(fit$num.trees), function(t) {
    info <- ranger::treeInfo(fit, t)
    as.character(info$prediction[info$terminal])
  })
}
# API end-------------------------------------------------------------




# experiment ---------------------------------------------------------

# email
counts <- table(email.data[[3002]])
print("count total member of different class in email dataset")
print(counts)
counts <- prop.table(counts)
print("count distribution of different class in email dataset")
print(round(counts,3))
res <- leaf_label_counts(email_forest500)
print("total number of label in email forest")
print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
print(sum(res$counts))
print("average leaf counts per tree")
print( sum(res$counts) / 500 )
print("propotion of label in email forest")
print(round(res$proportion, 3))



# fetal
counts <- table(fetal.data$fetal_health)
print("count total member of different class in fetal dataset")
print(counts)
counts <- prop.table(counts)
print("count distribution of different class in fetal dataset")
print(round(counts,3))
res <- leaf_label_counts(fetal_forest)
print("total number of label in fetal forest")
print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
print(sum(res$counts))
print("average leaf counts per tree")
print( sum(res$counts) / 500 )
print("propotion of label in fetal forest")
print(round(res$proportion, 3))

# gene
counts <- table(gene_y)
print("count total member of different class in gene dataset")
print(counts)
counts <- prop.table(counts)
print("count distribution of different class in gene dataset")
print(round(counts,3))
res <- leaf_label_counts(gene_forest500)
print("total number of label in gene forest")
print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
print(sum(res$counts))
print("average leaf counts per tree")
print( sum(res$counts) / 500 )
print("propotion of label in gene forest")
print(round(res$proportion, 3))

# mush
counts <- table(mush.data$class)
print("count total member of different class in mush dataset")
print(counts)
counts <- prop.table(counts)
print("count distribution of different class in mush dataset")
print(round(counts,3))
res <- leaf_label_counts(mush_forest500)
print("total number of label in mush forest")
print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
print(sum(res$counts))
print("average leaf counts per tree")
print( sum(res$counts) / 500 )
print("propotion of label in mush forest")
print(round(res$proportion, 3))

# loan
counts <- table(loan_data$loan_status)
print("count total member of different class in loan dataset")
print(counts)
counts <- prop.table(counts)
print("count distribution of different class in loan dataset")
print(round(counts,3))
res <- leaf_label_counts(loan_forest500)
print("total number of label in loan forest")
print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
print(sum(res$counts))
print("average leaf counts per tree")
print( sum(res$counts) / 500 )
print("propotion of label in loan forest")
print(round(res$proportion, 3))

# dataset provide from X-TIME
# gesture
counts <- table(gesture_data[[33]])
print("count total member of different class in gesture dataset")
print(counts)
counts <- prop.table(counts)
print("count distribution of different class in gesture dataset")
print(round(counts,3))
res <- leaf_label_counts(gesture_forest500)  
print("total number of label in gesture forest")
print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
print(sum(res$counts))
print("average leaf counts per tree")
print( sum(res$counts) / 500 )
print("propotion of label in gesture forest")
print(round(res$proportion, 3))

# cov_type
#counts <- table(cov_data$V55)
#print("count total member of different class in cov_type dataset")
#print(counts)
#counts <- prop.table(counts)
#print("count distribution of different class in cov_type dataset")
print(round(counts,3))
#res <- leaf_label_counts(cov_forest) 
#print("total number of label in cov_type forest")
#print(res$counts)       # 各類別在整個森林中被標成葉子標籤的次數
#print("average leaf counts per tree")
#print( sum(res$counts) / 500 )
#print("propotion of label in cov_type forest")
#print(round(res$proportion, 3))

# 如果想看每棵樹的葉標籤
#per_tree <- leaf_labels_per_tree(fetal_forest)
# 例如查看第1棵樹：
# per_tree[[1]]
