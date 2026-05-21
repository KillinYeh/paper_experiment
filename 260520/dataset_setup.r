############################################################
## Edge Tabular Dataset Loader + Ranger Training Script
## Datasets:
## 1. Energy Efficiency (ENB2012_data.xlsx)
## 2. Breast Cancer Wisconsin Diagnostic (wdbc.data)
## 3. Heart Disease / Cleveland (cleveland.data or processed.cleveland.data)
## 4. Parkinsons Telemonitoring (parkinsons_updrs.data)
############################################################

## ---------------------------------------------------------
## 0. Packages
## ---------------------------------------------------------

library(ranger)
library(readxl)

## ---------------------------------------------------------
## 1. Global setting
## ---------------------------------------------------------

set.seed(42)

## Put the four files in the same folder as this R script,
## or change data_dir to your dataset folder.
data_dir <- "./paper/ranger/R/ranger_proj_test/Dataset_Embedded"

energy_file     <- file.path(data_dir, "ENB2012_data.xlsx")
breast_file     <- file.path(data_dir, "wdbc.data")
heart_file      <- file.path(data_dir, "cleveland.data")
parkinsons_file <- file.path(data_dir, "parkinsons_updrs.data")

## UCI does not provide fixed train/test files for these four uploaded datasets.
## We therefore use reproducible 80/20 splits.
## For Parkinsons Telemonitoring, we split by subject# to avoid putting the
## same patient in both training and testing sets.

train_ratio <- 0.80


## ---------------------------------------------------------
## 2. Split helpers
## ---------------------------------------------------------

split_regression <- function(df, train_ratio = 0.80, seed = 42) {
    set.seed(seed)
    n <- nrow(df)
    train_idx <- sample(seq_len(n), size = floor(train_ratio * n), replace = FALSE)
    list(
        train = df[train_idx, , drop = FALSE],
        test  = df[-train_idx, , drop = FALSE]
    )
}


split_classification_stratified <- function(df, target_col, train_ratio = 0.80, seed = 42) {
    set.seed(seed)
    y <- df[[target_col]]
    train_idx <- unlist(lapply(split(seq_len(nrow(df)), y), function(idx) {
        sample(idx, size = floor(train_ratio * length(idx)), replace = FALSE)
    }))
    
    train_idx <- sort(as.integer(train_idx))
    
    list(
        train = df[train_idx, , drop = FALSE],
        test  = df[-train_idx, , drop = FALSE]
    )
}


split_by_group <- function(df, group_col, train_ratio = 0.80, seed = 42) {
    set.seed(seed)
    
    groups <- unique(df[[group_col]])
    train_groups <- sample(groups, size = floor(train_ratio * length(groups)), replace = FALSE)
    
    train_idx <- df[[group_col]] %in% train_groups
    
    list(
        train = df[train_idx, , drop = FALSE],
        test  = df[!train_idx, , drop = FALSE],
        train_groups = train_groups
    )
}


## ---------------------------------------------------------
## 3. Ranger training helper
## ---------------------------------------------------------

train_ranger_xy <- function(train_df, target_col, task = c("classification", "regression"),
                            num.trees = 100) {
    task <- match.arg(task)
    
    train_df <- as.data.frame(train_df)
    train_df <- na.omit(train_df)
    
    if (task == "classification") {
        train_df[[target_col]] <- as.factor(train_df[[target_col]])
    } else {
        train_df[[target_col]] <- as.numeric(train_df[[target_col]])
    }
    
    feature_cols <- setdiff(names(train_df), target_col)
    
    dataset.x <- train_df[, feature_cols, drop = FALSE]
    dataset.y <- train_df[[target_col]]
    
    fit <- ranger(
        x = dataset.x,
        y = dataset.y,
        importance = "permutation",
        scale.permutation.importance = TRUE,
        local.importance = TRUE,
        verbose = TRUE,
        num.trees = num.trees
    )
    
    list(
        model = fit,
        x = dataset.x,
        y = dataset.y,
        feature_cols = feature_cols,
        target_col = target_col,
        task = task
    )
}


## ---------------------------------------------------------
## 4. Dataset loaders
## ---------------------------------------------------------

load_energy_efficiency <- function(path) {
    df <- readxl::read_excel(path)
    df <- as.data.frame(df)
    
    ## Official targets:
    ## Y1 = Heating Load
    ## Y2 = Cooling Load
    ## Features:
    ## X1~X8
    names(df) <- trimws(names(df))
    
    df
}


load_breast_wdbc <- function(path) {
    wdbc_cols <- c(
        "ID", "Diagnosis",
        "radius_mean", "texture_mean", "perimeter_mean", "area_mean", "smoothness_mean",
        "compactness_mean", "concavity_mean", "concave_points_mean", "symmetry_mean", "fractal_dimension_mean",
        "radius_se", "texture_se", "perimeter_se", "area_se", "smoothness_se",
        "compactness_se", "concavity_se", "concave_points_se", "symmetry_se", "fractal_dimension_se",
        "radius_worst", "texture_worst", "perimeter_worst", "area_worst", "smoothness_worst",
        "compactness_worst", "concavity_worst", "concave_points_worst", "symmetry_worst", "fractal_dimension_worst"
    )
    
    df <- read.csv(path, header = FALSE, col.names = wdbc_cols, stringsAsFactors = FALSE)
    
    ## Official target:
    ## Diagnosis = M / B
    ## ID is an identifier and should not be used as a feature.
    df$ID <- NULL
    df$Diagnosis <- as.factor(df$Diagnosis)
    
    df
}


read_processed_cleveland <- function(path) {
    heart_cols <- c(
        "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
        "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num"
    )
    
    df <- read.csv(path, header = FALSE, col.names = heart_cols,
                   na.strings = c("?", "-9"), stringsAsFactors = FALSE)
    
    for (nm in names(df)) {
        df[[nm]] <- as.numeric(df[[nm]])
    }
    
    df <- na.omit(df)
    
    ## Official target:
    ## num = 0 means absence; num = 1,2,3,4 means presence.
    df$target <- factor(ifelse(df$num == 0, "absence", "presence"))
    df$num <- NULL
    
    df
}


read_raw_cleveland_76attr <- function(path) {
    ## The uploaded cleveland.data is the original 76-attribute wrapped format,
    ## where one patient record is spread across multiple lines and usually ends
    ## with the token "name".
    ##
    ## We extract the 14 commonly used attributes:
    ## age, sex, cp, trestbps, chol, fbs, restecg, thalach,
    ## exang, oldpeak, slope, ca, thal, num.
    ##
    ## If you have processed.cleveland.data, prefer that file instead.
    raw_lines <- readLines(path, warn = FALSE, encoding = "latin1")
    tokens <- unlist(strsplit(paste(raw_lines, collapse = " "), "\\s+"))
    tokens <- tokens[nchar(tokens) > 0]
    
    records <- list()
    current <- character(0)
    
    for (tok in tokens) {
        current <- c(current, tok)
        if (tok == "name") {
            if (length(current) == 76L) {
                records[[length(records) + 1L]] <- current
            }
            current <- character(0)
        }
    }
    
    if (length(records) == 0L) {
        stop("No valid Cleveland records were parsed. Please use processed.cleveland.data if available.")
    }
    
    ## 1-based positions in the original 76-attribute Cleveland format.
    pos <- c(
        age = 3, sex = 4, cp = 9, trestbps = 10, chol = 12,
        fbs = 16, restecg = 19, thalach = 32, exang = 38,
        oldpeak = 40, slope = 41, ca = 44, thal = 51, num = 58
    )
    
    mat <- do.call(rbind, lapply(records, function(rec) rec[as.integer(pos)]))
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
    names(df) <- names(pos)
    
    df[df == "-9"] <- NA
    df[df == "?"] <- NA
    
    for (nm in names(df)) {
        df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    }
    
    df <- na.omit(df)
    
    df$target <- factor(ifelse(df$num == 0, "absence", "presence"))
    df$num <- NULL
    
    if (nrow(df) < 300) {
        warning(sprintf(
            "Parsed Cleveland raw file has %d complete rows. The uploaded cleveland.data may be incomplete/corrupted. Prefer processed.cleveland.data if possible.",
            nrow(df)
        ))
    }
    
    df
}


load_heart_cleveland <- function(path) {
    ## If the file already looks like processed.cleveland.data, use the processed reader.
    ## Otherwise parse the original 76-attribute wrapped format.
    first_line <- readLines(path, n = 1, warn = FALSE, encoding = "latin1")
    first_fields_comma <- strsplit(first_line, ",")[[1]]
    first_fields_space <- strsplit(first_line, "\\s+")[[1]]
    
    if (length(first_fields_comma) >= 14L) {
        read_processed_cleveland(path)
    } else if (length(first_fields_space) == 14L) {
        read.table(path, header = FALSE, na.strings = c("?", "-9"),
                   col.names = c(
                       "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
                       "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num"
                   )) |>
            transform(target = factor(ifelse(num == 0, "absence", "presence"))) |>
            subset(select = -num) |>
            na.omit()
    } else {
        read_raw_cleveland_76attr(path)
    }
}


load_parkinsons_updrs <- function(path) {
    df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE)
    
    ## Official targets:
    ## motor_UPDRS and total_UPDRS
    ##
    ## subject. is an ID-like grouping variable. We use it for subject-level splitting,
    ## but remove it from model features.
    names(df) <- make.names(names(df))
    
    df
}


## ---------------------------------------------------------
## 5. Load datasets
## ---------------------------------------------------------

energy_df <- load_energy_efficiency(energy_file)
breast_df <- load_breast_wdbc(breast_file)
heart_df <- load_heart_cleveland(heart_file)
parkinsons_df <- load_parkinsons_updrs(parkinsons_file)


## ---------------------------------------------------------
## 6. Train/test split
## ---------------------------------------------------------

## Energy Efficiency: no official train/test split in the uploaded UCI file.
energy_split <- split_regression(energy_df, train_ratio = train_ratio, seed = 42)
energy_train <- energy_split$train
energy_test  <- energy_split$test

## Breast Cancer Wisconsin Diagnostic: no official train/test split in the uploaded UCI file.
breast_split <- split_classification_stratified(breast_df, "Diagnosis", train_ratio = train_ratio, seed = 42)
breast_train <- breast_split$train
breast_test  <- breast_split$test

## Cleveland Heart Disease: no official train/test split in the uploaded UCI file.
heart_split <- split_classification_stratified(heart_df, "target", train_ratio = train_ratio, seed = 42)
heart_train <- heart_split$train
heart_test  <- heart_split$test

## Parkinsons Telemonitoring:
## no official train/test split in the uploaded UCI file.
## Because there are repeated recordings for each subject, use subject-level split.
parkinsons_split <- split_by_group(parkinsons_df, "subject.", train_ratio = train_ratio, seed = 42)
parkinsons_train_full <- parkinsons_split$train
parkinsons_test_full  <- parkinsons_split$test

## Remove subject ID from features after splitting.
parkinsons_train <- parkinsons_train_full
parkinsons_test  <- parkinsons_test_full
parkinsons_train$subject. <- NULL
parkinsons_test$subject.  <- NULL


## ---------------------------------------------------------
## 7. Ranger training
## ---------------------------------------------------------

## Energy Efficiency: train two regression models.
energy_y1_fit <- train_ranger_xy(
    train_df = energy_train[, setdiff(names(energy_train), "Y2"), drop = FALSE],
    target_col = "Y1",
    task = "regression",
    num.trees = 100
)

energy_y2_fit <- train_ranger_xy(
    train_df = energy_train[, setdiff(names(energy_train), "Y1"), drop = FALSE],
    target_col = "Y2",
    task = "regression",
    num.trees = 100
)

energy_y1_forest100 <- energy_y1_fit$model
energy_y2_forest100 <- energy_y2_fit$model


## Breast Cancer Wisconsin Diagnostic: binary classification.
breast_fit <- train_ranger_xy(
    train_df = breast_train,
    target_col = "Diagnosis",
    task = "classification",
    num.trees = 100
)

breast_forest100 <- breast_fit$model


## Cleveland Heart Disease: binary classification.
heart_fit <- train_ranger_xy(
    train_df = heart_train,
    target_col = "target",
    task = "classification",
    num.trees = 100
)

heart_forest100 <- heart_fit$model


## Parkinsons Telemonitoring: train two regression models.
## For motor_UPDRS model, remove total_UPDRS from features to avoid target leakage.
parkinsons_motor_fit <- train_ranger_xy(
    train_df = parkinsons_train[, setdiff(names(parkinsons_train), "total_UPDRS"), drop = FALSE],
    target_col = "motor_UPDRS",
    task = "regression",
    num.trees = 100
)

## For total_UPDRS model, remove motor_UPDRS from features to avoid target leakage.
parkinsons_total_fit <- train_ranger_xy(
    train_df = parkinsons_train[, setdiff(names(parkinsons_train), "motor_UPDRS"), drop = FALSE],
    target_col = "total_UPDRS",
    task = "regression",
    num.trees = 100
)

parkinsons_motor_forest100 <- parkinsons_motor_fit$model
parkinsons_total_forest100 <- parkinsons_total_fit$model


## ---------------------------------------------------------
## 8. Collect x/y objects for later experiments
## ---------------------------------------------------------

energy_y1_train_x <- energy_y1_fit$x
energy_y1_train_y <- energy_y1_fit$y

breast_train_x <- breast_fit$x
breast_train_y <- breast_fit$y

heart_train_x <- heart_fit$x
heart_train_y <- heart_fit$y

parkinsons_motor_train_x <- parkinsons_motor_fit$x
parkinsons_motor_train_y <- parkinsons_motor_fit$y
parkinsons_total_train_x <- parkinsons_total_fit$x
parkinsons_total_train_y <- parkinsons_total_fit$y


## ---------------------------------------------------------
## 9. Quick report
## ---------------------------------------------------------

cat("\n================ Dataset Split Summary ================\n")
cat(sprintf("Energy Efficiency: train=%d, test=%d\n", nrow(energy_train), nrow(energy_test)))
cat(sprintf("Breast WDBC      : train=%d, test=%d\n", nrow(breast_train), nrow(breast_test)))
cat(sprintf("Heart Cleveland  : train=%d, test=%d\n", nrow(heart_train), nrow(heart_test)))
cat(sprintf("Parkinsons UPDRS : train=%d, test=%d, train_subjects=%d, test_subjects=%d\n",
            nrow(parkinsons_train), nrow(parkinsons_test),
            length(unique(parkinsons_train_full$subject.)),
            length(unique(parkinsons_test_full$subject.))))
cat("=======================================================\n")

cat("\nModels created:\n")
cat("  energy_y1_forest100\n")
cat("  energy_y2_forest100\n")
cat("  breast_forest100\n")
cat("  heart_forest100\n")
cat("  parkinsons_motor_forest100\n")
cat("  parkinsons_total_forest100\n")
