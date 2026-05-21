############################################################
## Total Energy Comparison:
## Read/Evaluate Energy + Write Energy
############################################################

calculate_total_energy_from_result <- function(result, dataset_name, testing_dataset_size) {
  if (is.null(result$energy$method_summary)) {
    stop(paste0("Input error: ", dataset_name, " does not contain result$energy$method_summary"))
  }

  df <- result$energy$method_summary

  required_cols <- c(
    "method",
    "fake_match_total_fj",
    "rewrite_total_fj",
    "fake_match_write_pj",
    "rewrite_write_pj"
  )

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(paste0(
      "Input error: ", dataset_name,
      " missing columns: ",
      paste(missing_cols, collapse = ", ")
    ))
  }

  out <- df[, required_cols]
  out$dataset <- dataset_name
  out$testing_dataset_size <- testing_dataset_size

  # Read / evaluate energy:
  # method_summary 裡的 read energy 是 fJ / inference，
  # 所以先除以 1000 轉成 pJ，再乘上 testing dataset size。
  out$fake_match_read_pj <-
    (out$fake_match_total_fj / 1000) * testing_dataset_size

  out$rewrite_read_pj <-
    (out$rewrite_total_fj / 1000) * testing_dataset_size

  # Total energy = testing-set read energy + one-time model write energy
  out$fake_match_total_energy_pj <-
    out$fake_match_read_pj + out$fake_match_write_pj

  out$rewrite_total_energy_pj <-
    out$rewrite_read_pj + out$rewrite_write_pj

  # Fake Match / Rewrite total energy ratio
  out$fake_over_rewrite_total_energy_ratio <- ifelse(
    out$rewrite_total_energy_pj == 0,
    NA_real_,
    out$fake_match_total_energy_pj / out$rewrite_total_energy_pj
  )

  # Positive means Fake Match saves energy.
  out$fake_total_energy_reduction_ratio <- ifelse(
    out$rewrite_total_energy_pj == 0,
    NA_real_,
    (out$rewrite_total_energy_pj - out$fake_match_total_energy_pj) /
      out$rewrite_total_energy_pj
  )

  out$fake_total_energy_saved_pj <-
    out$rewrite_total_energy_pj - out$fake_match_total_energy_pj

  out$winner <- ifelse(
    out$fake_match_total_energy_pj < out$rewrite_total_energy_pj,
    "Fake Match",
    "Rewrite"
  )

  out <- out[, c(
    "dataset",
    "method",
    "testing_dataset_size",
    "fake_match_read_pj",
    "rewrite_read_pj",
    "fake_match_write_pj",
    "rewrite_write_pj",
    "fake_match_total_energy_pj",
    "rewrite_total_energy_pj",
    "fake_total_energy_saved_pj",
    "fake_over_rewrite_total_energy_ratio",
    "fake_total_energy_reduction_ratio",
    "winner"
  )]

  rownames(out) <- NULL
  out
}

############################################################
## Dataset testing sizes
############################################################

total_energy_energy <- calculate_total_energy_from_result(
  result = sim_energy,
  dataset_name = "Energy",
  testing_dataset_size = 154
)

total_energy_breast <- calculate_total_energy_from_result(
  result = sim_breast,
  dataset_name = "Breast",
  testing_dataset_size = 115
)

total_energy_heart <- calculate_total_energy_from_result(
  result = sim_heart,
  dataset_name = "Heart",
  testing_dataset_size = 56
)

total_energy_parkinsons <- calculate_total_energy_from_result(
  result = sim_parkinsons,
  dataset_name = "Parkinsons",
  testing_dataset_size = 1201
)

############################################################
## Merge all results
############################################################

total_energy_comparison <- rbind(
  total_energy_energy,
  total_energy_breast,
  total_energy_heart,
  total_energy_parkinsons
)

############################################################
## Print and save
############################################################

print(total_energy_comparison)

write.csv(
  total_energy_comparison,
  file = "total_energy_comparison_embedded_datasets.csv",
  row.names = FALSE
)

saveRDS(
  total_energy_comparison,
  file = "total_energy_comparison_embedded_datasets.rds"
)
