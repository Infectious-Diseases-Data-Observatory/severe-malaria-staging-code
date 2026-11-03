# CREATE BINNED DATASETS FOR EACH MULTIPLE IMPUTATION
# Uses the same candidate thresholds as the complete-case GOSDT analysis

library(data.table)
options(scipen = 999)

# ---- Paths ----
DATA_FILE <- file.path("Staging/Imputed/outputs/CSV/dat_all_hierarchical.csv")
OUT_DIR <- file.path("GOSDT/Python/Imputed_dataset_binned")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Load MI hierarchical dataset ----
my_table <- fread(DATA_FILE)

cat("\nMI hierarchical dataset\n")
cat("Rows: ", nrow(my_table), "\n", sep = "")
cat("Unique patients: ", uniqueN(my_table$ID), "\n", sep = "")
cat("Imputations: ", uniqueN(my_table$.imp), "\n", sep = "")

# ---- Trial indicator ----
# Match the complete-case analysis:
# Trial = 1 for SEAQUAMAT and 0 for AQUAMAT
my_table[, Trial := as.integer(trial == "SEAQUAMAT")]

# ---- Threshold definitions ----
# Use exactly the same candidate thresholds as the complete-case GOSDT analysis
reps_of_each <- c(1, 1, 3, 1, 1, 3, 1, 1, 4, 4, 3, 4, 4)
threshold_names <- rep(c("coma", "aconvulsions", "hb", "aanaemia", "ashockcomp", "temp",
                         "ahypoglyc", "Trial", "be", "bun", "parasitaemia", "mu_age",
                         "tring_obs"), times = reps_of_each)
threshold_symbols <- rep(c("=", "=", "<", "=", "=", "<", "=", "=", "<", ">", ">", ">", "<"),
                         times = reps_of_each)
threshold_values <- c(1, 1, 5, 8, 11, 1, 1, 37, 39, 40, 1, 1, -11, -9, -6, -3, 8, 10, 15,
                      20, 100000, 250000, 500000, 5, 10, 15, 20, 30, 40, 55, 75)
num_thresholds <- length(threshold_names)

# ---- Column names ----
threshold_full <- character(num_thresholds)

for (i in seq_len(num_thresholds)) {
  threshold_full[i] <- sprintf("%s:%s%s", threshold_names[i], threshold_symbols[i],
                               threshold_values[i])
}


# ---- Verify candidate threshold definitions ----
threshold_check <- data.table(Variable = threshold_names, Operator = threshold_symbols,
                              Threshold = threshold_values, Feature = threshold_full)

cat("\nCandidate threshold definitions:\n")
print(threshold_check)

# Check expected number of predictors
stopifnot(num_thresholds == 31)

# Check for duplicated feature names
stopifnot(!anyDuplicated(threshold_full))

# Check the eight staging features
expected_staging <- c("mu_age:>5", "mu_age:>10", "mu_age:>15", "mu_age:>20",
                      "tring_obs:<30", "tring_obs:<40", "tring_obs:<55", "tring_obs:<75")

stopifnot(all(expected_staging %in% threshold_full))
cat("\nAll candidate threshold checks passed.\n")

# ---- Function to bin one imputation ----
bin_one_imputation <- function(dt, imp_number) {

  # Select one imputation
  d <- copy(dt[.imp == imp_number])

  # Create output table with the correct number of rows
  # Retain ID as metadata for the common train/test patient split
  out <- data.table(.row = seq_len(nrow(d)), ID = d$ID)

  for (i in seq_len(num_thresholds)) {
    x <- d[[threshold_names[i]]]

    if (threshold_symbols[i] == "<") {
      v <- as.integer(x < threshold_values[i])

    } else if (threshold_symbols[i] == ">") {
      v <- as.integer(x > threshold_values[i])

    } else {
      v <- as.integer(x == threshold_values[i])
    }
    v[is.na(x)] <- NA_integer_
    out[[threshold_full[i]]] <- v
  }

  # Remove temporary row index
  out[, .row := NULL]

  # Mortality outcome is the target variable
  out[, outcome := d$outcome]
  out
}

# ---- Bin and save all 20 imputations ----
imp_ids <- sort(unique(my_table$.imp))

for (imp in imp_ids) {
  binned_imp <- bin_one_imputation(my_table, imp)
  out_file <- file.path(OUT_DIR, sprintf("Imputed_dataset_binned_%02d.csv", imp))
  fwrite(binned_imp, out_file)
  cat("Saved imputation ", imp, ": ", nrow(binned_imp), " rows\n", sep = "")
}

# ---- Final summary ----
cat("\nBinning complete.\n")
cat("Number of imputations: ", length(imp_ids), "\n", sep = "")
cat("Candidate predictors: ", num_thresholds, "\n", sep = "")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")