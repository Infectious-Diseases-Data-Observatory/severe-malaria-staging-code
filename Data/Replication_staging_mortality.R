# REPLICATION OF STAGING-MORTALITY ASSOCIATIONS

library(data.table)

# ---- State directory paths ----
BASE_DIR <- "Data"
DATA_FILE <- file.path(BASE_DIR, "dat_all.csv")
OUT_DIR <- file.path(BASE_DIR, "Replication_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Read observed patient data ----
dat <- fread(DATA_FILE)

# State variables required for the analysis
stage_cols <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")
required_cols <- c("ID", "Trial", "outcome", stage_cols)

# Check that all required variables are present
stopifnot(all(required_cols %in% names(dat)))

# Retain required columns
dat <- copy(dat[, ..required_cols])

# Keep AQUAMAT and SEAQUAMAT
dat <- dat[Trial %in% c("AQUAMAT", "SEAQUAMAT")]

# Remove patients with missing staging or outcome values
dat <- dat[complete.cases(dat)]

# Check mortality coding and patient identifiers
stopifnot(all(dat$outcome %in% c(0, 1)))
stopifnot(!anyDuplicated(dat[, .(Trial, ID)]))

# ---- Apply staging inclusion criteria ----
# Calculate total number of staged parasites
dat[, total_stage := rowSums(.SD), .SDcols = stage_cols]

# Retain patients whose stage percentages sum to 93-107%
dat <- dat[total_stage >= 93 & total_stage <= 107]

# Convert each patient's observed stage counts to percentages
for (s in stage_cols) {
  dat[, (paste0(s, "_pct")) := 100 * get(s) / total_stage]
}

# ---- Calculate cohort mortality ----
cohort_summary <- dat[, .(N = .N, Deaths = sum(outcome == 1),
                          Survivors = sum(outcome == 0),
                          Mortality_pct = round(100 * mean(outcome == 1), 2)), by = Trial]
print(cohort_summary)

write.csv(cohort_summary, file.path(OUT_DIR, "Cohort_mortality_summary.csv"),
          row.names = FALSE)

# ---- Calculate observed stage distributions by outcome ----
stage_pct_cols <- paste0(stage_cols, "_pct")

# Calculate mean within-patient stage proportions
stage_summary <- dat[, lapply(.SD, function(x) round(mean(x), 3)),
                     by = .(Trial, outcome), .SDcols = stage_pct_cols]

# Label mortality outcomes
stage_summary[, Outcome := ifelse(outcome == 1, "Deaths", "Survivors")]

# Reorder columns
setcolorder(stage_summary, c("Trial", "Outcome", "outcome", stage_pct_cols))
print(stage_summary)

write.csv(stage_summary, file.path(OUT_DIR, "Observed_stage_distribution_by_outcome.csv"),
          row.names = FALSE)

# ---- Calculate historical staging measures ----
# Silamut and White (1993) defined mature stages as pigmented trophozoites and schizonts
# (approximately 26-48 h).
# These correspond to early, mid and late trophozoites plus schizonts in the current
# seven-stage classification.
dat[, tiny_rings_pct := tring_pct]
dat[, mature_stages_pct := etroph_pct + mtroph_pct + ltroph_pct + schizont_pct]

# Categorise patients using the historical thresholds
dat[, tiny_rings_gt50 := tiny_rings_pct > 50]
dat[, mature_stages_gt20 := mature_stages_pct > 20]

# ---- Function 1: Calculate mortality risk ratios and 95% confidence intervals ----
calculate_rr_ci <- function(data, threshold_col, threshold_label) {

  # Calculate mortality risks within each trial
  data[, {

    # Identify patients above and at or below the threshold
    exposed <- get(threshold_col)

    # Count deaths and patients above the threshold
    a <- sum(outcome[exposed] == 1)
    n1 <- sum(exposed)

    # Count deaths and patients at or below the threshold
    c <- sum(outcome[!exposed] == 1)
    n0 <- sum(!exposed)

    # Calculate mortality risks
    risk1 <- if (n1 > 0) a / n1 else NA_real_
    risk0 <- if (n0 > 0) c / n0 else NA_real_

    # Calculate risk ratio
    rr <- if (!is.na(risk0) && risk0 > 0) {
      risk1 / risk0
    } else {
      NA_real_
    }

    # Calculate Katz standard error of the log risk ratio
    if (a > 0 && c > 0 && n1 > 0 && n0 > 0) {
      se_log_rr <- sqrt(1/a - 1/n1 + 1/c - 1/n0)

      # Calculate 95% confidence interval
      ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
      ci_upper <- exp(log(rr) + 1.96 * se_log_rr)

    } else {
      ci_lower <- NA_real_
      ci_upper <- NA_real_
    }

    # Return mortality counts, risks and confidence intervals
    .(Threshold = threshold_label,
      N_above = n1,
      Deaths_above = a,
      N_not_above = n0,
      Deaths_not_above = c,
      Mortality_above_pct = round(100 * risk1, 2),
      Mortality_not_above_pct = round(100 * risk0, 2),
      Risk_ratio = round(rr, 3),
      CI_lower = round(ci_lower, 3),
      CI_upper = round(ci_upper, 3))

  }, by = Trial]
}

# ---- Calculate risk ratios for both historical thresholds ----
# Tiny rings >50%
rr_tiny <- calculate_rr_ci(dat, "tiny_rings_gt50", "Tiny rings >50%")

# Mature trophozoites and schizonts >20%
rr_mature <- calculate_rr_ci(dat, "mature_stages_gt20",
                             "Mature trophozoites and schizonts >20%")

# Combine threshold results
results <- rbindlist(list(rr_tiny, rr_mature))

# Create formatted risk ratios with 95% confidence intervals
results[, RR_95CI := sprintf("%.2f (%.2f-%.2f)", Risk_ratio, CI_lower, CI_upper)]

# Order results by trial and threshold
setorder(results, Trial, Threshold)

# Print and save results
print(results)

write.csv(results, file.path(OUT_DIR, "Historical_staging_risk_ratios_95CI_corrected.csv"),
          row.names = FALSE)

# ---- Function 2: Create mortality summaries for each threshold group ----
create_threshold_summary <- function(data) {

  # Create one row for each group above the threshold
  above <- data[, .(Trial, Threshold, Threshold_group = "Above threshold", N = N_above,
                    Deaths = Deaths_above, Survivors = N_above - Deaths_above,
                    Mortality_pct = Mortality_above_pct)]

  # Create one row for each group at or below the threshold
  below <- data[, .(Trial, Threshold, Threshold_group = "At or below threshold",
                    N = N_not_above, Deaths = Deaths_not_above,
                    Survivors = N_not_above - Deaths_not_above,
                    Mortality_pct = Mortality_not_above_pct)]

  # Combine the two groups
  result <- rbindlist(list(above, below))

  # Order results
  setorder(result, Trial, Threshold, Threshold_group)

  return(result)
}

# ---- Calculate threshold group summaries ----
threshold_results <- create_threshold_summary(results)
print(threshold_results)

write.csv(threshold_results,
          file.path(OUT_DIR, "Historical_staging_thresholds_corrected.csv"),
          row.names = FALSE)

# ---- Final summary ----
cat("\nReplication analysis completed.\n")
cat("Total patients:", nrow(dat), "\n")
cat("AQUAMAT:", sum(dat$Trial == "AQUAMAT"), "\n")
cat("SEAQUAMAT:", sum(dat$Trial == "SEAQUAMAT"), "\n")

cat("\nResults saved to:\n")
cat(OUT_DIR, "\n")