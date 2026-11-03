# CONSTRUCT SEQUESTRATION WEIGHTS FROM OBSERVED STAGING

library(data.table)

options(scipen = 999)

# ---- State directory paths ----
DATA_DIR <- file.path("Data")
DATA_FILE <- file.path(DATA_DIR, "dat_all.csv")
OUT_DIR <- file.path(DATA_DIR, "Sequestration_weights")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Define microscopy stages ----
stage_cols <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")

stage_names <- c("Tiny rings", "Small rings", "Large rings", "Early trophozoites",
                 "Mid trophozoites", "Late trophozoites", "Schizonts")

stage_breaks <- c(0, 6, 16, 26, 30, 34, 38, 48)
stage_duration <- diff(stage_breaks)

# Reference proportions under uniform underlying parasite age
reference_proportion <- stage_duration / 48

# Final circulation weights used in Stan
final_weights <- c(1.000, 0.500, 0.250, 0.200, 0.035, 0.027, 0.005)

stopifnot(final_weights[1] == 1)
stopifnot(all(diff(final_weights) <= 0))

# ---- Load observed staging dataset ----
dat <- fread(DATA_FILE)

if ("Trial" %in% names(dat) && !"trial" %in% names(dat)) {
  setnames(dat, "Trial", "trial")
}

# ---- Prepare observed staging ----
prepare_staging <- function(data) {
  required_cols <- c("trial", "ID", stage_cols)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  d <- copy(data[, ..required_cols])

  # Keep AQUAMAT and SEAQUAMAT
  d <- d[trial %in% c("AQUAMAT", "SEAQUAMAT")]

  # Remove missing staging
  d <- d[complete.cases(d)]

  # Ensure one row per patient
  stopifnot(!anyDuplicated(d, by = c("trial", "ID")))

  # Calculate total staged cells
  d[, total_stage := rowSums(.SD), .SDcols = stage_cols]

  # Apply established staging inclusion criteria
  d <- d[total_stage >= 93 & total_stage <= 107]

  # Convert stage counts to within-patient proportions
  for (s in stage_cols) {
    d[, (s) := get(s) / total_stage]
  }
  return(d)
}

dat <- prepare_staging(dat)

# ---- Calculate raw circulation weights ----
calculate_weights <- function(data, trial_name) {

  if (trial_name == "Combined") {
    d <- copy(data)
  } else {
    d <- copy(data[trial == trial_name])
  }

  # Mean observed peripheral stage proportions
  observed <- unlist(d[, lapply(.SD, mean), .SDcols = stage_cols], use.names = FALSE)

  # Compare observed proportions with uniform latent proportions
  raw_ratio <- observed / reference_proportion

  # Normalise relative to tiny rings
  raw_weights <- raw_ratio / raw_ratio[1]

  # Apply final weights to uniform underlying stage proportions
  predicted <- reference_proportion * final_weights
  predicted <- predicted / sum(predicted)

  # Return results
  data.table(Trial = trial_name, N = nrow(d), Stage = stage_names,
             Duration = stage_duration, Reference_pct = 100 * reference_proportion,
             Observed_pct = 100 * observed, Raw_weight = raw_weights,
             Final_weight = final_weights, Weight_difference = final_weights - raw_weights,
             Predicted_pct = 100 * predicted,
             Prediction_difference_pp = 100 * (predicted - observed))
}

# ---- Calculate trial-specific and combined weights ----
results <- rbindlist(list(calculate_weights(dat, "AQUAMAT"),
                          calculate_weights(dat, "SEAQUAMAT"),
                          calculate_weights(dat, "Combined")))

# ---- Check monotonicity ----
monotonicity <- results[, .(N = first(N), Tiny_weight = first(Raw_weight),
                            Non_increasing = all(diff(Raw_weight) <= 0),
                            Increasing_transitions = sum(diff(Raw_weight) > 0)), by = Trial]

# ---- Compare raw and final combined weights ----
comparison <- results[Trial == "Combined", .(Stage, Raw_weight, Final_weight,
                                             Weight_difference, Prediction_difference_pp)]

comparison[, Relative_difference_pct := 100 * (Final_weight / Raw_weight - 1)]

# ---- Print results ----
cat("\nCohort sizes\n")
print(results[, .(N = first(N)), by = Trial])

cat("\nCirculation weights\n")
print(results[, .(Trial, Stage, Reference_pct, Observed_pct, Raw_weight, Final_weight,
                  Weight_difference)], digits = 6)

cat("\nMonotonicity checks\n")
print(monotonicity)

cat("\nRaw versus final combined weights\n")
print(comparison, digits = 6)

# ---- Save outputs ----
fwrite(results, file.path(OUT_DIR, "Sequestration_weights.csv"))
fwrite(monotonicity, file.path(OUT_DIR, "Sequestration_monotonicity.csv"))
fwrite(comparison, file.path(OUT_DIR, "Raw_vs_final_sequestration_weights.csv"))

cat("\nAnalysis completed.\n")
cat("Results saved to:\n")
cat(OUT_DIR, "\n")