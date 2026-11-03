library(data.table)
library(dplyr)
library(mgcv)

# State directory paths and file names
DATA_DIR <- file.path("Staging/Complete case/outputs/CSV")
SOURCE_FILE <- file.path(DATA_DIR, "dat_all_hierarchical.csv")
OUT_DIR <- file.path("Data")

# Create the output folder if needed
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read the source dataset
Data <- fread(SOURCE_FILE) %>% as.data.frame()

# Drop the index column if present
if ("V1" %in% names(Data)) {
  Data$`V1` <- NULL
}

# ---- Function 1: Standardise the binary variables ----
coerce_binary01 <- function(x) {
  x_chr <- trimws(as.character(x))
  as.integer(x_chr %in% c("1", "Yes", "YES", "yes", "TRUE", "True", "true", "T"))
}

# Call Function 1 (standardise outcome and coma for source dataset)
Data <- Data %>%
  mutate(outcome = coerce_binary01(outcome),
         coma = factor(coerce_binary01(coma), levels = c(0, 1)))

# Ensure the outcome column has integers (0s and 1s)
Data <- Data %>%
  mutate(outcome = as.integer(outcome),
         site = factor(site),
         Trial = factor(Trial))

# State variables included in the analysis
analysis_vars <- c("mu_age", "sigma_age", "tring_obs", "mts_obs", "age", "hb", "temp",
                   "be", "coma", "site", "Trial", "outcome")

# Keep only complete cases (patients with no missing values)
Data <- Data %>%
  filter(complete.cases(across(all_of(analysis_vars))))

# State specifications of variables
marker_specs <- tibble::tribble(~marker, ~label,
                                "mu_age", "Mean age",
                                "sigma_age", "SD age",
                                "tring_obs", "Tiny rings",
                                "mts_obs", "MTS",
                                "age", "Patient age",
                                "hb", "Hb",
                                "temp", "Temperature",
                                "be", "Base excess")

# State the observed range where each threshold will be found
# Keep the biologically defined 0-48 h range for mu_age
threshold_windows <- list(mu_age = c(0, 48),
                          sigma_age = range(Data$sigma_age, na.rm = TRUE),
                          tring_obs = range(Data$tring_obs, na.rm = TRUE),
                          mts_obs = range(Data$mts_obs, na.rm = TRUE),
                          age = range(Data$age, na.rm = TRUE),
                          hb = range(Data$hb, na.rm = TRUE),
                          temp = range(Data$temp, na.rm = TRUE),
                          be = range(Data$be, na.rm = TRUE))

# ---- Function 3: Get the reference (site, trial) levels ----
get_ref_levels <- function(model, data) {
  mf <- tryCatch(model.frame(model), error = function(e) NULL)

  # Get the site levels from the GAM (Generalised Additive Model) if possible
  site_ref <- if (!is.null(mf) && "site" %in% names(mf)) {
    levels(droplevels(mf$site))[1]

  # Otherwise get them from the data
  } else {
    levels(data$site)[1]
  }

  # Get the trial levels from the GAM if possible
  trial_ref <- if (!is.null(mf) && "Trial" %in% names(mf)) {
    levels(droplevels(mf$Trial))[1]

  # Otherwise get them from the data
  } else {
    levels(data$Trial)[1]
  }

  # Return the site and trial levels together
  list(site_ref = site_ref, trial_ref = trial_ref)
}

# ---- Function 4: Return the search range for a variable ----
marker_range <- function(marker) {
  threshold_windows[[marker]]
}

# ---- Function 5: Build a small dataset for making predictions ----
build_newdata <- function(data, marker, x, coma = 0, site_ref, trial_ref, age_cat = NULL) {

  # Set the patient characteristics used for prediction
  nd <- data.frame(coma = factor(coma, levels = levels(data$coma)),
                   site = factor(site_ref, levels = levels(data$site)),
                   Trial = factor(trial_ref, levels = levels(data$Trial)))

  # Add the mu_age categories if the current variable is not mu_age
  if (!is.null(age_cat)) {
    nd$mu_age_cat <- factor(age_cat, levels = levels(data$mu_age_cat))
  }

  # Add the value for each variable
  # If a variable has more than one value then duplicate the row
  # And add the new value to the new row
  nd <- nd[rep(1, length(x)), , drop = FALSE]
  nd[[marker]] <- x
  nd
}

# ---- Function 6: Build a mortality curve for mu_age ----
pooled_curve_mu_age <- function(model, data, marker, coma = 0, n = 5000) {

  # Call Function 3 (get the reference levels)
  refs <- get_ref_levels(model, data)

  # Call Function 4 (get the variable search range)
  rng <- marker_range(marker)

  # Make a fine grid of values
  x <- seq(rng[1], rng[2], length.out = n)

  # Call Function 5 (build a small prediction dataset)
  nd <- build_newdata(data, marker, x, coma = coma, site_ref = refs$site_ref,
                      trial_ref = refs$trial_ref)

  # Predict mortality risk
  risk <- as.numeric(predict(model, newdata = nd, type = "response"))

  # Return the mortality curve
  data.frame(marker = marker, x = x, coma = coma, risk = risk)
}

# ---- Function 7: Find where the mortality curve reaches 5% ----
find_crossings <- function(x, p, target = 0.05) {

  # Compare the predicted mortality with the target value
  s <- p - target

  roots <- numeric(0)

  # Keep any values that are already within a small tolerance of the target
  exact_hits <- which(abs(s) < 1e-8)
  if (length(exact_hits) > 0) {
    roots <- c(roots, x[exact_hits])
  }

  # Find where the mortality curve passes from one side of the target to the other
  idx <- which(s[-length(s)] * s[-1] < 0)

  # Estimate the associated roots
  if (length(idx) > 0) {
    roots <- c(roots, mapply(function(i) {
      x[i] + (target - p[i]) * (x[i + 1] - x[i]) / (p[i + 1] - p[i])}, idx))
  }

  # Return the roots in order
  sort(unique(round(roots, 8)))
}

# ---- Function 8: Fit one GAM for one variable ----
fit_marker_gam <- function(data, marker, use_age_cat = FALSE) {

  # Allow the relationship to differ between the two mu_age groups
  if (use_age_cat) {
    fml <- as.formula(paste0("outcome ~ mu_age_cat + s(", marker,
                             ", bs = 'ts', k = 20, by = mu_age_cat) + ",
                             "s(site, bs = 're') + s(Trial, bs = 're') + coma"))

  # Otherwise fit a single relationship for the whole dataset
  } else {
    fml <- as.formula(paste0("outcome ~ s(", marker, ", bs = 'ts', k = 20) + ",
                             "s(site, bs = 're') + s(Trial, bs = 're') + coma"))
  }

  # Fit the model
  gam(fml, family = binomial(link = "logit"), data = data, select = TRUE, method = "REML",
      na.action = na.exclude)
}

# ---- Function 9: Find the 5% threshold for mu_age ----
find_5pct_threshold_mu_age <- function(model, data, marker, coma = 0, target = 0.05) {

  # Call Function 6 (build a mortality curve for mu_age)
  curve <- pooled_curve_mu_age(model, data, marker, coma = coma, n = 5000)

  # Call Function 7 (find where the mortality curve crosses the target)
  roots <- find_crossings(curve$x, curve$risk, target = target)

  # Use the first crossing when one exists
  if (length(roots) > 0) {
    threshold <- roots[1]

  # Otherwise no 5% mortality threshold exists
  } else {
    threshold <- NA_real_
  }

  # Return the result
  data.frame(marker = marker, label = "Mean age", coma = coma, threshold_5pct = threshold)
}

# 1) Find the threshold for the variable mu_age first
# Call Function 8 (fit GAM for mu_age)
mu_age_model <- fit_marker_gam(Data, "mu_age", use_age_cat = FALSE)

# Call Function 9 (find the 5% threshold for mu_age)
mu_age_threshold <- find_5pct_threshold_mu_age(mu_age_model, Data, "mu_age")

# Categorise mu_age into two levels, using the threshold to separate them
Data <- Data %>%
  mutate(mu_age_cat = factor(ifelse(mu_age < mu_age_threshold$threshold_5pct, "below",
                                    "above"),
                             levels = c("below", "above")))

# 2) Fit the remaining thresholds
other_markers <- c("sigma_age", "tring_obs", "mts_obs", "age", "hb", "temp", "be")

other_models <- list()

# Loop through the remaining variables and fit one model at a time
for (m in other_markers) {

  # Keep only the rows needed for this variable
  vars_needed <- c("outcome", "coma", "site", "Trial", "mu_age_cat", m)
  cc <- complete.cases(Data[, vars_needed])

  # Call Function 8 (fit GAM for each variable)
  other_models[[m]] <- fit_marker_gam(Data[cc, ], m, use_age_cat = TRUE)
}

# Keep the models in the same order as in the variable list
other_models <- other_models[other_markers]

# ---- Function 10: Calculate weights for the two mu_age groups ----
mu_age_cat_weights <- function(data, coma = 0) {
  w <- prop.table(table(data$mu_age_cat[data$coma == coma]))
  w <- w[levels(data$mu_age_cat)]
  w[is.na(w)] <- 0
  as.numeric(w)
}

# ---- Function 11: Build a mortality curve for each of the other variables ----
pooled_curve <- function(model, data, marker, n = 5000) {

  # Call Function 3 (get the reference levels)
  refs <- get_ref_levels(model, data)

  # Call Function 4 (get the variable search range)
  rng <- marker_range(marker)

  # Make a fine grid of values
  x <- seq(rng[1], rng[2], length.out = n)

  # Get the two mu_age groups
  mu_age_cat_levels <- levels(data$mu_age_cat)

  # Call Function 10 (calculate weights for the two mu_age groups)
  w <- mu_age_cat_weights(data, coma = 0)

  # Predict the risk for each mu_age group
  risk_mat <- sapply(mu_age_cat_levels, function(a) {

    # Call Function 5 (build a small prediction dataset)
    nd <- build_newdata(data, marker = marker, x = x, coma = 0,
                        site_ref = refs$site_ref, trial_ref = refs$trial_ref, age_cat = a)

    # Predict mortality risk
    as.numeric(predict(model, newdata = nd, type = "response"))
  })

  # Keep the result as a matrix when only one mu_age group is present
  if (is.null(dim(risk_mat))) {
    risk_mat <- matrix(risk_mat, ncol = 1)
  }

  # Average the two mu_age groups using the weights
  pooled_risk <- as.numeric(risk_mat %*% w / sum(w))
  data.frame(marker = marker, x = x, risk = pooled_risk)
}

# ---- Function 12: Find the 5% threshold for each of the other variables ----
find_5pct_threshold <- function(model, data, marker, target = 0.05) {

  # Call Function 11 (build a mortality curve)
  curve <- pooled_curve(model, data, marker, n = 5000)

  # Call Function 7 (find where the mortality curve crosses the target)
  roots <- find_crossings(curve$x, curve$risk, target = target)

  # Use the first crossing when one exists
  if (length(roots) > 0) {
    threshold <- roots[1]

  # Otherwise mark the threshold as missing
  } else {
    threshold <- NA_real_
  }

  # Return the result
  data.frame(marker = marker, label = marker_specs$label[match(marker,
                                                               marker_specs$marker)],
             coma = 0, threshold_5pct = threshold)
}

# Find thresholds for all remaining variables
other_thresholds <- bind_rows(lapply(other_markers, function(m) {

  # Call Function 12
  find_5pct_threshold(other_models[[m]], Data, m)})) %>%

  # Keep the remaining variables in order
  select(marker, label, coma, threshold_5pct) %>%
  arrange(match(marker, c("mu_age", other_markers)))

# Combine the mu_age threshold with the other thresholds
threshold_summary <- bind_rows(mu_age_threshold, other_thresholds)

# Save the threshold table
write.csv(threshold_summary, file.path(OUT_DIR, "Thresholds_coma0.csv"), row.names = FALSE)

# Save the fitted models and results for later use
saveRDS(list(data = Data, marker_specs = marker_specs, mu_age_model = mu_age_model,
             other_models = other_models, mu_age_threshold = mu_age_threshold,
             threshold_summary = threshold_summary),
        file.path(OUT_DIR, "Thresholds_for_mortality.rds"))

threshold_summary