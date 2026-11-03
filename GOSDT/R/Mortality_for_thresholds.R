library(data.table)
library(dplyr)
library(mgcv)

# State directory paths
OUT_DIR <- "Data"

# Read an external RDS file
# It contains Generalised Additive Models (GAMs) that capture non-linear relationships
# between predictors (stated below) and an outcome (mortality)
obj <- readRDS(file.path(OUT_DIR, "Thresholds_for_mortality.rds"))

# Extract stored objects from the RDS file
Data <- obj$data
marker_specs <- obj$marker_specs
mu_age_model <- obj$mu_age_model
other_models <- obj$other_models
mu_age_threshold <- obj$mu_age_threshold
threshold_summary <- obj$threshold_summary

# State threshold values of the predictors at which the mortality risk will be predicted
threshold_grid <- tibble::tribble(~marker, ~label, ~value,
                                  "mu_age", "Mean age", 5,
                                  "mu_age", "Mean age", 10,
                                  "mu_age", "Mean age", 15,
                                  "mu_age", "Mean age", 20,
                                  "tring_obs", "Tiny rings", 30,
                                  "tring_obs", "Tiny rings", 40,
                                  "tring_obs", "Tiny rings", 55,
                                  "tring_obs", "Tiny rings", 75,
                                  "hb", "Hb", 5,
                                  "hb", "Hb", 8,
                                  "hb", "Hb", 11,
                                  "temp", "Temperature", 37,
                                  "temp", "Temperature", 39,
                                  "temp", "Temperature", 40,
                                  "be", "Base excess", -11,
                                  "be", "Base excess", -9,
                                  "be", "Base excess", -6,
                                  "be", "Base excess", -3)

# Categorise mu_age into two levels, using the threshold to separate them
Data <- Data %>%
  mutate(mu_age_cat = factor(ifelse(mu_age < mu_age_threshold$threshold_5pct, "below",
                                    "above"),
                             levels = c("below", "above")))

# ---- Function 1: Get the reference (site, trial) levels ----
get_ref_levels <- function(model, data) {
  mf <- tryCatch(model.frame(model), error = function(e) NULL)

  # Get the site levels from the GAM if possible
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

# ---- Function 2: Calculate weights for the two mu_age groups ----
mu_age_cat_weights <- function(data, coma = 0) {
  w <- prop.table(table(data$mu_age_cat[data$coma == coma]))
  w <- w[levels(data$mu_age_cat)]
  w[is.na(w)] <- 0
  as.numeric(w)
}

# ---- Function 3: Build a small dataset for making predictions ----
build_newdata <- function(data, marker, value, coma, site_ref, trial_ref, age_cat = NULL) {

  # Set the patient characteristics used for prediction
  nd <- data.frame(coma = factor(coma, levels = levels(data$coma)),
                   site = factor(site_ref, levels = levels(data$site)),
                   Trial = factor(trial_ref, levels = levels(data$Trial)))

  # Add the mu_age categories if the current variable is not mu_age
  if (!is.null(age_cat)) {
    nd$mu_age_cat <- factor(age_cat, levels = levels(data$mu_age_cat))
  }

  # Add each value for each variable
  nd[[marker]] <- value
  nd
}

# ---- Function 4: Estimate mortality risk for a threshold ----
pooled_risk_at_threshold <- function(model, data, marker, value, coma = 0) {

  # Call Function 1 (get the reference levels)
  refs <- get_ref_levels(model, data)

  # Get the two mu_age groups
  mu_age_cat_levels <- levels(data$mu_age_cat)

  # Call Function 2 (calculate weights for the two mu_age groups)
  w <- mu_age_cat_weights(data, coma)

  # Predict the risk for each mu_age group
  risks <- sapply(mu_age_cat_levels, function(a) {

    # Call Function 3 (build a small prediction dataset)
    nd <- build_newdata(data, marker, value, coma, refs$site_ref, refs$trial_ref,
                        age_cat = a)

    # Predict mortality risk
    as.numeric(predict(model, newdata = nd, type = "response"))
  })

  # Average the two mu_age groups using the weights
  as.numeric(sum(risks * w) / sum(w))
}

# Loop through each threshold value and predict the mortality risk
result_table <- bind_rows(lapply(seq_len(nrow(threshold_grid)), function(i) {

  # Extract the variable name, threshold value, and display label
  m <- threshold_grid$marker[i]
  v <- threshold_grid$value[i]
  lab <- threshold_grid$label[i]

  # If the current variable is mu_age, predict directly from the mu_age model
  if (m == "mu_age") {
    r <- as.numeric(predict(mu_age_model,
                            newdata = data.frame(coma = factor(0,
                                                               levels = levels(Data$coma)),
                                                 site = factor(levels(Data$site)[1],
                                                               levels = levels(Data$site)),
                                                 Trial = factor(levels(Data$Trial)[1],
                                                                levels = levels(Data$Trial)),
                                                 mu_age = v), type = "response"))

  # Otherwise, call Function 4 (estimate mortality risk for a threshold)
  } else {
    r <- pooled_risk_at_threshold(other_models[[m]], Data, m, v, coma = 0)
  }

  # Store the predicted mortality risk for the current threshold
  data.frame(marker = m, label = lab, threshold_value = v, coma = 0,
             mortality_risk_pct = round(r * 100, 1))})) %>%

  # Order the results by variable and threshold value
  arrange(marker, threshold_value)

# Write the result table to a CSV file
write.csv(result_table, file.path(OUT_DIR, "Mortality_for_thresholds.csv"),
          row.names = FALSE)

result_table