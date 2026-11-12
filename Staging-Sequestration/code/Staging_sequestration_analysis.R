# STAGING-SEQUESTRATION ANALYSIS
# Parasite biomass is estimated using Hendriksen et al. (2012)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Install mice: install.packages('mice')")
}

# Input and output locations
ANALYSIS_DIR <- file.path("Staging-Sequestration")

CSV_OUT <- file.path(ANALYSIS_DIR, "outputs", "CSV")
FIG_OUT <- file.path(ANALYSIS_DIR, "outputs", "figures")
dir.create(CSV_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

read_input <- function(...) {
  fread(file = file.path(...))
}

assert <- function(ok, msg) {
  if (!isTRUE(ok)) {
    stop(msg, call. = FALSE)
  }
}

# Check that each patient has a unique identifier
unique_keys <- function(d, keys, label) {
  assert(all(keys %in% names(d)), paste(label, "missing key column"))
  assert(!anyNA(d[, ..keys]) && anyDuplicated(d, by = keys) == 0L,
         paste(label, "has missing or duplicate keys"))
}

# Published biomass conversion factor and input units
PFHRP2_UNIT <- "ng/mL"
PARASITAEMIA_UNIT <- "parasites/uL"

BIOMASS_K <- 7.3e7

# Load and prepare AQUAMAT staging data
STAGES <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")

# Use six stage proportions, with tring as the reference stage
STAGE_PREDICTORS <- setdiff(STAGES, "tring")

obs <- read_input("Data", "dat_all.csv")[Trial == "AQUAMAT"]
obs[, ID := as.character(ID)]
unique_keys(obs, "ID", "Observed AQUAMAT")
obs[, (STAGES) := lapply(.SD, as.numeric), .SDcols = STAGES]

# Calculate the total percentage across the seven parasite stages
obs[, stage_sum := rowSums(.SD), .SDcols = STAGES]
obs[, stage_complete := rowSums(!is.finite(as.matrix(.SD))) == 0L, .SDcols = STAGES]

# Retain patients with stage percentages summing to 93-107%
obs[, stage_eligible := stage_complete & is.finite(stage_sum) &
      stage_sum >= 93 & stage_sum <= 107]

# Convert observed stage percentages into proportions
for (s in STAGES) obs[, (paste0("p_", s)) := fifelse(stage_eligible, get(s) / stage_sum,
                                                     NA_real_)]
STAGE_PREDICTORS <- paste0("p_", STAGE_PREDICTORS)

hrp <- obs[is.finite(hrp2_plasma) & hrp2_plasma >= 0,
           c("ID", "hrp2_plasma", "stage_eligible", STAGE_PREDICTORS), with = FALSE]

imp <- read_input("Data", "Imputed_dataset_long.csv")
imp[, ID := as.character(ID)]
unique_keys(imp, c("ID", ".imp"), "Imputed clinical data")
assert(uniqueN(imp$.imp) == 20L, "Expected exactly 20 imputations")
stage_imp <- read_input("Staging", "Imputed", "outputs", "CSV", "stage_proportions.csv")
stage_imp[, ID := as.character(ID)]
unique_keys(stage_imp, "ID", "Imputed staging export")

stage_imp <- stage_imp[trial == "AQUAMAT", .(ID, mu_age, sigma_age, model)]

# Match observed PfHRP2 measurements to the imputed staging data
primary <- merge(merge(imp, hrp, by = "ID"), stage_imp, by = "ID")

if ("trial" %in% names(primary)) {
  primary <- primary[trial == "AQUAMAT"]
}

unique_keys(primary, c("ID", ".imp"), "Primary cohort")

assert(all(primary[, .(n = uniqueN(.imp)), by = ID]$n == 20L),
       "Incomplete primary imputations")

assert(all(is.finite(primary$hct) & primary$hct > 0 & primary$hct < 100 &
           is.finite(primary$parasitaemia) & primary$parasitaemia > 0 &
           is.finite(primary$mu_age) & is.finite(primary$sigma_age)),
       "Invalid primary biomass/staging inputs")

assert(!any(primary[, .(vary = uniqueN(mu_age) > 1 | uniqueN(sigma_age) > 1 |
                               uniqueN(hrp2_plasma) > 1 |
                               uniqueN(model) > 1), by = ID]$vary),
       "Observed measurements vary across imputations")

stage_cc <- read_input("Staging", "Complete case", "outputs", "CSV",
                       "stage_proportions.csv")

stage_cc[, ID := as.character(ID)]
unique_keys(stage_cc, "ID", "Complete case staging export")
stage_cc <- stage_cc[Trial == "AQUAMAT", .(ID, mu_age, sigma_age, model)]

# Match observed clinical measurements to complete case staging estimates
cc <- merge(obs, stage_cc, by = "ID")

cc <- cc[is.finite(hrp2_plasma) & hrp2_plasma >= 0 &
         is.finite(hct) & hct > 0 & hct < 100 &
         is.finite(parasitaemia) & parasitaemia > 0 &
         is.finite(mu_age) & is.finite(sigma_age)]

unique_keys(cc, "ID", "Complete case cohort")
cc[, .imp := 1L]

# Check staging eligibility once per patient
stage_audit_one <- function(d, population) {
  assert(all(STAGES %in% names(d)), paste(population, "missing original stage columns"))

  x <- unique(d, by = "ID")
  vals <- as.matrix(x[, ..STAGES])
  complete <- rowSums(is.finite(vals)) == length(STAGES)
  sums <- rowSums(vals)
  eligible <- complete & is.finite(sums) & sums >= 93 & sums <= 107
  valid_sums <- sums[complete & is.finite(sums)]

  if (length(valid_sums)) {
    median_sum <- median(valid_sums)
    min_sum <- min(valid_sums)
    max_sum <- max(valid_sums)

  } else {
    median_sum <- NA_real_
    min_sum <- NA_real_
    max_sum <- NA_real_
  }

  data.table(population = population, patients = nrow(x),
             complete_stage_profiles = sum(complete),
             incomplete_stage_profiles = sum(!complete),
             within_93_to_107 = sum(eligible),
             below_93 = sum(complete & sums < 93),
             above_107 = sum(complete & sums > 107),
             median_stage_sum = median_sum,
             min_stage_sum = min_sum,
             max_stage_sum = max_sum)
}

stage_audit <- rbindlist(list(stage_audit_one(obs, "Original AQUAMAT"),
                              stage_audit_one(primary, "Primary matched"),
                              stage_audit_one(cc, "Complete case matched")))

assert(stage_audit[population == "Primary matched",
                   within_93_to_107] == uniqueN(primary$ID),
       "Primary cohort violates 93-107% staging eligibility")

assert(stage_audit[population == "Complete case matched",
                   within_93_to_107] == uniqueN(cc$ID),
       "Complete case cohort violates 93-107% staging eligibility")

fwrite(stage_audit, file.path(CSV_OUT, "stage_composition_audit.csv"))

# Calculate parasite biomass
biomass <- function(d, K) {
  z <- copy(d)

  # Estimate total parasite biomass from plasma PfHRP2
  z[, total_per_kg := K * hrp2_plasma * (1 - hct / 100)]

  # Estimate circulating parasite biomass from peripheral parasitaemia
  z[, circulating_per_kg := 80000 * parasitaemia]

  # Estimate sequestered biomass as total minus circulating biomass
  z[, sequestered_raw_per_kg := total_per_kg - circulating_per_kg]

  # Calculate the sequestered fraction where total biomass is positive
  z[, fraction_raw := fifelse(total_per_kg > 0,
                              sequestered_raw_per_kg / total_per_kg, NA_real_)]

  # Transform biomass estimates for regression analysis
  z[, total_log10 := fifelse(total_per_kg > 0,
                             log10(pmax(total_per_kg, .Machine$double.xmin)), NA_real_)]

  z[, circulating_log10 := log10(circulating_per_kg)]
  z[, sequestered_asinh := asinh(sequestered_raw_per_kg / 1e11)]

  # Represent mean parasite age over the circular 48-hour cycle
  z[, `:=`(age_sin = sin(2*pi*mu_age/48), age_cos = cos(2*pi*mu_age/48),
           sigma_5h = sigma_age/5)]

  assert(all(is.finite(z$circulating_log10) & is.finite(z$sequestered_asinh)),
         "Nonfinite biomass transforms")

  assert(all(is.na(z[hrp2_plasma == 0]$fraction_raw) &
               is.na(z[hrp2_plasma == 0]$total_log10)),
         "Zero PfHRP2 boundary handling failed")
  z
}

# Compare the full and reduced models across imputations
extract_d1 <- function(fits_full, fits_reduced) {
  z <- mice::D1(mice::as.mira(fits_full), mice::as.mira(fits_reduced))$result
  required <- c("F.value", "df1", "df2", "P(>F)")

  assert(is.matrix(z) && nrow(z) == 1L && all(required %in% colnames(z)),
         "Unexpected mice::D1 structure")

  # Extract the test statistic, degrees of freedom and p-value
  data.table(test = "MI D1", statistic = as.numeric(z[1, "F.value"]),
             df1 = as.numeric(z[1, "df1"]), df2 = as.numeric(z[1, "df2"]),
             p_value = as.numeric(z[1, "P(>F)"]))
}

# Compare observed microscopy with parasite-age estimates
fit_comparison <- function(d, outcome_name, population, coefficient, selection,
                           comparison) {

  assert(outcome_name %in% names(d), paste("Missing outcome:", outcome_name))

  # Define the reduced and full regression models
  predictors <- switch(comparison,

                       # A to B: add mean parasite age to observed microscopy
                       "A_to_B" = list(reduced = STAGE_PREDICTORS,
                                       full = c(STAGE_PREDICTORS, "age_sin", "age_cos")),

                       # B to C: additionally include parasite age standard deviation
                       "B_to_C" = list(reduced = c(STAGE_PREDICTORS, "age_sin", "age_cos"),
                                       full = c(STAGE_PREDICTORS, "age_sin", "age_cos",
                                                "sigma_5h")),

                       stop("Unknown comparison"))

  # Identify the measurements required for each model
  needed <- c(outcome_name, predictors$full)
  eligible <- which(rowSums(is.finite(as.matrix(d[, ..needed]))) == length(needed))
  z <- d[eligible]

  if (nrow(z) == 0L) {
    return(NULL)
  }

  # Separate the data into individual imputations
  pieces <- split(z, by = ".imp", keep.by = TRUE)
  f0 <- reformulate(predictors$reduced, response = outcome_name)
  f1 <- reformulate(predictors$full, response = outcome_name)

  # Fit both models separately within each imputation
  reduced <- lapply(pieces, function(part) lm(f0, data = part))
  full <- lapply(pieces, function(part) lm(f1, data = part))

  assert(all(vapply(c(reduced, full), function(mod) !anyNA(coef(mod)), logical(1))),
         paste("Singular microscopy model:", comparison, outcome_name,
               "- inspect stage distributions and consider a simpler baseline"))

  m <- length(full)

  # Pool estimates and compare models for multiply imputed data
  if (m > 1L) {
    tab <- as.data.table(summary(mice::pool(mice::as.mira(full)), conf.int = TRUE))
    tst <- extract_d1(full, reduced)

  } else {

    # Extract regression results for the complete case analysis
    co <- summary(full[[1]])$coefficients
    ci <- confint(full[[1]])

    tab <- data.table(term = rownames(co), estimate = co[, 1], std.error = co[, 2],
                      statistic = co[, 3], p.value = co[, 4], conf.low = ci[, 1],
                      conf.high = ci[, 2])

    a <- anova(reduced[[1]], full[[1]])

    tst <- data.table(test = "Complete case partial F", statistic = a$F[2],
                      df1 = a$Df[2], df2 = a$Res.Df[2], p_value = a$`Pr(>F)`[2])
  }

  # Record the population, outcome and model comparison
  prefix <- data.table(population = population, coefficient = coefficient,
                       selection = selection, comparison = comparison,
                       outcome = outcome_name, patients = uniqueN(z$ID), m = m)

  # Normalise CI column names across mice::pool and complete case lm
  if ("conf.low" %in% names(tab) && "2.5 %" %in% names(tab)) {
    tab[, `2.5 %` := NULL]
  }

  if ("conf.high" %in% names(tab) && "97.5 %" %in% names(tab)) {
    tab[, `97.5 %` := NULL]
  }

  if ("2.5 %" %in% names(tab)) {
    setnames(tab, "2.5 %", "conf.low")
  }

  if ("97.5 %" %in% names(tab)) {
    setnames(tab, "97.5 %", "conf.high")
  }

  assert(all(c("conf.low", "conf.high") %in% names(tab)) &&
           !anyDuplicated(names(tab)), "Missing or duplicate coefficient confidence limits")

  # Return the regression coefficients and model comparison test
  list(coef = cbind(prefix[rep(1L, nrow(tab))], tab), test = cbind(prefix, tst))
}

# Calculate pooled predictions and 95% confidence intervals
mi_curve_prediction <- function(fits, X, dfcom) {
  m <- length(fits)
  assert(m >= 2L, "MI curve requires at least two imputations")
  nm <- names(coef(fits[[1]]))

  assert(all(vapply(fits, function(f) identical(names(coef(f)), nm), logical(1))),
         "Inconsistent curve coefficients across imputations")

  assert(all(nm %in% colnames(X)), "Prediction design matrix is missing terms")

  # Match prediction variables to the fitted model
  X <- X[, nm, drop = FALSE]
  Q <- do.call(rbind, lapply(fits, function(f) as.numeric(coef(f)[nm])))

  # Calculate prediction uncertainty within each imputation
  U <- vapply(fits, function(f) {
    V <- vcov(f)[nm, nm, drop = FALSE]
    pmax(0, rowSums((X %*% V) * X))
  }, numeric(nrow(X)))

  if (is.null(dim(U))) {
    U <- matrix(U, ncol = m)
  }

  # Combine predictions across imputations
  predictions <- X %*% t(Q)
  qbar <- rowMeans(predictions)
  ubar <- rowMeans(U)
  b <- apply(predictions, 1L, var)

  # Combine within- and between-imputation uncertainty
  total_var <- ubar + (1 + 1/m) * b
  se <- sqrt(total_var)

  # Adjust confidence intervals for uncertainty from missing data
  lambda <- ifelse(total_var > 0, (1 + 1/m) * b / total_var, 0)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- ((dfcom + 1) / (dfcom + 3)) * dfcom * (1 - lambda)

  df <- ifelse(is.infinite(df_old), df_obs,
               1 / (1 / df_old + 1 / df_obs))

  # Calculate the critical value for the 95% confidence intervals
  critical <- qt(0.975, df = df)

  data.table(estimate = as.numeric(qbar), lower_95 = as.numeric(qbar - critical * se),
             upper_95 = as.numeric(qbar + critical * se), df = as.numeric(df))
}

# Generate biomass curves across mean parasite age
fit_curve <- function(d, outcome_name) {

  # Keep observations with a biomass outcome
  good <- which(is.finite(d[[outcome_name]]))
  pieces <- split(d[good], by = ".imp", keep.by = TRUE)

  # Fit the regression model in each imputation
  fits <- lapply(pieces, function(part)
    lm(reformulate(c("age_sin", "age_cos", "sigma_5h"), outcome_name), data = part))

  # Predict across the full 48-hour parasite cycle
  ages <- seq(0, 48, by = 0.25)

  # Fix parasite age standard deviation at its median
  sigma_ref <- median(unique(d[.imp == min(.imp), .(ID, sigma_age)], by = "ID")$sigma_age)

  X <- cbind(`(Intercept)` = 1, age_sin = sin(2*pi*ages/48),
             age_cos = cos(2*pi*ages/48), sigma_5h = sigma_ref/5)

  dfcom <- min(vapply(fits, df.residual, numeric(1)))

  # Pool curves and their confidence intervals
  ans <- mi_curve_prediction(fits, X, dfcom)
  cbind(data.table(outcome = outcome_name, age_hours = ages), ans)
}

# Analyse total, circulating and sequestered parasite biomass
OUTCOMES <- c("total_log10", "circulating_log10", "sequestered_asinh")
coef_tables <- list(); test_tables <- list(); count_tables <- list()
summary_tables <- list()

# Analyse the imputed and complete case cohorts separately
for (population in c("Primary_MI", "Complete_case")) {
  if (population == "Primary_MI") {
    source_data <- primary
  } else {
    source_data <- cc
  }

  # Calculate biomass using the Hendriksen coefficient
  coefficient <- "Hendriksen"
  full_data <- biomass(source_data, BIOMASS_K)

  # Compare analyses including and excluding zero PfHRP2 values
  for (selection in c("Nonnegative_PfHRP2", "Positive_only")) {
    if (selection == "Positive_only") {
      z <- full_data[hrp2_plasma > 0]
    } else {
      z <- full_data
    }

    # Summarise biomass estimates within each imputation
    perimp <- z[, .(zeros = sum(hrp2_plasma == 0),
                    negative_sequestered = sum(sequestered_raw_per_kg < 0),
                    median_total = median(total_per_kg),
                    median_circulating = median(circulating_per_kg),
                    median_sequestered = median(sequestered_raw_per_kg)), by = .imp]

    summary_tables[[length(summary_tables) + 1L]] <-
      data.table(population = population, coefficient = coefficient, selection = selection,
                 patients = uniqueN(z$ID), imputations = uniqueN(z$.imp),
                 zeros_per_imputation = mean(perimp$zeros),
                 negative_sequestered_per_imputation = mean(perimp$negative_sequestered),
                 mean_median_total_per_kg = mean(perimp$median_total),
                 mean_median_circulating_per_kg = mean(perimp$median_circulating),
                 mean_median_sequestered_per_kg = mean(perimp$median_sequestered))

    # Analyse each of the three biomass outcomes
    for (outcome_name in OUTCOMES) {
      valid <- z[which(is.finite(z[[outcome_name]]))]

      count_tables[[length(count_tables) + 1L]] <-
        data.table(population = population, coefficient = coefficient,
                   selection = selection, outcome = outcome_name,
                   patients = uniqueN(valid$ID),
                   microscopy_complete_patients = uniqueN(valid$ID),
                   zero_pfhrp2_in_outcome = uniqueN(valid[hrp2_plasma == 0]$ID))

      # Test the additional contribution of inferred parasite age parameters
      for (comparison in c("A_to_B", "B_to_C")) {
        ans <- fit_comparison(z, outcome_name, population, coefficient, selection,
                              comparison)

        assert(!is.null(ans), paste("No model data:", population, outcome_name))
        coef_tables[[length(coef_tables) + 1L]] <- ans$coef
        test_tables[[length(test_tables) + 1L]] <- ans$test
      }
    }
  }
}

assert(length(test_tables) == 24L && length(coef_tables) == 24L,
       "Incomplete regression results")

assert(all(is.finite(rbindlist(test_tables)$p_value)), "Nonfinite nested-model p-values")

# Export biomass summaries, sample counts and regression results
fwrite(rbindlist(summary_tables, fill = TRUE), file.path(CSV_OUT, "biomass_summary.csv"))
fwrite(rbindlist(count_tables, fill = TRUE), file.path(CSV_OUT, "outcome_counts.csv"))
coef_export <- rbindlist(coef_tables, fill = TRUE)
fwrite(coef_export, file.path(CSV_OUT, "pooled_coefficients.csv"))

fwrite(rbindlist(test_tables, fill = TRUE),
       file.path(CSV_OUT, "incremental_model_tests.csv"))

# Summarise patient numbers and zero PfHRP2 measurements
flow <- data.table(population = c("Original AQUAMAT", "Observed nonnegative PfHRP2",
                                  "Primary matched", "Complete case matched"),

                   patients = c(nrow(obs), nrow(hrp), uniqueN(primary$ID), uniqueN(cc$ID)),

                   zero_pfhrp2 = c(sum(obs$hrp2_plasma == 0, na.rm = TRUE),
                                   sum(hrp$hrp2_plasma == 0),
                                   uniqueN(primary[hrp2_plasma == 0]$ID),
                                   uniqueN(cc[hrp2_plasma == 0]$ID)))

fwrite(flow, file.path(CSV_OUT, "cohort_flow.csv"))

# Track patient exclusions at each stage of cohort selection
sequential_flow <- function(population) {
  base <- copy(obs)

  # Identify patients with hierarchical staging estimates
  if (population == "Primary_MI") {
    stage_ids <- stage_imp$ID
  } else {
    stage_ids <- stage_cc$ID
  }

  # Identify patients with valid clinical measurements
  if (population == "Primary_MI") {
    clinical_ids <- imp[is.finite(hct) & hct > 0 & hct < 100 &
                        is.finite(parasitaemia) & parasitaemia > 0,
                        .(n = uniqueN(.imp)), by = ID][n == 20L, ID]
  } else {
    clinical_ids <- obs[is.finite(hct) & hct > 0 & hct < 100 &
                        is.finite(parasitaemia) & parasitaemia > 0, ID]
  }

  # Apply the eligibility criteria in sequence
  ids <- base$ID
  steps <- list(c("All AQUAMAT", length(ids)))
  ids <- intersect(ids, base[is.finite(hrp2_plasma) & hrp2_plasma >= 0, ID])
  steps[[2]] <- c("Nonnegative observed PfHRP2", length(ids))
  ids <- intersect(ids, base[stage_eligible == TRUE, ID])
  steps[[3]] <- c("Stage sum 93-107%", length(ids))
  ids <- intersect(ids, stage_ids)

  if (population == "Primary_MI") {
    stage_label <- "Included in imputed hierarchical staging export"
  } else {
    stage_label <- "Included in complete case hierarchical staging export"
  }

  steps[[4]] <- c(stage_label, length(ids))

  ids <- intersect(ids, clinical_ids)
  steps[[5]] <- c("Valid biomass inputs", length(ids))

  if (population == "Primary_MI") {
    final_n <- uniqueN(primary$ID)
  } else {
    final_n <- uniqueN(cc$ID)
  }

  # Confirm that the final count matches the analytical cohort
  assert(length(ids) == final_n, paste("Sequential flow does not reproduce", population))

  data.table(population = population, step = vapply(steps, `[`, "", 1L),
             patients = as.integer(vapply(steps, function(z) z[[2]], "")),

             excluded_at_step = c(NA_integer_,
                                  -diff(as.integer(vapply(steps, function(z) z[[2]],
                                                          "")))))
}

fwrite(rbindlist(lapply(c("Primary_MI", "Complete_case"), sequential_flow)),
       file.path(CSV_OUT, "sequential_cohort_flow.csv"))

# Compare biomass estimates at selected parasite ages
mi_linear_contrast <- function(fits, delta, dfcom) {
  nm <- names(coef(fits[[1]]))
  delta <- delta[nm]

  # Calculate the difference within each imputation
  q <- vapply(fits, function(f) sum(delta * coef(f)[nm]), numeric(1))

  # Calculate uncertainty in each difference
  u <- vapply(fits, function(f) {
    as.numeric(t(delta) %*% vcov(f)[nm, nm, drop = FALSE] %*% delta)
  }, numeric(1))

  # Combine estimates and uncertainty across imputations
  m <- length(fits)
  within <- mean(u); between <- var(q)
  total <- within + (1 + 1/m) * between

  if (total > 0) {
    lambda <- (1 + 1/m) * between / total
  } else {
    lambda <- 0
  }

  if (lambda > 0) {
    dfold <- (m - 1) / lambda^2
  } else {
    dfold <- Inf
  }

  dfobs <- ((dfcom + 1) / (dfcom + 3)) * dfcom * (1 - lambda)

  if (is.infinite(dfold)) {
    df <- dfobs
  } else {
    df <- 1/(1/dfold + 1/dfobs)
  }

  assert(is.finite(df) && df > 0 && is.finite(total) && total >= 0,
         "Invalid contrast uncertainty")

  # Calculate the pooled difference and 95% confidence interval
  est <- mean(q); half <- qt(.975, df) * sqrt(total)
  c(estimate = est, conf.low = est - half, conf.high = est + half, df = df)
}

age_once <- unique(primary[.imp == min(.imp), .(ID, mu_age, sigma_age)], by = "ID")

# Calculate the mean parasite age while accounting for its circular scale
circular_mean_48 <- function(x) {
  theta <- 2*pi*x/48
  phase <- atan2(mean(sin(theta)), mean(cos(theta)))
  (phase %% (2*pi))*48/(2*pi)
}

contrast_rows <- list()

# Calculate illustrative contrasts for each biomass outcome
for (outcome_name in OUTCOMES) {
  d <- biomass(primary, BIOMASS_K)
  pieces <- split(d[which(is.finite(d[[outcome_name]]))], by = ".imp")

  fits <- lapply(pieces, function(z) lm(reformulate(c("age_sin", "age_cos", "sigma_5h"),
                                                    outcome_name), data = z))

  dfcom <- min(vapply(fits, df.residual, numeric(1)))

  # 30 h versus 6 h, both within observed mean age coverage
  age_delta <- c("(Intercept)" = 0,
                 age_sin = sin(2*pi*30/48) - sin(2*pi*6/48),
                 age_cos = cos(2*pi*30/48) - cos(2*pi*6/48), sigma_5h = 0)

  # 8 h versus 3 h, within the central SD distribution
  sd_delta <- c("(Intercept)" = 0, age_sin = 0, age_cos = 0, sigma_5h = (8-3)/5)

  for (comparison in c("Mean age: 30 h versus 6 h", "Age SD: 8 h versus 3 h")) {
    if (startsWith(comparison, "Mean")) {
      delta <- age_delta
    } else {
      delta <- sd_delta
    }

    z <- mi_linear_contrast(fits, delta, dfcom)
    is_log <- outcome_name %in% c("total_log10", "circulating_log10")

    # Convert difference on the log scale into biomass ratios
    if (is_log) {
      ratio <- 10^z[["estimate"]]
      ratio_low <- 10^z[["conf.low"]]
      ratio_high <- 10^z[["conf.high"]]

    } else {
      ratio <- NA_real_
      ratio_low <- NA_real_
      ratio_high <- NA_real_
    }

    if (is_log) {
      contrast_scale <- "log10 difference"
    } else {
      contrast_scale <- "asinh difference"
    }

    contrast_rows[[length(contrast_rows) + 1L]] <- data.table(
      outcome = outcome_name, comparison = comparison, scale = contrast_scale,
      estimate = z[["estimate"]], conf.low = z[["conf.low"]],
      conf.high = z[["conf.high"]], df = z[["df"]],
      ratio = ratio, ratio_low = ratio_low, ratio_high = ratio_high
    )
  }
}

fwrite(rbindlist(contrast_rows), file.path(CSV_OUT, "effect_size_contrasts.csv"))

# Flexible observed microscopy sensitivity
# Use one observed microscopy profile per patient
microscopy_once <- unique(primary[.imp == min(.imp), c("ID", STAGE_PREDICTORS),
                                  with = FALSE], by = "ID")

P <- as.matrix(microscopy_once[, ..STAGE_PREDICTORS])

# Identify the main patterns in observed parasite stage composition
pca <- prcomp(P, center = TRUE, scale. = FALSE)

# Allow nonlinear relationships with observed microscopy
add_nonlinear_stage_terms <- function(d) {
  z <- copy(d)
  X <- sweep(as.matrix(z[, ..STAGE_PREDICTORS]), 2L, pca$center, "-")
  scores <- X %*% pca$rotation[, 1:2, drop = FALSE]

  # Add cubic terms for the first two microscopy components
  z[, `:=`(stage_pc1_cubic = scores[,1]^3,
           stage_pc2_cubic = scores[,2]^3)]
  z
}

# Repeat model comparisons with a more flexible microscopy baseline
fit_flexible_comparison <- function(d, outcome_name, population, coefficient, selection,
                                    comparison) {

  extra <- c("stage_pc1_cubic", "stage_pc2_cubic")

  # Define models with and without inferred parasite age parameters
  pred <- switch(comparison,

                 A_to_B = list(reduced = c(STAGE_PREDICTORS, extra),
                               full = c(STAGE_PREDICTORS, extra, "age_sin", "age_cos")),

                 B_to_C = list(reduced = c(STAGE_PREDICTORS, extra, "age_sin", "age_cos"),
                               full = c(STAGE_PREDICTORS, extra, "age_sin", "age_cos",
                                        "sigma_5h")))

  need <- c(outcome_name, pred$full)
  z <- d[which(rowSums(is.finite(as.matrix(d[, ..need]))) == length(need))]
  parts <- split(z, by = ".imp")
  f0 <- reformulate(pred$reduced, response = outcome_name)
  f1 <- reformulate(pred$full, response = outcome_name)

  # Fit the reduced and full models within each imputation
  reduced <- lapply(parts, function(x) lm(f0, data = x))
  full <- lapply(parts, function(x) lm(f1, data = x))

  if (any(vapply(c(reduced, full), function(f) anyNA(coef(f)), logical(1)))) {
    return(data.table(population = population, coefficient = coefficient,
                      selection = selection, comparison = comparison,
                      outcome = outcome_name, patients = uniqueN(z$ID),
                      status = "Singular flexible model", statistic = NA_real_,
                      df1 = NA_real_, df2 = NA_real_, p_value = NA_real_))
  }

  # Pool model comparisons for multiply imputed data
  if (length(full) > 1L) {
    test <- extract_d1(full, reduced)

  } else {
    a <- anova(reduced[[1]], full[[1]])
    test <- data.table(statistic = a$F[2], df1 = a$Df[2],
                       df2 = a$Res.Df[2], p_value = a$`Pr(>F)`[2])
  }

  data.table(population = population, coefficient = coefficient,
             selection = selection, comparison = comparison, outcome = outcome_name,
             patients = uniqueN(z$ID), status = "OK", statistic = test$statistic,
             df1 = test$df1, df2 = test$df2, p_value = test$p_value)
}

# Repeat the sensitivity analysis across all populations and outcomes
flexible_results <- list()
for (population in c("Primary_MI", "Complete_case")) {
  if (population == "Primary_MI") {
    src <- primary
  } else {
    src <- cc
  }

  # Repeat the sensitivity analysis using the Hendriksen coefficient
  coefficient <- "Hendriksen"
  d <- add_nonlinear_stage_terms(biomass(src, BIOMASS_K))

  for (selection in c("Nonnegative_PfHRP2", "Positive_only")) {
    if (selection == "Positive_only") {
      z <- d[hrp2_plasma > 0]
    } else {
      z <- d
    }

    for (outcome_name in OUTCOMES) {
      for (comparison in c("A_to_B", "B_to_C")) {
        flexible_results[[length(flexible_results) + 1L]] <-
          fit_flexible_comparison(z, outcome_name, population, coefficient, selection,
                                  comparison)
      }
    }
  }
}

# Combine and export the sensitivity analysis results
flex_dt <- rbindlist(flexible_results, fill = TRUE)
assert(nrow(flex_dt) == 24L, "Incomplete flexible sensitivity grid")
fwrite(flex_dt, file.path(CSV_OUT, "flexible_microscopy_sensitivity.csv"))

# Sensitivity analyses: excluding negative sequestered biomass estimates
seq_mi <- biomass(primary, BIOMASS_K)
seq_cc <- biomass(cc, BIOMASS_K)

# Summarise estimator disagreement across imputations
seq_mi[, circ_total_ratio := fifelse(total_per_kg > 0, circulating_per_kg / total_per_kg,
                                     NA_real_)]

seq_cc[, circ_total_ratio := fifelse(total_per_kg > 0, circulating_per_kg / total_per_kg,
                                     NA_real_)]

# Report negative and zero estimates per imputation
negative_audit <- seq_mi[, .(patients = .N, negative = sum(sequestered_raw_per_kg < 0),
                             zero = sum(sequestered_raw_per_kg == 0)), by = .imp]

negative_summary <- negative_audit[, .(patients = mean(patients),
                                       mean_negative = mean(negative),
                                       percent_negative = 100 * mean(negative / patients),
                                       mean_zero = mean(zero))]

# Keep the same patients across all 20 imputations
mi_flags <- seq_mi[, .(n_imputations = .N,
                       any_negative = any(sequestered_raw_per_kg < 0),
                       median_ratio = if (all(is.na(circ_total_ratio))) {NA_real_}
                       else {median(circ_total_ratio, na.rm = TRUE)}), by = ID]

assert(all(mi_flags$n_imputations == 20L), "Incomplete sensitivity imputation coverage")

cc_flags <- seq_cc[, .(any_negative = sequestered_raw_per_kg < 0,
                       median_ratio = circ_total_ratio), by = ID]

# Compare exclusions based on estimator disagreement
# Ratio thresholds are exploratory, not validated cut-offs
# Zero PfHRP2 patients remain in the analyses using these ratios
sensitivity_rules <- c("Exclude_any_negative", "Exclude_ratio_ge_2", "Exclude_ratio_ge_5")

sensitivity_counts <- list()
sensitivity_tests <- list()
sensitivity_coefficients <- list()

for (population in c("Primary_MI", "Complete_case")) {
  if (population == "Primary_MI") {
    source_data <- seq_mi
    flags <- mi_flags
  } else {
    source_data <- seq_cc
    flags <- cc_flags
  }

  for (selection in sensitivity_rules) {
    if (selection == "Exclude_any_negative") {
      excluded <- flags[any_negative == TRUE, ID]
    } else if (selection == "Exclude_ratio_ge_2") {
      excluded <- flags[!is.na(median_ratio) & median_ratio >= 2, ID]
    } else {
      excluded <- flags[!is.na(median_ratio) & median_ratio >= 5, ID]
    }

    z <- source_data[!ID %in% excluded]

    if (population == "Primary_MI") {
      assert(all(z[, .(n = uniqueN(.imp)), by = ID]$n == 20L),
             "Incomplete sensitivity cohort imputations")
    }

    sensitivity_counts[[length(sensitivity_counts) + 1L]] <-
      data.table(population = population, selection = selection,
                 original_patients = uniqueN(source_data$ID),
                 excluded_patients = length(excluded), retained_patients = uniqueN(z$ID))

    # Standard nested comparisons
    for (comparison in c("A_to_B", "B_to_C")) {
      ans <- fit_comparison(z, outcome_name = "sequestered_asinh", population = population,
                            coefficient = "Hendriksen", selection = selection,
                            comparison = comparison)

      assert(!is.null(ans), "Sensitivity regression failed")

      sensitivity_tests[[length(sensitivity_tests) + 1L]] <-
        ans$test[, model_type := "Standard"][]

      # Retain pooled coefficients and confidence intervals
      sensitivity_coefficients[[length(sensitivity_coefficients) + 1L]] <- ans$coef
    }

    # Nonlinear microscopy comparisons
    z_flexible <- add_nonlinear_stage_terms(z)

    for (comparison in c("A_to_B", "B_to_C")) {
      ans <- fit_flexible_comparison(z_flexible, outcome_name = "sequestered_asinh",
                                     population = population, coefficient = "Hendriksen",
                                     selection = selection, comparison = comparison)

      sensitivity_tests[[length(sensitivity_tests) + 1L]] <-
        ans[, model_type := "Nonlinear_microscopy"][]
    }
  }
}

sensitivity_counts_dt <- rbindlist(sensitivity_counts)
sensitivity_tests_dt <- rbindlist(sensitivity_tests, fill = TRUE)
sensitivity_coef_dt <- rbindlist(sensitivity_coefficients, fill = TRUE)

assert(nrow(sensitivity_counts_dt) == 6L, "Incomplete sensitivity cohort summary")
assert(nrow(sensitivity_tests_dt) == 24L, "Incomplete sensitivity model comparisons")
assert(nrow(sensitivity_coef_dt) == 114L, "Incomplete sensitivity coefficients")

# Export aggregate results without patient identifiers
fwrite(negative_summary, file.path(CSV_OUT, "sequestered_negative_summary.csv"))
fwrite(sensitivity_counts_dt, file.path(CSV_OUT, "sequestered_sensitivity_counts.csv"))
fwrite(sensitivity_tests_dt, file.path(CSV_OUT, "sequestered_sensitivity_tests.csv"))
fwrite(sensitivity_coef_dt, file.path(CSV_OUT, "sequestered_sensitivity_coefficients.csv"))

cat("\nNegative sequestered biomass summary:\n")
print(negative_summary)

cat("\nSequestered biomass sensitivity cohorts:\n")
print(sensitivity_counts_dt)

cat("\nSequestered biomass sensitivity comparisons:\n")
print(sensitivity_tests_dt[, .(population, selection, model_type, comparison, patients,
                               statistic, p_value)])

# Record published input units and conversion factors
coefficient_audit <-
  data.table(estimator = c("Hendriksen (2012)", "Circulating biomass"),
             coefficient_per_kg = c(BIOMASS_K, 80000),
             input_unit = c(PFHRP2_UNIT, PARASITAEMIA_UNIT),
             formula_per_kg = c("7.3e7 * PfHRP2 * (1 - hct/100)", "80000 * parasitaemia"),
             aquamat_literature_units_supported = TRUE, unit_metadata_confirmed = FALSE,
             literature_reference = c(paste0("Hendriksen et al. (2012), ",
                                             "doi:10.1371/journal.pmed.1001297"),
                                      paste0("Hendriksen et al. (2012), ",
                                             "doi:10.1371/journal.pmed.1001297")),
             note = paste0("Units supported by AQUAMAT publications; ",
                           "original files have no unit metadata, ",
                           "so preprocessing provenance remains unverified"))

fwrite(coefficient_audit, file.path(CSV_OUT, "biomass_coefficient_audit.csv"))

# Generate six individual figures for publication
main <- biomass(primary, BIOMASS_K)

# Summarise the estimated sequestered fraction across imputations
fraction_audit <- main[is.finite(fraction_raw),
                       .(patients = .N, negative = sum(fraction_raw < 0),
                         below_minus_one = sum(fraction_raw < -1),
                         within_zero_one = sum(fraction_raw >= 0 & fraction_raw <= 1),
                         median_fraction = median(fraction_raw),
                         p01 = as.numeric(quantile(fraction_raw, 0.01)),
                         p99 = as.numeric(quantile(fraction_raw, 0.99))), by = .imp]

# Export aggregate diagnostics without patient identifiers
fraction_summary <- fraction_audit[, lapply(.SD, mean),
                                   .SDcols = setdiff(names(fraction_audit), ".imp")]

fraction_summary[, `:=`(population = "Primary_MI", coefficient = "Hendriksen",
                        selection = "Positive_only")]

fwrite(fraction_summary, file.path(CSV_OUT, "sequestered_fraction_audit.csv"))

curves <- rbindlist(lapply(OUTCOMES, function(nm) fit_curve(main, nm)))

# Audit non-positive sequestered biomass estimates across imputations
seq_audit <- main[, .(negative = sum(sequestered_raw_per_kg < 0),
                      zero = sum(sequestered_raw_per_kg == 0),
                      nonpositive = sum(sequestered_raw_per_kg <= 0),
                      patients = .N), by = .imp]

seq_summary <- seq_audit[, .(mean_negative = mean(negative), mean_zero = mean(zero),
                             mean_nonpositive = mean(nonpositive),
                             percent_nonpositive = 100 * mean(nonpositive / patients))]

fwrite(seq_summary, file.path(CSV_OUT, "sequestered_nonpositive_audit.csv"))

# Back-transform pooled descriptive predictions at 6 and 30 hours
seq_curve <- curves[outcome == "sequestered_asinh" & age_hours %in% c(6, 30),
                    .(age_hours, estimate)]

seq_curve[, biomass_per_kg := 1e11 * sinh(estimate)]

seq_6 <- seq_curve[age_hours == 6, biomass_per_kg]
seq_30 <- seq_curve[age_hours == 30, biomass_per_kg]

seq_contrast <- data.table(age_6_biomass_per_kg = seq_6, age_30_biomass_per_kg = seq_30,
                           difference_per_kg = seq_30 - seq_6, ratio = seq_30 / seq_6,
                           percent_increase = 100 * (seq_30 / seq_6 - 1))

fwrite(seq_contrast, file.path(CSV_OUT, "sequestered_original_scale_contrast.csv"))

AQ_LINE <- "#C9343D"  # muted AQUAMAT red
AQ_BAND <- "#EC767D"
AQ_RUG <- "#C9343D"

# Exactly one row per patient; never count all 20 imputation copies
ages_once <- unique(primary[.imp == min(.imp), .(ID, mu_age)], by = "ID")
sd_once <- unique(primary[.imp == min(.imp), .(ID, sigma_age)], by = "ID")

age_ref <- circular_mean_48(ages_once$mu_age)

sigma_curves <- list()

for (outcome_name in OUTCOMES) {
  pieces <- split(main[which(is.finite(main[[outcome_name]]))], by = ".imp")

  fits <- lapply(pieces, function(z)
    lm(reformulate(c("age_sin", "age_cos", "sigma_5h"), outcome_name), data = z))

  assert(length(fits) == 20L, "Expected 20 fitted MI models for SD curves")

  sd_grid <- seq(quantile(sd_once$sigma_age, .05),
                 quantile(sd_once$sigma_age, .95), length.out = 100)

  X <- cbind(`(Intercept)` = 1, age_sin = sin(2*pi*age_ref/48),
             age_cos = cos(2*pi*age_ref/48), sigma_5h = sd_grid/5)

  dfcom <- min(vapply(fits, df.residual, numeric(1)))
  ans <- mi_curve_prediction(fits, X, dfcom)

  sigma_curves[[outcome_name]] <- cbind(data.table(outcome = outcome_name,
                                                   sigma_hours = sd_grid), ans)
}

sigdt <- rbindlist(sigma_curves)

# Compare estimated sequestered biomass at age SDs of 3 and 8 hours
seq_sd_curve <- sigma_curves[["sequestered_asinh"]]

seq_sd_3 <- approx(seq_sd_curve$sigma_hours, seq_sd_curve$estimate, xout = 3)$y
seq_sd_8 <- approx(seq_sd_curve$sigma_hours, seq_sd_curve$estimate, xout = 8)$y

seq_sd_contrast <- data.table(sd_3_biomass_per_kg = 1e11 * sinh(seq_sd_3),
                              sd_8_biomass_per_kg = 1e11 * sinh(seq_sd_8))

seq_sd_contrast[, percent_increase :=
                  100 * (sd_8_biomass_per_kg / sd_3_biomass_per_kg - 1)]

fwrite(seq_sd_contrast, file.path(CSV_OUT, "sequestered_sd_original_scale_contrast.csv"))

# Use matching vertical limits for each biomass outcome
shared_y_limits <- lapply(OUTCOMES, function(nm) {
  z <- rbind(curves[outcome == nm, .(lower_95, upper_95)],
             sigdt[outcome == nm, .(lower_95, upper_95)])
  
  y <- range(c(z$lower_95, z$upper_95), na.rm = TRUE)
  pad <- diff(y) * 0.05
  y + c(-pad, pad)
})

names(shared_y_limits) <- OUTCOMES

# Display original biomass units while retaining transformed regression scales
figure_ylabels <- list(total_log10 = expression(10^{10}*" parasites/kg"),
  circulating_log10 = expression(10^9*" parasites/kg"),
  sequestered_asinh = expression(10^{10}*" parasites/kg"))

# Convert transformed y-axis positions to original biomass units
biomass_axis <- function(outcome_name) {
  if (outcome_name == "total_log10") {
    return(scale_y_continuous(breaks = seq(8, 13, by = 0.1),
                              labels = function(y) format(round(10^y / 1e10, 1),
                                                          trim = TRUE,
                                                          scientific = FALSE)))
  }

  if (outcome_name == "circulating_log10") {
    return(scale_y_continuous(breaks = seq(7, 12, by = 0.1),
                              labels = function(y) format(round(10^y / 1e9, 1),
                                                          trim = TRUE,
                                                          scientific = FALSE)))
  }


  if (outcome_name == "sequestered_asinh") {
    return(scale_y_continuous(
      breaks = function(limits) seq(limits[1], limits[2], length.out = 4),
      labels = function(y) format(round(sinh(y) * 10, 1), trim = TRUE,
                                  scientific = FALSE)))
  }
}

# Format figures for publication
journal_theme <- theme_bw(base_size = 13) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        axis.title = element_text(size = 12, colour = "grey20"),
        axis.text = element_text(size = 11, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.35),
        plot.margin = margin(10, 12, 10, 10))

# Rug marks show patient concentration, not parasite density
plot_one <- function(dat, outcome_name, predictor = c("age", "sd")) {
  predictor <- match.arg(predictor)

  if (predictor == "age") {
    x_col <- "age_hours"
    rug_col <- "mu_age"
    rug_data <- ages_once
    x_label <- "Mean parasite age (h; 48-hour cycle)"

  } else {
    x_col <- "sigma_hours"
    rug_col <- "sigma_age"
    rug_data <- sd_once
    x_label <- "Parasite age standard deviation (h)"
  }

  subset_data <- dat[dat$outcome == outcome_name]
  assert(nrow(subset_data) > 0L, paste("No curve for", outcome_name))

  if (predictor == "age") {
    x_scale <- scale_x_continuous(breaks = seq(0, 48, 6), limits = c(0, 48),
                                  expand = c(0, 0))
  } else {
    x_scale <- scale_x_continuous(limits = range(subset_data[[x_col]]),
                                  breaks = seq(0, 18, 2),
                                  expand = expansion(mult = c(0.02, 0.02)))
  }

  p <- ggplot(subset_data, aes(x = .data[[x_col]], y = estimate)) +

    geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
                fill = AQ_BAND, alpha = 0.24) +

    geom_line(colour = AQ_LINE, linewidth = 0.65) +

    geom_rug(data = rug_data, aes(x = .data[[rug_col]]),
             inherit.aes = FALSE, sides = "b", colour = AQ_RUG,
             alpha = 0.095, linewidth = 0.25) +

    labs(x = x_label, y = figure_ylabels[[outcome_name]]) +
    journal_theme + x_scale + biomass_axis(outcome_name) +
    coord_cartesian(ylim = shared_y_limits[[outcome_name]])

  # Allow vertical space for the long rotated sequestered-biomass label
  if (outcome_name == "sequestered_asinh") {
    p <- p + theme(
      plot.margin = margin(t = 24, r = 12, b = 18, l = 24, unit = "pt")
    )
  }

  p
}

# Export six separate vector PDFs, one per outcome and predictor
for (predictor in c("age", "sd")) {
  if (predictor == "age") {
    d <- curves
  } else {
    d <- sigdt
  }

  for (nm in OUTCOMES) {
    plt <- plot_one(d, nm, predictor)

    stem <- paste0(predictor, "_", c(total_log10 = "total_biomass",
                                     circulating_log10 = "circulating_biomass",
                                     sequestered_asinh = "sequestered_biomass")[[nm]])

    ggsave(file.path(FIG_OUT, paste0(stem, ".pdf")), plt,
           width = 6.5, height = 2.8, units = "in", device = "pdf", bg = "white")
  }
}

# Check that all six figures were saved
figure_names <- c("age_total_biomass.pdf", "age_circulating_biomass.pdf",
                  "age_sequestered_biomass.pdf", "sd_total_biomass.pdf",
                  "sd_circulating_biomass.pdf", "sd_sequestered_biomass.pdf")

figure_paths <- file.path(FIG_OUT, figure_names)
assert(all(file.exists(figure_paths)), "One or more expected figures are missing")
assert(all(file.info(figure_paths)$size > 0), "One or more figures are empty")

cat(sprintf("SD figure reference mean parasite age: %.2f h (circular mean).\n", age_ref))

# Report the completed analysis and cohort sizes
cat("Final validation candidate completed. Primary:", uniqueN(primary$ID),
    "patients; Complete case:", uniqueN(cc$ID), "patients.\n")

cat("Model tests:", length(test_tables), "; coefficient tables:", length(coef_tables), "\n")
cat("Stage composition audit and outputs:", CSV_OUT, "\n")
