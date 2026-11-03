# COMPLETE-CASE BIMODALITY SCREEN

library(data.table)
options(scipen = 999)

# ---- Paths ----
DATA_FILE <- file.path("Data/dat_all.csv")
PROJECT_DIR <- file.path("Staging/Complete case")

OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
CSV_DIR <- file.path(OUTPUT_DIR, "CSV")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")

dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Screening settings ----
SIGMA_MIN_VALUES <- c(1, 2, 3)
SIGMA_MAX <- 24

# Fixed stage-specific sequestration profile.
# Values represent the relative contribution of each latent parasite stage
# to the peripheral-blood stage distribution. Progressive reduction with
# parasite maturity represents sequestration of infected red blood cells
# from the peripheral circulation.
SEQUESTRATION_PROFILE <- c(1.000, 0.500, 0.250, 0.200, 0.035, 0.027, 0.005)

stage_names <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")
stage_boundaries_mix <- c(0, 6, 16, 26, 30, 34, 38, 48)

# ---- Output files ----
# Full bimodality screening results
SCREEN_FILE <- file.path(RDS_DIR, "hierarchical_bimodality_screen.rds")

# Final classification used by hierarchical_main.R
CLASSIFICATION_FILE <- file.path(RDS_DIR, "hierarchical_bimodality_classification.rds")

# Human-readable screening outputs
CANDIDATES_CSV <- file.path(CSV_DIR, "hierarchical_bimodal_candidates.csv")
CLASSIFICATION_CSV <- file.path(CSV_DIR, "hierarchical_bimodality_classification.csv")

# ---- Load complete-case cohort ----
dat_all_mix <- fread(DATA_FILE)
dat_all_mix <- dat_all_mix[, .(ID, age, coma, aconvulsions, hb, hct, aanaemia, ashockcomp,
                               temp, ahypoglyc, outcome, site, Trial, be, bun,
                               parasitaemia, tring, sring, lring, etroph, mtroph, ltroph,
                               schizont)]

dat_all_mix <- dat_all_mix[complete.cases(dat_all_mix)]
dat_all_mix[, total_stage := tring + sring + lring + etroph + mtroph + ltroph + schizont]
dat_all_mix <- dat_all_mix[total_stage >= 93 & total_stage <= 107]
dat_all_mix[, total_stage := NULL]

stage_matrix_mix <- as.matrix(dat_all_mix[, ..stage_names])
n_patients_mix <- nrow(stage_matrix_mix)

metadata <- dat_all_mix[, .(patient = seq_len(.N), ID, age, coma, aconvulsions, hb, hct,
                            aanaemia, ashockcomp, temp, ahypoglyc, outcome, site, Trial,
                            be, bun, parasitaemia, tring, sring, lring, etroph, mtroph,
                            ltroph, schizont)]

# ---- Peripheral stage probabilities under sequestration ----
# Calculate latent stage probabilities over m = -1:1
# Apply the fixed stage-specific sequestration profile to obtain
# the corresponding peripheral stage probabilities
# M_eff = 1, so no additional generation weighting is applied
stage_probabilities_mix <- function(mu_age, sigma_age) {
  stage_probs <- numeric(7)

  for (j in 1:7) {
    for (m in -1:1) {

      upper <- (stage_boundaries_mix[j + 1] + m*48 - mu_age) / sigma_age
      lower <- (stage_boundaries_mix[j] + m*48 - mu_age) / sigma_age

      stage_probs[j] <- stage_probs[j] + pt(upper, df = 10) - pt(lower, df = 10)
    }
  }

  # Apply the fixed stage-specific sequestration profile before
  # normalising to the predicted peripheral stage distribution
  stage_probs <- stage_probs * SEQUESTRATION_PROFILE

  total <- sum(stage_probs)
  stage_probs / total
}

# ---- Single-component log-likelihood ----
single_loglik_mix <- function(params, counts) {
  mu <- params[1]
  sigma <- exp(params[2])
  probs <- stage_probabilities_mix(mu, sigma)
  sum(counts * log(probs))
}

# ---- Find single-component starting values ----
mu_grid_initial <- seq(0, 48, by = 2)
sigma_grid_initial <- c(1, 2, 3, 4, 6, 8, 10, 12, 16, 20)

find_single_start_mix <- function(counts) {
  best_ll <- -Inf
  best_mu <- NA_real_
  best_sigma <- NA_real_

  for (mu in mu_grid_initial) {
    for (sigma in sigma_grid_initial) {

      ll <- single_loglik_mix(c(mu, log(sigma)), counts)

      if (is.finite(ll) && ll > best_ll) {
        best_ll <- ll
        best_mu <- mu
        best_sigma <- sigma
      }
    }
  }
  list(mu = best_mu, sigma = best_sigma, loglik = best_ll)
}

# ---- Fit single-component model ----
fit_single_patient_mix <- function(counts) {
  start <- find_single_start_mix(counts)

  fit <- optim(par = c(start$mu, log(start$sigma)),
               fn = function(x) {-single_loglik_mix(x, counts)},
               method = "Nelder-Mead", control = list(maxit = 3000))

  list(mu = fit$par[1] %% 48, sigma = exp(fit$par[2]), loglik = -fit$value,
       convergence = fit$convergence)
}

# ---- Circular distance ----
circular_distance_48_mix <- function(x, y) {
  d <- abs((x %% 48) - (y %% 48))
  min(d, 48 - d)
}

# ---- Starting pairs for two-component model ----
mixture_mu_pairs <- list(c(2, 18), c(2, 26), c(2, 40), c(10, 24), c(10, 34), c(18, 34),
                         c(22, 38))

# ---- Bounded SD transformation ----
sigma_from_raw <- function(raw, sigma_min) {
  sigma_min + (SIGMA_MAX - sigma_min) * plogis(raw)
}

raw_from_sigma <- function(sigma, sigma_min) {
  sigma <- min(max(as.numeric(sigma), sigma_min + 1e-6), SIGMA_MAX - 1e-6)
  p <- (sigma - sigma_min) / (SIGMA_MAX - sigma_min)
  qlogis(p)
}

# ---- Two-component mixture log-likelihood ----
mixture_loglik_mix <- function(params, counts, sigma_min) {
  mu1 <- params[1]
  mu2 <- params[2]
  sigma1 <- sigma_from_raw(params[3], sigma_min)
  sigma2 <- sigma_from_raw(params[4], sigma_min)
  weight1 <- plogis(params[5])

  if (!all(is.finite(c(mu1, mu2, sigma1, sigma2, weight1)))) {
    return(-Inf)
  }

  p1 <- stage_probabilities_mix(mu1, sigma1)
  p2 <- stage_probabilities_mix(mu2, sigma2)

  probs <- weight1*p1 + (1 - weight1)*p2

  if (any(!is.finite(probs)) || any(probs <= 0)) {
    return(-Inf)
  }
  sum(counts * log(probs))
}

# ---- Fit constrained two-component model ----
fit_mixture_patient_constrained <- function(counts, sigma_min) {
  single_start <- find_single_start_mix(counts)

  best_fit <- NULL
  best_value <- Inf

  for (pair in mixture_mu_pairs) {
    initial_sigma <- min(max(single_start$sigma / 2, sigma_min + 0.5), SIGMA_MAX - 0.5)
    initial_sigma_parameter <- raw_from_sigma(initial_sigma, sigma_min)

    initial_values <- c(pair[1], pair[2], initial_sigma_parameter, initial_sigma_parameter,
                        qlogis(0.5))

    fit <- try(optim(par = initial_values,
                     fn = function(x) {-mixture_loglik_mix(x, counts, sigma_min)},
                     method = "Nelder-Mead", control = list(maxit = 4000)), silent = TRUE)

    if (inherits(fit, "try-error")) {
      next
    }

    if (is.finite(fit$value) && fit$value < best_value) {
      best_fit <- fit
      best_value <- fit$value
    }
  }

  if (is.null(best_fit)) {
    return(NULL)
  }

  mu1 <- best_fit$par[1] %% 48
  mu2 <- best_fit$par[2] %% 48

  sigma1 <- sigma_from_raw(best_fit$par[3], sigma_min)
  sigma2 <- sigma_from_raw(best_fit$par[4], sigma_min)

  weight1 <- plogis(best_fit$par[5])

  if (mu1 > mu2) {
    tmp_mu <- mu1
    mu1 <- mu2
    mu2 <- tmp_mu
    tmp_sigma <- sigma1
    sigma1 <- sigma2
    sigma2 <- tmp_sigma
    weight1 <- 1 - weight1
  }

  list(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, weight1 = weight1,
       loglik = -best_fit$value, convergence = best_fit$convergence)
}

# ---- Run each constrained bimodality screen ----
all_results_constrained <- list()

for (sigma_min in SIGMA_MIN_VALUES) {

  cat("\nMinimum component sigma = ", sigma_min, " hours; maximum = ", SIGMA_MAX, " hours\n",
      sep = "")

  results <- data.frame(patient = seq_len(n_patients_mix), ID = dat_all_mix$ID,
                        n_cells = rowSums(stage_matrix_mix), single_loglik = NA_real_,
                        mixture_loglik = NA_real_, delta_loglik = NA_real_,
                        single_BIC = NA_real_, mixture_BIC = NA_real_,
                        delta_BIC = NA_real_, mu_single = NA_real_, sigma_single = NA_real_,
                        mu1 = NA_real_, sigma1 = NA_real_, mu2 = NA_real_, sigma2 = NA_real_,
                        weight1 = NA_real_, circular_separation = NA_real_,
                        mixture_convergence = NA_integer_, bimodal_candidate =FALSE)

  for (patient in seq_len(n_patients_mix)) {
    counts <- stage_matrix_mix[patient, ]

    single_fit <- fit_single_patient_mix(counts)
    mixture_fit <- fit_mixture_patient_constrained(counts, sigma_min)

    if (is.null(mixture_fit)) {
      next
    }

    n_cells <- sum(counts)

    single_BIC <- -2*single_fit$loglik + 2*log(n_cells)
    mixture_BIC <- -2*mixture_fit$loglik + 5*log(n_cells)
    delta_BIC <- single_BIC - mixture_BIC

    separation <- circular_distance_48_mix(mixture_fit$mu1, mixture_fit$mu2)

    results[patient, "single_loglik"] <- single_fit$loglik
    results[patient, "mixture_loglik"] <- mixture_fit$loglik
    results[patient, "delta_loglik"] <- mixture_fit$loglik - single_fit$loglik

    results[patient, "single_BIC"] <- single_BIC
    results[patient, "mixture_BIC"] <- mixture_BIC
    results[patient, "delta_BIC"] <- delta_BIC

    results[patient, "mu_single"] <- single_fit$mu
    results[patient, "sigma_single"] <- single_fit$sigma
    results[patient, "mu1"] <- mixture_fit$mu1
    results[patient, "sigma1"] <- mixture_fit$sigma1
    results[patient, "mu2"] <- mixture_fit$mu2
    results[patient, "sigma2"] <- mixture_fit$sigma2

    results[patient, "weight1"] <- mixture_fit$weight1
    results[patient, "circular_separation"] <- separation
    results[patient, "mixture_convergence"] <- mixture_fit$convergence

    results[patient, "bimodal_candidate"] <- is.finite(delta_BIC) && delta_BIC > 10 &&
                                             mixture_fit$weight1 >= 0.10 &&
                                             mixture_fit$weight1 <= 0.90 &&
                                             mixture_fit$sigma1 >= sigma_min &&
                                             mixture_fit$sigma1 <= SIGMA_MAX &&
                                             mixture_fit$sigma2 >= sigma_min &&
                                             mixture_fit$sigma2 <= SIGMA_MAX &&
                                             separation >= 6 && mixture_fit$convergence == 0

    if (patient %% 250 == 0) {
      cat("Completed ", patient, " / ", n_patients_mix, "\n", sep = "")
    }
  }

  all_results_constrained[[paste0("sigma_min_", sigma_min)]] <- results
  n_candidate <- sum(results$bimodal_candidate, na.rm = TRUE)

  cat("\nMinimum sigma = ", sigma_min, " h: ", n_candidate, " / ", n_patients_mix,
      " candidates (", round(100 * n_candidate / n_patients_mix, 2), "%)\n", sep = "")
}

# ---- Identify robust bimodal patients ----
candidate_sets <- lapply(SIGMA_MIN_VALUES, function(x) {
  r <- all_results_constrained[[paste0("sigma_min_", x)]]
  r$patient[r$bimodal_candidate]})

names(candidate_sets) <- paste0("sigma_min_", SIGMA_MIN_VALUES)
stable_candidates <- Reduce(intersect, candidate_sets)

cat("\nRobust bimodal patients\n\n")
cat("Number stable across all three thresholds: ", length(stable_candidates), "\n", sep = "")

# ---- Build classification table ----
r3 <- all_results_constrained[["sigma_min_3"]]

bimodality_classification <- data.frame(patient = seq_len(n_patients_mix),
                                        ID = dat_all_mix$ID,
                                        pass_sigma1 = seq_len(n_patients_mix) %in%
                                          candidate_sets$sigma_min_1,
                                        pass_sigma2 = seq_len(n_patients_mix) %in%
                                          candidate_sets$sigma_min_2,
                                        pass_sigma3 = seq_len(n_patients_mix) %in%
                                          candidate_sets$sigma_min_3)

bimodality_classification$n_screens_passed <-
  rowSums(bimodality_classification[, c("pass_sigma1", "pass_sigma2", "pass_sigma3")])

bimodality_classification$category <-
  ifelse(bimodality_classification$n_screens_passed == 3, "Robust bimodal",
    ifelse(bimodality_classification$n_screens_passed == 2, "Probable bimodal",
      ifelse(bimodality_classification$n_screens_passed == 1, "Possible bimodal",
             "No strong bimodal evidence")))

# Retain sigma = 3 fitted values for downstream initialisation
bimodality_classification$delta_BIC_sigma3 <- r3$delta_BIC
bimodality_classification$mu1_sigma3 <- r3$mu1
bimodality_classification$sigma1_sigma3 <- r3$sigma1
bimodality_classification$mu2_sigma3 <- r3$mu2
bimodality_classification$sigma2_sigma3 <- r3$sigma2
bimodality_classification$weight1_sigma3 <- r3$weight1
bimodality_classification$weight2_sigma3 <- 1 - r3$weight1
bimodality_classification$circular_separation_sigma3 <- r3$circular_separation
bimodality_classification$mixture_convergence_sigma3 <- r3$mixture_convergence

# ---- Classification summary ----
category_counts <- as.data.frame(table(bimodality_classification$category))
names(category_counts) <- c("category", "n_patients")

category_counts$percentage <- 100 * category_counts$n_patients /
  nrow(bimodality_classification)

cat("\nFixed-sequestration bimodality classification\n\n")
print(category_counts, row.names = FALSE)
cat("\nSensitivity-screen candidate counts\n\n")

print(sapply(all_results_constrained, function(x) {
  sum(x$bimodal_candidate, na.rm = TRUE)}))

# ---- Save outputs ----
bimodal_candidates <- bimodality_classification[
  bimodality_classification$n_screens_passed > 0, ]

bimodal_candidates <- merge(metadata, bimodal_candidates, by = c("patient", "ID"),
                            sort = FALSE)

bimodality_classification_csv <- merge(metadata, bimodality_classification,
                                       by = c("patient", "ID"), sort = FALSE)

fwrite(bimodal_candidates, CANDIDATES_CSV)
fwrite(bimodality_classification_csv, CLASSIFICATION_CSV)

saveRDS(all_results_constrained, SCREEN_FILE)
saveRDS(bimodality_classification, CLASSIFICATION_FILE)

cat("\nSaved fixed-sequestration bimodality screening results.\n")
cat("Robust bimodal patients: ", length(stable_candidates), "\n", sep = "")
cat("\nOutputs:\n")

print(c(SCREEN_FILE, CLASSIFICATION_FILE, CANDIDATES_CSV, CLASSIFICATION_CSV),
      row.names = FALSE)