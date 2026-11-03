# HIERARCHICAL PARASITE-AGE MODELLING
# COMPLETE CASE

library(data.table)
library(rstan)

options(scipen = 999)
options(mc.cores = 4)

# Prevent compiled Stan-model RDS files being written beside the Stan source files
rstan_options(auto_write = FALSE)

# ---- Directories ----
PROJECT_DIR <- file.path("Staging/Complete case")

DATA_DIR <- "Data"
STAN_DIR <- "Staging/code/Stan"
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
CSV_DIR <- file.path(OUTPUT_DIR, "CSV")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")

dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Files ----
DATA_FILE <- file.path(DATA_DIR, "dat_all.csv")
UNIMODAL_MODEL_FILE <- file.path(STAN_DIR, "biomass_hierarchical_unimodal.stan")
BIMODAL_MODEL_FILE <- file.path(STAN_DIR, "biomass_hierarchical_bimodal.stan")
CLASSIFICATION_FILE <- file.path(RDS_DIR, "hierarchical_bimodality_classification.rds")

UNIMODAL_FIT_FILE <- file.path(RDS_DIR, "fit_hierarchical_unimodal.rds")
UNIMODAL_SUMMARY_FILE <- file.path(RDS_DIR, "hierarchical_unimodal_patient_summary.rds")
BIMODAL_FIT_FILE <- file.path(RDS_DIR, "fit_hierarchical_bimodal.rds")
BIMODAL_MEMBERSHIP_FILE <- file.path(RDS_DIR, "hierarchical_bimodal_membership.rds")

BIMODAL_COMPONENT_SUMMARY_FILE <- file.path(RDS_DIR,
                                            "hierarchical_bimodal_component_summary.rds")

HIERARCHICAL_OUTPUT_FILE <- file.path(CSV_DIR, "dat_all_hierarchical.csv")

# ---- Model constraints ----
AA <- c(0, 6, 16, 26, 30, 34, 38, 48)
STAGE_NAMES <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")
SEQUESTRATION_PROFILE <- c(1.000, 0.500, 0.250, 0.200, 0.035, 0.027, 0.005)

RHO_MU_AGE <- 0.5
N_CHAINS <- 4
MAX_TREEDEPTH <- 12

# ---- Functions ----
circular_mean_48 <- function(x) {
  theta <- 2 * pi * ((x %% 48) / 48)
  phase <- atan2(mean(sin(theta)), mean(cos(theta)))
  (phase %% (2 * pi)) * 48 / (2 * pi)
}

mixture_circular_centre_48 <- function(mu1, mu2, weight1) {
  theta1 <- 2 * pi * ((mu1 %% 48) / 48)
  theta2 <- 2 * pi * ((mu2 %% 48) / 48)

  x <- weight1 * cos(theta1) + (1 - weight1) * cos(theta2)
  y <- weight1 * sin(theta1) + (1 - weight1) * sin(theta2)

  phase <- atan2(y, x)
  (phase %% (2 * pi)) * 48 / (2 * pi)
}

mixture_sd_draw <- function(mu1, delta, sigma1, sigma2, weight1) {
  mu2 <- mu1 + delta
  mu_mix <- weight1 * mu1 + (1 - weight1) * mu2

  variance_mix <- weight1 * (sigma1^2 + (mu1 - mu_mix)^2) +
    (1 - weight1) * (sigma2^2 + (mu2 - mu_mix)^2)

  sqrt(pmax(variance_mix, 0))
}

# ---- Prepare canonical cohort ----
cat("\nPreparing hierarchical-model cohort...\n")
dat_all <- fread(DATA_FILE)

required_columns <- c("ID", "age", "coma","aconvulsions", "hb", "hct", "aanaemia",
                      "ashockcomp", "temp", "ahypoglyc", "outcome", "site", "Trial", "be",
                      "bun", "parasitaemia", STAGE_NAMES)

dat_all <- dat_all[, ..required_columns]
dat_all <- dat_all[complete.cases(dat_all)]
dat_all[, total_stage := tring + sring + lring + etroph + mtroph + ltroph + schizont]
dat_all <- dat_all[total_stage >= 93 & total_stage <= 107]

N_TOTAL <- nrow(dat_all)

# ---- Patient partition ----
bimodality_classification <- readRDS(CLASSIFICATION_FILE)

robust_bimodal <- bimodality_classification[
  as.character(bimodality_classification$category) == "Robust bimodal", , drop = FALSE]

robust_bimodal <- robust_bimodal[order(robust_bimodal$patient), , drop = FALSE]

robust_rows <- robust_bimodal$patient
unimodal_rows <- setdiff(seq_len(N_TOTAL), robust_rows)

N_BIMODAL <- length(robust_rows)
N_UNIMODAL <- length(unimodal_rows)

cat("Total patients: ", N_TOTAL, "\n", sep = "")
cat("Unimodal patients: ", N_UNIMODAL, "\n", sep = "")
cat("Bimodal patients: ", N_BIMODAL, "\n", sep = "")

# ---- Shared data objects ----
dat_stage <- dat_all[, ..STAGE_NAMES]

dat_meta <- dat_all[, .(ID, age, coma, aconvulsions, hb, hct, aanaemia, ashockcomp, temp,
                        ahypoglyc, outcome, site, Trial, be, bun, parasitaemia)]

# ---- Unimodal model ----
cat("\nHierarchical unimodal model\n")

unimodal_stage <- dat_stage[unimodal_rows, , drop = FALSE]
unimodal_meta <- dat_meta[unimodal_rows, , drop = FALSE]

# ---- Unimodal starting values ----
mu_grid <- seq(0, 48, by = 1)
sigma_grid <- seq(1, 12, by = 1)

unimodal_stage_probabilities <- function(mu_age, sigma_age) {
  stage_props <- numeric(7)

  for (j in 1:7) {
    for (m in -1:1) {
      upper <- (AA[j + 1] + m*48 - mu_age) / sigma_age
      lower <- (AA[j] + m*48 - mu_age) / sigma_age

      stage_props[j] <- stage_props[j] + pt(upper, df = 10) - pt(lower, df = 10)
    }
  }
  stage_props <- stage_props * SEQUESTRATION_PROFILE
  stage_props / sum(stage_props)
}

mu_start <- numeric(N_UNIMODAL)
sigma_start <- numeric(N_UNIMODAL)

cat("\nCalculating unimodal starting values...\n")
start_initialisation <- Sys.time()

for (i in seq_len(N_UNIMODAL)) {
  counts_i <- as.numeric(unimodal_stage[i, ])
  positive_counts <- counts_i > 0

  best_value <- -Inf
  best_mu <- NA_real_
  best_sigma <- NA_real_

  for (sigma_i in sigma_grid) {
    for (mu_i in mu_grid) {
      p <- unimodal_stage_probabilities(mu_i, sigma_i)
      ll <- sum(counts_i[positive_counts] * log(p[positive_counts]))

      if (is.finite(ll) && ll > best_value) {
        best_value <- ll
        best_mu <- mu_i
        best_sigma <- sigma_i
      }
    }
  }
  mu_start[i] <- best_mu
  sigma_start[i] <- best_sigma
}

cat("Initialisation time: ",
    round(as.numeric(difftime(Sys.time(), start_initialisation, units = "mins")), 2),
    " minutes\n", sep = "")

# ---- Unimodal initialisation ----
make_unimodal_init <- function(chain_id) {
  set.seed(1684 + chain_id)

  mu_chain <- mu_start + rnorm(N_UNIMODAL, 0, 0.20)
  sigma_chain <- pmax(sigma_start + rnorm(N_UNIMODAL, 0, 0.10), 0.5)
  log_sigma_chain <- log(sigma_chain)

  mu_mean_start <- mean(mu_chain)
  mu_sd_start <- max(sd(mu_chain), 1)

  log_sigma_mean_start <- mean(log_sigma_chain)
  log_sigma_sd_start <- max(sd(log_sigma_chain), 0.25)

  z_mu_age_start <- (mu_chain - mu_mean_start) / (RHO_MU_AGE * mu_sd_start)

  list(mu_age = mu_chain, z_mu_age = z_mu_age_start,
       mu_age_mean = rnorm(1, mu_mean_start, 0.10),
       mu_age_sd = max(rnorm(1, mu_sd_start, 0.10), 0.5), log_sigma_age = log_sigma_chain,
       log_sigma_age_mean = rnorm(1, log_sigma_mean_start, 0.05),
       log_sigma_age_sd = max(rnorm(1, log_sigma_sd_start, 0.05), 0.15))
}

unimodal_init_list <- lapply(seq_len(N_CHAINS), make_unimodal_init)

# ---- Fit unimodal model ----
cat("\nCompiling unimodal model...\n")
unimodal_model <- stan_model(UNIMODAL_MODEL_FILE)

cat("Running unimodal fit...\n")
unimodal_start_time <- Sys.time()

fit_unimodal <- sampling(object = unimodal_model,
                         data = list(K_stages = 7, n = N_UNIMODAL,
                                     stage_obs = as.matrix(unimodal_stage), aa = AA,
                                     sequestration_profile = SEQUESTRATION_PROFILE,
                                     rho_mu_age = RHO_MU_AGE),
                         seed = 1684, warmup = 500, iter = 1500, chains = N_CHAINS,
                         init = unimodal_init_list, refresh = 100,
                         control = list(adapt_delta = 0.99, max_treedepth = MAX_TREEDEPTH))

cat("Unimodal runtime: ",
    round(as.numeric(difftime(Sys.time(), unimodal_start_time, units = "mins")), 2),
    " minutes\n", sep = "")

saveRDS(fit_unimodal, UNIMODAL_FIT_FILE)

# ---- Unimodal summary ----
unimodal_posterior <- rstan::extract(fit_unimodal, pars = c("mu_age", "sigma_age"),
                                     permuted = FALSE, inc_warmup = FALSE)

unimodal_mu_names <- grep("^mu_age\\[", dimnames(unimodal_posterior)[[3]], value = TRUE)
unimodal_sigma_names <- grep("^sigma_age\\[", dimnames(unimodal_posterior)[[3]],
                             value = TRUE)

unimodal_mu_age <- vapply(seq_len(N_UNIMODAL), function(i) {
  circular_mean_48(unimodal_posterior[, , unimodal_mu_names[i]])}, numeric(1))

unimodal_sigma_age <- vapply(seq_len(N_UNIMODAL), function(i) {
  mean(unimodal_posterior[, , unimodal_sigma_names[i]])}, numeric(1))

unimodal_patient_summary <- data.frame(ID = unimodal_meta$ID, mu_age = unimodal_mu_age,
                                       sigma_age = unimodal_sigma_age)

saveRDS(unimodal_patient_summary, UNIMODAL_SUMMARY_FILE)

# The full fit is already saved
# Remove large objects before starting the bimodal fit
rm(unimodal_posterior, unimodal_model, fit_unimodal)
gc()

# ---- Bimodal model ----
cat("\nHierarchical bimodal model\n")

bimodal_stage <- dat_stage[robust_rows, , drop = FALSE]
bimodal_meta <- dat_meta[robust_rows, , drop = FALSE]

# ---- Bimodal starting values ----
canonical_mu1 <- robust_bimodal$mu1_sigma3
canonical_mu2 <- robust_bimodal$mu2_sigma3
canonical_sigma1 <- robust_bimodal$sigma1_sigma3
canonical_sigma2 <- robust_bimodal$sigma2_sigma3
canonical_weight1 <- robust_bimodal$weight1_sigma3

raw_gap <- canonical_mu2 - canonical_mu1
wrap_crossing <- raw_gap > 24

if (any(wrap_crossing)) {

  tmp <- canonical_mu1[wrap_crossing]
  canonical_mu1[wrap_crossing] <- canonical_mu2[wrap_crossing] - 48
  canonical_mu2[wrap_crossing] <- tmp

  tmp <- canonical_sigma1[wrap_crossing]
  canonical_sigma1[wrap_crossing] <- canonical_sigma2[wrap_crossing]
  canonical_sigma2[wrap_crossing] <- tmp

  canonical_weight1[wrap_crossing] <- 1 - canonical_weight1[wrap_crossing]
}

canonical_delta <- canonical_mu2 - canonical_mu1

# ---- Bimodal population starts ----
mu_age_mean_init <- mean(canonical_mu1)
mu_age_sd_init <- max(sd(canonical_mu1), 1)

log_sigma_age_1_init <- log(pmax(canonical_sigma1, 0.5))
log_sigma_age_2_init <- log(pmax(canonical_sigma2, 0.5))

log_sigma_centre_init <- 0.5 * (log_sigma_age_1_init + log_sigma_age_2_init)
log_sigma_difference_init <- log_sigma_age_1_init - log_sigma_age_2_init

log_sigma_centre_mean_init <- mean(log_sigma_centre_init)
log_sigma_centre_sd_init <- max(sd(log_sigma_centre_init), 0.1)

log_sigma_difference_mean_init <- mean(log_sigma_difference_init)
log_sigma_difference_sd_init <- max(sd(log_sigma_difference_init), 0.1)

# ---- Bimodal initialisation ----
make_bimodal_init <- function(chain_id) {
  set.seed(51000 + chain_id)

  mu1_start <- canonical_mu1 + rnorm(N_BIMODAL, 0, 0.15)
  delta_start <- pmin(pmax(canonical_delta + rnorm(N_BIMODAL, 0, 0.10), 0.25), 23.75)
  sigma1_start <- pmin(pmax(canonical_sigma1 + rnorm(N_BIMODAL, 0, 0.05), 0.5), 23.75)
  sigma2_start <- pmin(pmax(canonical_sigma2 + rnorm(N_BIMODAL, 0, 0.05), 0.5), 23.75)
  weight1_start <- pmin(pmax(canonical_weight1 + rnorm(N_BIMODAL, 0, 0.005), 0.02), 0.98)

  z_mu_start <- (mu1_start - mu_age_mean_init) / (RHO_MU_AGE * mu_age_sd_init)
  z_mu_start[!is.finite(z_mu_start)] <- 0
  z_mu_start <- pmax(pmin(z_mu_start, 5), -5)

  log_sigma_centre_start <- 0.5 * (log(sigma1_start) + log(sigma2_start))
  log_sigma_difference_start <- log(sigma1_start) - log(sigma2_start)

  z_log_sigma_centre_start <- (log_sigma_centre_start - log_sigma_centre_mean_init) /
                               log_sigma_centre_sd_init

  z_log_sigma_difference_start <- (log_sigma_difference_start -
                                   log_sigma_difference_mean_init) /
                                   log_sigma_difference_sd_init

  z_log_sigma_centre_start[!is.finite(z_log_sigma_centre_start)] <- 0
  z_log_sigma_difference_start[!is.finite(z_log_sigma_difference_start)] <- 0

  z_log_sigma_centre_start <- pmax(pmin(z_log_sigma_centre_start, 5), -5)
  z_log_sigma_difference_start <- pmax(pmin(z_log_sigma_difference_start, 5), -5)

  list(mu_age_1 = mu1_start, z_mu_age = z_mu_start, mu_age_mean = mu_age_mean_init,
       mu_age_sd = mu_age_sd_init, delta_age = delta_start,
       z_log_sigma_centre = z_log_sigma_centre_start,
       z_log_sigma_difference = z_log_sigma_difference_start,
       log_sigma_centre_mean = log_sigma_centre_mean_init,
       log_sigma_centre_sd = log_sigma_centre_sd_init,
       log_sigma_difference_mean = log_sigma_difference_mean_init,
       log_sigma_difference_sd = log_sigma_difference_sd_init, weight1 = weight1_start)
}

bimodal_init_list <- lapply(seq_len(N_CHAINS), make_bimodal_init)

# ---- Fit bimodal model ----
cat("\nCompiling bimodal model...\n")
bimodal_model <- stan_model(BIMODAL_MODEL_FILE)

cat("Running bimodal fit...\n")
bimodal_start_time <- Sys.time()

fit_bimodal <- sampling(object = bimodal_model,
                        data = list(K_stages = 7, n = N_BIMODAL,
                                    stage_obs = as.matrix(bimodal_stage), aa = AA,
                                    sequestration_profile = SEQUESTRATION_PROFILE,
                                    rho_mu_age = RHO_MU_AGE),
                        seed = 51001, warmup = 1000, iter = 5000, chains = N_CHAINS,
                        init = bimodal_init_list, refresh = 100,
                        control = list(adapt_delta = 0.99, max_treedepth = MAX_TREEDEPTH))

cat("Bimodal runtime: ",
    round(as.numeric(difftime(Sys.time(), bimodal_start_time, units = "mins")), 2),
    " minutes\n", sep = "")

saveRDS(fit_bimodal, BIMODAL_FIT_FILE)

# ---- Bimodal summary ----
bimodal_posterior <- rstan::extract(fit_bimodal,
                                    pars = c("mu_age_1_mod48", "mu_age_2_mod48",
                                             "delta_age", "sigma_age_1", "sigma_age_2",
                                             "weight1"))

mu_age_1 <- numeric(N_BIMODAL)
mu_age_2 <- numeric(N_BIMODAL)
sigma_age_1 <- numeric(N_BIMODAL)
sigma_age_2 <- numeric(N_BIMODAL)
weight1_mean <- numeric(N_BIMODAL)
delta_age_mean <- numeric(N_BIMODAL)
bimodal_mu_age <- numeric(N_BIMODAL)
bimodal_sigma_age <- numeric(N_BIMODAL)

for (i in seq_len(N_BIMODAL)) {

  mu1_draw <- bimodal_posterior$mu_age_1_mod48[, i]
  mu2_draw <- bimodal_posterior$mu_age_2_mod48[, i]
  sigma1_draw <- bimodal_posterior$sigma_age_1[, i]
  sigma2_draw <- bimodal_posterior$sigma_age_2[, i]
  weight1_draw <- bimodal_posterior$weight1[, i]
  delta_draw <- bimodal_posterior$delta_age[, i]

  mu_age_1[i] <- circular_mean_48(mu1_draw)
  mu_age_2[i] <- circular_mean_48(mu2_draw)
  sigma_age_1[i] <- mean(sigma1_draw)
  sigma_age_2[i] <- mean(sigma2_draw)
  weight1_mean[i] <- mean(weight1_draw)
  delta_age_mean[i] <- mean(delta_draw)

  mixture_mu_draw <- mapply(mixture_circular_centre_48, mu1 = mu1_draw, mu2 = mu2_draw,
                            weight1 = weight1_draw)

  bimodal_mu_age[i] <- circular_mean_48(mixture_mu_draw)

  mixture_sigma_draw <- mapply(mixture_sd_draw, mu1 = mu1_draw, delta = delta_draw,
                               sigma1 = sigma1_draw, sigma2 = sigma2_draw,
                               weight1 = weight1_draw)

  bimodal_sigma_age[i] <- mean(mixture_sigma_draw)
}

bimodal_membership <- data.frame(patient = robust_bimodal$patient, ID = bimodal_meta$ID,
                                 mu_age = bimodal_mu_age, sigma_age = bimodal_sigma_age,
                                 mu_age_1 = mu_age_1, sigma_age_1 = sigma_age_1,
                                 mu_age_2 = mu_age_2, sigma_age_2 = sigma_age_2,
                                 weight1 = weight1_mean, delta_age = delta_age_mean)

bimodal_membership$weight2 <- 1 - bimodal_membership$weight1

bimodal_membership <- merge(bimodal_membership,
                            robust_bimodal[, c("ID", "delta_BIC_sigma3",
                                               "circular_separation_sigma3",
                                               "weight1_sigma3", "mu1_sigma3",
                                               "mu2_sigma3"), drop = FALSE],
                            by = "ID", all.x = TRUE, sort = FALSE)

bimodal_membership <- bimodal_membership[match(bimodal_meta$ID, bimodal_membership$ID), ,
                                         drop = FALSE]

rownames(bimodal_membership) <- NULL
saveRDS(bimodal_membership, BIMODAL_MEMBERSHIP_FILE)

bimodal_component_summary <- bimodal_membership[, c("ID", "mu_age", "sigma_age", "mu_age_1",
                                                    "sigma_age_1", "mu_age_2", "sigma_age_2",
                                                    "weight1", "weight2", "delta_age")]

saveRDS(bimodal_component_summary, BIMODAL_COMPONENT_SUMMARY_FILE)
rm(bimodal_posterior, bimodal_model)
gc()

# ---- Combine model outputs ----
cat("\nCombining hierarchical outputs...\n")

bimodal_out <- as.data.table(copy(bimodal_component_summary))
bimodal_out[, model := "bimodal"]

setnames(bimodal_out, old = c("weight1", "weight2"), new = c("weight_1", "weight_2"))

bimodal_out <- bimodal_out[, .(ID, model, mu_age, sigma_age, mu_age_1, sigma_age_1,
                               weight_1, mu_age_2, sigma_age_2, weight_2, delta_age)]

unimodal_out <- as.data.table(copy(unimodal_patient_summary))
unimodal_out[, model := "unimodal"]

unimodal_out[, `:=`(mu_age_1 = mu_age, sigma_age_1 = sigma_age, weight_1 = NA_real_,
                    mu_age_2 = NA_real_, sigma_age_2 = NA_real_, weight_2 = NA_real_,
                    delta_age = NA_real_)]

unimodal_out <- unimodal_out[, .(ID, model, mu_age, sigma_age, mu_age_1, sigma_age_1,
                                 weight_1, mu_age_2, sigma_age_2, weight_2, delta_age)]

hierarchical_model <- rbindlist(list(bimodal_out, unimodal_out),
                                use.names = TRUE, fill = TRUE)

# ---- Observed microscopy predictors ----
dat_all[, tring_obs := 100 * tring / total_stage]
dat_all[, mts_obs := 100 * (etroph + mtroph + ltroph + schizont) / total_stage]

metadata <- dat_all[, .(ID, age, coma, aconvulsions, hb, hct, aanaemia, ashockcomp, temp,
                        ahypoglyc, outcome, site, Trial, be, bun, parasitaemia, tring_obs,
                        mts_obs)]

# ---- Final dataset ----
hierarchical_output <- merge(metadata, hierarchical_model, by = "ID", all.x = TRUE,
                             all.y = FALSE, sort = FALSE)

hierarchical_output <- hierarchical_output[match(dat_all$ID, hierarchical_output$ID)]

# ---- Final column order ----
setcolorder(hierarchical_output, c("ID", "age", "coma", "aconvulsions", "hb", "hct",
                                   "aanaemia", "ashockcomp", "temp", "ahypoglyc", "outcome",
                                   "site", "Trial", "be", "bun", "parasitaemia", "model",
                                   "mu_age", "sigma_age", "mu_age_1", "sigma_age_1",
                                   "weight_1", "mu_age_2", "sigma_age_2", "weight_2",
                                   "delta_age", "tring_obs", "mts_obs"))

# ---- Save final output ----
fwrite(hierarchical_output, HIERARCHICAL_OUTPUT_FILE)

# ---- Final summary ----
cat("\nHierarchical modelling complete.\n")
cat("Total patients: ", nrow(hierarchical_output), "\n", sep = "")
cat("Unimodal patients: ", sum(hierarchical_output$model == "unimodal"), "\n", sep = "")
cat("Bimodal patients: ", sum(hierarchical_output$model == "bimodal"), "\n", sep = "")

cat("mu_age range: ", paste(round(range(hierarchical_output$mu_age), 3),
                            collapse = " - "), " h\n", sep = "")

cat("sigma_age range: ", paste(round(range(hierarchical_output$sigma_age), 3),
                               collapse = " - "), " h\n", sep = "")

cat("\nFinal dataset:\n", HIERARCHICAL_OUTPUT_FILE, "\n", sep = "")