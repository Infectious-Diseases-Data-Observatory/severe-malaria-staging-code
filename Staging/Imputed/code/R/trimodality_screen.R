# IMPUTED COHORT TRIMODALITY SCREEN

library(data.table)
options(scipen = 999)

# ---- Paths ----
DATA_DIR <- "Data"
PROJECT_DIR <- file.path("Staging/Imputed")

OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
CSV_DIR <- file.path(OUTPUT_DIR, "CSV")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")

dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Files ----
DATA_FILE <- file.path(DATA_DIR, "Imputed_dataset_long.csv")

# ---- Screening settings ----
SIGMA_MAX <- 24
DF_T <- 10
MIN_WEIGHT <- 0.10
MIN_SEPARATION <- 6
DELTA_BIC_THRESHOLD <- 10

# Fixed stage-specific sequestration profile.
# Values represent the relative contribution of each latent parasite stage
# to the peripheral-blood stage distribution. Progressive reduction with
# parasite maturity represents sequestration of infected red blood cells
# from the peripheral circulation.
SEQUESTRATION_PROFILE <- c(1.000, 0.500, 0.250, 0.200, 0.035, 0.027, 0.005)

stage_names <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")
aa <- c(0, 6, 16, 26, 30, 34, 38, 48)

# ---- Load multiply imputed dataset ----
cat("\nLoading multiply imputed dataset...\n")
dat_all_imputed <- fread(DATA_FILE)

# Sort by imputation
imp_ids <- sort(unique(dat_all_imputed$.imp))

# ---- Construct canonical staging cohort ----
# The staging variables are identical across all 20 imputations
# The bimodality screen is run once using imputation 1 as the canonical copy
dat_all_mix <- copy(dat_all_imputed[.imp == imp_ids[1]])

# Preserve the canonical order explicitly by .id
setorder(dat_all_mix, .id)

# The 93-107 restriction was applied before multiple imputation
# So all patients should already satisfy this restriction
dat_all_mix[, total_stage := tring + sring + lring + etroph + mtroph + ltroph + schizont]

dat_stage <- dat_all_mix[, ..stage_names]

metadata <- dat_all_mix[, .(patient = seq_len(.N), .id, ID, age, coma, aconvulsions, hb,
                            hct, aanaemia, ashockcomp, temp, ahypoglyc, outcome, site,
                            trial, be, bun, parasitaemia, tring, sring, lring, etroph,
                            mtroph, ltroph, schizont)]

# ---- Load 2-component screening results ----
two_component_file <- file.path(RDS_DIR, "hierarchical_bimodality_screen.rds")
existing_two <- readRDS(two_component_file)

# ---- Circular functions ----
mod48 <- function(x) {
  ((x %% 48) + 48) %% 48
}

circular_distance <- function(a, b) {
  d <- abs(mod48(a) - mod48(b))
  pmin(d, 48 - d)
}

# ---- Component peripheral stage probabilities under sequestration ----
component_stage_probabilities <- function(mu, sigma) {
  p <- numeric(7)

  for (j in seq_len(7)) {

    # One 48-hour cycle before and after the current cycle
    for (m in -1:1) {
      upper <- (aa[j + 1] + 48*m - mu) / sigma
      lower <- (aa[j] + 48*m - mu) / sigma

      p[j] <- p[j] + pt(upper, df = DF_T) - pt(lower, df = DF_T)
    }
  }

  # Apply the fixed stage-specific sequestration profile before
  # normalising to the predicted peripheral stage distribution
  p <- p * SEQUESTRATION_PROFILE

  total <- sum(p)
  p / total
}

# ---- Bounded SD transformations ----
sigma_from_raw <- function(raw, sigma_min) {
  sigma_min + (SIGMA_MAX - sigma_min) * plogis(raw)
}

sigma_to_raw <- function(sigma, sigma_min) {
  sigma <- min(max(as.numeric(sigma), sigma_min + 1e-6), SIGMA_MAX - 1e-6)

  p <- (sigma - sigma_min) / (SIGMA_MAX - sigma_min)
  qlogis(p)
}

# ---- Parameter transformation ----
unpack_three <- function(par, sigma_min) {
  mus <- as.numeric(par[1:3])
  sigmas <- sigma_from_raw(par[4:6], sigma_min)

  eta <- as.numeric(par[7:9])
  eta <- eta - max(eta)
  ew <- exp(eta)
  weights <- ew / sum(ew)

  list(mu = mus, sigma = sigmas, weight = weights)
}

# ---- Three-component negative log-likelihood ----
three_component_nll <- function(par, counts, sigma_min) {
  u <- unpack_three(par, sigma_min)
  if (is.null(u)) return(1e100)

  probs <- matrix(NA_real_, nrow = 3, ncol = 7)

  for (k in 1:3) {
    probs[k, ] <- component_stage_probabilities(u$mu[k], u$sigma[k])
  }

  if (any(!is.finite(probs))) return(1e100)

  mixture_p <- as.numeric(u$weight %*% probs)
  mixture_p <- pmax(mixture_p, 1e-300)
  mixture_p <- mixture_p / sum(mixture_p)

  ll <- sum(counts * log(mixture_p))
  if (!is.finite(ll)) return(1e100)

  -ll
}

# ---- Create starting values from 2-component fit ----
weights_to_eta <- function(w) {
  w <- pmax(as.numeric(w), 1e-4)
  w <- w / sum(w)
  log(w)
}

make_starts <- function(two_row, sigma_min) {
  mu1 <- as.numeric(two_row$mu1)
  mu2 <- as.numeric(two_row$mu2)

  s1 <- max(as.numeric(two_row$sigma1), sigma_min + 0.10)
  s2 <- max(as.numeric(two_row$sigma2), sigma_min + 0.10)

  s1 <- min(s1, SIGMA_MAX - 0.10)
  s2 <- min(s2, SIGMA_MAX - 0.10)

  w1 <- as.numeric(two_row$weight1)
  if (!is.finite(w1)) w1 <- 0.5
  w1 <- min(max(w1, 0.05), 0.95)
  w2 <- 1 - w1

  candidate_specs <- list(list(mu = c(mu1 - 2, mu1 + 2, mu2), sigma = c(s1, s1, s2),
                               weight = c(w1 / 2, w1 / 2, w2)),
                          list(mu = c(mu1, mu2 - 2, mu2 + 2), sigma = c(s1, s2, s2),
                               weight = c(w1, w2 / 2, w2 / 2)),
                          list(mu = c(mu1, mu2, mu1 + 0.5 * (mu2 - mu1)),
                               sigma = c(s1, s2, max((s1 + s2) / 2, sigma_min + 0.10)),
                               weight = c(0.45 * w1, 0.45 * w2, 0.55)),
                          list(mu = c(mu1, mu2, mu1 + 16), sigma = c(s1, s2, max(s1, s2)),
                               weight = c(0.45 * w1, 0.45 * w2, 0.55)),
                          list(mu = c(mu1, mu2, mu2 + 16), sigma = c(s1, s2, max(s1, s2)),
                               weight = c(0.45 * w1, 0.45 * w2, 0.55)))

  lapply(candidate_specs, function(z) {
    z$weight <- z$weight / sum(z$weight)

    c(as.numeric(z$mu), vapply(z$sigma, sigma_to_raw, numeric(1), sigma_min = sigma_min),
      weights_to_eta(z$weight))})
}

# ---- Fit one patient ----
fit_three_patient <- function(counts, two_row, sigma_min) {
  starts <- make_starts(two_row, sigma_min)

  best_fit <- NULL
  best_value <- Inf

  for (start in starts) {
    fit <- tryCatch(optim(par = start, fn = three_component_nll, counts = counts,
                          sigma_min = sigma_min, method = "BFGS",
                          control = list(maxit = 1500, reltol = 1e-10)),
                    error = function(e) NULL)

    if (!is.null(fit) && length(fit$value) == 1 &&
        is.finite(fit$value) && fit$value < best_value) {
      best_fit <- fit
      best_value <- fit$value
    }
  }

  if (is.null(best_fit)) return(list(success = FALSE))

  u <- unpack_three(best_fit$par, sigma_min)
  if (is.null(u)) return(list(success = FALSE))

  phases <- mod48(u$mu)
  ord <- order(phases)

  phases <- phases[ord]
  sigmas <- u$sigma[ord]
  weights <- u$weight[ord]

  sep12 <- circular_distance(phases[1], phases[2])
  sep13 <- circular_distance(phases[1], phases[3])
  sep23 <- circular_distance(phases[2], phases[3])

  loglik3 <- -best_fit$value
  k3 <- 8
  n_cells <- sum(counts)
  bic3 <- -2 * loglik3 + k3 * log(n_cells)

  list(success = TRUE, convergence = best_fit$convergence, loglik3 = loglik3, BIC_3 = bic3,
       mu1 = phases[1], sigma1 = sigmas[1], weight1 = weights[1],
       mu2 = phases[2], sigma2 = sigmas[2], weight2 = weights[2],
       mu3 = phases[3], sigma3 = sigmas[3], weight3 = weights[3],
       separation12 = sep12, separation13 = sep13, separation23 = sep23)
}

# ---- Run one trimodality screen ----
run_sigma_screen <- function(sigma_min, existing_two_table) {
  cat("\nTrimodal screen: sigma minimum = ", sigma_min,
      " hours; maximum = ", SIGMA_MAX, " hours\n\n", sep = "")

  results <- data.frame(patient = seq_len(nrow(dat_stage)), ID = dat_all_mix$ID,
                        n_cells = rowSums(as.matrix(dat_stage)), BIC_2 = NA_real_,
                        loglik3 = NA_real_, BIC_3 = NA_real_, delta_BIC_2_to_3 = NA_real_,
                        mu1 = NA_real_, sigma1 = NA_real_, weight1 = NA_real_,
                        mu2 = NA_real_, sigma2 = NA_real_, weight2 = NA_real_,
                        mu3 = NA_real_, sigma3 = NA_real_, weight3 = NA_real_,
                        separation12 = NA_real_, separation13 = NA_real_,
                        separation23 = NA_real_, convergence = NA_integer_,
                        trimodal_candidate = FALSE)

  match_index <- match(results$patient, existing_two_table$patient)
  results$BIC_2 <- as.numeric(existing_two_table$mixture_BIC[match_index])

  for (i in seq_len(nrow(dat_stage))) {
    counts <- as.numeric(unlist(dat_stage[i, ..stage_names], use.names = FALSE))

    two_row <- existing_two_table[match_index[i], , drop = FALSE]

    fit3 <- fit_three_patient(counts = counts, two_row = two_row, sigma_min = sigma_min)

    if (isTRUE(fit3$success)) {
      delta32 <- as.numeric(results$BIC_2[i]) - as.numeric(fit3$BIC_3)

      results$loglik3[i] <- fit3$loglik3
      results$BIC_3[i] <- fit3$BIC_3
      results$delta_BIC_2_to_3[i] <- delta32

      results$mu1[i] <- fit3$mu1
      results$sigma1[i] <- fit3$sigma1
      results$weight1[i] <- fit3$weight1
      results$mu2[i] <- fit3$mu2
      results$sigma2[i] <- fit3$sigma2
      results$weight2[i] <- fit3$weight2
      results$mu3[i] <- fit3$mu3
      results$sigma3[i] <- fit3$sigma3
      results$weight3[i] <- fit3$weight3

      results$separation12[i] <- fit3$separation12
      results$separation13[i] <- fit3$separation13
      results$separation23[i] <- fit3$separation23
      results$convergence[i] <- fit3$convergence

      results$trimodal_candidate[i] <-
        is.finite(delta32) &&
        delta32 > DELTA_BIC_THRESHOLD &&
        all(c(fit3$weight1, fit3$weight2, fit3$weight3) >= MIN_WEIGHT) &&
        all(c(fit3$separation12, fit3$separation13, fit3$separation23) >= MIN_SEPARATION) &&
        all(c(fit3$sigma1, fit3$sigma2, fit3$sigma3) >= sigma_min) &&
        all(c(fit3$sigma1, fit3$sigma2, fit3$sigma3) <= SIGMA_MAX) &&
        fit3$convergence == 0
    }

    if (i %% 100 == 0 || i == nrow(dat_stage)) {
      cat("Completed ", i, " / ", nrow(dat_stage), " patients; current candidates = ",
          sum(results$trimodal_candidate[seq_len(i)], na.rm = TRUE), "\n", sep = "")
    }
  }

  cat("\nSigma minimum ", sigma_min,
      " h: ", sum(results$trimodal_candidate, na.rm = TRUE),
      " / ", nrow(dat_stage), " trimodal candidates (",
      round(100 * mean(results$trimodal_candidate, na.rm = TRUE), 2), "%)\n", sep = "")

  results
}

# ---- Run trimodality screen for sigma = 1, 2, 3 ----
all_trimodal_results <- list()

for (sigma_min in c(1, 2, 3)) {
  existing_name <- paste0("sigma_min_", sigma_min)

  all_trimodal_results[[existing_name]] <- run_sigma_screen(sigma_min = sigma_min,
                                                            existing_two_table = existing_two[[existing_name]])
}

# ---- Robustness classification ----
r1 <- all_trimodal_results$sigma_min_1
r2 <- all_trimodal_results$sigma_min_2
r3 <- all_trimodal_results$sigma_min_3

trimodality_classification <- data.frame(patient = r1$patient, ID = r1$ID,
                                         pass_sigma1 = r1$trimodal_candidate,
                                         pass_sigma2 = r2$trimodal_candidate,
                                         pass_sigma3 = r3$trimodal_candidate,
                                         delta_BIC_sigma1 = r1$delta_BIC_2_to_3,
                                         delta_BIC_sigma2 = r2$delta_BIC_2_to_3,
                                         delta_BIC_sigma3 = r3$delta_BIC_2_to_3)

for (nm in c("pass_sigma1", "pass_sigma2", "pass_sigma3")) {
  trimodality_classification[[nm]][is.na(trimodality_classification[[nm]])] <- FALSE
}

trimodality_classification$n_screens_passed <- rowSums(
  trimodality_classification[, c("pass_sigma1", "pass_sigma2", "pass_sigma3")]
)

trimodality_classification$category <- ifelse(
  trimodality_classification$n_screens_passed == 3, "Robust trimodal",
  ifelse(trimodality_classification$n_screens_passed == 2, "Probable trimodal",
         ifelse(trimodality_classification$n_screens_passed == 1, "Possible trimodal",
                "No strong trimodal evidence")))

trimodality_classification$category <- factor(trimodality_classification$category,
                                              levels = c("Robust trimodal",
                                                         "Probable trimodal",
                                                         "Possible trimodal",
                                                         "No strong trimodal evidence"))

# ---- Summary ----
category_counts <- as.data.frame(table(trimodality_classification$category))
names(category_counts) <- c("category", "n_patients")
category_counts$percentage <- 100 * category_counts$n_patients /
  nrow(trimodality_classification)

cat("\nFinal trimodality classification\n\n")
print(category_counts, row.names = FALSE)

# ---- Exact pass patterns ----
pass_patterns <- aggregate(patient ~ pass_sigma1 + pass_sigma2 + pass_sigma3,
                           data = trimodality_classification, FUN = length)

names(pass_patterns)[names(pass_patterns) == "patient"] <- "n_patients"
pass_patterns <- pass_patterns[order(pass_patterns$n_patients, decreasing = TRUE), ]

cat("\nExact trimodal sensitivity patterns\n\n")
print(pass_patterns, row.names = FALSE)

# ---- Robust trimodal patients ----
robust_trimodal <- trimodality_classification[
  trimodality_classification$category == "Robust trimodal",
]

robust_trimodal <- robust_trimodal[order(robust_trimodal$delta_BIC_sigma3,
                                         decreasing = TRUE, na.last = TRUE), ]

cat("\nRobust trimodal candidates\n\n")
cat("Number robust trimodal: ", nrow(robust_trimodal), "\n\n", sep = "")

if (nrow(robust_trimodal) > 0) {
  print(robust_trimodal, row.names = FALSE)
} else {
  cat("None.\n")
}

# ---- Save outputs ----

TRIMODALITY_SCREEN_FILE <- file.path(RDS_DIR, "hierarchical_trimodality_screen.rds")

TRIMODALITY_CLASSIFICATION_FILE <- file.path(RDS_DIR,
                                             "hierarchical_trimodality_classification.rds")

TRIMODALITY_CLASSIFICATION_CSV <- file.path(CSV_DIR,
                                            "hierarchical_trimodality_classification.csv")

TRIMODAL_CANDIDATES_CSV <- file.path(CSV_DIR, "hierarchical_trimodal_candidates.csv")

saveRDS(all_trimodal_results, TRIMODALITY_SCREEN_FILE)
saveRDS(trimodality_classification, TRIMODALITY_CLASSIFICATION_FILE)
fwrite(trimodality_classification, TRIMODALITY_CLASSIFICATION_CSV)

trimodal_candidates <- trimodality_classification[trimodality_classification$pass_sigma1 |
                                                    trimodality_classification$pass_sigma2 |
                                                    trimodality_classification$pass_sigma3,]

trimodal_candidates <- merge(metadata, trimodal_candidates, by = c("patient", "ID"),
                             sort = FALSE)

fwrite(trimodal_candidates, TRIMODAL_CANDIDATES_CSV)

cat("\nSaved ", nrow(trimodal_candidates), " trimodal candidate patients to:\n",
    TRIMODAL_CANDIDATES_CSV, "\n", sep = "")

cat("\nSaved trimodality screening results.\n")
cat("Robust trimodal patients: ", nrow(robust_trimodal), "\n", sep = "")