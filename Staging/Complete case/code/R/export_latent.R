# EXPORT POSTERIOR-MEAN LATENT STAGE PROPORTIONS
# COMPLETE CASE

library(data.table)
library(rstan)
options(scipen = 999)

# ---- Paths ----
DATA_DIR <- "Data"
PROJECT_DIR <- file.path("Staging/Complete case")

OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
CSV_DIR <- file.path(OUTPUT_DIR, "CSV")

# Number of posterior draws used to estimate each posterior-mean latent stage proportion
N_POSTERIOR_DRAWS <- 250L

# ---- Input and output files ----
UNIMODAL_FIT_FILE <- file.path(RDS_DIR, "fit_hierarchical_unimodal.rds")
BIMODAL_FIT_FILE <- file.path(RDS_DIR, "fit_hierarchical_bimodal.rds")
CLASSIFICATION_FILE <- file.path(RDS_DIR, "hierarchical_bimodality_classification.rds")
COMBINED_FILE <- file.path(CSV_DIR, "dat_all_hierarchical.csv")
OUTPUT_FILE <- file.path(CSV_DIR, "stage_proportions.csv")

# ---- Reconstruct patient ordering used by Stan ----
dat_all <- fread(file.path(DATA_DIR, "dat_all.csv"))

required_cols <- c("ID", "age", "coma", "aconvulsions", "hb", "hct", "aanaemia",
                   "ashockcomp", "temp", "ahypoglyc", "outcome", "site", "Trial", "be",
                   "bun", "parasitaemia", "tring", "sring", "lring", "etroph", "mtroph",
                   "ltroph", "schizont")

dat_all <- dat_all[, ..required_cols]
dat_all <- dat_all[complete.cases(dat_all)]
dat_all[, total_stage := tring + sring + lring + etroph + mtroph + ltroph + schizont]
dat_all <- dat_all[total_stage >= 93 & total_stage <= 107]
dat_all[, total_stage := NULL]

classification <- readRDS(CLASSIFICATION_FILE)

robust_bimodal <- classification[as.character(classification$category) == "Robust bimodal",
                                 , drop = FALSE]

robust_bimodal <- robust_bimodal[order(robust_bimodal$patient), , drop = FALSE]

bimodal_rows <- robust_bimodal$patient
unimodal_rows <- setdiff(seq_len(nrow(dat_all)), bimodal_rows)

bimodal_ids <- as.character(dat_all$ID[bimodal_rows])
unimodal_ids <- as.character(dat_all$ID[unimodal_rows])

N_BIMODAL <- length(bimodal_ids)
N_UNIMODAL <- length(unimodal_ids)

# ---- Functions ----
DF_T <- 10
aa <- c(0, 6, 16, 26, 30, 34, 38, 48)

latent_names <- c("tring_latent", "sring_latent", "lring_latent", "etroph_latent",
                  "mtroph_latent", "ltroph_latent", "schizont_latent")

# Fixed stage-specific sequestration profile.
# Values represent the relative contribution of each latent parasite stage
# to the peripheral-blood stage distribution.
SEQUESTRATION_PROFILE <- c(1.000, 0.500, 0.250, 0.200, 0.035, 0.027, 0.005)

peripheral_names <- c("tring_peripheral", "sring_peripheral", "lring_peripheral",
                      "etroph_peripheral", "mtroph_peripheral", "ltroph_peripheral",
                      "schizont_peripheral")

select_draws <- function(n, target) {
  if (target >= n) {
    return(seq_len(n))
  }
  unique(as.integer(round(seq(1, n, length.out = target))))
}

latent_stage_probs <- function(mu, sigma) {
  p <- numeric(7)

  for (j in 1:7) {
    for (m in -1:1) {
      upper <- (aa[j + 1] + m*48 - mu) / sigma
      lower <- (aa[j] + m*48 - mu) / sigma

      p[j] <- p[j] + pt(upper, df = DF_T) - pt(lower, df = DF_T)
    }
  }
  p / sum(p)
}

# Convert latent stage probabilities to predicted peripheral stage
# probabilities under the fixed stage-specific sequestration assumption
peripheral_stage_probs <- function(latent_probs) {
  p <- latent_probs * SEQUESTRATION_PROFILE
  p / sum(p)
}

# ---- Unimodal latent and predicted peripheral stages ----
cat("\nLoading unimodal fit...\n")
fit_unimodal <- readRDS(UNIMODAL_FIT_FILE)

post_u <- rstan::extract(fit_unimodal, pars = c("mu_age", "sigma_age"))
draw_u <- select_draws(nrow(post_u$mu_age), N_POSTERIOR_DRAWS)
cat("Using ", length(draw_u), " unimodal posterior draws.\n", sep = "")

latent_u <- matrix(0, nrow = N_UNIMODAL, ncol = 7)
peripheral_u <- matrix(0, nrow = N_UNIMODAL, ncol = 7)

for (d in draw_u) {
  for (i in seq_len(N_UNIMODAL)) {
    p_latent <- latent_stage_probs(post_u$mu_age[d, i], post_u$sigma_age[d, i])

    # Apply the sequestration assumption to obtain the corresponding
    # predicted peripheral stage distribution
    p_peripheral <- peripheral_stage_probs(p_latent)

    latent_u[i, ] <- latent_u[i, ] + p_latent
    peripheral_u[i, ] <- peripheral_u[i, ] + p_peripheral
  }
}

latent_u <- latent_u / length(draw_u)
peripheral_u <- peripheral_u / length(draw_u)

colnames(latent_u) <- latent_names
colnames(peripheral_u) <- peripheral_names

latent_u <- data.table(ID = unimodal_ids, model = "unimodal", latent_u * 100,
                       peripheral_u * 100)

rm(fit_unimodal, post_u, draw_u)
gc()

# ---- Bimodal latent and predicted peripheral stages ----
cat("\nLoading bimodal fit...\n")
fit_bimodal <- readRDS(BIMODAL_FIT_FILE)

post_b <- rstan::extract(fit_bimodal, pars = c("mu_age_1_mod48", "mu_age_2_mod48",
                                               "sigma_age_1", "sigma_age_2", "weight1"))

draw_b <- select_draws(nrow(post_b$mu_age_1_mod48), N_POSTERIOR_DRAWS)
cat("Using ", length(draw_b), " bimodal posterior draws.\n", sep = "")

latent_b <- matrix(0, nrow = N_BIMODAL, ncol = 7)
peripheral_b <- matrix(0, nrow = N_BIMODAL, ncol = 7)

for (d in draw_b) {
  for (i in seq_len(N_BIMODAL)) {

    p1 <- latent_stage_probs(post_b$mu_age_1_mod48[d, i], post_b$sigma_age_1[d, i])
    p2 <- latent_stage_probs(post_b$mu_age_2_mod48[d, i], post_b$sigma_age_2[d, i])

    w1 <- post_b$weight1[d, i]

    # Combine the latent component distributions using the mixture weight
    p_latent <- w1*p1 + (1 - w1)*p2
    
    # Apply sequestration and normalise each component separately
    # before combining, consistent with the Stan mixture likelihood
    p_peripheral <- w1*peripheral_stage_probs(p1) + (1 - w1)*peripheral_stage_probs(p2)
    
    latent_b[i, ] <- latent_b[i, ] + p_latent
    peripheral_b[i, ] <- peripheral_b[i, ] + p_peripheral
  }
}

latent_b <- latent_b / length(draw_b)
peripheral_b <- peripheral_b / length(draw_b)

colnames(latent_b) <- latent_names
colnames(peripheral_b) <- peripheral_names

latent_b <- data.table(ID = bimodal_ids, model = "bimodal", latent_b * 100,
                       peripheral_b * 100)

rm(fit_bimodal, post_b)
gc()

# ---- Combine ----
latent <- rbindlist(list(latent_u, latent_b), use.names = TRUE)

# --- Attach final mu_age / sigma_age / Trial ----
combined <- fread(COMBINED_FILE)

latent <- merge(combined[, .(ID, Trial, mu_age, sigma_age)],
                latent, by = "ID", all.x = TRUE, sort = FALSE)

setcolorder(latent, c("ID", "Trial", "model", "mu_age", "sigma_age", latent_names,
                      peripheral_names))

# ---- Save ----
fwrite(latent, OUTPUT_FILE)
cat("\nSaved posterior-mean latent and model-predicted peripheral stage proportions to:\n",
    OUTPUT_FILE, "\n", sep = "")