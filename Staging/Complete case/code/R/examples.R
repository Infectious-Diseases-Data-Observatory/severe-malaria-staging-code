# STAGING MODEL PRESENTATION EXAMPLES
# COMPLETE CASE

library(data.table)
library(rstan)

# ---- Paths ----
DATA_DIR <- "Data"
PROJECT_DIR <- file.path("Staging/Complete case")

OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
CSV_DIR <- file.path(OUTPUT_DIR, "CSV")
EXAMPLES_DIR <- file.path(OUTPUT_DIR, "examples")
dir.create(EXAMPLES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Input files ----
DATA_FILE <- file.path(DATA_DIR, "dat_all.csv")
HIERARCHICAL_FILE <- file.path(CSV_DIR, "dat_all_hierarchical.csv")
BIMODAL_FIT_FILE <- file.path(OUTPUT_DIR, "RDS", "fit_hierarchical_bimodal.rds")

BIMODAL_SUMMARY_FILE <- file.path(OUTPUT_DIR, "RDS",
                                  "hierarchical_bimodal_component_summary.rds")

# ---- Load data ----
dat_all <- fread(DATA_FILE)
dat_hierarchical <- fread(HIERARCHICAL_FILE)

# ---- Define presentation examples ----
example_1_ID <- "AQ_I0015"
example_2_ID <- "AQ_Z0818"
example_3_ID <- "SQ_C0168"

# ---- Extract final hierarchical parameter estimates ----
example_1 <- dat_hierarchical[ID == example_1_ID]
example_2 <- dat_hierarchical[ID == example_2_ID]
example_3 <- dat_hierarchical[ID == example_3_ID]

cat("\nPresentation examples:\n")

print(rbindlist(list(example_1, example_2, example_3), fill = TRUE)[,
    .(ID, Trial, model, mu_age, sigma_age, mu_age_1, sigma_age_1, weight_1, mu_age_2,
      sigma_age_2, weight_2, delta_age)])

# ---- Extract observed peripheral staging ----
stage_names <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")
stage_breaks <- c(0, 6, 16, 26, 30, 34, 38, 48)

AQ_colour <- "#E31B23"
SQ_colour <- "#1481E7"

get_observed_staging <- function(patient_ID) {
  observed <- dat_all[ID == patient_ID, c("ID", stage_names), with = FALSE]
  stage_counts <- as.numeric(observed[1, ..stage_names])
  100 * stage_counts / sum(stage_counts)
}

observed_1 <- get_observed_staging(example_1_ID)
observed_2 <- get_observed_staging(example_2_ID)
observed_3 <- get_observed_staging(example_3_ID)

# ---- Define plotting grid ----
hours <- seq(0, 48, length.out = 2000)

x_ticks <- seq(0, 48, by = 4)
y_ticks <- seq(0, 0.30, by = 0.05)
y_ticks_bimodal <- seq(0, 0.35, by = 0.05)

DF_T <- 10

# ---- Calculate wrapped Student-t density ----
wrapped_student_t_density <- function(hours, mu_age, sigma_age) {
  density <- numeric(length(hours))

  for (m in -1:1) {
    density <- density + dt((hours + m * 48 - mu_age) / sigma_age, df = DF_T) / sigma_age
  }
  density
}

# Fixed stage boundaries and peripheral weights
SEQUESTRATION_PROFILE <- c(1, 0.5, 0.25, 0.2, 0.035, 0.027, 0.005)

# Calculate latent stage probabilities
latent_stage_probs <- function(mu, sigma) {
  p <- numeric(7)

  for (j in seq_len(7)) {
    for (m in -1:1) {
      p[j] <- p[j] + pt((stage_breaks[j + 1] + 48*m - mu) / sigma, df = DF_T) -
        pt((stage_breaks[j] + 48*m - mu) / sigma, df = DF_T)
    }
  }
  p / sum(p)
}

# Convert peripheral mixture weight to latent weight
latent_mixture_weight <- function(p1, p2, w) {
  z1 <- sum(p1 * SEQUESTRATION_PROFILE)
  z2 <- sum(p2 * SEQUESTRATION_PROFILE)

  stopifnot(z1 > 0, z2 > 0)
  (w / z1) / (w / z1 + (1 - w) / z2)
}

# ---- Plot observed peripheral staging ----
plot_observed_staging <- function(stage_percentages, patient_ID, trial) {

  if (trial == "AQUAMAT") {
    my_colour <- AQ_colour
  } else {
    my_colour <- SQ_colour
  }

  plot(NA, xlim = c(0, 48), ylim = c(0, 100), xaxs = "i", yaxs = "i",
       xlab = "Parasite age (h)", ylab = "Observed stage proportion (%)",
       main = paste0("Observed peripheral staging: ", patient_ID),
       xaxt = "n", yaxt = "n", bty = "l")

  rect(xleft = stage_breaks[-length(stage_breaks)], ybottom = 0, xright = stage_breaks[-1],
       ytop = stage_percentages, col = my_colour, border = "black")

  axis(1, at = stage_breaks)
  axis(2, at = seq(0, 100, 20), las = 1)
}

# ---- Plot first unimodal example ----
mean_1 <- example_1$mu_age
sd_1 <- example_1$sigma_age
density_1 <- wrapped_student_t_density(hours, mean_1, sd_1)

png(file.path(EXAMPLES_DIR, paste0(example_1_ID, "_latent_distribution.png")),
    width = 2400, height = 1800, res = 300)

plot(hours, density_1, type = "l", lwd = 1, xlim = c(0, 48), ylim = c(0, 0.30),
     xaxs = "i", yaxs = "i", xlab = "Parasite age (h)", ylab = "Density",
     main = paste0("Latent parasite age distribution: ", example_1_ID, "\n",
                   "mu_age = ", round(mean_1, 2), " h; sigma_age = ", round(sd_1, 2), " h"),
     xaxt = "n", yaxt = "n")

axis(1, at = x_ticks)
axis(2, at = y_ticks, labels = sprintf("%.2f", y_ticks), las = 3)
dev.off()

png(file.path(EXAMPLES_DIR, paste0(example_1_ID, "_observed_staging.png")),
    width = 2400, height = 1800, res = 300)

plot_observed_staging(observed_1, example_1_ID, example_1$Trial)
dev.off()

# ---- Plot second unimodal example ----
mean_2 <- example_2$mu_age
sd_2 <- example_2$sigma_age
density_2 <- wrapped_student_t_density(hours, mean_2, sd_2)

png(file.path(EXAMPLES_DIR, paste0(example_2_ID, "_latent_distribution.png")),
    width = 2400, height = 1800, res = 300)

plot(hours, density_2, type = "l", lwd = 1, xlim = c(0, 48), ylim = c(0, 0.30),
     xaxs = "i", yaxs = "i", xlab = "Parasite age (h)", ylab = "Density",
     main = paste0("Latent parasite age distribution: ", example_2_ID, "\n",
                   "mu_age = ", round(mean_2, 2), " h; sigma_age = ", round(sd_2, 2), " h"),
     xaxt = "n", yaxt = "n")

axis(1, at = x_ticks)
axis(2, at = y_ticks, labels = sprintf("%.2f", y_ticks), las = 3)
dev.off()

png(file.path(EXAMPLES_DIR, paste0(example_2_ID, "_observed_staging.png")),
    width = 2400, height = 1800, res = 300)

plot_observed_staging(observed_2, example_2_ID, example_2$Trial)
dev.off()

# ---- Plot bimodal example ----
# Recover the saved Stan patient ordering
bimodal_summary <- readRDS(BIMODAL_SUMMARY_FILE)

patient_index <- match(example_3_ID, as.character(bimodal_summary$ID))
stopifnot(!is.na(patient_index))

# Load the existing bimodal posterior fit
fit_bimodal <- readRDS(BIMODAL_FIT_FILE)

post_b <- rstan::extract(fit_bimodal, pars = c("mu_age_1_mod48", "mu_age_2_mod48",
                                               "sigma_age_1", "sigma_age_2", "weight1"))

# Use the same posterior draws as the heatmaps
select_draws <- function(n, target) {
  if (target >= n) return(seq_len(n))
  unique(as.integer(round(seq(1, n, length.out = target))))
}

draw_b <- select_draws(nrow(post_b$weight1), 250L)
density_3 <- numeric(length(hours))

for (r in seq_along(draw_b)) {
  d <- draw_b[r]

  mu1 <- post_b$mu_age_1_mod48[d, patient_index]
  mu2 <- post_b$mu_age_2_mod48[d, patient_index]
  s1 <- post_b$sigma_age_1[d, patient_index]
  s2 <- post_b$sigma_age_2[d, patient_index]
  w <- post_b$weight1[d, patient_index]

  p1 <- latent_stage_probs(mu1, s1)
  p2 <- latent_stage_probs(mu2, s2)

  w_latent <- latent_mixture_weight(p1, p2, w)

  density_3 <- density_3 + w_latent * wrapped_student_t_density(hours, mu1, s1) +
    (1 - w_latent) * wrapped_student_t_density(hours, mu2, s2)
}

density_3 <- density_3 / length(draw_b)
rm(fit_bimodal, post_b)
gc()

png(file.path(EXAMPLES_DIR, paste0(example_3_ID, "_latent_distribution.png")),
    width = 2400, height = 1800, res = 300)

plot(hours, density_3, type = "l", lwd = 1, xlim = c(0, 48), ylim = c(0, 0.3),
     xaxs = "i", yaxs = "i", xlab = "Parasite age (h)", ylab = "Density",
     main = paste0("Latent bimodal parasite age distribution: ", example_3_ID),
     xaxt = "n", yaxt = "n")

axis(1, at = x_ticks)
axis(2, at = y_ticks_bimodal, labels = sprintf("%.2f", y_ticks_bimodal), las = 3)
dev.off()

png(file.path(EXAMPLES_DIR, paste0(example_3_ID, "_observed_staging.png")),
    width = 2400, height = 1800, res = 300)

plot_observed_staging(observed_3, example_3_ID, example_3$Trial)
dev.off()