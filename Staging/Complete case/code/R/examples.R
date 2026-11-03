# STAGING MODEL PRESENTATION EXAMPLES
# COMPLETE CASE

library(data.table)

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

# ---- Load data ----
dat_all <- fread(DATA_FILE)
dat_hierarchical <- fread(HIERARCHICAL_FILE)

# Specify example patient IDs from the restricted dataset before running this script
example_1_ID <- NA_character_
example_2_ID <- NA_character_
example_3_ID <- NA_character_

stopifnot(!anyNA(c(example_1_ID, example_2_ID, example_3_ID)))

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
     main = paste0("Latent parasite-age distribution: ", example_1_ID, "\n",
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
     main = paste0("Latent parasite-age distribution: ", example_2_ID, "\n",
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
mean_3a <- example_3$mu_age_1
sd_3a <- example_3$sigma_age_1
weight_3a <- example_3$weight_1

mean_3b <- example_3$mu_age_2
sd_3b <- example_3$sigma_age_2
weight_3b <- example_3$weight_2

density_3a <- wrapped_student_t_density(hours, mean_3a, sd_3a)
density_3b <- wrapped_student_t_density(hours, mean_3b, sd_3b)
density_3 <- weight_3a * density_3a + weight_3b * density_3b

png(file.path(EXAMPLES_DIR, paste0(example_3_ID, "_latent_distribution.png")),
    width = 2400, height = 1800, res = 300)

plot(hours, density_3, type = "l", lwd = 1, xlim = c(0, 48), ylim = c(0, 0.35),
     xaxs = "i", yaxs = "i", xlab = "Parasite age (h)", ylab = "Density",
     main = paste0("Latent bimodal parasite-age distribution: ", example_3_ID, "\n",
                   "Component 1: mu = ", round(mean_3a, 2), " h; sigma = ",
                   round(sd_3a, 2), " h; weight = ", round(weight_3a, 2), "\n",
                   "Component 2: mu = ", round(mean_3b, 2), " h; sigma = ",
                   round(sd_3b, 2), " h; weight = ", round(weight_3b, 2)),
     xaxt = "n", yaxt = "n", cex.main = 0.9)

axis(1, at = x_ticks)
axis(2, at = y_ticks_bimodal, labels = sprintf("%.2f", y_ticks_bimodal), las = 3)
dev.off()

png(file.path(EXAMPLES_DIR, paste0(example_3_ID, "_observed_staging.png")),
    width = 2400, height = 1800, res = 300)

plot_observed_staging(observed_3, example_3_ID, example_3$Trial)
dev.off()