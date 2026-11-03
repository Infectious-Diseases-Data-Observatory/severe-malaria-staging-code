# MULTIPLE IMPUTATION - MAR
# 20 imputed datasets
# 20 iterations per dataset
#
# MICE (Multivariate Imputation by Chained Equations) methods:
# - pmm = predictive mean matching, used for continuous/numeric variables
# - logreg = logistic regression, used for binary variables
# - polyreg = polytomous logistic regression, used for unordered categorical variables
#             with >2 levels
# - polr = proportional odds logistic regression, used for ordered categorical variables

# Load data.table package to read large CSV file
library(data.table)

# Load MICE package to implement the above methods
library(mice)

# Turn off scientific notation
options(scipen = 999)

INFILE <- file.path("Data/dat_all.csv")

# Choose whether to remove patients with RBC (red blood cell) total < 93 or > 107
apply_rbc_restriction <- TRUE

# Choose whether to keep hrp2_plasma in the multiple imputation or not
keep_hrp2_plasma <- FALSE

# Read large CSV file
my_table <- fread(INFILE, na.strings = c("", "NA", "NaN"))

# Remove:
# - import/export artefacts (Unnamed: 0, V1)
# - identifier variables (ID)
# - redundant variables (country)
# - variables not used in the analysis (gam)
cols_to_drop <- intersect(names(my_table), c("Unnamed: 0", "V1", "country", "gam"))

if (length(cols_to_drop) > 0) {
  my_table[, (cols_to_drop) := NULL]
}

# Remove patients with RBC-total < 93 or > 107 if requested
if (apply_rbc_restriction) {
  rbc_vars <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont")

  if (!all(rbc_vars %in% names(my_table))) {
    stop("Not all RBC stage variables were found in the dataset.")
  }

  my_table[, rbc_total := tring + sring + lring + etroph + mtroph + ltroph + schizont]
  my_table <- my_table[between(rbc_total, 93, 107)]
  my_table[, rbc_total := NULL]
}

# Gather binary variables that exist in the dataset
binary_vars <- intersect(c("outcome", "coma", "aconvulsions", "aanaemia", "ashockcomp",
                           "aacidosis", "ahypoglyc"), names(my_table))

# Remove unused categories from each binary variable
for (v in binary_vars) {
  my_table[[v]] <- droplevels(factor(my_table[[v]]))
}

# Gather categorical variables that exist in the dataset
cat_vars <- intersect(c("site", "drug", "Trial"), names(my_table))

# Remove unused categories from each categorical variable
for (v in cat_vars) {
  my_table[[v]] <- droplevels(factor(my_table[[v]]))
}

# Gather numeric variables that exist in the dataset
numeric_vars <- intersect(c("age", "hb", "hct", "temp", "be", "bun", "parasitaemia",
                            "hrp2_plasma", "tring", "sring", "lring", "etroph", "mtroph",
                            "ltroph", "schizont"), names(my_table))

# Remove unused categories from each numeric variable
for (v in numeric_vars) {
  my_table[[v]] <- as.numeric(my_table[[v]])
}

# Function 1: Create a table displaying a summary of missingness in the dataset
create_missingness_summary <- function(dt, out_file) {

  miss_sum <- data.table(variable = names(dt), n_missing = as.integer(colSums(is.na(dt))),
                         pct_missing = round(100 * colSums(is.na(dt)) /
                                               nrow(dt), 1))[order(-n_missing)]
  fwrite(miss_sum, out_file)
  invisible(miss_sum)
}

# Function 2: Create a table displaying missingness by site
create_missingness_by_site <- function(dt, vars, out_file) {

  # Keep only variables that are present in the dataset
  vars <- intersect(vars, names(dt))

  # If none remain, or if site is unavailable, skip the site-specific missingness summary
  if (length(vars) == 0 || !("site" %in% names(dt))) return(invisible(NULL))

  # Calculate the percentage of missing values for each variable within each site
  miss_site <- dt[, lapply(.SD, function(x) round(mean(is.na(x)) * 100, 1)),
                  by = site, .SDcols = vars]

  fwrite(miss_site, out_file)
  invisible(miss_site)
}

# Function 3: Choose a reasonable MICE method
set_mice_method <- function(x) {
  if (is.ordered(x)) return("polr")
  if (is.factor(x)) {
    if (nlevels(x) == 2) return("logreg")
    return("polyreg")
  }
  return("pmm")
}

# Do not impute:
# - outcome (but remains a predictor)
# - site (but remains a predictor)
# - drug (remove from the predictor list)
# - hrp2_plasma (it is structurally missing in some sites)
dont_impute <- intersect(c("outcome", "site", "drug", "hrp2_plasma"), names(my_table))

# Do impute this list of variables
imputed_vars <- intersect(c("hb", "hct", "temp", "be", "bun", "parasitaemia"),
                          names(my_table))

# Impute AQUAMAT and SEAQUAMAT separately
trial_levels <- intersect(c("AQUAMAT", "SEAQUAMAT"), unique(as.character(my_table$Trial)))

if (length(trial_levels) == 0) {
  stop("No AQUAMAT / SEAQUAMAT rows found in Trial column.")
}

# Create a folder for all plots
plot_root <- file.path("Data/Imputation_plots")
dir.create(plot_root, showWarnings = FALSE, recursive = TRUE)

all_imputed_long <- list()
imp_objects <- list()

# Create labels for the trace plots
trace_name_map <- c(hb = "Hb [g/dL]", hct = "HCT [L/L]", bun = "BUN [mmol/L]",
                    be = "Base Excess [mmol/L]", temp = "Temperature [C]",
                    parasitaemia = "Parasitaemia [/uL]")

# Create labels for the density plots
density_name_map <- c(hb = "Hb [g/dL]", hct = "HCT [L/L]", bun = "BUN [mmol/L]",
                      be = "Base Excess [mmol/L]", temp = "Temperature [C]",
                      parasitaemia = "Parasitaemia [/uL]")

# Create nicer tick marks for the density plots
density_ticks <- list(hb = c(0, 5, 10, 15, 20), hct = c(0, 10, 20, 30, 40, 50, 60, 70),
                      bun = c(0, 50, 100, 150, 200, 250, 300),
                      be = c(-40, -30, -20, -10, 0, 10, 20), temp = c(34, 36, 38, 40, 42),
                      parasitaemia = c(0, 500000, 1000000, 1500000, 2000000))

# Function 4: Create custom trace plots for a MIDS (Multiply Imputed Data Set)
plot_mids_custom <- function(x, y = NULL, theme = mice.theme(), layout = c(2, 3), type = "l",
                             col = 1:10, lty = 1, var_labels = NULL, ylab = "Value") {

  # Function 5: Create custom strip labels for the trace plots
  strip.combined <- function(which.given, which.panel, factor.levels, ...) {
    lab <- as.character(factor.levels[which.panel[which.given]])

    # Relabel the variable-name strip on the right
    if (!is.null(var_labels) && which.given == 2 && lab %in% names(var_labels)) {
      lab <- unname(var_labels[lab])
    }

    # Draw the left strip (mean / sd)
    if (which.given == 1) {
      lattice::panel.rect(0, 0, 1, 1, col = theme$strip.background$col, border = 1)
      lattice::panel.text(x = 0, y = 0.5, pos = 4, lab = lab)
    }

    # Draw the right strip (variable name)
    if (which.given == 2) {
      lattice::panel.text(x = 1, y = 0.5, pos = 2, lab = lab)
    }
  }

  # Check that the supplied object is a MIDS object
  if (!is.mids(x)) {
    stop("argument 'x' must be a 'mids' object", call. = FALSE)
  }

  # Check that convergence diagnostics are available
  if (is.null(x$chainMean)) {
    stop("no convergence diagnostics found", call. = FALSE)
  }

  # Extract chain means and standard deviations
  mn <- x$chainMean
  sm <- sqrt(x$chainVar)

  # Select subset of non-missing entries
  obs <- apply(!(is.nan(mn) | is.na(mn)), 1, all)
  varlist <- names(obs)[obs]

  # Create formula if not given in y
  if (missing(y)) {
    formula <- as.formula(paste0(paste0(varlist, collapse = "+"), "~.it|.ms"))
  } else {
    formula <- NULL

    # Plot all variables if y is NULL
    if (is.null(y)) {
      formula <- as.formula(paste0(paste0(varlist, collapse = "+"), "~.it|.ms"))
    }

    # Plot variables supplied by name
    if (is.character(y)) {
      formula <- if (length(y) == 1) {
        as.formula(paste0(y, "~.it|.ms"))
      } else {
        as.formula(paste0(paste0(y, collapse = "+"), "~.it|.ms"))
      }
    }

    # Plot variables supplied by position
    if (is.integer(y) || is.logical(y)) {
      vars <- varlist[y]
      formula <- if (length(vars) == 1) {
        as.formula(paste0(vars, "~.it|.ms"))
      } else {
        as.formula(paste0(paste0(vars, collapse = "+"), "~.it|.ms"))
      }
    }

    # Accept a user-supplied plotting formula
    if (is.null(formula)) {
      formula <- as.formula(y)
    }
  }

  # Extract the number of imputations and iterations
  m <- x$m
  it <- x$iteration

  # Rearrange chain summaries into plotting format
  mn <- matrix(aperm(mn[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it)
  sm <- matrix(aperm(sm[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it)

  # Create plotting dataset for all iterations and imputations
  adm <- expand.grid(seq_len(it), seq_len(m), c("mean", "sd"))
  data <- cbind(adm, rbind(mn, sm))
  colnames(data) <- c(".it", ".m", ".ms", varlist)

  # Remove local object so lattice uses the data column
  .m <- NULL
  rm(.m)

  # Function 5: Keep readable y-axis ticks, but keep them between 2 and 4 per panel
  yscale_components_min2 <- function(lim, packet.number = 0, packet.list = NULL,
                                     right = TRUE, ...) {

    # Create default y-axis components
    ans <- lattice::yscale.components.default(lim, packet.number = packet.number,
                                              packet.list = packet.list, right = right, ...)

    # Determine the plotting range
    lim <- range(lim, finite = TRUE)

    # Return the default axis if the plotting range is invalid
    if (!all(is.finite(lim))) {
      return(ans)
    }

    # Create a small range if all values are identical
    if (diff(lim) == 0) {
      at <- lim[1] + c(-0.5, 0, 0.5)
    } else {
      at <- pretty(lim, n = 4)
      at <- at[at >= lim[1] & at <= lim[2]]

      # Reduce the number of ticks if there are too many
      if (length(at) > 4) {
        at <- pretty(lim, n = 3)
        at <- at[at >= lim[1] & at <= lim[2]]
      }

      # Further reduce the number of tick marks if required
      if (length(at) > 4) {
        at <- pretty(lim, n = 2)
        at <- at[at >= lim[1] & at <= lim[2]]
      }

      # Create additional tick marks if there are too few
      if (length(at) < 2) {
        at <- seq(lim[1], lim[2], length.out = 3)
      }

      # Limit the final axis to at most four tick marks
      if (length(at) > 4) {
        at <- at[round(seq(1, length(at), length.out = 4))]
        at <- unique(at)
      }
    }

    # Update the y-axis tick positions and labels
    ans$left$ticks$at <- at
    ans$left$labels$at <- at
    ans$left$labels$labels <- format(at, scientific = FALSE, trim = TRUE)

    ans
  }

  # Create trace plot
  tp <- lattice::xyplot(x = formula, data = data, groups = .m, type = type, lty = lty,
                        col = col, layout = layout, between = list(x = 0.9, y = 0.4),
                        scales = list(y = list(relation = "free", alternating = FALSE),
                                      x = list(alternating = FALSE)),
                        as.table = TRUE, xlab = "Iteration", ylab = ylab,

                        # Call Function 5 (create custom strip labels)
                        strip = strip.combined, par.strip.text = list(lines = 0.5),
                        yscale.components = yscale_components_min2)

  # Apply the MICE plotting theme
  update(tp, par.settings = theme)
}

# Run the multiple imputation separately for each trial
for (trial_name in trial_levels) {
  message("Imputing trial: ", trial_name)

  dt <- as.data.table(copy(my_table[Trial == trial_name]))

  # Preserve the original patient IDs for this trial
  trial_ID <- dt$ID

  # Remove identifiers before multiple imputation
  dt[, c("ID", "Trial") := NULL]

  # Gather binary and categorical variables
  for (v in c("outcome", "coma", "aconvulsions", "aanaemia", "ashockcomp", "aacidosis",
              "ahypoglyc", "site", "drug")) {

    # Remove unused categories from each variable
    if (v %in% names(dt)) {
      dt[[v]] <- droplevels(factor(dt[[v]]))
    }
  }

  # Remove hrp2_plasma from the multiple imputation if requested
  if (!keep_hrp2_plasma && "hrp2_plasma" %in% names(dt)) {
    dt[, hrp2_plasma := NULL]
  }

  # Create a sub-folder for each trial
  trial_dir <- file.path(plot_root, trial_name)
  dir.create(trial_dir, showWarnings = FALSE, recursive = TRUE)

  # Call Function 1 (create summary of missingness)
  create_missingness_summary(dt, file.path(paste0("Data/Missingness_", trial_name, ".csv")))

  # Call Function 2 (create summary of missingness by site)
  create_missingness_by_site(dt, c("hb", "hct", "temp", "be", "bun", "parasitaemia"),
                             file.path(paste0("Data/Missingness_site_", trial_name, ".csv")))

  # Gather variables with at least one missing value
  missing_vars <- names(which(colSums(is.na(dt)) > 0))

  # Do not impute outcome / site / drug / hrp2_plasma
  impute_vars <- setdiff(missing_vars, dont_impute)

  # Build a vector of MICE methods
  method_vec <- make.method(dt)
  method_vec[] <- ""

  # Call Function 3 (choose reasonable MICE method)
  for (v in impute_vars) {
    method_vec[v] <- set_mice_method(dt[[v]])
  }

  # Set prediction matrix to be as broad as possible given prior restrictions
  pred_matrix <- make.predictorMatrix(dt)
  pred_matrix[,] <- 1

  # Other than self-prediction
  diag(pred_matrix) <- 0

  # Do not use drug as a predictor
  if ("drug" %in% colnames(pred_matrix)) {
    pred_matrix[, "drug"] <- 0
  }

  # Fetch the variables not to be imputed
  target_not_impute <- intersect(dont_impute, names(dt))

  # Add this information to the predictor matrix
  if (length(target_not_impute) > 0) {
    pred_matrix[target_not_impute, ] <- 0
  }

  # Outcome is not imputed but used as a predictor
  if ("outcome" %in% names(dt)) {
    pred_matrix["outcome", ] <- 0
  }

  message("Starting MICE for ", trial_name)

  # Run the MICE algorithm with the internal bookkeeping count suppressed
  imp <- suppressWarnings(suppressMessages(mice(dt, m = 20, maxit = 20, method = method_vec,
                                                predictorMatrix = pred_matrix, donors = 10,
                                                seed = 123, printFlag = FALSE)))

  message("Finished MICE for ", trial_name)

  # Create trace plots for variables that were imputed
  imputed_vars_this_trial <- impute_vars
  trace_sets <- list(intersect(c("hb", "hct", "bun"), imputed_vars_this_trial),
                     intersect(c("be", "temp", "parasitaemia"), imputed_vars_this_trial))

  for (i in seq_along(trace_sets)) {
    vars <- trace_sets[[i]]

    if (length(vars) > 0) {

      # Save each trace plot
      jpeg(file.path(trial_dir, paste0("trace_plot_", i, ".jpg")),
           width = 1800, height = 1400, res = 200)

      # Call Function 4 (create each custom trace plot)
      print(plot_mids_custom(imp, vars, var_labels = trace_name_map))
      dev.off()
    }
  }

  # Create density plots for variables that were imputed
  for (v in imputed_vars) {
    if (!(v %in% imputed_vars_this_trial)) next

    # Use the observed values to decide whether a density plot is possible
    x_obs <- dt[[v]][!is.na(dt[[v]])]

    if (is.numeric(dt[[v]]) && length(x_obs) >= 2 && length(unique(x_obs)) >= 2) {
      jpeg(file.path(trial_dir, paste0("density_", v, ".jpg")),
           width = 1800, height = 1400, res = 200)

      # Plot observed density (blue line) vs imputed densities (red lines)
      # Skip variables that cannot support a density estimate
      tryCatch({print(densityplot(imp, as.formula(paste("~", v)),
                                  scales = list(x = list(cex = 1.6, at = density_ticks[[v]],
                                                         labels = format(density_ticks[[v]],
                                                                         big.mark = ",",
                                                                         scientific = FALSE,
                                                                         trim = TRUE)),
                                                y = list(cex = 1.6)),
                                  xlab = list(label = density_name_map[[v]], cex = 2.0),
                                  ylab = list(label = "Density", cex = 2.0)))},
               error = function(e) {message("Skipping density plot for ", v, ": ",
                                            e$message)})

      dev.off()
    }
  }

  # Create long-form imputed datasets
  imputed_long <- complete(imp, action = "long", include = TRUE)
  imputed_long <- as.data.table(imputed_long)
  imputed_long <- imputed_long[imputed_long$.imp != 0, ]

  # Restore the original patient ID
  imputed_long[, ID := trial_ID[as.integer(.id)]]

  imputed_long[, trial := trial_name]

  all_imputed_long[[trial_name]] <- imputed_long
  message("Stored long dataset for ", trial_name, " with ", nrow(imputed_long), " rows")
}

# Combine and save all imputed datasets
message("Combining imputed datasets")
combined_imputed_long <- rbindlist(all_imputed_long, fill = TRUE)
message("Combined dataset dimensions:")
print(dim(combined_imputed_long))

# Check patient identifiers
stopifnot(!anyNA(combined_imputed_long$ID), !anyNA(combined_imputed_long$.id),
          !anyNA(combined_imputed_long$.imp))

# Check that every imputation contains 4326 unique patients
patients_per_imputation <- combined_imputed_long[, uniqueN(ID), by = .imp]

stopifnot(nrow(patients_per_imputation) == 20, all(patients_per_imputation$V1 == 4326))

# Check that every patient occurs exactly once in each imputation
patient_imputation_counts <- combined_imputed_long[, .N, by = .(ID, .imp)]

stopifnot(all(patient_imputation_counts$N == 1))

# Check that the same 4326 patients occur in all 20 imputations
stopifnot(uniqueN(combined_imputed_long$ID) == 4326)

message("Patient ID checks passed")

fwrite(combined_imputed_long, file.path("Data/Imputed_dataset_long.csv"))
message("Finished writing Imputed_dataset_long.csv")
message("All imputations and exports completed successfully.")