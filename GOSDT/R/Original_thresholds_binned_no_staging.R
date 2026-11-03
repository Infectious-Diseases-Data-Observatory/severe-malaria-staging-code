# CREATE A TABLE CONTAINING A BINARY COLUMN FOR EACH SPLIT VALUE

# Load data.table package to read large CSV file
library(data.table)

# Turn off scientific notation
options(scipen = 999)

DATA_DIR <- file.path("Staging/Complete case/outputs/CSV")
my_table <- fread(file.path(DATA_DIR, "dat_all_hierarchical.csv"))
OUT_DIR <- file.path("GOSDT/Python")

# Create a binary indicator for Trial
my_table[, Trial := as.integer(Trial == "SEAQUAMAT")]

# Provide threshold names, symbols, and values
reps_of_each <- c(1, 1, 3, 1, 1, 3, 1, 1, 4, 4, 3)
threshold_names <- rep(c("coma", "aconvulsions", "hb", "aanaemia", "ashockcomp", "temp",
                         "ahypoglyc", "Trial", "be", "bun", "parasitaemia"),
                       times = reps_of_each)
threshold_symbols <- rep(c("=", "=", "<", "=", "=", "<", "=", "=", "<", ">", ">"),
                         times = reps_of_each)
threshold_values <- c(1, 1, 5, 8, 11, 1, 1, 37, 39, 40, 1, 1, -11, -9, -6, -3, 8, 10, 15, 20,
                      100000, 250000, 500000)
num_thresholds <- length(threshold_names)

# Write full names of threshold columns
threshold_full <- 1:num_thresholds

for (i in 1:num_thresholds) {
  threshold_full[i] <- sprintf("%s:%s%s", threshold_names[i], threshold_symbols[i],
                               threshold_values[i])
}

# Create new data frame to store columns of binary data
Original_dataset_binned <- data.frame(matrix(NA, nrow = nrow(my_table),
                                             ncol = num_thresholds))

names(Original_dataset_binned) <- threshold_full

# Create function which increments itself
inc <- function(x) {
  eval.parent(substitute(x <- x + 1))
}

# Fill data frame with columns of 1s & 0s to describe presence/absence of each feature
for (i in 1:num_thresholds) {
  my_condition <- 1:nrow(my_table)
  j <- 0

  # Loop through elements of a feature of interest
  # Want to split them in two (0s & 1s) based on a condition
  for (x in my_table[[threshold_names[i]]]) {
    inc(j)

    # Handle missing values
    if (is.na(x)) {
      my_condition[j] = NA

      # Separate conditions according to symbol
    } else if (threshold_symbols[i] == "<") {

      # If condition true then set current element of vector to be 1
      if (x < threshold_values[i]) {
        my_condition[j] = 1

        # Else set current element of vector to be 0
      } else {
        my_condition[j] = 0
      }
    } else if (threshold_symbols[i] == ">") {
      if (x > threshold_values[i]) {
        my_condition[j] = 1
      } else {
        my_condition[j] = 0
      }
    } else {
      if (x == threshold_values[i]) {
        my_condition[j] = 1
      } else {
        my_condition[j] = 0
      }
    }
  }

  # Set the vector to be a column of new data frame
  Original_dataset_binned[[threshold_full[i]]] <- my_condition
  Original_dataset_binned
}

# Append mortality outcome as the final column
Original_dataset_binned$outcome <- my_table$outcome

# State location where data frame will be saved
my_output <- file.path(OUT_DIR, "Original_dataset_binned_no_staging.csv")

# Save data frame as CSV without row names
write.csv(Original_dataset_binned, file = my_output, row.names = FALSE)