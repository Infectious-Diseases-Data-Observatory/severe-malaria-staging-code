import pandas as pd
from gosdt import GOSDTClassifier
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import train_test_split
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             confusion_matrix, precision_score, recall_score,
                             f1_score)

# Hyperparameters for the model
REGULARIZATION = 1e-3
SIMILAR_SUPPORT = False
DEPTH_BUDGET = 4
TIME_LIMIT = 60
VERBOSE = False
ALLOW_SMALL_REG = True
BALANCE = True

# Validation settings
TEST_SIZE = 0.20
RANDOM_STATES = range(123, 143)

N_IMPUTATIONS = 20

# Staging predictors excluded from the no-staging analysis
STAGING_FEATURES = ["mu_age:>5", "mu_age:>10", "mu_age:>15", "mu_age:>20",
                    "tring_obs:<30", "tring_obs:<40", "tring_obs:<55",
                    "tring_obs:<75"]

# Main function: Run repeated validation of the GOSDT algorithm
def gosdt_repeated_validation():

    # Read imputation 1 as the reference dataset
    # Patient IDs and outcomes are identical across all 20 imputations
    reference_file = "GOSDT/Python/Imputed_dataset_binned/Imputed_dataset_binned_01.csv"
    reference_data = pd.read_csv(reference_file, sep = ",")

    reference_ids = reference_data["ID"]
    reference_target = reference_data["outcome"]

    # Function 1: Extract dictionary from node
    def extract_dict(node):

        # If node is empty then return an empty dictionary
        if node is None:
            return {}

        # If node is already a dictionary then return it
        if isinstance(node, dict):
            return node

        # If node is an object with a dictionary attribute then return this dictionary
        if hasattr(node, "__dict__"):
            return vars(node)

        return {}

    # Function 2: Read one value from a node
    def read_node_value(node, *names, default = None):

        # Call Function 1 (extract dictionary from node)
        data = extract_dict(node)

        # Try to match input name with each dictionary name
        for name in names:

            # If a matching name is found then return its value
            if name in data:
                return data[name]

            # If input name exists as attribute of node then return its value
            if hasattr(node, name):
                return getattr(node, name)

        # If no matching value is found, return the default value
        return default

    # Function 3: Check whether this node is a leaf
    def is_leaf(node):

        # Call Function 2 (get prediction and children of node)
        has_prediction = read_node_value(node, "prediction") is not None
        has_left_child = read_node_value(node, "left_child") is not None
        has_right_child = read_node_value(node, "right_child") is not None

        # The node is a leaf if it has a prediction and no children
        return has_prediction or (not has_left_child and not has_right_child)

    # Function 4: Predict the outcome of one patient using the fitted tree
    def predict_one(node, row):

        # Start at the root of the tree
        current_node = node

        # Continue until a leaf is reached
        while not is_leaf(current_node):

            # Call Function 2 (get feature and children of node)
            feature = int(read_node_value(current_node, "feature"))
            left_child = read_node_value(current_node, "left_child")
            right_child = read_node_value(current_node, "right_child")

            # Get the patient's value for this binary feature
            feature_value = row.iloc[feature]

            # Follow the appropriate branch
            # This direction reproduces the loss reported internally by GOSDT
            if feature_value <= 0.5:
                current_node = right_child
            else:
                current_node = left_child

        # Return the prediction at the leaf
        return int(read_node_value(current_node, "prediction"))

    # Function 5: Predict the outcomes of all patients in a dataset
    def predict_dataset(node, features):

        # Call Function 4 for every patient
        predictions = [predict_one(node, row) for _, row in features.iterrows()]

        return predictions

    # Function 6: Calculate model performance
    def evaluate_model(node, features, target, dataset_name, imp, repeat, random_state):

        # Call Function 5 (predict outcomes)
        predictions = predict_dataset(node, features)

        # Calculate the confusion matrix
        tn, fp, fn, tp = confusion_matrix(target, predictions,
                                           labels = [0, 1]).ravel()

        # Return one row containing the validation results
        return pd.DataFrame({"Imputation": [imp], "Repeat": [repeat],
                             "RandomState": [random_state],
                             "Dataset": [dataset_name], "N": [len(target)],
                             "Accuracy": [accuracy_score(target, predictions)],
                             "BalancedAccuracy":
                                 [balanced_accuracy_score(target, predictions)],
                             "Sensitivity":
                                 [recall_score(target, predictions, pos_label = 1,
                                               zero_division = 0)],
                             "Specificity":
                                 [recall_score(target, predictions, pos_label = 0,
                                               zero_division = 0)],
                             "Precision":
                                 [precision_score(target, predictions, pos_label = 1,
                                                  zero_division = 0)],
                             "F1":
                                 [f1_score(target, predictions, pos_label = 1,
                                           zero_division = 0)],
                             "TN": [tn], "FP": [fp], "FN": [fn], "TP": [tp]})

    # Create an empty list for results from all repetitions and imputations
    all_results = []

    # Repeat the stratified 80:20 patient-level validation
    for repeat, random_state in enumerate(RANDOM_STATES, start = 1):
        print(f"\nRepetition {repeat} of {len(RANDOM_STATES)} "
              f"(random state = {random_state})\n")

        # Create ONE patient-level split for this repetition
        # The same patients are then used for training/test in all 20 imputations
        train_ids, test_ids = train_test_split(reference_ids, test_size = TEST_SIZE,
                                               random_state = random_state,
                                               stratify = reference_target)
        train_ids = set(train_ids)
        test_ids = set(test_ids)

        print(f"Training patients: {len(train_ids)}")
        print(f"Test patients: {len(test_ids)}")

        # Fit one independent tree for each imputation
        for imp in range(1, N_IMPUTATIONS + 1):
            print(f"Imputation {imp:02d}/{N_IMPUTATIONS}")

            input_file = f"GOSDT/Python/Imputed_dataset_binned/Imputed_dataset_binned_{imp:02d}.csv"
            data = pd.read_csv(input_file, sep = ",")

            # Apply this repetition's common patient split
            train_data = data[data["ID"].isin(train_ids)].copy()
            test_data = data[data["ID"].isin(test_ids)].copy()

            # ID is metadata and must never be supplied to GOSDT
            # Remove the eight staging predictors for the no-staging analysis
            X_train = train_data.drop(columns = ["ID", "outcome"] + STAGING_FEATURES)
            y_train = train_data["outcome"]
            X_test = test_data.drop(columns = ["ID", "outcome"] + STAGING_FEATURES)
            y_test = test_data["outcome"]

            if imp == 1:
                print(f"Candidate predictors: {X_train.shape[1]}")

            # Initialise the GOSDT model
            gosdt_model = GOSDTClassifier(regularization = REGULARIZATION,
                                          similar_support = SIMILAR_SUPPORT,
                                          time_limit = TIME_LIMIT,
                                          depth_budget = DEPTH_BUDGET,
                                          verbose = VERBOSE,
                                          allow_small_reg = ALLOW_SMALL_REG,
                                          balance = BALANCE)

            # Train on this imputation's training patients only
            gosdt_model.fit(X_train, y_train)

            # Extract the first fitted tree
            tree0 = gosdt_model.trees_[0].tree

            # Evaluate training dataset
            training_results = evaluate_model(tree0, X_train, y_train, "Training",
                                              imp, repeat, random_state)

            # Evaluate held-out test dataset
            test_results = evaluate_model(tree0, X_test, y_test, "Test", imp,
                                          repeat, random_state)

            all_results.append(training_results)
            all_results.append(test_results)

    # Combine the results from all repetitions
    validation_results = pd.concat(all_results, ignore_index = True)

    # Save all repeated validation results as a CSV
    output_csv = "GOSDT/Python/Imputed_dataset_repeated_validation_no_staging.csv"
    validation_results.to_csv(output_csv, index = False)
    print(f"\nSaved '{output_csv}'")

    # Keep held-out test results
    test_results = validation_results[validation_results["Dataset"] == "Test"].copy()

    # Validation statistics to summarise
    metrics = ["Accuracy", "BalancedAccuracy", "Sensitivity", "Specificity",
               "Precision", "F1"]

    # For each train/test repetition, first average performance
    # across the 20 imputations
    repeat_results = (test_results.groupby(["Repeat", "RandomState"],
                                           as_index = False)[metrics].mean())

    # Summarise variation across the 20 independent
    # train/test partitions
    summary_results = (repeat_results[metrics].agg(["mean", "std", "median",
                                                    "min", "max"]).T.reset_index())

    summary_results.columns = ["Metric", "Mean", "SD", "Median", "Minimum", "Maximum"]

    summary_results["Q25"] = [repeat_results[metric].quantile(0.25)
                              for metric in metrics]

    summary_results["Q75"] = [repeat_results[metric].quantile(0.75)
                              for metric in metrics]

    summary_results = summary_results[["Metric", "Mean", "SD", "Median", "Q25",
                                       "Q75", "Minimum", "Maximum"]]

    repeat_csv = "GOSDT/Python/Imputed_dataset_repeated_validation_by_repeat_no_staging.csv"
    repeat_results.to_csv(repeat_csv, index = False)
    print(f"Saved '{repeat_csv}'")

    summary_csv = "GOSDT/Python/Imputed_dataset_repeated_validation_summary_no_staging.csv"
    summary_results.to_csv(summary_csv, index = False)
    print(f"Saved '{summary_csv}'")

    print("\nRepeated held-out validation summary:")
    print(summary_results.to_string(index = False))

# Run the main function
gosdt_repeated_validation()