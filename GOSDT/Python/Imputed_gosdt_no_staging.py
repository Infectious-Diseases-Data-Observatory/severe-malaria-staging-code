import pandas as pd
from gosdt import GOSDTClassifier
import warnings
warnings.filterwarnings("ignore")
import re
import os
from sklearn.model_selection import train_test_split

# Hyperparameters for the model
REGULARIZATION = 1e-3
SIMILAR_SUPPORT = False
DEPTH_BUDGET = 4
TIME_LIMIT = 60
VERBOSE = True
ALLOW_SMALL_REG = True
BALANCE = True

N_IMPUTATIONS = 20

# Staging predictors excluded from the no-staging analysis
STAGING_FEATURES = ["mu_age:>5", "mu_age:>10", "mu_age:>15", "mu_age:>20",
                    "tring_obs:<30", "tring_obs:<40", "tring_obs:<55",
                    "tring_obs:<75"]

# Main function: Run the GOSDT algorithm
def gosdt_function(input_csv, output_csv, train_ids, test_ids):

    # Read the data from the CSV file
    data = pd.read_csv(input_csv, sep = ",")

    # Apply the same patient-level training/test split to every imputation
    train_data = data[data["ID"].isin(train_ids)].copy()
    test_data = data[data["ID"].isin(test_ids)].copy()

    # ID identifies patients but is never supplied to GOSDT.
    # Remove the eight staging predictors for this analysis.
    X_train = train_data.drop(columns = ["ID", "outcome"] + STAGING_FEATURES)
    y_train = train_data["outcome"]
    X_test = test_data.drop(columns = ["ID", "outcome"] + STAGING_FEATURES)
    y_test = test_data["outcome"]

    print(f"Training rows: {len(X_train)}")
    print(f"Test rows: {len(X_test)}")
    print(f"Candidate predictors: {X_train.shape[1]}")

    # Initialise the GOSDT model
    gosdt_model = GOSDTClassifier(regularization = REGULARIZATION,
                                  similar_support = SIMILAR_SUPPORT,
                                  time_limit = TIME_LIMIT,
                                  depth_budget = DEPTH_BUDGET, verbose = VERBOSE,
                                  allow_small_reg = ALLOW_SMALL_REG, balance = BALANCE)

    # Train the GOSDT model on the training dataset only
    gosdt_model.fit(X_train, y_train)

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

    # Function 4: Split a feature label into predictor name, comparison symbol,
    # and threshold value
    def extract_split_details(label):

        # Convert label to text and remove extra spaces
        label = str(label).strip()

        # Clean up label by removing the ending "<= 0.5"
        label = re.sub(r"\s*<=\s*0\.5\s*$", "", label)

        # Search for predictor name, comparison symbol, and threshold value
        m = re.search(r"(.+?)(?:[:\s]*)(<=|>=|<|>|==|=)\s*([-\d\.]+)$", label)
        if not m:
            return label.strip(), "", ""

        # Return the predictor name, comparison symbol, and threshold value separately
        predictor = m.group(1).strip()
        symbol = m.group(2).strip()
        value = m.group(3).strip()
        return predictor, symbol, value

    # Function 5: Get the left and right children of a node
    def get_children(node):

        # Call Function 2 (get children of node)
        return (read_node_value(node, "left_child", default = None),
                read_node_value(node, "right_child", default = None),)

    # Function 6: Walk through the tree and fill one table row per node/leaf with
    # information
    def tree_to_table_rows(node, path = ""):
        rows = []

        # If node is empty then return an empty table
        if node is None:
            return rows

        # Call Function 3 (if a leaf is reached)
        if is_leaf(node):

            # Call Function 2 (get its path and prediction)
            rows.append({"Path": path, "NodeType": "leaf", "Predictor": "",
                         "Symbol": "", "Value": "",
                         "Prediction": int(read_node_value(node, "prediction")),})
            return rows

        # Call Function 2 (get feature of node)
        feat = read_node_value(node, "feature")

        # Set colname to be feature name whether empty or not
        if feat is not None:
            colname = str(X_train.columns[int(feat)])
        else:
            colname = ""

        # Call Function 4 (split feature label)
        predictor, symbol, value = extract_split_details(colname)

        # Save the path, type, predictor name, comparison symbol, and threshold value
        rows.append({"Path": path, "NodeType": "node", "Predictor": predictor,
                     "Symbol": symbol, "Value": value, "Prediction": "",})

        # Call Function 5 (get left and right children of node)
        left_child, right_child = get_children(node)

        # Call Function 6 recursively (fill table rows with tree information)
        rows.extend(tree_to_table_rows(left_child, path + "L"))
        rows.extend(tree_to_table_rows(right_child, path + "R"))
        return rows

    # Extract the first tree from the fitted GOSDT model
    tree0 = gosdt_model.trees_[0].tree

    # Call Function 6 (fill table rows with information of first tree)
    rows = tree_to_table_rows(tree0)

    # Save the tree structure as a CSV
    data_out = pd.DataFrame(rows, columns = ["Path", "NodeType", "Predictor",
                                             "Symbol", "Value", "Prediction"])
    data_out.to_csv(output_csv, index = False)
    print(f"Saved '{output_csv}'")

# Create the same common training/test patient split as Imputed_gosdt.py
reference_file = "GOSDT/Python/Imputed_dataset_binned/Imputed_dataset_binned_01.csv"
reference_data = pd.read_csv(reference_file)

train_ids, test_ids = train_test_split(reference_data["ID"],
                                       test_size = 0.20,
                                       random_state = 123,
                                       stratify = reference_data["outcome"])

train_ids = set(train_ids)
test_ids = set(test_ids)

print(f"Training patients: {len(train_ids)}")
print(f"Test patients: {len(test_ids)}")

# Create output directory
directory = "GOSDT/Python/Imputed_dataset_splits"
os.makedirs(directory, exist_ok = True)

# Fit one independent no-staging GOSDT tree per imputation
for imp in range(1, N_IMPUTATIONS + 1):
    print(f"\nImputation {imp:02d}/{N_IMPUTATIONS}\n")

    input_csv = f"GOSDT/Python/Imputed_dataset_binned/Imputed_dataset_binned_{imp:02d}.csv"
    output_csv = f"GOSDT/Python/Imputed_dataset_splits/Imputed_dataset_splits_no_staging_{imp:02d}.csv"

    gosdt_function(input_csv, output_csv, train_ids, test_ids)