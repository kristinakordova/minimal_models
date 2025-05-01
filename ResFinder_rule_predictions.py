import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from imblearn.metrics import sensitivity_score, specificity_score

# File paths
annotations_path = "ResFinder_annotations.csv"
phenotype_path = "Patric_AMR_data"
kleborate_path = "concatinated_Kleborate.txt"
antibiotic_classes_path = "Drugs_per_class.csv"
genes_to_class_path = "Kleborate_class_to_columns.csv"

# Load data
res_finer_automatic_predictions = pd.read_csv("phenotype_table.txt", sep="\t")
AMR = pd.read_csv(phenotype_path)
annotations = pd.read_csv(annotations_path, low_memory=False)
antibiotic_classes = pd.read_csv(antibiotic_classes_path)
genes_to_class = pd.read_csv(genes_to_class_path)

# Process Kleborate data if applicable
if kleborate_path:
    Kleb = pd.read_csv(kleborate_path, delimiter="\t")
    Kleb = Kleb[Kleb["species"] == "Klebsiella pneumoniae"]
    Kleb["strain"] = Kleb["strain"].astype(str)
    AMR["Genome ID"] = AMR["Genome ID"].astype(str)
    AMR = AMR[AMR["Genome ID"].isin(Kleb["strain"])]
else:
    AMR["Genome ID"] = AMR["Genome ID"].astype(str)

# Prepare antibiotics list
antibiotics = list(np.unique(AMR["Antibiotic"]))
antibiotics.remove("ampicillin")  # Remove specific antibiotic if needed

# Replace '+' with '/' in Antibiotic column
res_finer_automatic_predictions["Antibiotic"] = res_finer_automatic_predictions["Antibiotic"].str.replace("+", "/")

# Handle trimethoprim/sulfamethoxazole formatting
res_finer_automatic_predictions_trip = res_finer_automatic_predictions[
    res_finer_automatic_predictions["Antibiotic"] == "trimethoprim"
]
res_finer_automatic_predictions_sul = res_finer_automatic_predictions[
    res_finer_automatic_predictions["Antibiotic"] == "sulfamethoxazole"
]
res_finer_automatic_predictions_trip_sul = pd.merge(
    res_finer_automatic_predictions_trip,
    res_finer_automatic_predictions_sul,
    on="Genome.ID",
    how="inner"
)
res_finer_automatic_predictions_trip_sul["Phenotype_Res_Finder"] = np.where(
    res_finer_automatic_predictions_trip_sul["Phenotype_Res_Finder_x"] == res_finer_automatic_predictions_trip_sul["Phenotype_Res_Finder_y"],
    res_finer_automatic_predictions_trip_sul["Phenotype_Res_Finder_x"],
    "No resistance"
)
res_finer_automatic_predictions_trip_sul["Antibiotic"] = "trimethoprim/sulfamethoxazole"
res_finer_automatic_predictions_trip_sul = res_finer_automatic_predictions_trip_sul[
    ["Genome.ID", "Antibiotic", "Class_x", "Phenotype_Res_Finder"]
]
res_finer_automatic_predictions_trip_sul.columns = res_finer_automatic_predictions.columns

# Concatenate the new data
res_finer_automatic_predictions = pd.concat(
    [res_finer_automatic_predictions, res_finer_automatic_predictions_trip_sul],
    ignore_index=True
)

# Filter antibiotics to those present in both datasets
antibiotics = set(antibiotics).intersection(set(res_finer_automatic_predictions["Antibiotic"]))

# Evaluate performance for each antibiotic
results = []
for antibiotic in antibiotics:
    # Format AMR data for the current antibiotic
    AMR_subset = AMR[AMR["Antibiotic"] == antibiotic]
    AMR_subset.loc[:, "Genome ID"] = AMR_subset["Genome ID"].astype(str)
    AMR_subset = AMR_subset[["Genome ID", "Resistant Phenotype"]]
    AMR_subset.rename(columns={"Genome ID": "Genome.ID"}, inplace=True)
    # Format ResFinder predictions for the current antibiotic
    res_finer_automatic_predictions_subset = res_finer_automatic_predictions[
        res_finer_automatic_predictions["Antibiotic"] == antibiotic
    ]
    res_finer_automatic_predictions_subset.loc[:, "Genome.ID"] = res_finer_automatic_predictions_subset["Genome.ID"].astype(str)
    res_finer_automatic_predictions_subset = res_finer_automatic_predictions_subset[
        ["Genome.ID", "Phenotype_Res_Finder"]
    ]
    # Convert phenotypes to binary
    AMR_subset.loc[:, "Resistant Phenotype"] = np.where(AMR_subset["Resistant Phenotype"] == 'Susceptible', 0, 1)
    res_finer_automatic_predictions_subset.loc[:, "Phenotype_Res_Finder"] = np.where(
        res_finer_automatic_predictions_subset["Phenotype_Res_Finder"] == 'No resistance', 0, 1
    )
    AMR_subset["Genome.ID"] = AMR_subset["Genome.ID"].astype(str) + ".fna"
    # Merge datasets and calculate metrics
    merged = pd.merge(AMR_subset, res_finer_automatic_predictions_subset, on="Genome.ID", how="inner")
    merged["Resistant Phenotype"] = pd.to_numeric(merged["Resistant Phenotype"], errors="coerce")
    merged["Phenotype_Res_Finder"] = pd.to_numeric(merged["Phenotype_Res_Finder"], errors="coerce")
    AUC = roc_auc_score(merged["Resistant Phenotype"], merged["Phenotype_Res_Finder"])
    sensitivity = sensitivity_score(merged["Resistant Phenotype"], merged["Phenotype_Res_Finder"])
    specificity = specificity_score(merged["Resistant Phenotype"], merged["Phenotype_Res_Finder"])
    # Store results
    results.append((antibiotic, AUC, sensitivity, specificity))
    print(f"{antibiotic}, Res_finder_rule, {AUC}, {sensitivity}, {specificity}")

