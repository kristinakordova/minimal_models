import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from imblearn.metrics import sensitivity_score, specificity_score
from xgboost import XGBClassifier
from Annotation_formatting import *

def main():
    parser = argparse.ArgumentParser(description='Predict phenotype based on annotations.')
    parser.add_argument('--annotations', required=True, help='Path to the annotations file')
    parser.add_argument('--phenotype', required=True, help='Phenotype table containing phenotypes in the same format as the BV-BRC database')
    parser.add_argument('--Kleborate', required=False, help='Path to the Kleborate annotations if applicable')
    parser.add_argument('--model', required=True, help='Which model to use - XGboost or ElasticNet')
    parser.add_argument('--annotation_tool', required=True, help='Which annotation tool are we using for prediction')
    parser.add_argument('--antibipotic_classes', required=True, help='Antibiotic classes to consider')
    parser.add_argument('--genes_to_class', required=True, help='Genes to class mapping file')

    args = parser.parse_args()

    predict_phenotype(
        annotations=args.annotations,
        phenotype=args.phenotype,
        Kleborate=args.Kleborate,
        model=args.model,
        annotation_tool=args.annotation_tool,
        antibipotic_classes=args.antibipotic_classes,
        genes_to_class=args.genes_to_class
    )

if __name__ == '__main__':
    main()
