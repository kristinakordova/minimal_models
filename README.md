# minimal_models
Minimal machine learning models for predicting antimicrobial resistance phenotype.

Use:

python predict_phenotype.py --annotations RGI_concat.csv --phenotype Patric_AMR_data --model "ElasticNet" --annotation_tool "RGI" --output_path Minimal_model --antibipotic_classes Drugs_per_class.csv --genes_to_class Kleborate_class_to_columns.csv

Data: https://zenodo.org/records/14126594?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjNiZTcxYzBmLTEwMWMtNGMxOS1hZWE1LTZlYjExNDI2NGZlMyIsImRhdGEiOnt9LCJyYW5kb20iOiJjZGFiZTVjYWJmOWQ2NDhkMTgwYmQ3ZjhhNGY2MjJhMyJ9.UbacJcicqSqhfTQ7gMir8uwCQpMJdl8UlKtBt9ifkVsHGaOf4zzorRXR6Bw2IVXvsLYhlg_7IKuas0Nl3V_SkA
