def predict_phenotype(
    annotations: str,
    phenotype: str, 
    Kleborate: str,
    model: str, 
    annotation_tool: str, 
    antibiotic_classes: str, 
    genes_to_class: str):
    """
    Predicts phenotype based on annotations using minimal models.

    Args:
        annotations (str): Path to the annotations CSV file.
        phenotype (str): Path to the phenotype CSV file.
        Kleborate (str): Path to Kleborate annotations (optional).
        model (str): 'XGBoost' or 'ElasticNet'.
        annotation_tool (str): Tool used for annotation.
        antibiotic_classes (str): File path for antibiotic classes.
        genes_to_class (str): File path for gene to class mappings.
    """
    # Load data
    AMR = pd.read_csv(phenotype)
    annotations = pd.read_csv(annotations, low_memory=False)
    antibiotic_classes = pd.read_csv(antibiotic_classes)
    if Kleborate == True: 
        Kleb = pd.read_csv(Kleborate, delimiter = "\t")
        Kleb = Kleb[Kleb["species"] == "Klebsiella pneumoniae"]
        Kleb["strain"] = Kleb["strain"].astype(str)
        AMR["Genome ID"] = AMR["Genome ID"].astype(str)
        AMR = AMR[AMR["Genome ID"].isin(Kleb["strain"])]
    else:
        AMR["Genome ID"] = AMR["Genome ID"].astype(str)
    antibiotics = list(np.unique(AMR["Antibiotic"]))
    antibiotics = list(np.unique(AMR["Antibiotic"]))
    genes_to_class = pd.read_csv(genes_to_class)

    for antibiotic in antibiotics:
        #format amr file so it matches the corresponding antibiotic
        AMR_subset = AMR[AMR["Antibiotic"]==antibiotic]
        AMR_subset.loc[:, "Genome ID"] = AMR_subset["Genome ID"].astype(str)
        AMR_subset = AMR_subset[["Genome ID", "Resistant Phenotype"]]
        AMR_subset["Resistant Phenotype"] = np.where(AMR_subset["Resistant Phenotype"] == 'Susceptible', 0, 1)

        # Extract relevant annotation genes based on annotation tool
        if annotation_tool == "RGI":
            df = RGI_dataframe(annotations, antibiotic, AMR_subset)

        elif annotation_tool == "AMRFinderPlus":
            df = AMRFinder_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        elif annotation_tool == "ResFinder":
            df = ResFinder_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        elif annotation_tool == "DeepARG":
            df = DeepARG_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        elif annotation_tool == "Abricate":
            df = Abricate_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        elif annotation_tool == "StarAMR":
            df = StarAMR_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        elif annotation_tool == "SraX":
            df = SraX_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        elif annotation_tool == "Kleborate":
            df = Kleborate_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset)

        
        AUC_performance, sensitivity_performances, specificity_performances = predict(df, model)
        summarise_results(AUC_performance, sensitivity_performances, specificity_performances, df, antibiotic)
