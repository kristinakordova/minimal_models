import argparse

def predict_phenotype(annotations, phenotype, Kleborate, model, annotation_tool, output_path, antibiotic_classes, genes_to_class):
    """
    Preducts phenotype based on annotations - minimal models
    
    Parameters:
    annotations: Path to the annotations file
    phenotype: Phenotype table containing phenotypes in the same format as the BV-BRC database
    Kleborate: Path to the Kleborate annotations to subset species in analysing Klebsiella pneumoniae
    model: Which model to use - XGboost or ElasticNet
    annotation_tool: Which annotation tool are we using for prediction - See supported tools in the README
    output_path: Folder to save the output
    antibipotic_classes: Antibiotic to class or subclass mapping file
    genes_to_class: Genes to class mapping file for Kleborate
    """
    # Replace this with actual prediction code
    print("Preducts phenotype based on annotations - minimal models") 

    import pandas as pd
    from tqdm import tqdm
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    import matplotlib.pyplot as plt
    from sklearn.linear_model import LogisticRegression
    from imblearn.metrics import sensitivity_score, specificity_score
    from sklearn.metrics import roc_auc_score, recall_score, confusion_matrix
    from sklearn.model_selection import train_test_split
    from sklearn.model_selection import KFold, cross_val_score
    from xgboost import XGBClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier  # Assuming you're using RandomForestClassifier
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.model_selection import StratifiedKFold
    import pandas as pd
    from sklearn.model_selection import StratifiedKFold
    import matplotlib.pyplot as plt
    from imblearn.metrics import sensitivity_score, specificity_score
    from sklearn.metrics import roc_auc_score, recall_score, confusion_matrix

    tool = f"{annotation_tool}_{model}"
    AMR = pd.read_csv(phenotype)
    annotations = pd.read_csv(annotations, low_memory=False)
    file = open(f'{output_path}/AUC_across_tools.txt', 'a')
    file_2 = open(f'{output_path}/number_of_genes_across_tools.txt', 'a')
    file_3 = open(f'{output_path}/performances_across_tools_sensitivity.txt', 'a')
    file_4 = open(f'{output_path}/performances_across_tools_specificity.txt', 'a')
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
    average_across_antibiotics = []
    average_across_antibiotics_sensitivity = []
    average_across_antibiotics_specificity = []
    antibiotics = list(np.unique(AMR["Antibiotic"]))
    number_of_genes = []
    genes_to_class = pd.read_csv(genes_to_class)

    for antibiotic in antibiotics:
        #format amr file so it matches the corresponding antibiotic
        AMR_subset = AMR[AMR["Antibiotic"]==antibiotic]
        AMR_subset.loc[:, "Genome ID"] = AMR_subset["Genome ID"].astype(str)
        AMR_subset = AMR_subset[["Genome ID", "Resistant Phenotype"]]
        AMR_subset["Resistant Phenotype"] = np.where(AMR_subset["Resistant Phenotype"] == 'Susceptible', 0, 1)
        #extract class 
        if annotation_tool == "RGI":
            annotations["Genome.ID"] = annotations["Genome.ID"].str.replace('.fna$', '', regex=True)
            #extract class 
            if '/' in antibiotic:
                antibiotic = antibiotic.split('/')
                genes1 = list(set(annotations[annotations["Antibiotic"].str.contains(antibiotic[0], case=False, na=False)]["Best_Hit_ARO"]))
                genes2 = list(set(annotations[annotations["Antibiotic"].str.contains(antibiotic[1], case=False, na=False)]["Best_Hit_ARO"]))
                genes = genes1 + genes2
            else:
                genes = list(set(annotations[annotations["Antibiotic"].str.contains(antibiotic, case=False, na=False)]["Best_Hit_ARO"]))
            RGI_subset = annotations[["Best_Hit_ARO","Genome.ID","Best_Identities"]]
            RGI_subset = RGI_subset.drop_duplicates(subset=["Genome.ID", "Best_Hit_ARO"])
            RGI_subset = RGI_subset.dropna()
            df_pivot = RGI_subset.pivot(index="Genome.ID", columns="Best_Hit_ARO", values="Best_Identities")
            df_pivot = df_pivot.notna().astype(int)
            df_pivot = df_pivot[genes]
            merged = pd.merge(df_pivot, AMR_subset, left_on= df_pivot.index, right_on = "Genome ID", how = "inner")
        
        elif annotation_tool == "AMRFinderPlus":
            antibiotic_classes_subset = antibipotic_classes[antibipotic_classes["Drugs"]==antibiotic]
            antibiotic_classes_subset = list(set(antibiotic_classes_subset["Subclass"]))
            annotations.loc[:, "Contig id"] = annotations["Contig id"].astype(str)
            amr_finder_drug = annotations[annotations["Contig id"].isin(AMR_subset["Genome ID"])]
            genes = list(set(amr_finder_drug[amr_finder_drug["Subclass"].isin(antibiotic_classes_subset)]["Gene symbol"]))
            amr_finder_drug = amr_finder_drug[["Contig id","Gene symbol","% Identity to reference sequence"]]
            amr_finder_drug_change = amr_finder_drug.drop_duplicates(subset=["Contig id", "Gene symbol"])
            df_pivot = amr_finder_drug_change.pivot(index="Contig id", columns="Gene symbol", values="% Identity to reference sequence")
            df_pivot = df_pivot.notna().astype(int)
            df_pivot = df_pivot[genes]
            merged = pd.merge(df_pivot, AMR_subset, left_on= df_pivot.index, right_on = "Genome ID", how = "inner")
            print(merged)
        
        elif annotation_tool == annotation_tool == "ResFinder":
            antibiotic_classes_subset = antibipotic_classes[antibipotic_classes["Drugs"] == antibiotic]
            antibiotic_classes_subset = list(set(antibiotic_classes_subset["Subclass"]))
            antibiotic_classes_subset = [x.lower() for x in antibiotic_classes_subset]
            annotations["Genome.ID"] = annotations["Genome.ID"].str.replace('.fna$', '', regex=True)
            annotations.loc[:, "Genome.ID"] = annotations["Genome.ID"].astype(str)
            ResFinder_drug = annotations[annotations["Genome.ID"].isin(AMR_subset["Genome ID"])]
            if '/' in antibiotic:
                antibiotic = antibiotic.split('/')
                genes1 = list(set(ResFinder_drug[ResFinder_drug["Phenotype"].str.contains(antibiotic[0], case=False, na=False)]["Resistance gene"]))
                genes2 = list(set(ResFinder_drug[ResFinder_drug["Phenotype"].str.contains(antibiotic[1], case=False, na=False)]["Resistance gene"]))
                genes = genes1 + genes2
            else:
                genes = list(set(ResFinder_drug[ResFinder_drug["Phenotype"].str.contains(antibiotic, case=False, na=False)]["Resistance gene"]))
            ResFinder_drug = ResFinder_drug[["Genome.ID","Resistance gene","Identity"]]
            ResFinder_drug = ResFinder_drug.drop_duplicates(subset=["Genome.ID", "Resistance gene"])
            df_pivot_res_finder = ResFinder_drug.pivot(index="Genome.ID", columns="Resistance gene", values="Identity")
            df_pivot_res_finder = df_pivot_res_finder.notna().astype(int)
            df_pivot_res_finder = df_pivot_res_finder[genes]
            merged = pd.merge(df_pivot_res_finder, AMR_subset, left_on= df_pivot_res_finder.index, right_on = "Genome ID", how = "right")
            print(merged)

        elif annotation_tool == "DeepARG":
            antibiotic_classes_subset = antibipotic_classes[antibipotic_classes["Drugs"] == antibiotic]
            antibiotic_classes_subset = list(set(antibiotic_classes_subset["Subclass"]))
            antibiotic_classes_subset = [x.lower() for x in antibiotic_classes_subset]
            annotations["Genome.ID"] = annotations["Genome.ID"].str.replace('.fna$', '', regex=True)
            annotations.loc[:, "Genome.ID"] = annotations["Genome.ID"].astype(str)
            DeepARG_drug = annotations[annotations["Genome.ID"].isin(AMR_subset["Genome ID"])]
            genes = list(set(DeepARG_drug[DeepARG_drug["predicted_ARG-class"].isin(antibiotic_classes_subset)]["best-hit"]))
            genes_multiclass = list(DeepARG_drug[DeepARG_drug["predicted_ARG-class"]=="multidrug"]["best-hit"])
            genes = genes + genes_multiclass
            genes = list(set([x.split('|')[-1] for x in genes]))
            DeepARG_drug = DeepARG_drug[["Genome.ID","best-hit","identity"]]
            DeepARG_drug["best-hit"] = DeepARG_drug["best-hit"].str.split('|').str[-1]
            DeepARG_drug = DeepARG_drug.drop_duplicates(subset=["Genome.ID", "best-hit"])
            df_pivot = DeepARG_drug.pivot(index="Genome.ID", columns="best-hit", values="identity")
            df_pivot = df_pivot.notna().astype(int)
            df_pivot = df_pivot[genes]
            merged = pd.merge(df_pivot, AMR_subset, left_on= df_pivot.index, right_on = "Genome ID", how = "inner")

        elif annotation_tool == "Abricate":
            antibiotic_classes_subset = antibiotic_classes[antibiotic_classes["Drugs"] == antibiotic]
            antibiotic_classes_subset["Subclass"] = antibiotic_classes_subset["Subclass"].str.lower()
            annotations["#FILE"] = annotations["#FILE"].str.rstrip(".fna")
            annotations["RESISTANCE"]= annotations["RESISTANCE"].str.lower()
            annotations["#FILE"] = annotations["#FILE"].astype(str)
            Abricate_subset = annotations[annotations["#FILE"].isin(AMR_subset["Genome ID"])]
            mask = Abricate_subset["RESISTANCE"].apply(lambda x: any(subclass in x for subclass in antibiotic_classes_subset["Subclass"]))
            genes = Abricate_subset[mask]["GENE"]
            Abricate_subset = Abricate_subset[["#FILE","%COVERAGE", "GENE"]]
            Abricate_subset = Abricate_subset.drop_duplicates(subset=["#FILE", "GENE"])
            Abricate_subset = Abricate_subset.pivot(index="#FILE", columns="GENE", values="%COVERAGE")
            Abricate_subset = Abricate_subset.notna().astype(int)
            Abricate_subset = Abricate_subset[list(set(genes))]
            merged = pd.merge(Abricate_subset,AMR_subset, left_on =  Abricate_subset.index,right_on = "Genome ID", how = "inner")

        elif annotation_tool == "StarAMR":
            annotations["Isolate ID"] = annotations["Isolate ID"].str.replace('.fna$', '', regex=True)
            annotations.loc[:, "Isolate ID"] = annotations["Isolate ID"].astype(str)
            if '/' in antibiotic:
                antibiotic = antibiotic.split('/')
                genes1 = list(set(annotations[annotations["CGE Predicted Phenotype"].str.contains(antibiotic[0], case=False, na=False)]["Data"]))
                genes2 = list(set(annotations[annotations["CGE Predicted Phenotype"].str.contains(antibiotic[1], case=False, na=False)]["Data"]))
                genes = genes1 + genes2
            else:
                genes = list(set(annotations[annotations["CGE Predicted Phenotype"].str.contains(antibiotic, case=False, na=False)]["Data"]))
            StarAMR_subset = annotations[["Data","Isolate ID","%Identity"]]
            StarAMR_subset = StarAMR_subset.drop_duplicates(subset=["Isolate ID", "Data"])
            df_pivot = StarAMR_subset.pivot(index="Isolate ID", columns="Data", values="%Identity")
            df_pivot = df_pivot.notna().astype(int)
            df_pivot = df_pivot[genes]
            merged = pd.merge(df_pivot, AMR_subset, left_on= df_pivot.index, right_on = "Genome ID", how = "inner")
        
        elif annotation_tool == "SraX":
            annotations["Genome.ID"] = annotations["Genome.ID"].str.replace('.fna$', '', regex=True)
            annotations.loc[:, "Genome.ID"] = annotations["Genome.ID"].astype(str)
            antibiotic_classes_subset = antibipotic_classes[antibipotic_classes["Drugs"] == antibiotic]
            antibiotic_classes_subset["Subclass"] = antibiotic_classes_subset["Subclass"].str.lower()
            annotations["Drug class"]= annotations["Drug class"].str.lower()
            SraX_subset = annotations[annotations["Genome.ID"].isin(AMR_subset["Genome ID"])]
            mask = SraX_subset["Drug class"].isin(antibiotic_classes_subset["Subclass"])
            genes = SraX_subset[mask]["ARG"]
            SraX_subset = SraX_subset[["Genome.ID","Identity (%)", "ARG"]]
            SraX_subset = SraX_subset.drop_duplicates(subset=["Genome.ID", "ARG"])
            SraX_subset = SraX_subset.pivot(index="Genome.ID", columns="ARG", values="Identity (%)")
            SraX_subset = SraX_subset.notna().astype(int)
            SraX_subset = SraX_subset[list(set(genes))]
            merged = pd.merge(SraX_subset,AMR_subset, left_on =  SraX_subset.index,right_on = "Genome ID", how = "inner")

        if annotation_tool == "Kleborate":
            print(antibiotic)
            merge_antibiotic_classes = pd.merge(antibiotic_classes, genes_to_class, left_on = "Subclass", right_on = "Antibiotic", how = "inner")
            resistances = annotations.set_index('strain')
            arr = resistances.values
            flat_list = arr.flatten().tolist()
            split_values = [item.split(';') for item in flat_list]
            split_values = [item for sublist in split_values for item in sublist]
            split_values = list(set(split_values))
            merge_antibiotic_class = list(set(merge_antibiotic_classes[merge_antibiotic_classes["Drugs"] == antibiotic]["Kleborate_genes"]))
            set_class = resistances[merge_antibiotic_class]
            array_to_add = set_class.values
            flat_list = array_to_add.flatten().tolist()
            split_values = [item.split(';') for item in flat_list]
            split_values = [item for sublist in split_values for item in sublist]
            split_values = list(set(split_values))
            # Assuming resistances is a DataFrame and resistances.index gives the index
            df_class = pd.DataFrame(index=resistances.index, columns=split_values)
            for index, row in tqdm(resistances.iterrows(), total=resistances.shape[0]):
                row_values = [item.split(';') for item in row.values]
                row_values = [item for sublist in row_values for item in sublist]
                for col in df_class.columns:
                    if col in row_values:
                        df_class.loc[index,col] = 1
                    else: 
                        df_class.loc[index,col] = 0
            AMR_subset["Genome ID"] = AMR_subset["Genome ID"].astype(str)
            df_class.index = df_class.index.astype(str)
            AMR_subset = AMR[AMR["Antibiotic"]==antibiotic]
            AMR_subset = AMR_subset[AMR_subset["Genome ID"].isin(df_class.index)]
            AMR_subset = AMR_subset[["Genome ID", "Resistant Phenotype"]]
            AMR_subset["Resistant Phenotype"] = np.where(AMR_subset["Resistant Phenotype"] == 'Susceptible', 0, 1)
            Kleborate_antibiotic = df_class[df_class.index.isin(AMR_subset["Genome ID"])]
            merged = pd.merge(Kleborate_antibiotic,AMR_subset, left_on =  Kleborate_antibiotic.index,right_on = "Genome ID", how = "inner")
              
        antibiotic_performance = []
        sensitivity_performances = []
        specificity_performances = []
        #start predictions
        if len(merged.columns) >2:
            X = merged.drop(['Resistant Phenotype','Genome ID'], axis=1)
            y = merged['Resistant Phenotype']
            skf = StratifiedKFold(n_splits=5)
            for i, (train_index, test_index) in enumerate(skf.split(X, y)):
                X_train, X_test = X.iloc[train_index], X.iloc[test_index]
                y_train, y_test = y.iloc[train_index], y.iloc[test_index]
                columns = X.columns
                scaler = StandardScaler()
                X_train = scaler.fit_transform(X_train)
                X_test = scaler.transform(X_test)
                y_train = np.array(y_train)
                y_test = np.array(y_test)
                if model == "XGBoost":
                    model = XGBClassifier(random_state=42, max_depth = 3, n_estimators = 30, eval_metric='logloss',seed = 1024)
                elif model == "ElasticNet":
                    model = LogisticRegression(C=0.1, penalty='elasticnet', solver='saga',l1_ratio=0.5, random_state=10)
                model.fit(X_train,y_train)
                y_pred = model.predict(X_test)
                AUC = roc_auc_score(y_test, y_pred)
                sensitivity = sensitivity_score(y_test, y_pred)
                specificity = specificity_score(y_test, y_pred)
                antibiotic_performance.append(AUC)
                sensitivity_performances.append(sensitivity)
                specificity_performances.append(specificity)

        #summarise results
        if len(antibiotic_performance) > 0:
            average = sum(antibiotic_performance) / len(antibiotic_performance)
            average = (f"{average:.3f}")
            sensitivity = sum(sensitivity_performances) / len(sensitivity_performances)
            sensitivity = (f"{sensitivity:.3f}")
            specificity = sum(specificity_performances) / len(specificity_performances)
            specificity = (f"{specificity:.3f}")
            average_across_antibiotics.append(average)
            average_across_antibiotics_sensitivity.append(sensitivity)
            average_across_antibiotics_specificity.append(specificity)
            number_of_genes.append(len(X.columns))
        elif len(antibiotic_performance) == 0:
            average_across_antibiotics.append(0)
            average_across_antibiotics_sensitivity.append(0)
            average_across_antibiotics_specificity.append(0)
            number_of_genes.append(0)
    antibiotics_to_write = antibiotics.copy()
    antibiotics_to_write.insert(0, "tool")
    average_across_antibiotics.insert(0, tool)
    average_across_antibiotics_sensitivity.insert(0, tool)
    average_across_antibiotics_specificity.insert(0, tool)
    number_of_genes.insert(0, tool)
    antibiotics_to_write = ', '.join(antibiotics_to_write).replace("'", "").replace("[", "")
    average_across_antibiotics = ', '.join(average_across_antibiotics).replace("'", "").replace("[", "")
    average_across_antibiotics_sensitivity = ', '.join(average_across_antibiotics_sensitivity).replace("'", "").replace("[", "")
    average_across_antibiotics_specificity = ', '.join(average_across_antibiotics_specificity).replace("'", "").replace("[", "")
    
    if file.tell() == 0:
        file.write(f'{antibiotics_to_write}.\n')
    if file_2.tell() == 0:
        file_2.write(f'{antibiotics_to_write}.\n')
    if file_3.tell() == 0:
        file_3.write(f'{antibiotics_to_write}.\n')
    if file_4.tell() == 0:
        file_4.write(f'{antibiotics_to_write}.\n')

    file.write(f'{average_across_antibiotics}\n')
    file.close()

    file_3.write(f'{average_across_antibiotics_sensitivity}\n')
    file_3.close()

    file_4.write(f'{average_across_antibiotics_specificity}\n')
    file_4.close()

    number_of_genes.insert(0, tool)
    file_2.write(f'{number_of_genes}\n')
    file_2.close()
    pass

def main():
    parser = argparse.ArgumentParser(description='Predict phenotype based on annotations.')
    parser.add_argument('--annotations', required=True, help='Path to the annotations file')
    parser.add_argument('--phenotype', required=True, help='Phenotype table containing phenotypes in the same format as the BV-BRC database')
    parser.add_argument('--Kleborate', required=False, help='Path to the Kleborate annotations if applicable')
    parser.add_argument('--model', required=True, help='Which model to use - XGboost or ElasticNet')
    parser.add_argument('--annotation_tool', required=True, help='Which annotation tool are we using for prediction')
    parser.add_argument('--output_path', required=True, help='Folder to save the output')
    parser.add_argument('--antibipotic_classes', required=True, help='Antibiotic classes to consider')
    parser.add_argument('--genes_to_class', required=True, help='Genes to class mapping file')

    args = parser.parse_args()

    predict_phenotype(
        annotations=args.annotations,
        phenotype=args.phenotype,
        Kleborate=args.Kleborate,
        model=args.model,
        annotation_tool=args.annotation_tool,
        output_path=args.output_path,
        antibipotic_classes=args.antibipotic_classes,
        genes_to_class=args.genes_to_class
    )

if __name__ == '__main__':
    main()
