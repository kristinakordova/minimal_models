def RGI_dataframe(annotations, antibiotic, AMR_subset):
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
    return merged

def AMRFinder_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset):
    antibiotic_classes_subset = antibiotic_classes[antibiotic_classes["Drugs"]==antibiotic]
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
    return merged

def ResFinder_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset):
    antibiotic_classes_subset = antibiotic_classes[antibiotic_classes["Drugs"] == antibiotic]
    antibiotic_classes_subset = list(set(antibiotic_classes_subset["Subclass"]))
    antibiotic_classes_subset = [x.lower() for x in antibiotic_classes_subset]
    annotations["Genome.ID"] = annotations["Genome.ID"].str.replace('.fna$', '', regex=True)
    annotations.loc[:, "Genome.ID"] = annotations["Genome.ID"].astype(str)
    ResFinder_drug = annotations[annotations["Genome.ID"].isin(AMR_subset["Genome ID"])]
    AMR_subset = AMR_subset[AMR_subset["Genome ID"].isin(ResFinder_drug["Genome.ID"])]
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
    return merged

def DeepARG_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset):
    antibiotic_classes_subset = antibiotic_classes[antibiotic_classes["Drugs"] == antibiotic]
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
    return merged

def Abricate_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset):
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
    return merged

def StarAMR_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset):
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
    return merged

def SraX_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset):
    annotations["Genome.ID"] = annotations["Genome.ID"].str.replace('.fna$', '', regex=True)
    annotations.loc[:, "Genome.ID"] = annotations["Genome.ID"].astype(str)
    antibiotic_classes_subset = antibiotic_classes[antibiotic_classes["Drugs"] == antibiotic]
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
    return merged

def Kleborate_dataframe(annotations, antibiotic_classes, antibiotic, AMR_subset,genes_to_class):
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
    return merged

def predict(merged, model):
    AUC_performance = []
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
            AUC_performance.append(AUC)
            sensitivity_performances.append(sensitivity)
            specificity_performances.append(specificity)
            return AUC_performance, sensitivity_performances, specificity_performances

      def summarise_results(AUC_performance, sensitivity_performances, specificity_performances, merged, antibiotic):
    #summarise results

    if len(AUC_performance) > 0:
        AUC = sum(AUC_performance) / len(AUC_performance)
        AUC = (f"{AUC:.3f}")
        sensitivity = sum(sensitivity_performances) / len(sensitivity_performances)
        sensitivity = (f"{sensitivity:.3f}")
        specificity = sum(specificity_performances) / len(specificity_performances)
        specificity = (f"{specificity:.3f}")
        X = merged.drop(['Resistant Phenotype','Genome ID'], axis=1)
        number_of_genes = len(X.columns)
    elif len(AUC_performance) == 0:
        AUC = 0
        sensitivity = 0
        specificity = 0
        number_of_genes = 0
    print(f"{antibiotic}, {AUC}, {sensitivity}, {specificity}, {number_of_genes}")

