import os
import pandas as pd
final_dataframe = pd.DataFrame()
with open("complete.txt", 'r') as file:
    for file_name in file:
        file_name = file_name.strip()
        file_name = file_name.replace("/", "")
        print(file_name)
        # Define the file path
        file_path = f"{file_name}/pheno_table.txt"
        # Check if the file exists
        if not os.path.exists(file_path):
            print(f"{file_path} is missing")
            continue
        data = []
        with open(file_path, 'r') as phenotype_file:
            for line in phenotype_file:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    line = line[:3]
                    line.insert(0, file_name)
                    data.append(line)
        # Create a DataFrame from the structured data
        columns = ["Genome.ID", "Antibiotic", "Class", "Phenotype_Res_Finder"]
        df = pd.DataFrame(data, columns=columns)
        # Concatenate the DataFrame
        final_dataframe = pd.concat([final_dataframe, df], ignore_index=True)
