import os
import pandas as pd
import bambi as bmb
import arviz as az
import numpy as np

# Set this to the parent directory that contains genotype folders (e.g., "WT", "KO", etc.)
parent_dir = "D:\E209K-PRG Current Clamp + Minis p21-24\minis\exctitatory"

# List to collect all dataframes
all_data = []

# Iterate over each subfolder (assumed to be a genotype folder)
for genotype_folder in os.listdir(parent_dir):
    genotype_path = os.path.join(parent_dir, genotype_folder)

    # Ensure it is a directory
    if os.path.isdir(genotype_path):
        genotype = genotype_folder  # Folder name is the genotype label

        # Look through files in this folder
        for filename in os.listdir(genotype_path):
            if filename.endswith(".xlsx") or filename.endswith(".xls"):
                file_path = os.path.join(genotype_path, filename)

                # Extract IDs
                animal_id = filename[:8]
                cell_id = os.path.splitext(filename)[0]

                try:
                    # Load Excel file (first/default sheet)
                    full_df = pd.read_excel(file_path)

                    # Add metadata columns
                    full_df['genotype'] = genotype
                    full_df['animal_id'] = animal_id
                    full_df['cell_id'] = cell_id

                    all_data.append(full_df)

                except Exception as e:
                    print(f"⚠️ Error reading {file_path}: {e}")

# Combine all loaded data
if all_data:
    full_df = pd.concat(all_data, ignore_index=True)
    print(full_df.head())

    # Optional: Save to CSV
    # full_df.to_csv("combined_mEPSC_data.csv", index=False)
else:
    print("No data found.")

full_df = full_df.rename(columns={"IEI (ms)": "iei"})
full_df = full_df[full_df['iei'] > 0]  # remove zero or negative intervals if any

# Calculate count per event and log(IEI) for offset
full_df['count'] = 1  # each row is one event
full_df['log_iei'] = np.log(full_df['iei'])  # offset (exposure to time)

# Ensure identifiers are categorical
full_df['genotype'] = full_df['genotype'].astype('category')
full_df['animal_id'] = full_df['animal_id'].astype('category')
full_df['cell_id'] = full_df['cell_id'].astype('category')

# Fit hierarchical Poisson regression
# Nesting: cells within animals, animals within genotypes (via grouping)
model = bmb.Model(
    formula="count ~ genotype + (1|animal_id) + (1|cell_id) + offset(log_iei)",
    data=full_df,
    family="poisson"
)

# Fit the model using MCMC
results = model.fit(draws=1000, tune=1000, chains=4, cores=4)

# Summarize posterior estimates

az.summary(results, hdi_prob=0.95)