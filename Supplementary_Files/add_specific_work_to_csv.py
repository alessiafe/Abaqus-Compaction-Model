import pandas as pd
import os

def process_csvs(folder_path):
    # Loop over all csv file in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path)
            # Check if columns exist
            if {'work', 'tissue_volume', 'cavity_volume'}.issubset(df.columns):
                # Extract firs value of tissue_volume and cavity_volume
                tissue_vol_0 = df['tissue_volume'].iloc[0]
                cavity_vol_0 = df['cavity_volume'].iloc[0]
                vol_tot = tissue_vol_0 + cavity_vol_0
                # Add new columns
                df['work/vol_tot'] = df['work'] / vol_tot
                df['work/vol_cw'] = df['work'] / tissue_vol_0
                # Save file
                df.to_csv(file_path, index=False)
                print(f"File updated: {filename}")
            else:
                print(f"Missing columns in {filename}, Skip.")
    return

folder_path = '/home/aferrara/Desktop/abaqus-compaction-model/MD/LW3'
process_csvs(folder_path)
