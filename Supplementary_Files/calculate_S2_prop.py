# Imports
from math import *
import numpy as np
import pandas as pd
import importlib
import os
import Micromaterial_Calculation as micromat
importlib.reload(micromat)

# Set path
runpath = os.path.dirname(os.path.abspath(__file__))+'/' # main path

# Set material parameters
moisture_arr = [0, 0.12, 0.3] #np.arange(0.,0.31,0.02) # moisture content (from 0 to 1)
ligred_arr = [1.] #np.arange(0.1,1.05,0.1) # lignin stiffness reduction factor (from 0 -> no stiffness to 1 -> full stiffness)

# Calculate material properties
materialEngConst = []
for k in range(len(moisture_arr)):
    for j in range(len(ligred_arr)):
        moisture = moisture_arr[k] # moisture content (from 0 to 0.3)
        ligred = ligred_arr[j] # lignin stiffness reduction factor (from 0 -> no stiffness to 1 -> full stiffness)
        materialStiffMat = micromat.layer_mat_prop(moisture,runpath,'mht_full',ligred)[0][4] # [GPa] local stiffness matrices of S2 layer
        # Calculate engineering constants: [E1, E2, E3, nu23, nu13, nu12, G23, G13, G12]
        # with:
        # - 1 = fibrils direction
        # - 2 = normal to fibrils direction
        # - 3 = normal material plane
        materialEngConst.append(np.insert(micromat.convert_mat_prop(materialStiffMat, 'o', 'c', 'e')[0], 0, [moisture, ligred])) # [GPa] local engineering constants of S2 layer

# Create dataframe
column_names = ['mc', 'lig', 'E1', 'E2', 'E3', 'nu23', 'nu13', 'nu12', 'G23', 'G13', 'G12']
df = pd.DataFrame(materialEngConst, columns=column_names)

# Save dataframe as csv
df.to_csv(runpath+'S2_eng_const.csv', index=False)
print(df)
print("Done!")

