import numpy as np
import pandas as pd
import os
from collections import OrderedDict

# Set paths
structure = 'EW3'
scriptpath = os.path.abspath(__file__) # absolute path of current script
mainpath = os.path.dirname(scriptpath) # directory of current script
savepath = mainpath + '/' + structure + '/' + structure + '_stiffness_coeff.npz'

# Read npz file and convert to csv
data = np.load(savepath, allow_pickle=True)
column_names = ['structure', 'mc', 'lig', 'C11', 'C22', 'C33', 'C12', 'C13', 'C23', 'C44', 'C55', 'C66'] 
data = OrderedDict((col, data[col].tolist()) for col in column_names)
for key, value in data.items():
    if isinstance(value[0], bytes):  # Check if the first element is a byte string
        data[key] = [v.decode('utf-8') for v in value]

df = pd.DataFrame(data=OrderedDict((key, data[key]) for key in data))
df = df.sort_values(by=['mc','lig'])
savepath = mainpath + '/' + structure + '/' + structure + '_stiffness_coeff.csv'
df.to_csv(savepath, index=False)

print('Done!')
