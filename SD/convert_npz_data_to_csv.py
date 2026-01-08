import numpy as np
import pandas as pd
import os
from collections import OrderedDict

# Set paths
structure = 'EW3'
scriptpath = os.path.abspath(__file__) # absolute path of current script
mainpath = os.path.dirname(scriptpath) # directory of current script
savepath = mainpath + '/' + structure + '/'

# Get file names
files = os.listdir(mainpath)
jobName = [os.path.splitext(file)[0] for file in files if file.endswith('.odb') and file.startswith(structure)]

# Read npz file and convert to csv
for i in range(len(jobName)):   
    data = np.load(savepath+jobName[i]+'.npz', allow_pickle=True)
    df = pd.DataFrame(data=OrderedDict((key, data[key]) for key in data))
    df.to_csv(savepath+jobName[i]+'.csv', index=False)

print('Done!')