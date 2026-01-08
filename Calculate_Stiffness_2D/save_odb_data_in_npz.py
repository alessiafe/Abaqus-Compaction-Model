import imp
import os
import glob
import numpy as np
from collections import OrderedDict

# Import Compaction_Class
import Compaction_Class
imp.reload(Compaction_Class)
from Compaction_Class import Compaction_Class
compaction = Compaction_Class() # create compaction class instance

# Set paths
structure = 'EW3'
scriptpath = os.path.abspath('save_odb_data_in_npz_py') # absolute path of current script
mainpath = os.path.dirname(scriptpath) # directory of current script
structpath = mainpath + '/Structure_mat/'
savepath = mainpath + '/' + structure + '/' + structure + '_stiffness_coeff.npz'
if not os.path.exists(mainpath + '/' + structure): os.makedirs(mainpath + '/' + structure)

# Read job names    
files = os.listdir(mainpath)
jobName = list({os.path.splitext(file)[0].rsplit('-', 1)[0] for file in files if file.endswith('.odb') and file.startswith(structure)})
# Open cae file
filepath = ', '.join(glob.glob(structpath + '*' + structure + '_mstruct.mat'))
caepath = ', '.join(glob.glob(structpath + '*' + structure + '.cae'))
modelName = compaction.open_model(caepath, filepath)

# Store results
C = []
for i in range(len(jobName)):
    print(jobName[i])
    stiffMat = compaction.calculate_2D_stiffness_matrix(jobName[i])
    info = compaction.extract_info(jobName[i])
    C.append(info + [stiffMat[0][0], stiffMat[1][1], stiffMat[0][1], stiffMat[2][2]])

# Save results as dictionary
column_names = ['structure', 'mc', 'lig', 'C11', 'C22', 'C12', 'C66']   
C_dict = OrderedDict((col, [row[i] for row in C]) for i, col in enumerate(column_names))
column_names = list(C_dict.keys())
np.savez(savepath,**C_dict)

print('Done!')
