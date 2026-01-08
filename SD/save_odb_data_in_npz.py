import imp
import os
import glob

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
savepath = mainpath + '/' + structure + '/'
if not os.path.exists(mainpath + '/' + structure): os.makedirs(mainpath + '/' + structure)

# Read job names    
files = os.listdir(mainpath)
jobName = [os.path.splitext(file)[0] for file in files if file.endswith('.odb') and file.startswith(structure)]

# Open cae file
filepath = ', '.join(glob.glob(structpath + '*' + structure + '_mstruct.mat'))
caepath = ', '.join(glob.glob(structpath + '*' + structure + '.cae'))
modelName = compaction.open_model(caepath, filepath)
# Save results
for i in range(len(jobName)):
    compaction.save_results(jobName[i], savepath)
    
print('Done!')
