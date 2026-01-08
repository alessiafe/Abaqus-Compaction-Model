import imp
import os
import glob
import time
from mpmath import *
mp.dps = 150 # increase precision to stabilize matrix operations (default = 15)

# Import Compaction_Class
import Compaction_Class
imp.reload(Compaction_Class)
from Compaction_Class import Compaction_Class
compaction = Compaction_Class() # create compaction class instance

# Set paths
mainpath = os.path.abspath('compaction_simulator.py')
mainpath = str(os.path.dirname(mainpath)) + '/'
structpath = mainpath + 'Structure_mat/'
structname = 'testspruce_TW3'
structfile = structname + '_mstruct'
filepath = structpath + structfile
savename = structname + '.cae'
caepath = structpath + savename

# Material data   
matpath = structpath + structname + '_materials.pkl' # path to material data files     
matrix_density = 1500. # [kg/m3]
fiber_density = 1500. # [kg/m3]

# Mesh parameters
meshSizeFiber = 2.5 # fiber mesh global size
minSizeFiber = 0.8 # fiber mesh size factor
meshSizeMatrix = 1.5 # matrix mesh global size
minSizeMatrix = 0.8 # matrix mesh size factor
lengthRatio = 0.1 # distortion aspect ratio
massScale = 1e12 # mass scaling factor

# Simulation parameters
ligred_arr = [1.] #np.arange(0.8,1.1,0.2) # lignin reduction factor (from 0=no stiffness to 1=full stiffness)
moisture_arr = [0.12] #np.arange(0.,0.31,0.06) # moisture content (from 0 to 0.3)
jobName = []

# Create model from scratch and save
#modelName = compaction.setup_geometry(filepath) 
#compaction.save_model(caepath)

for k in range(len(ligred_arr)):
    for i in range(len(moisture_arr)):

        # Open file
        modelName = compaction.open_model(caepath, filepath) # the model has been already created

        # Set material properties
        moisture = moisture_arr[i] # moisture content
        ligred = ligred_arr[k] # lignin stiffness reduction factor
        compaction.setup_materials(matpath, moisture, matrix_density, fiber_density, ligred=ligred)

        # Set mesh
        compaction.create_mesh(meshSizeFiber, meshSizeMatrix, minSizeFiber, minSizeMatrix, lengthRatio)
        
        # Create step
        stepName = 'Step-1'
        stepTime = 1 # step time
        
        # Set mechanical densification booundary conditions
        compactionLevel = 0.5 # degree of compaction as fraction of the total height
        numIntervals = 1000 # number of intervals to store data
        compaction.mechanical_densification_bc(modelName, stepName, stepTime, compactionLevel, numIntervals, massScale=massScale)

        # Run simulation
        jobName.append(structname[-3:]+'-MD'+'-mc{}'.format(int(moisture*100))+'-lig{}'.format(int(ligred*100)))
        compaction.submit_job(jobName[-1], modelName)

time.sleep(10)
for j in range(len(jobName)):
    # Remove simulation files except odb
    path_all = glob.glob(os.getcwd() + '/' + jobName[j] + '.*')
    extension = ('.odb') # extension to keep
    path_remove = [filename for filename in path_all if not filename.endswith(extension)]
    [os.remove(filePath) for filePath in path_remove]

