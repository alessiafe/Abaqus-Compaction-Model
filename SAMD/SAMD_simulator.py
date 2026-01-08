import imp
import os
import glob
import time
import math
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
structname = 'testspruce_EW3'
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
minSizeFiber = 1. # fiber mesh size factor
meshSizeMatrix = 1. # matrix mesh global size
minSizeMatrix = 0.9 # matrix mesh size factor
lengthRatio = 0.1 # distortion aspect ratio
massScale = 1e12 # mass scaling factor

# Simulation parameters
rps_arr = [2.] #np.arange(0.5, 5.1, 0.5) # frequency [Hz]
amp_arr = [5.] #np.arange(0.5, 10.1, 0.5) # amplitude [microm]
jobName = []

# Create model from scratch and save
#modelName = compaction.setup_geometry(filepath) 
#compaction.save_model(caepath)

for i in range(len(amp_arr)):
    for k in range(len(rps_arr)):

        # Open file
        modelName = compaction.open_model(caepath, filepath) # the model has been already created

        # Set material properties
        moisture = 0.24 # moisture content (from 0 to 1)
        ligred = 0.8 # lignin stiffness reduction factor (from 0 -> no stiffness to 1 -> full stiffness)
        compaction.setup_materials(matpath, moisture, matrix_density, fiber_density, ligred=ligred)
        
        # Set mesh
        compaction.create_mesh(meshSizeFiber, meshSizeMatrix, minSizeFiber, minSizeMatrix, lengthRatio)
        
        # Create step
        stepName = 'Step-1'
        stepTime = 4 # step time
        
        # Set shear-assisted mechanical densification boundary conditions
        rps = rps_arr[k] # frequency [Hz]
        shearFreq = rps * 2*math.pi # angular frequency [rad/s]
        shearAmp = amp_arr[i] # amplitude [microm]
        numIntervals = 1000 # # number of intervals to store data
        compactionLevel = 0.59 # degree of compaction as fraction of the total height
        compaction.shear_mechanical_densification_bc(modelName, stepName, stepTime, compactionLevel, shearAmp, shearFreq, numIntervals)
        
        # Run simulation
        jobName.append(structname[-3:]+'-SAMD'+'-rps{}'.format(int(rps))+'-amp{}'.format(int(10*shearAmp)))
        compaction.submit_job(jobName[-1], modelName)

time.sleep(10)
for j in range(len(jobName)):
    # Remove simulation files except odb
    path_all = glob.glob(os.getcwd() + '/' + jobName[j] + '.*')
    extension = ('.odb') # extension to keep
    path_remove = [filename for filename in path_all if not filename.endswith(extension)]
    [os.remove(filePath) for filePath in path_remove]


