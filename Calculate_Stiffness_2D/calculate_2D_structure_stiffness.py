import imp
import os
import glob
import time

# Set paths
mainpath = os.path.abspath('calculate_2D_structure_stiffness.py')
mainpath = str(os.path.dirname(mainpath)) + '/'
structpath = mainpath + 'Structure_mat/'
savepath = mainpath

# Import Compaction_Class
import Compaction_Class
imp.reload(Compaction_Class)
from Compaction_Class import Compaction_Class
compaction = Compaction_Class() # create compaction class instance

# Material data   
matrix_density = 1500. # matrix density [kg/m3]
fiber_density = 1500. # fiber density [kg/m3]

# Mesh parameters
struct_arr = ['testspruce_LW3'] # set structure(s) to test: testspruce_EW3, testspruce_TW3, testspruce_LW3
meshSizeFiber = 4. # fiber mesh global size
minSizeFiber = 0.9 # fiber mesh size factor
meshSizeMatrix = 1.5 # matrix mesh global size
minSizeMatrix = 0.9 # matrix mesh size factor
lengthRatio = 0.1 # distortion aspect ratio

# Simulation parameters
ligred_arr = [1.] #np.arange(0.2,1.1,0.2) # lignin reduction factor (from 0=no stiffness to 1=full stiffness)
moisture_arr = [0.12] #np.arange(0.12,0.31,0.06) # moisture content (from 0 to 0.3)

jobName = []
for k in range(len(struct_arr)):
    structname = struct_arr[k]
    structfile = structname + '_mstruct'
    filepath = structpath + structfile
    caename = structname + '.cae'
    caepath = structpath + caename
    matpath = structpath + structname + '_materials.pkl'

    # Create model from scratch and save
    #modelName = compaction.setup_geometry(filepath) 
    #compaction.save_model(caepath)

    for i in range(len(moisture_arr)):
        for j in range(len(ligred_arr)):

            # Open file
            modelName = compaction.open_model(caepath, filepath) # the model has been already created

            # Set material properties
            moisture = moisture_arr[i] # moisture content 
            ligred = ligred_arr[j] # lignin stiffness reduction factor 
            compaction.setup_materials(matpath, moisture, matrix_density, fiber_density, ligred=ligred)
            compaction.create_surface_partition(1., sides=True) # create surfaces to apply boundary conditions

            # Set mesh       
            compaction.create_mesh_std(meshSizeFiber, meshSizeMatrix, minSizeFiber, minSizeMatrix, lengthRatio)

            # Set general constraints
            stepName = 'Step-1'
            compaction.set_general_constraints(modelName, stepName)
            
            # Calculate stiffness matrix
            eps = 0.02 # strain [-]
            jobName.append(structname[-3:]+'-mc{}'.format(int(moisture*100))+'-lig{}'.format(int(ligred*100)))
            compaction.compute_2D_stiffness_matrix(eps, stepName, jobName[-1], savepath)
          
time.sleep(10)
for j in range(len(jobName)):
    # Remove simulation files except odb
    path_all = glob.glob(os.getcwd() + '/' + jobName[j] + '*')
    extension = ('.odb') # extension to keep
    path_remove = [filename for filename in path_all if not filename.endswith(extension)]
    [os.remove(filePath) for filePath in path_remove]         