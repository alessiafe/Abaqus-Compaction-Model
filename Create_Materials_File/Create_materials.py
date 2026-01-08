import numpy as np
from scipy.io import loadmat
import importlib as imp
import os
from copy import copy
import Micromaterial_Calculation as micromat
imp.reload(micromat)
import pickle
from collections import OrderedDict
from mpmath import *
mp.dps = 150 # increase precision to stabilize matrix operations (default = 15)

# Functions
def extract_elements(matrix):
    elements = OrderedDict([
        ('C11', matrix[0, 0]),
        ('C22', matrix[1, 1]),
        ('C33', matrix[2, 2]),
        ('C12', matrix[0, 1]),
        ('C13', matrix[0, 2]),
        ('C23', matrix[1, 2]),
        ('C44', matrix[3, 3]),
        ('C55', matrix[4, 4]),
        ('C66', matrix[5, 5]),
    ])
    return elements

def make_isotropic(matrix):
    diag_mean1 = np.mean([matrix[0, 0], matrix[1, 1], matrix[2, 2]])
    diag_mean2 = np.mean([matrix[3, 3], matrix[4, 4], matrix[5, 5]])
    matrix_temp = matrix[:3,:3]
    off_diag_mean = np.mean(matrix_temp[np.triu_indices(3, k=1)])
    isotropic_matrix = copy(matrix)
    isotropic_matrix[:3, :3] = np.full_like(matrix[:3, :3], off_diag_mean)  # Fill with off-diagonal mean
    np.fill_diagonal(isotropic_matrix[:3, :3], diag_mean1)  # Replace first 3 diagonal elements
    np.fill_diagonal(isotropic_matrix[3:, 3:], diag_mean2)
    return isotropic_matrix

# Set structure to process
structname = 'LW3' # user input: EW3, TW3, LW3

# Set paths
mainpath = '/home/aferrara/Desktop/abaqus-compaction-model/Calculate_Stiffness_3D/'
structpath = mainpath + 'Structure_mat/'
structfile = 'testspruce_' + structname + '_mstruct'
filepath = structpath + structfile

script_path = os.path.abspath(__file__)
matpath = os.path.dirname(script_path) + '/'

# Load structure data
data = loadmat(filepath, squeeze_me=True, struct_as_record=False)
mstruct = data['mstruct']

# Calculate properties
moisture_arr = np.arange(0., 0.31, 0.01)
ligred_arr = np.arange(0.01, 1.01, 0.01)
alpha=[0., 0., 60., -60., 15., 75.]

# Initialize a dictionary to hold all matrices
all_matrices = {}

# Loop over moisture and lignin content
for moisture in moisture_arr:
    for ligred in ligred_arr:
        ########## CREATE LAYERS MATERIALS ##########
        materialStiffMat = micromat.layer_mat_prop(moisture, matpath, 'mht_full', ligred)[0]
        
        ########## STORE MATRIX MATERIAL ##########
        C = [make_isotropic(micromat.full_laminate_matrix(alpha[:2], mstruct.fibers[0].lthicks[:2], materialStiffMat[:2]))]  # Initialize C for each combination
        ########## STORE FIBERS MATERIAL ##########       
        for i in range(len(mstruct.fibers)):
            f1 = mstruct.fibers[i]  # fiber data
            lthicks = f1.lthicks  # layers thickness
            C.append(micromat.full_laminate_matrix(alpha[2:], lthicks[2:], materialStiffMat[2:]))  # add the laminate matrix

        # Store matrices in dictionary with a key of moisture-lignin combination
        key = f"{int(moisture * 100)}-{int(ligred * 100)}"
        extracted_data = [extract_elements(matrix) for matrix in C]
        all_matrices[key] = np.array(extracted_data)  # Convert the list to a numpy array

# Save all matrices as a .npz file
npz_file_path = structpath + 'testspruce_' + structname + '_materials.pkl'
with open(npz_file_path, 'wb') as f:
    # Use pickle to dump the data with protocol 2
    pickle.dump(all_matrices, f, protocol=2)
#np.savez(npz_file_path, **all_matrices)

print(f"Data saved to {npz_file_path}")