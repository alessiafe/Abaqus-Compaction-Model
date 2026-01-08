# Micromaterial calculation for equivalent single-layer cell wall model

This folder contains all scripts and files used for the **micromaterial calculations of the equivalent single-layer cell wall model**, composed of an isotropic compound middle lamella (CML), combining the middle lamella and primary wall, and an equivalent single-layer representing all S layers (S1, S2, S3).

## Scripts and files

- **`Create_materials.py`**  
  Script to calculate the stiffness matrix of the CML (matrix) and the equivalent single layer (fiber) for different moisture contents (0 to 0.3) and lignin reduction factors (0 excluded to 1) for a selected tissue (user input). This calculation is performed in Python, outside the Abaqus simulations, to reduce approximation errors in the matrix computations.
  The stiffness coefficients in GPa are saved as a `.pkl` file in the specified folder (user input).

- **`Micromaterial_Calculation.py`**  
  Library of functions to compute the hygroelastic properties of the wood cell wall layers, as well as the equivalent properties of planar and 3D laminate composites.

- **`basic_components.npy`**  
  NumPy array storing the stiffness matrices and hygroexpansion coefficients of the basic wood cell wall components (cellulose, hemicellulose, and lignin).  

- **`materials_composite.npy`**  
  NumPy array storing data required to compute the composite materials (rule, fiber, fiber volume fraction, matrix, matrix volume fraction, and geometry coefficients) derived from the basic components (rule, engineering constants, hygroexpansion coefficients).
