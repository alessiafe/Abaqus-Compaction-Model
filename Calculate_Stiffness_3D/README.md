# Calculation of stiffness of 3D tissue models

This folder contains all scripts and files used for the **calculation of the equivalent stiffness of 3D tissue models**, obtained by applying 8 elementary strain cases: two biaxial, three triaxial, and three shear.

It is recommended to apply strain in both the radial and tangential directions. Fixing one of them to 0 can lead to over-constrained conditions, which may distort the results.

The coordinate system follows the anatomical directions of wood tissues:  
**x** → tangential, **y** → radial, **z** → longitudinal.

## Folder structure

- **`Structure_mat`**  
  Contains:  
  - all `CAE models` (and corresponding `.jnl` files) **generated from the `.mat` files** also included in this folder,  
  - the `.mat` files defining the coordinates of the points for each fiber and the outer contour of the matrix (previously **produced via MATLAB**), and  
  - the `.pkl` files containing the **pre-calculated material parameters** for matrix and fiber at various moisture contents and lignin reduction factors.  

  The available tissue structures are one for each type: EW (`testspruce_EW3`), TW (`testspruce_TW3`), and LW (`testspruce_LW3`).

  The CAE files containing only the geometry are typically **generated and saved before running the simulations**, as their creation is time-consuming. When creating the geometry, the **thickness in the longitudinal (z) direction** must be specified as a user input. These pre-generated CAE models are then opened and used to perform the simulations.

  There are also CAE files ending in `_full` that contain the complete model ready to run.
  
- **`calculate_2D_structure_stiffness.py`**  
  Main script to be **run from Abaqus**. Performs simulations of the four strain cases and computes the equivalent stiffness tensor of the 2D tissue model. The tissue to be simulated must be specified by the user (`testspruce_EW3`, `testspruce_TW3`, or `testspruce_LW3`), along with:
  - *moisture content* from 0 to 0.3 (array)
  - *lignin reduction factor* from 0 (excluded) to 1 (array)
  - *mesh parameters* (global size and size factors for both fiber and matrix).  
    The values used are:
    | Parameter             | EW  | TW  | LW  |
    |-----------------------|-----|-----|-----|
    | Fiber mesh size       | 2.5 | 3.  | 4.  |
    | Fiber size factor     | 0.9 | 0.9 | 0.3 |
    | Matrix mesh size      | 1.  | 1.  | 1.  |
    | Matrix size factor    | 0.9 | 1.  | 1.  |

- **`run_simulation_from_terminal.py`**  
  Script to **run the Abaqus simulation (`calculate_2D_structure_stiffness.py`) directly from the terminal**.

- **`Compaction_Class.py`**  
  Library of functions used to define and manage simulation objects, boundary conditions, and data extraction.

- **`save_odb_data_in_npz.py`**  
  Script to be **run from Abaqus**. Extracts simulation results from the `.odb` file and saves them in the corresponding tissue folder (created automatically if it does not exist) in NumPy `.npz` format. The tissue to process must be specified by the user (`EW3`, `TW3`, or `LW3`).

- **`convert_npz_data_to_csv.py`**  
  Converts the saved `.npz` data files into `.csv` format for post-processing and visualization. The tissue to process must be specified by the user (`EW3`, `TW3`, or `LW3`).
  
- **'Calculate_engineering_constants.ipynb'**
  Reads the stored stiffness coefficients, converts them into engineering constants, and saves the result. The tissue to process must be specified by the user.
   
## NOTE:
In Abaqus, the following units are adopted and used in the data stored in the CSV files:  
- Length, area, and volume: µm, µm², µm³  
- Time: s  
- Pressure and stress: GPa  
- Force: mN  
- Mass: 10³ kg = ton  
- Energy: nJ


