# Simulation of Self-Densification (SD)

This folder contains all scripts and files used for the **simulation of the self-densification** of wood through shrinking hydrogel fillings that collapse tracheids.  
The phenomenon is modeled using *fluid cavity conditions* with a **compressible fluid**, while all external surfaces are free to move.

The global coordinate system corresponds to the anatomical directions of wood:  
**x** → tangential, **y** → radial, **z** → longitudinal.

## Folder structure

- **`Structure_mat`**  
  Contains:  
  - all `CAE models` (and corresponding `.jnl` files) **generated from the `.mat` files** also included in this folder,  
  - the `.mat` files defining the coordinates of the points for each fiber and the outer contour of the matrix (previously **produced via MATLAB**), and  
  - the `.pkl` files containing the **pre-calculated material parameters** for matrix and fiber at various moisture contents and lignin reduction factors.  

  The available tissue structures are one for each type: EW (`testspruce_EW3`), TW (`testspruce_TW3`), and LW (`testspruce_LW3`).

  The CAE files containing only the geometry are typically **generated and saved before running the simulations**, as their creation is time-consuming. These pre-generated CAE models are then opened and used to perform the simulations.

  There are also CAE files ending in `_full` that contain the complete model ready to run.

- **`SD_simulator.py`**  
  Main script to be **run from Abaqus**. Performs the self-densification simulations. The tissue to be simulated must be specified by the user (`testspruce_EW3`, `testspruce_TW3`, or `testspruce_LW3`), along with:
  - *moisture content* from 0 to 0.3 (array)
  - *lignin reduction factor* from 0 (excluded) to 1 (array)
  - *mesh parameters* (global size and size factors for both fiber and matrix).  
    The values used are:
    | Parameter             | EW  | TW  | LW  |
    |-----------------------|-----|-----|-----|
    | Fiber mesh size       | 2.  | 2.5 | 4.  |
    | Fiber size factor     | 1.  | 0.8 | 1.  |
    | Matrix mesh size      | 1.  | 1.5 | 1.  |
    | Matrix size factor    | 0.9 | 0.8 | 1.  |

  - *fluid properties* (density, expansion coefficient, bulk modulus)
  - *temperature variation* (array). The values used are 15 for EW, 36 for TW, and 84 for LW.

- **`run_simulation_from_terminal.py`**  
  Script to **run the Abaqus simulation (`SD_simulator.py`) directly from the terminal**.

- **`Compaction_Class.py`**  
  Library of functions used to define and manage simulation objects, boundary conditions, and data extraction.

- **`save_odb_data_in_npz.py`**  
  Script to be **run from Abaqus**. Extracts simulation results from the `.odb` file and saves them in the corresponding tissue folder (created automatically if it does not exist) in NumPy `.npz` format. The tissue to process must be specified by the user (`EW3`, `TW3`, or `LW3`).

- **`convert_npz_data_to_csv.py`**  
  Converts the saved `.npz` data files into `.csv` format for post-processing and visualization. The tissue to process must be specified by the user (`EW3`, `TW3`, or `LW3`).
  
## NOTE: 
In Abaqus, the following units are adopted and used in the data stored in the CSV files:  
- Length, area, and volume: µm, µm², µm³  
- Time: s  
- Pressure and stress: GPa  
- Force: mN  
- Mass: 10³ kg = ton  
- Energy: nJ

The names of the result files (`.odb`, `.npz`, `.csv`) indicate the tissue type, the densification protocol, the corresponding moisture content and lignin in percentage, and the temperature variation.
