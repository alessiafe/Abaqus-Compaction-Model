# Simulation of Shear-Assisted Mechanical Densification (SAMD)

This folder contains all scripts and files used for the **simulation of mechanical densification** of wood in the radial direction with superimposed **transverse oscillatory excitation**.  
The phenomenon is modeled by applying opposite quasi-static vertical displacements to reference points kinematically coupled to the top and bottom plates in contact with the tissue structure. Additionally, a harmonic horizontal oscillatory excitation is applied.

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

- **`SAMD_simulator.py`**  
  Main script to be **run from Abaqus**. Performs the mechanical densification simulations. The tissue to be simulated must be specified by the user (`testspruce_EW3`, `testspruce_TW3`, or `testspruce_LW3`), along with:
  - *moisture content* from 0 to 0.3
  - *lignin reduction factor* from 0 (excluded) to 1
  - *mesh parameters* (global size and size factors for both fiber and matrix).  
    The values used are:
    | Parameter             | EW  |
    |-----------------------|-----|
    | Fiber mesh size       | 2.5 |
    | Fiber size factor     | 1.  |
    | Matrix mesh size      | 1.  |
    | Matrix size factor    | 0.9 |

  - *amplitude* and *frequency* in Hz (arrays) of the harmonic oscillation. The maximum values that produce results with EW3 are 10 µm and 10 Hz, respectively.

- **`run_simulation_from_terminal.py`**  
  Script to **run the Abaqus simulation (`SAMD_simulator.py`) directly from the terminal**.

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

The names of the result files (`.odb`, `.npz`, `.csv`) indicate the tissue type, the densification protocol, and the frequency and amplitude (multiplied by 10 to avoid decimals). Moisture content and lignin are not included, as they are fixed values and should be known in advance.
