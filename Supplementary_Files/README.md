# Supplementary files

This folder contains the following supplementary scripts:

- **`add_specific_work_to_csv.py`**  
  Loads all CSV files in a user-specified folder that contain results of densification simulations, calculates the **specific work** in nJ/µm³, and adds it as a new column in each file.

- **`calculate_S2_prop.py`**  
  Calculates the **engineering constants** in GPa of the S2 layer for a given array of *moisture content* (0 to 0.3) and *lignin reduction factor* (0 excluded to 1).  
  **Note:** To run this script, the class `Micromaterial_Calculation` is required, along with the files `basic_components.npy` and `materials_composite.npy`, as provided in the folder **Create_Materials File**.
