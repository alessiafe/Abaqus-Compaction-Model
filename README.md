# Simulations of wood densification

This folder contains all scripts and data used for the **multiscale modeling of wood tissues**, which serves as the foundation for simulating various **densification protocols** and also calculating 2D and 3D stiffness.

## Folder structure

- **`Create_Materials_Files`**  
  Contains scripts for **micromaterial calculations of the equivalent single-layer cell wall model**, composed of an isotropic compound middle lamella (CML), combining the middle lamella and primary wall, and an equivalent single-layer representing all S layers (S1, S2, S3). This calculation is performed in Python, outside the Abaqus simulations, to reduce approximation errors in the matrix computations.

- **`Create_Tissue_Structures`**  
  Contains scripts for **generating cellular structures** from microscopic wood images from microscopic images of wood sections and use the geometry datasets to generate a model.

- **`Calculate_Stiffness_2D`**  
  Contains scripts for the **calculation of the equivalent stiffness of 2D tissue models**, obtained by applying four elementary strain cases: two uniaxial, one biaxial, and one shear.

- **`Calculate_Stiffness_3D`**  
  Contains scripts for the **calculation of the equivalent stiffness of 3D tissue models**, obtained by applying 8 elementary strain cases: two biaxial, three triaxial, and three shear.

- **`MD`**  
  Contains scripts for the **simulation of mechanical densification** of wood in the radial direction.

- **`SAMD`**  
  Contains scripts for the **simulation of mechanical densification** of wood in the radial direction with superimposed **transverse oscillatory excitation**.

- **`SD`**  
  Contains scripts for the **simulation of the self-densification** of wood through shrinking hydrogel fillings that collapse tracheids.

- **`Supplementary_Files`** 
  Contains two supplementary scripts: one to calculate the **specific work** of densification and the other to calculate the **engineering constants** of the S2 layer.

The simulations in **`Calculate_Stiffness_2D`**, **`Calculate_Stiffness_3D`**, **`MD`**, **`SAMD`**, and **`SD`** can be carried out independently, as each folder contains the necessary scripts and files. However, they can only be run after generating the material files (**`Create_Materials_Files`**) and the geometry data (**`Create_Tissue_Structures`**), which are already created and present in each folder in the current state.

### Note:
Each folder contains its own set of functions, all of which are called in the same way through `Compaction_Class.py`. However, it is recommended **not to mix these files**, as there may be differences between them. Each set is tailored to a specific densification protocol, and while they may share the same function names, they may not have been updated for the boundary conditions of other protocols. Regarding stiffness calculations, `Compaction_Class.py` includes specific functions dedicated to this task.
