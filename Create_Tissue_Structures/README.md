# Generation of cellular structure of tissues

This folder contains all scripts and files used for the **generation of cellular structures** from microscopic images of wood sections:

- **`testspruce_EW3.bmp`**, **`testspruce_TW3.bmp`**, **`testspruce_LW3.bmp`**  
  Images in `.bmp` format obtained from microscopic images of wood sections.

- **`Create_tissue_struct_EW3.m`**, **`Create_tissue_struct_TW3.m`**, **`Create_tissue_struct_LW3.m`**  
  Main scripts to create a `.mat` file for a tissue based on the pre-generated `.bmp` image (user input). The `.mat` file contains the coordinates of all fibers (outer and inner splines) and the outer matrix, which can be used to generate a model.  
  Each tissue has its own script, as some steps are customized:  
  - Specific manipulation to remove intersections in some specific splines  
  - Scaling factor to exclude CML from the rest of the cell wall  
  - Coordinates of the points used to flatten the top and bottom of the mask of the matrix  

  The following functions are used:
  - **`scale_spline.m`** to scale a spline by a given thickness (outward -, inward +), keeping the same centroid, so that the new spline lies at a constant distance from the original.  
  - **`calculate_centroid.m`** to calculate the centroid coordinates of a spline.  
  - **`close_spline.m`** to remove duplicate points and close a spline (first point = last point).  
  - **`resample_spline.m`** to re-sample the spline to increase point density for a more uniform distribution. The resolution of the output spline is controlled in two steps: first, the spline is evaluated at a high number of points (2000 in the current setup) to get a smooth curve; then, a subset of points (currently 500) are then selected.

- **`Calculate_mean_thickness_and_width.m`**  
  Calculates the mean **cell wall thickness** (S layers excluding CML), **S2 layer thickness**, and **radial and tangential width** of all tracheids for a given tissue (`testspruce_EW3`, `testspruce_TW3`, `testspruce_LW3`).

- **`Calculate_mean_tangential_size.m`**  
  Calculates the mean **tangential size** across all tissues.

- **`Read_wall_angle.m`**  
  Opens a given `.bmp` image and allows the user to draw lines by clicking two points. It calculates the angle of each drawn line, collects all values until the figure is closed, and returns the average angle.
