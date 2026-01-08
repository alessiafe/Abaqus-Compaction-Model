from abaqus import *
from abaqusConstants import *
from caeModules import *
from part import *
from mesh import *
from visualization import *
from copy import copy
from scipy.io import loadmat
import numpy as np
import time
from re import *
from mpmath import *
mp.dps = 150 # increase precision to stabilize matrix operations (default = 15)

class Compaction_Class:

    def __init__(self):
    # Constructor of Compaction_Class
    # Args: none

        # Model name       
        self.modelName = 'Model-1'

        # Geometry parameters initialization
        self.lthick = 1. # longitudinal thickness
        self.matrixName = 'Matrix' # matrix part name
        self.fibersName = [] # names of each fiber part
        self.fibersTotName  = 'Fibers' # name of merged fibers part
        self.tissueName = 'Tissue' # name of tissue part (= matrix + fibers)
        self.fibersSetsName = [] # names of fibers set
        self.lumenSurfName = [] # names of fibers inner surfaces
        self.fibersExtSurfName = [] # names of fibers outer surfaces
        self.fibersTotExtSurfName = ['ExtSurf-Fibers', 'LeftSurf-Fibers', 'RightSurf-Fibers'] # names of merged fibers external surfaces
        self.fibersTotInnerSurfName = 'InnerSurf-Fibers' # name of fibers total inner surface
        self.matrixInnerSurfName = 'InnerSurf-Matrix' # name of matrix inner surface
        self.matrixExtSurfName = ['ExtSurf-Matrix', 'BottomSurf-Matrix',
                                  'TopSurf-Matrix', 'LeftSurf-Matrix', 'RightSurf-Matrix'] # names of matrix external surfaces
        self.generalContactSurfName = ['ExtSurf-Tissue', 'LateralSurf-Tissue']#'LeftSurf-Tissue', 'RightSurf-Tissue'] # names of lateral surfaces for general contact interaction
        
        # Material parameters initialization
        self.matrix_eng_const = [] # array of matrix engineering constants [E1, E2, E3, nu23, nu13, nu12, G23, G13, G12]
        self.fiber_eng_const = [] # array of fiber engineering constants [E1, E2, E3, nu23, nu13, nu12, G23, G13, G12]
        self.matrix_density = [] # matrix density [ton/microm3]
        self.fiber_density = [] # fibers density [ton/microm3]

    	# Element features initialization
        self.meshSize = 1. # global mesh size
        self.minMeshSize = 0.1 # minimum size factor

        # Output request
        self.addRPsRequest = False # request displacement and force output (only for RPs compation control)

        return
    # end: _init_

    def totuple(self, a, scale=1):
    # Turn point list in tuple of tuples
    # Args:
    # - a = point list
    # - scale = scale factor (=1 by default)
    # Returns:
    # - a = tuple of tuples
        try:
            aL=[]
            for i in range(a[0].size):
                aL.append((a[0,i]*scale,a[1,i]*scale))
            return tuple(aL)
        except TypeError:
            return a
    # end: totuple

    def setup_geometry(self, structpath, scalefib=1):
    # Setup geometry of the tissue structure: create parts, sets and surfaces
    # Args:
    # - structpath = path of matlab file containing the tissue structure data
    # - scalefib = scale factor of fibers
    # Returns:
    # - modelName = name of the model

        self.mdb = Mdb() # start model
        myModel=self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        ############### LOAD DATA ###############
        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab structure
        self.mstruct = data['mstruct'] # get data

        ############### CREATE FIBERS ###############
        start = time.time() # start timer
        print("Creating fibers...\n")
        fiber_out = [] # initialize outer points array
        fiber_in = [] # initialize inner points array
        for i in range(self.mstruct.fibers.size): # repeat for each fiber
            print("    {}/{}\n".format(i+1,self.mstruct.fibers.size))
            # Get fiber data
            f1 = self.mstruct.fibers[i] # get fiber data
            a_in = f1.innerpts # internal points
            a_out = f1.outerpts # external points
            fiber_out.append(a_out) # append external points
            fiber_in.append(a_in) # append inner points
            s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
            s.setPrimaryObject(option=STANDALONE) # set sketch object in viewport
            a_out = self.totuple(a_out)
            a_in = self.totuple(a_in)
            s.Spline(points=a_out) # draw external fiber
            s.Spline(points=a_in) # draf inner fiber
            
            self.fibersName.append('fiber-'+str(i+1)) # append fiber name
            myPart = myModel.Part(name=self.fibersName[-1], dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY) # create part
            myPart.BaseShell(sketch=s) # extrude fiber
            s.unsetPrimaryObject() # unset sketch object in viewport
            del myModel.sketches['__profile__'] # close sketch

            # Create set of fiber
            self.fibersSetsName.append('Set-'+self.fibersName[-1]) # append set name
            myPart = myModel.parts[self.fibersName[-1]]
            myPart.Set(faces=myPart.faces, name=self.fibersSetsName[-1]) # create set
            
            # Create inner surface of fiber
            #point_in = self.totuple(a_in) # tuple of inner points coordinates
            surface = myPart.edges.findAt(((a_in[0][0], a_in[0][1], 0.), )) #[myPart.edges.findAt(((point[0], point[1], 0.), )) for point in a_in]  # find inner face
            self.lumenSurfName.append('Lumen-'+self.fibersName[-1]) # append surface name
            myPart.Surface(side1Edges=surface, name=self.lumenSurfName[-1]) # create surface
            myPart.Set(edges=surface, name=self.lumenSurfName[-1]) # create set
            
            # Create outer surface of fiber
            #point_out = self.totuple(a_out) # tuple of outer points coordinates
            surface = myPart.edges.findAt(((a_out[0][0], a_out[0][1], 0.), )) #[myPart.edges.findAt(((point[0], point[1], 0.), )) for point in a_out] # find inner face
            self.fibersExtSurfName.append('ExtFiber-'+self.fibersName[-1]) # append surface name
            myPart.Surface(side1Edges=surface, name=self.fibersExtSurfName[-1]) # create surface
            myPart.Set(edges=surface, name=self.fibersExtSurfName[-1]) # create set
            
            # Scale fibers
            if scalefib != 1:    
                s = myPart.features['Shell planar-1'].sketch
                myModel.ConstrainedSketch(name='__edit__', objectToCopy=s)
                s1 = myModel.sketches['__edit__']
                s1.setPrimaryObject(option=SUPERIMPOSE)
                # External points
                fiber_temp = self.totuple(fiber_out[i]) # get external coordinates
                fiber_x = [] # x-coordinates of fibers
                fiber_y = [] # y-coordinates of fibers
                for k in range(len(fiber_temp)): # for each point
                    out_temp = fiber_temp[k] # point coordinates
                    fiber_x.append(out_temp[0]) # append x
                    fiber_y.append(out_temp[1]) # append y
                g = s1.geometry
                myPart.projectReferencesOntoSketch(sketch=s1, 
                    upToFeature=myPart.features['Shell planar-1'], filter=COPLANAR_EDGES)
                s1.scale(scaleValue=scalefib, scaleCenter=(min(fiber_x)+(max(fiber_x)-min(fiber_x))/2., min(fiber_y)+(max(fiber_y)-min(fiber_y))/2.),
                    objectList=(g[2],))
                # Inner points
                fiber_temp = self.totuple(fiber_in[i]) # get external coordinates
                fiber_x = [] # x-coordinates of fibers
                fiber_y = [] # y-coordinates of fibers
                for k in range(len(fiber_temp)): # for each point
                    in_temp = fiber_temp[k] # point coordinates
                    fiber_x.append(in_temp[0]) # append x
                    fiber_y.append(in_temp[1]) # append y
                s1.scale(scaleValue=scalefib, scaleCenter=(min(fiber_x)+(max(fiber_x)-min(fiber_x))/2., min(fiber_y)+(max(fiber_y)-min(fiber_y))/2.),
                    objectList=(g[3],))
                s1.unsetPrimaryObject()
                myPart.features['Shell planar-1'].setValues(sketch=s1)
                del myModel.sketches['__edit__']
                myPart.regenerate()
            
            # Create fiber instances
            myAssembly.Instance(name=self.fibersName[-1], part=myPart, dependent=ON)
        
        # Merge fibers in one single part
        myAssembly.InstanceFromBooleanMerge(name=self.fibersTotName , instances=([myAssembly.instances[self.fibersName[i]]
            for i in range(len(self.fibersName))] ), keepIntersections=ON, domain=GEOMETRY, originalInstances=DELETE)
        myAssembly.features.changeKey(fromName=self.fibersTotName +'-1', toName=self.fibersTotName )
        # Merge outer surfaces of fibers
        myPart = myModel.parts[self.fibersTotName ]
        myPart.SurfaceByBoolean(name=self.fibersTotExtSurfName[0], surfaces=([myPart.surfaces[self.fibersExtSurfName[i]]
            for i in range(len(self.fibersExtSurfName))]))
        # Merge inner surfaces of fibers
        myPart.SurfaceByBoolean(name=self.fibersTotInnerSurfName, surfaces=([myPart.surfaces[self.lumenSurfName[i]]
            for i in range(len(self.lumenSurfName))]))
        # Create fibers assembly set 
        myModel.parts[self.fibersTotName ].Set(faces=myModel.parts[self.fibersTotName ].faces, name='Set-'+self.fibersTotName) # create set
        
        ############### CREATE MATRIX ###############
        print("Creating matrix...\n")
        mask = self.mstruct.mask # get mask data
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
        s.setPrimaryObject(option=STANDALONE) # set sketch objetc in viewport
        mask = self.totuple(mask)
        s.Spline(points=mask) # draw external fiber
        #s.rectangle(point1=(min(fiber_x)-deltax,min(fiber_y)-deltay), point2=(max(fiber_x)+deltax,max(fiber_y)+deltay)) # draw rectangle
        for i in range(self.mstruct.fibers.size): # for each fiber 
            f2 = self.mstruct.fibers[i] # get fiber data
            a_in = self.totuple(f2.innerpts) # inner points
            s.Spline(points=a_in) # mask inner points
        myPart = myModel.Part(name=self.matrixName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY) # create part
        myPart.BaseShell(sketch=s) # extrude mask
        s.unsetPrimaryObject() # unset sketch objetc in viewport
        del myModel.sketches['__profile__'] # close sketch
        
        # Create mask assembly
        myPart = myModel.parts[self.matrixName]
        myAssembly.Instance(name=self.matrixName, part=myPart, dependent=ON) # create instance

        # Create outer surface
        myPart = myModel.parts[self.matrixName]
        surface = myPart.edges #myPart.edges.findAt(((mask[0][0], mask[1][0], 0.), ))
        myPart.Surface(side1Edges=surface, name=self.matrixExtSurfName[0]) # create outer surface
        #
        # Cut fibers from mask
        tempName = 'mask-temp'
        myAssembly.InstanceFromBooleanCut(cuttingInstances=(myAssembly.instances[self.fibersTotName ], ), 
            instanceToBeCut=myAssembly.instances[self.matrixName], 
            name=tempName, originalInstances=SUPPRESS)
        myAssembly.features[self.fibersTotName ].resume()
        del myModel.parts[self.matrixName]
        del myAssembly.features[self.matrixName]
        myModel.parts.changeKey(fromName=tempName, toName=self.matrixName)
        myAssembly.features.changeKey(fromName=tempName+'-1', toName=self.matrixName)
     
        myAssembly.regenerate()
        # Create inner surface of matrix
        myPart = myModel.parts[self.matrixName]
        surface = myPart.edges # find inner face
        myPart.Surface(side1Edges=surface, name=self.matrixInnerSurfName) # create surface
        myPart.SurfaceByBoolean(name=self.matrixInnerSurfName, operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixInnerSurfName], myPart.surfaces[self.matrixExtSurfName[0]], ))
        myPart.Set(edges=surface, name=self.matrixInnerSurfName) # create set
        # Create matrix set
        myModel.parts[self.matrixName].Set(faces=myModel.parts[self.matrixName].faces, name='Set-'+self.matrixName) # create mask set
        
        myAssembly.regenerate() # regenerate assembly
        #
        end = time.time() # stop timer
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))

        return self.modelName
    # end: setup_geometry

    def create_material(self, materialname, elasticprop, density):
    # Create material
    # Args:
    # - materialname = material name
    # - elasticprop = array of elastic engineering constants (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23) [GPa]
    # - density = material density [ton/microm3]

        myModel = self.mdb.models[self.modelName] # model
        myModel.Material(name=materialname) # create material
        myMaterial = myModel.materials[materialname]
        #myModel.materials[materialname].Elastic(type=ENGINEERING_CONSTANTS, table=(elasticprop, )) # set engineering constants
        myMaterial.Elastic(type=ORTHOTROPIC, table=((elasticprop[2], elasticprop[5], elasticprop[1],
                                                     elasticprop[4], elasticprop[3], elasticprop[0],
                                                     elasticprop[6], elasticprop[7], elasticprop[8]),))       
        myModel.materials[materialname].Density(table=((density, ), )) # set density

        return
    # end: create_material

    def assign_section(self, partname, setname, sectionname, materialname):
    # Create and assign section
    # Args:
    # - partname = part name
    # - sectionname = name of section to create
    # - materialname = name of material to assign

        myModel = self.mdb.models[self.modelName] # model
        myPart = myModel.parts[partname] # part
        mySet = myPart.sets['Set-'+setname] # part set
        # Create section
        myModel.HomogeneousSolidSection(name=sectionname, material=materialname, thickness=self.lthick)
        # Assign section to part
        myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
            region=mySet, sectionName=sectionname, thicknessAssignment=FROM_SECTION)

        return
    # end: assignMatrixSection
 
    def setup_materials(self, matpath, moisture, matrix_dens, fiber_dens, ligred=1., lthick=1., scale=1e-21):
    # Setup material properties
    # Args:
    # - structpath = path of matlab file containing the tissue structure data
    # - matrix_elprop = array of matrix engineering constants (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23) [GPa]
    # - fiber_elprop = array of fiber engineering constants (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23) [GPa]
    # - matrix_dens = matrix density [kg/m3]
    # - fiber_dens = fiber density [kg/m3]
    # - alpha = array of ply microfibril orientations (MFA) (=[0., 0., 60., -60., 15., 75.])
    # - lthick = longitudinal thickness of tissue [microm3]
    # - scale = density scale factor for unit consistency (=1e-21 by default: kg/m3 -> ton/microm3)

        self.lthick = copy(lthick) # store longitudinal thickness
        self.matrix_density = matrix_dens*scale # matrix density [ton/microm3]
        self.fiber_density = fiber_dens*scale # fibers density [ton/microm3]
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        print("Setting materials...\n")
        start = time.time() # start timer

        ### CALCULATE LOCAL STIFFNESS OF EACH CELL WALL LAYER
        # Micromaterial calc. = 
        #  - 1-direction is the longitudinal direction
        #  - 2-direction is tangential to the cell wall
        #  - 3-direction is normal to the cell wall
        # Cell wall coord. sys. in Abaqus =
        #  - 1-direction is normal to the cell wall
        #  - 2-direction is tangential to the cell wall
        #  - 3-direction is the longitudinal direction
        
        # Import materials data
        import pickle
        with open(matpath, 'rb') as f:
            loaded_data = pickle.load(f)
        #loaded_data = np.load(matpath, allow_pickle=True)
        key = str(int(moisture*100))+ '-'+ str(int(ligred*100))
        C = loaded_data[key]

        ### SET MATRIX MATERIAL
        self.matrix_stiff_coeff = C[0].values()
        materialname = self.matrixName # material name
        sectionname = 'Section-'+materialname # section name
        self.create_material(materialname, self.matrix_stiff_coeff, self.matrix_density) # create material
        self.assign_section(self.matrixName, self.matrixName, sectionname, materialname) # create and assign section
        myPart = myModel.parts[self.matrixName] # matrix part
        myRegion = myPart.sets['Set-'+self.matrixName] # matrix region
        # Set discrete material orientation = global coord. sys.
        myPart.MaterialOrientation(region=myRegion, 
            orientationType=GLOBAL, axis=AXIS_1, additionalRotationType=ROTATION_NONE, 
            localCsys=None, fieldName='', stackDirection=STACK_3)
        
        ### SET FIBERS MATERIAL
        myPart = myModel.parts[self.fibersTotName] # merged fibers part
        for i in range(len(self.fibersName)): # repeat for each fiber
            # Get stiffness coefficients
            self.fiber_stiff_coeff = C[i+1].values()
            materialname = self.fibersName[i] # material name
            sectionname = 'Section-'+materialname # section name
            self.create_material(materialname, self.fiber_stiff_coeff, self.fiber_density) # create material
            self.assign_section(self.fibersTotName, self.fibersName[i], sectionname, materialname) # create and assign section
            # Set discrete material orientation
            myRegion = myPart.sets[self.fibersSetsName[i]] # fiber region (set)
            primaryAxisRegion = myPart.sets[self.lumenSurfName[i]] # primary axis normal surface           
            myPart.MaterialOrientation(region=myRegion, 
                orientationType=DISCRETE, axis=AXIS_3, normalAxisDefinition=VECTOR, 
                normalAxisVector=(0.0, 0.0, 1.0), flipNormalDirection=False, 
                normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
                primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_2, 
                flipPrimaryDirection=False, additionalRotationType=ROTATION_NONE, 
                angle=0.0, additionalRotationField='', stackDirection=STACK_3)

        end = time.time() # stop timer
        # Print elapsed time to create and set materials
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))

        myAssembly.regenerate() # regenerate assembly
        
        return
    # end: setup_materials_new

    def create_mesh(self, meshSizeFiber=2., meshSizeMatrix=1., minSizeFiber=0.5, minSizeMatrix=0.9, lengthRatio=0.1):
    # Create part mesh
    # Args:
    # - meshSize = global element size of fibers mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the specified global element size
    # - scale = meshSize scale factor to get global element size of matrix mesh (1/5 by default)

        self.meshSizeFiber = meshSizeFiber # global element size of the mesh
        self.meshSizeMatrix = meshSizeMatrix # global element size of the mesh
        self.minSizeFiber = minSizeFiber # fraction of minimum mesh size
        self.minSizeMatrix = minSizeMatrix # fraction of minimum mesh size
        #self.minMeshSize = self.meshSize * self.minSizeFactor # minimum mesh size
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        print("Create mesh...\n")
        start = time.time() # start timer
        # Element type
        elemType1 = ElemType(elemCode=UNKNOWN_QUAD, elemLibrary=EXPLICIT)
        elemType2 = ElemType(elemCode=CPE6M, elemLibrary=EXPLICIT,
            secondOrderAccuracy=ON, hourglassControl=DEFAULT, distortionControl=ON, lengthRatio=0.1)
        
        ## SET FIBERS MESH
        myPart = myModel.parts[self.fibersTotName] # merged fibers part
        # Set mesh control
        myPart.setMeshControls(elemShape=TRI, regions=myPart.faces, allowMapped=True)
        # Set element type
        myRegion = myPart.sets['Set-'+self.fibersTotName ] # fibers region (set)
        myPart.setElementType(elemTypes=(elemType1, elemType2), regions=myRegion)
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeFiber, size=self.meshSizeFiber)
        # Generate mesh
        myPart.generateMesh()

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set mesh control
        myPart.setMeshControls(elemShape=TRI, regions=myPart.faces, allowMapped=True)
        # Set element type
        elemType2 = ElemType(elemCode=CPE6M, elemLibrary=EXPLICIT,
            secondOrderAccuracy=ON, hourglassControl=DEFAULT, distortionControl=ON, lengthRatio=lengthRatio)
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2), regions=myRegion)
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeMatrix, size=self.meshSizeMatrix)
        # Generate mesh
        myPart.generateMesh()

        myAssembly.regenerate() # regenerate assembly

        end = time.time() # stop timer
        # Print elapsed time to mesh the model
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))  

        return
    # end: create_mesh

    def create_single_mesh(self, meshSize=2., minSize=0.5):
    # Create part mesh
    # Args:
    # - meshSize = global element size of fibers mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the specified global element size
    # - scale = meshSize scale factor to get global element size of matrix mesh (1/5 by default)

        self.meshSize = meshSize # global element size of the mesh
        self.minSize = minSize # fraction of minimum mesh size
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        # Element type
        elemType1 = ElemType(elemCode=UNKNOWN_QUAD, elemLibrary=EXPLICIT)
        elemType2 = ElemType(elemCode=CPE6M, elemLibrary=EXPLICIT,
            secondOrderAccuracy=ON, hourglassControl=DEFAULT, distortionControl=ON, lengthRatio=0.1)

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set element type
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2), regions=myRegion)
        # Set mesh control
        myPart.setMeshControls(elemShape=TRI, regions=myPart.faces, allowMapped=True)
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSize, size=self.meshSize)
        # Generate mesh
        myPart.generateMesh()

        myAssembly.regenerate() # regenerate assembly
        
        return
    # end: create_single_mesh

    def find_points_within_tolerance(self, mask, tolerance):
    #    
        mask_array = np.array(mask) # convert mask to npArray
        mean_index = np.mean(mask_array[:, 0]) # calculate mean index for the x-coordinates
        close_points = [] # initialize a list to store points within the specified tolerance
        
        # Iterate through the points in the mask
        for point in mask_array:
            distance = abs(point[0] - mean_index) # calculate the distance from the mean index
            # Check if the distance is within the tolerance
            if distance <= tolerance:
                close_points.append(point) # append the point to the list if it's within tolerance
        
        return np.array(close_points)
    # end: find_points_within_tolerance

    def find_max_min_y(self, close_points):
        # Extract y-coordinates (position j) from close_points
        if close_points.size == 0:
            return None, None  # Handle the case with no close points        
        y_coordinates = close_points[:, 1]
        
        # Find maximum and minimum y values
        max_y = np.max(y_coordinates)
        min_y = np.min(y_coordinates)
        
        return max_y, min_y
    # end: find_max_min_y

    def create_surface_partition(self, meshSize,sides=False):
    #
        myModel = self.mdb.models[self.modelName] # model
        myPart = myModel.parts[self.matrixName] # matrix part
        mask = self.mstruct.mask # get mask data
        mask = self.totuple(mask)
        delta = meshSize

        # Find points on surfaces
        tolerance = 25.  # Define your tolerance
        close_points = self.find_points_within_tolerance(mask, tolerance)
        ymax, ymin = self.find_max_min_y(close_points)
        xmin = min(row[0] for row in mask)
        xmax = max(row[0] for row in mask)

        # Partition top surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymax-delta/3.)
        myPart.PartitionEdgeByDatumPlane(datumPlane=myPart.datums[7], edges=myPart.edges)
        myPart.Surface(name=self.matrixExtSurfName[2], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymax-delta, 0., xmax+delta, ymax+delta))
        #myPart.Set(name=self.matrixExtSurfName[2], edges=myPart.edges.getByBoundingBox(xmin-delta, ymax-delta, 0., xmax+delta, ymax+delta))
        # Partition bottom surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymin+delta/3.)
        myPart.PartitionEdgeByDatumPlane(datumPlane=myPart.datums[10], edges=myPart.edges)
        myPart.Surface(name=self.matrixExtSurfName[1], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymin-delta, 0., xmax+delta, ymin+delta))
        #myPart.Set(name=self.matrixExtSurfName[1], edges=myPart.edges.getByBoundingBox(xmin-delta, ymin-delta, 0., xmax+delta, ymin+delta))
        
        # Partition sides surfaces
        if sides:
            surfName = 'SidesSurf-Matrix'
            myPart.SurfaceByBoolean(name=surfName, operation=DIFFERENCE, surfaces=(myPart.surfaces[self.matrixExtSurfName[0]], 
                myPart.surfaces[self.matrixExtSurfName[1]],  myPart.surfaces[self.matrixExtSurfName[2]], ))        
            faces = myPart.surfaces[surfName].edges
            xmax = max(mask, key=lambda x: x[0])[0]
            xmin = min(mask, key=lambda x: x[0])[0]
            leftSide = ()
            rightSide = ()
            for i in range(len(faces)):
                if faces[i].pointOn[0][0] < (xmin+xmax)/2.:
                    leftSide += faces[i],
                else:
                    rightSide += faces[i],
            myPart.Surface(name=self.matrixExtSurfName[4], side1Edges=EdgeArray(rightSide))
            myPart.Surface(name=self.matrixExtSurfName[3], side1Edges=EdgeArray(leftSide))
     
        return
    # end:create_surface_partition

    def create_surfaces(self):
    #
        myModel = self.mdb.models[self.modelName] # model
        myPart = myModel.parts[self.matrixName] # matrix part
        mask = self.mstruct.mask # get mask data
        mask = self.totuple(mask)
        delta = 1e-3
        depth = 0.

        # Find points on surfaces
        ymax = max(row[1] for row in mask)
        ymin = min(row[1] for row in mask)
        xmin = min(row[0] for row in mask)
        xmax = max(row[0] for row in mask)

        # Top surface
        myPart.Surface(name=self.matrixExtSurfName[2], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymax-delta, 0.-delta, xmax+delta, ymax+delta, depth+delta))
        # Bottom surface
        myPart.Surface(name=self.matrixExtSurfName[1], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymin-delta, 0.-delta, xmax+delta, ymin+delta, depth+delta))
        # Left surface
        myPart.Surface(name=self.matrixExtSurfName[3], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymin-delta, 0.-delta, xmin+delta, ymax+delta, depth+delta))
        # Right surface
        myPart.Surface(name=self.matrixExtSurfName[4], side1Edges=myPart.edges.getByBoundingBox(xmax-delta, ymin-delta, 0.-delta, xmax+delta, ymax+delta, depth+delta))

        return
    # end:create_surfaces

    def create_model_copy(self, numCopy):
    # Create copies of the first model
    # Args:
    # - numCopy = number of the copies of the model to create
    # Returns:
    # - modelName = array of model names
        
        if numCopy==0:
            # if 0 has been passed as number of copies, warn the user and return a 1-element array with the name of the first and only model
            modelName = [self.modelName] # array of model name
            print('You passed 0 as number of model copies to create. No copy has been created.')
        else:
            newModelName = ['']*numCopy # initialize array of model names
            for k in range(numCopy): # repeat numCopy times
                newModelName[k] = 'Model-' + str(k+2) # new model name
                self.mdb.Model(name=newModelName[k], objectToCopy=self.mdb.models[self.modelName]) # copy model
            modelName = [self.modelName] + newModelName # append model name to array

        return modelName
    # end: create_model_copy

    def submit_job(self, jobName, modelName):
    # Create and submit job
    # Args:
    # - jobName = name of the job
    # - modelName = name of the model

        import os
        fileExt = '.lck'
        filePath = jobName + fileExt
        if os.path.exists(filePath):
            os.remove(filePath)
        
        self.mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
            resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=1, 
            activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1) # create job

        myJob = self.mdb.jobs[jobName] # job
        start = time.time() # start timer
        myJob.submit(consistencyChecking=OFF) # submit job
        myJob.waitForCompletion() # wait for job completion

        end = time.time() # stop timer
        # Print elapsed time to create and set materials
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))
        
        return
    # end: submit_job

    def save_model(self, savepath):
    # Save .cae file
    # Args:
    # - savepath = file path to save the .cae file
    # Returns:
    # - path = full path to the file

        self.mdb.saveAs(pathName=savepath) # save file
    
        return
    # end: save_model

    def open_model(self, modelpath, structpath):
    # Open .cae file with structure and materials already created.
    # Args:
    # - savepath = folder path where the .cae file is saved
    # - filename = name of the .cae file to open
    # Returns:
    # - path = full path to the file

        self.mdb = openMdb(pathName=modelpath) # open file
        self.modelName = mdb.models.keys()[0]
        myModel=self.mdb.models[self.modelName]
        if myModel.sections.values():
            self.lthick = myModel.sections.values()[0].thickness
        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab data structure
        self.mstruct = data['mstruct'] # get points data
        
        self.fibersName = []
        self.fibersSetsName = []
        self.lumenSurfName = []
        self.fibersExtSurfName = []
        # Initialize name arrays (setup_geometry)
        for i in range(len(myModel.parts)-2): # repeat for each fiber
            self.fibersName.append('fiber-'+str(i+1)) # append fiber name
            self.fibersSetsName.append('Set-'+self.fibersName[-1]) # append set name
            self.lumenSurfName.append('Lumen-'+self.fibersName[-1]) # append surface name
            self.fibersExtSurfName.append('ExtFiber-'+self.fibersName[-1]) # append surface name

        return self.modelName
    # end: open_model

    def create_mesh_std(self, meshSizeFiber=2., meshSizeMatrix=1., minSizeFiber=0.5, minSizeMatrix=0.9, lengthRatio=0.1):
    # Create part mesh
    # Args:
    # - meshSize = global element size of fibers mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the specified global element size
    # - scale = meshSize scale factor to get global element size of matrix mesh (1/5 by default)

        self.meshSizeFiber = meshSizeFiber # global element size of the mesh
        self.meshSizeMatrix = meshSizeMatrix # global element size of the mesh
        self.minSizeFiber = minSizeFiber # fraction of minimum mesh size
        self.minSizeMatrix = minSizeMatrix # fraction of minimum mesh size
        #self.minMeshSize = self.meshSize * self.minSizeFactor # minimum mesh size
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        print("Create mesh...\n")
        start = time.time() # start timer
        # Element type
        elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD, 
            secondOrderAccuracy=ON, distortionControl=DEFAULT)
                
        ## SET FIBERS MESH
        myPart = myModel.parts[self.fibersTotName] # merged fibers part
        # Set mesh control
        myPart.setMeshControls(elemShape=TRI, regions=myPart.faces, allowMapped=True)
        # Set element type
        myRegion = myPart.sets['Set-'+self.fibersTotName ] # fibers region (set)
        myPart.setElementType(elemTypes=(elemType1, elemType2), regions=myRegion)
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeFiber, size=self.meshSizeFiber)
        # Generate mesh
        myPart.generateMesh()

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set mesh control
        myPart.setMeshControls(elemShape=TRI, regions=myPart.faces, allowMapped=True)
        # Set element type
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2), regions=myRegion)
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeMatrix, size=self.meshSizeMatrix)
        # Generate mesh
        myPart.generateMesh()

        myAssembly.regenerate() # regenerate assembly

        end = time.time() # stop timer
        # Print elapsed time to mesh the model
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))  

        return
    # end: create_mesh_std

    def set_general_constraints(self, modelName, stepName, stepTime=1):
    # Set periodic boundary consitions (PBC)
    # Args:
    # - strains = array of strains to apply
    # - stepName = name of the step
        
        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        myInstance = myAssembly.instances[self.matrixName]
        myPart = myModel.parts[self.matrixName]

        ## CREATE STEP
        self.stepTime = stepTime        
        myModel.StaticStep(name=stepName, previous='Initial', timePeriod=self.stepTime, maxNumInc=10000, initialInc=1e-5, minInc=1e-10, nlgeom=ON)        
        #myModel.ExplicitDynamicsStep(name=stepName, timePeriod=self.stepTime, previous='Initial', improvedDtMethod=ON) # create step
        myStep = myModel.steps[stepName] # step
        #myStep.setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 1e12, 0.0, None, 0, 0, 0.0, 0.0, 0, None), )) # set mass scaling

        ## CREATE INTERACTION CONSTRAINTS
        # Tie matrix (inner) and fibers (outer)
        contIntName = 'FiberMatrix-Tie' # interaction name
        region1 = myAssembly.instances[self.fibersTotName ].surfaces[self.fibersTotExtSurfName[0]] # region 1 = fibers external surface
        region2 = myInstance.surfaces[self.matrixInnerSurfName] # region 2 = matrix inner surface
        myModel.Tie(name=contIntName, master=region1, slave=region2, constraintEnforcement=SURFACE_TO_SURFACE, 
            positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, thickness=ON) # create tie interaction
    
        myAssembly.regenerate() # regenerate assembly

        ## CREATE RPs
        # Find bounding box
        myPart.Set(nodes=myPart.nodes, name=self.matrixName+'-nodes') # create matrix nodes set
        nodeSet = myPart.sets[self.matrixName+'-nodes'] # matrix nodes set
        xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
        yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
        self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord)] # store matrix bounding box coordinates
        self.refPointName = []
        t = 10
        # Create bottom RP
        self.refPointName.append('RP-Bottom') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[2]-t, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create top RP
        self.refPointName.append('RP-Top') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[3]+t, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create left RP
        self.refPointName.append('RP-Left') # append RP name        
        myAssembly.ReferencePoint(point=(self.tissueBox[0]-t, (self.tissueBox[3]+self.tissueBox[2])/2., 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create right RP
        self.refPointName.append('RP-Right') # append RP name        
        myAssembly.ReferencePoint(point=(self.tissueBox[1]+t, (self.tissueBox[3]+self.tissueBox[2])/2., 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
               
        # Create node sets
        setName = ['Bottom-nodes', 'Top-nodes', 'Left-nodes', 'Right-nodes']
        for i in range(len(setName)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName[i], nodes=myPart.surfaces[self.matrixExtSurfName[1+i]].nodes)
        myPart.SetByBoolean(name=setName[2], operation=DIFFERENCE, sets=(
            myPart.sets[setName[2]], myPart.sets[setName[0]], myPart.sets[setName[1]],))
        myPart.SetByBoolean(name=setName[3], operation=DIFFERENCE, sets=(
            myPart.sets[setName[3]], myPart.sets[setName[0]], myPart.sets[setName[1]],))
        
        return
    #end: set_general_constraints

    def set_single_constraints(self, modelName, stepName, stepTime=1):
    # Set periodic boundary consitions (PBC)
    # Args:
    # - strains = array of strains to apply
    # - stepName = name of the step
        
        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        myPart = myModel.parts[self.matrixName]

        ## CREATE STEP
        self.stepTime = stepTime        
        myModel.StaticStep(name=stepName, previous='Initial', timePeriod=self.stepTime, maxNumInc=10000, initialInc=1e-5, minInc=1e-10, nlgeom=ON)        

        ## CREATE RPs
        # Find bounding box
        myPart.Set(nodes=myPart.nodes, name=self.matrixName+'-nodes') # create matrix nodes set
        nodeSet = myPart.sets[self.matrixName+'-nodes'] # matrix nodes set
        xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
        yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
        self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord)] # store matrix bounding box coordinates
        self.refPointName = []
        t = 10
        # Create bottom RP
        self.refPointName.append('RP-Bottom') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[2]-t, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create top RP
        self.refPointName.append('RP-Top') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[3]+t, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create left RP
        self.refPointName.append('RP-Left') # append RP name        
        myAssembly.ReferencePoint(point=(self.tissueBox[0]-t, (self.tissueBox[3]+self.tissueBox[2])/2., 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create right RP
        self.refPointName.append('RP-Right') # append RP name        
        myAssembly.ReferencePoint(point=(self.tissueBox[1]+t, (self.tissueBox[3]+self.tissueBox[2])/2., 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
               
        # Create sides node sets
        self.create_surfaces()
        setName1 = ['Bottom-nodes', 'Top-nodes', 'Left-nodes', 'Right-nodes']
        for i in range(len(setName1)):
            myPart.Set(name=setName1[i], nodes=myPart.surfaces[self.matrixExtSurfName[1+i]].nodes)
        myPart.SetByBoolean(name=setName1[2], operation=DIFFERENCE, sets=(
            myPart.sets[setName1[2]], myPart.sets[setName1[0]], myPart.sets[setName1[1]],))
        myPart.SetByBoolean(name=setName1[3], operation=DIFFERENCE, sets=(
            myPart.sets[setName1[3]], myPart.sets[setName1[0]], myPart.sets[setName1[1]],))
        
        return
    #end: set_general_constraints

    def set_periodic_boundary_conditions(self, modelName, stepName, strain, depth=1.0):
    # Set periodic boundary conditions (PBC)
    # Args:
    # - strains: array of strains to apply
    # - stepName: name of the step

        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        myInstance = myAssembly.instances[self.matrixName]
        
        # Define lattice vectors
        self.LatticeVector = [
            [self.tissueBox[1] - self.tissueBox[0], 0.0, 0.0],
            [0.0, self.tissueBox[3] - self.tissueBox[2], 0.0],
        ]
        
        # Calculate displacements
        # Array of strains
        disp_x = strain[0]*self.LatticeVector[0][0]/2.
        disp_y = strain[1]*self.LatticeVector[1][1]/2.
        #################################################
        disp_xy_y = strain[2]*self.LatticeVector[0][0]/2.
        disp_xy_x = strain[2]*self.LatticeVector[1][1]/2.
       
        # Set displacements and coupling DOF
        if strain[0] != 0 and all(x == 0.0 for x in strain[1:]):  # Case 1: x-strain
            strains = [[UNSET, 0, UNSET], [UNSET, 0, UNSET],      # bottom, top
                       [-disp_x, 0, UNSET], [disp_x, 0, UNSET]]   # left, right
            
        elif strain[1] != 0 and all(x == 0.0 for x in [strain[0]] + strain[2:]):  # Case 2: y-strain
            strains = [[0, -disp_y, UNSET], [0, disp_y, UNSET],   # bottom, top
                       [0, UNSET, UNSET], [0, UNSET, UNSET]]      # left, right
            
        ##############################################################################################################
            
        elif strain[2] != 0 and all(x == 0.0 for x in strain[:2]):  # Case 3: xy-strain
            strains = [[disp_xy_x, UNSET, UNSET], [-disp_xy_x, UNSET, UNSET], # bottom, top
                       [UNSET, disp_xy_y, UNSET], [UNSET, -disp_xy_y, UNSET]] # left, right 

        ##############################################################################################################

        elif strain[0] != 0 and strain[1] != 0 and all(x == 0.0 for x in strain[2:]):  # Case 4: x+y-strain
            strains = [[UNSET, -disp_y, UNSET], [UNSET, disp_y, UNSET],   # bottom, top
                       [-disp_x, UNSET, UNSET], [disp_x, UNSET, UNSET]]   # left, right
            
        else:
            raise ValueError("The strain case you provided is not implemented: ", strain)
        
        coup = [[ON if value != UNSET else OFF for value in row] for row in strains]
        strains = [[value if value != UNSET else 0 for value in row] for row in strains]

        # Define reference and set names
        constName = ['Couple-Bottom', 'Couple-Top', 'Couple-Left', 'Couple-Right']
        setName = ['Bottom-nodes', 'Top-nodes', 'Left-nodes', 'Right-nodes']

        # Apply Coupling constraints
        for i, (const, set_name) in enumerate(zip(constName, setName)):
            region1 = myAssembly.sets[self.refPointName[i]]
            region2 = myInstance.sets[set_name]
            myModel.Coupling(name=const, controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
                localCsys=None, u1=coup[i][0], u2=coup[i][1], ur3=coup[i][2])

        # Create amplitude for displacements
        ampName = 'Amp-Disp'
        myModel.EquallySpacedAmplitude( name=ampName, timeSpan=STEP,
            smooth=SOLVER_DEFAULT, fixedInterval=1.0, begin=0.0, data=(0.0, 1.0))

        # Apply boundary conditions and history output requests
        for j, refName in enumerate(self.refPointName):
            myModel.DisplacementBC(name='BC-'+refName, createStepName=stepName, region=myAssembly.sets[refName],
                u1=strains[j][0], u2=strains[j][1], ur3=strains[j][2],
                amplitude=ampName, distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None)

            myModel.HistoryOutputRequest(name='HO-'+refName,
                createStepName=stepName, variables=('U1', 'U2', 'RF1', 'RF2'),
                region=myAssembly.sets[refName], sectionPoints=DEFAULT, rebar=EXCLUDE, numIntervals=10)

        # Update the number of intervals for history output
        myModel.historyOutputRequests['H-Output-1'].setValues(numIntervals=10)
        myAssembly.regenerate()

        return
    # end: def set_periodic_boundary_conditio

    def calculate_stress_strain_from_RP(self, jobName, inc=-1):
    # Calculate stress and strain from the reaction forces and displacements of the reference points.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - stress = array of stress components
    # - strain = array of strain components

        numElCase = 4
        numComp = 3

        # Initialize arrays
        RF = []
        U = []
        strain = [0.]*numComp
        stress = [0.]*numComp
        # Open odb file
        odb = openOdb(jobName + ".odb")
        stepName = odb.steps.keys()[0]
        myStep = odb.steps[stepName]
        myAssembly = odb.rootAssembly
        
        # Store displacements and forces
        for i in range(numElCase):
            myRegion = 'Node ASSEMBLY.{}'.format(i+1)
            myHistory = myStep.historyRegions[myRegion]
            RF.append([myHistory.historyOutputs['RF1'].data[inc][1], myHistory.historyOutputs['RF2'].data[inc][1]])
            U.append([myHistory.historyOutputs['U1'].data[inc][1], myHistory.historyOutputs['U2'].data[inc][1]])
        
        # Find bounding box
        nodeSet = myAssembly.instances['MATRIX'].nodeSets['SET-MATRIX']
        xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
        yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
        self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord)] # store matrix bounding box coordinates
        self.LatticeVector = [[self.tissueBox[1]-self.tissueBox[0], 0., 0.], [0., self.tissueBox[3]-self.tissueBox[2], 0.]]
        
        # Calculate stresses from reaction forces
        stress[0] = (abs(RF[2][0])+abs(RF[3][0]))/(2*self.LatticeVector[1][1])
        stress[1] = (abs(RF[0][1])+abs(RF[1][1]))/(2*self.LatticeVector[0][0]) 
        stress[2] = ((abs(RF[0][0])+abs(RF[1][0]))/self.LatticeVector[0][0]+
                     (abs(RF[2][1])+abs(RF[3][1]))/self.LatticeVector[1][1])/4. # xy

        # Calculate strains from displacements
        strain[0] = (abs(U[2][0])+abs(U[3][0]))/self.LatticeVector[0][0]
        strain[1] = (abs(U[0][1])+abs(U[1][1]))/self.LatticeVector[1][1]
        strain[2] = ((abs(U[0][0])+abs(U[1][0]))/self.LatticeVector[1][1]+
                     (abs(U[2][1])+abs(U[3][1]))/self.LatticeVector[0][0]) # engineering xy
        
        odb.close()
        return stress, strain
    # end: def calculate_stress_strain_from_RP

    def calculate_2D_stiffness_matrix(self, jobName, inc=-1):
    # Calculate 2D stiffness matrix from 2 elementary cases results.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - StiffMat = stiffness matrix

        numElCase = 4 # Number of elementary cases
        # Initialize arrays
        stress_mat = [[0.]*3 for i in range(numElCase)]
        strain_mat = [[0.]*3 for i in range(numElCase)]
        self.StiffMat = [[0.]*3 for i in range(3)]
        # Get stresses and strains from odb files
        for i in range(numElCase):
            stress_mat[i], strain_mat[i] = self.calculate_stress_strain_from_RP(jobName=jobName+("-{}".format(i+1)), inc=inc)
        # Calculate stiffness matrix coefficients
        # Diagonal elements
        for i in range(3):
            stress = stress_mat[i]
            strain = strain_mat[i]
            self.StiffMat[i][i] = stress[i]/strain[i]
        # Out-of-diagonal elements
        k = [3]
        i = [0]
        j = [1]
        for n in range(1):
            stress = stress_mat[k[n]]
            strain = strain_mat[k[n]]
            Cij = (stress[i[n]]-self.StiffMat[i[n]][i[n]]*strain[i[n]])/strain[j[n]]
            Cji = (stress[j[n]]-self.StiffMat[j[n]][j[n]]*strain[j[n]])/strain[i[n]]
            self.StiffMat[i[n]][j[n]] = (Cij+Cji)/2.
            self.StiffMat[j[n]][i[n]] = self.StiffMat[i[n]][j[n]]

        self.CompMat = np.linalg.inv(self.StiffMat)

        self.calculate_engineering_constants()
        
        return self.StiffMat
    # end: def calculate_2D_stiffness_matrix

    def compute_2D_stiffness_matrix(self, eps, stepName, jobName, savepath):
    # Compute 2D stiffness matrix by running 4 elementary cases.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - StiffMat = stiffness matrix

        # Set elementary cases strains
        eps_x = [eps, 0., 0., eps]
        eps_y = [0., eps, 0., eps]
        eps_xy = [0., 0., eps, 0.]

        # Copy model
        numElCases = 4
        modelName = self.create_model_copy(numElCases-1)
    
        # Run elementary cases
        for k in range(numElCases):
            strains = [eps_x[k], eps_y[k], eps_xy[k]]
            # Set boundary conditions
            self.set_periodic_boundary_conditions(modelName[k], stepName, strains)
            # Run job
            self.submit_job(jobName+'-{}'.format(k+1), modelName[k])

        return #StiffMat
    # end: def compute_2D_stiffness_matrix

    def extract_info(self, jobName):
    # 
        #Split the string by '-' to get the parts
        parts = jobName.split('-')

        # Extract the structure, moisture and lignin content
        structure = parts[0]
        mc_match = search(r'mc(\d+)', jobName)
        lig_match = search(r'lig(\d+)', jobName)           
        if mc_match:
            mc = int(mc_match.group(1))/100.
        if lig_match:
            lig = int(lig_match.group(1))/100.
                
        return [structure, mc, lig]
    # end: extract_info

    def calculate_engineering_constants(self):
        # Extract engineering constants from the compliance matrix
        E_T = 1 / self.CompMat[0, 0]
        E_R = 1 / self.CompMat[1, 1]

        nu_TR = -self.CompMat[0, 1] / self.CompMat[0, 0]
        nu_RT = -self.CompMat[1, 0] / self.CompMat[1, 1]
        
        G_TR = 1 / self.CompMat[2, 2]

        # Print the results using the % operator
        print("\nE_T = {:.3f}, E_R = {:.3f}".format(E_T, E_R))
        print("nu_TR = {:.3f}".format(nu_TR))
        print("nu_RT = {:.3f}".format(nu_RT))
        print("G_RT = {:.3f}".format(G_TR))
    
        return
    # end: def calculate_engineering_constants


