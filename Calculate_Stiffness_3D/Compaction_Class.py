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
import imp
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
        self.matrixFaceSurfName = ['Face-0-Matrix', 'Face-1-Matrix']
        self.fibersFaceSurfName = ['Face-0-Fibers', 'Face-1-Fibers']
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

    def create_single_part(self, structpath, depth=1.):
    # Setup geometry of the tissue structure: create parts, sets and surfaces
    # Args:
    # - structpath = path of matlab file containing the tissue structure data
    # - scalefib = scale factor of fibers
    # Returns:
    # - modelName = name of the model

        self.mdb = Mdb() # start model
        myModel=self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        self.lthick = depth

        ############### LOAD DATA ###############
        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab structure
        self.mstruct = data['mstruct'] # get data

        ############### CREATE MATRIX ###############
        print("Creating matrix...\n")
        mask = self.totuple(self.mstruct.mask) # get mask data
        min_x = min(point[0] for point in mask)
        max_x = max(point[0] for point in mask)
        min_y = min(point[1] for point in mask)
        max_y = max(point[1] for point in mask)
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
        s.setPrimaryObject(option=STANDALONE) # set sketch objetc in viewport
        s.rectangle(point1=(min_x, min_y), point2=(max_x, max_y))
        myPart = myModel.Part(name=self.matrixName, dimensionality=THREE_D, type=DEFORMABLE_BODY) # create part
        myPart.BaseSolidExtrude(sketch=s, depth=depth) # extrude mask
        s.unsetPrimaryObject() # unset sketch objetc in viewport
        del myModel.sketches['__profile__'] # close sketch        
        
        # Create mask assembly
        myPart = myModel.parts[self.matrixName]
        myAssembly.Instance(name=self.matrixName, part=myPart, dependent=ON) # create instance
        
        # Create front and back surface of matrix
        min_x = min(point[0] for point in mask)
        max_x = max(point[0] for point in mask)
        min_y = min(point[1] for point in mask)
        max_y = max(point[1] for point in mask)
        delta = 1             
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, -1e-5, max_x+delta, max_y+delta, 1e-5)
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[0]) # create back surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[0]) # create set
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, depth-1e-5, max_x+delta, max_y+delta, depth+1e-5)      
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[1]) # create front surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[1]) # create set

        # Create surface of matrix
        surface = myPart.faces
        myPart.Surface(side1Faces=surface, name=self.matrixExtSurfName[0]) # create outer surface
           
        # Create external surface of matrix
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[0], operation=DIFFERENCE, surfaces=(myPart.surfaces[self.matrixExtSurfName[0]],
            myPart.surfaces[self.matrixFaceSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[1]],))
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[0]].faces, name=self.matrixExtSurfName[0]) # create set
        # Create matrix set
        myModel.parts[self.matrixName].Set(cells=myModel.parts[self.matrixName].cells, name='Set-'+self.matrixName) # create mask set        
        myAssembly.regenerate() # regenerate assembly

        return self.modelName
    # end: create_single_part

    def setup_geometry(self, structpath, depth=1., scalefib=1):
    # Setup geometry of the tissue structure: create parts, sets and surfaces
    # Args:
    # - structpath = path of matlab file containing the tissue structure data
    # - scalefib = scale factor of fibers
    # Returns:
    # - modelName = name of the model

        self.mdb = Mdb() # start model
        myModel=self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        self.lthick = depth

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
            myPart = myModel.Part(name=self.fibersName[-1], dimensionality=THREE_D, type=DEFORMABLE_BODY) # create part
            myPart.BaseSolidExtrude(sketch=s, depth=depth) # extrude fiber
            s.unsetPrimaryObject() # unset sketch object in viewport
            del myModel.sketches['__profile__'] # close sketch
            # Create set of fiber
            self.fibersSetsName.append('Set-'+self.fibersName[-1]) # append set name
            myPart = myModel.parts[self.fibersName[-1]]
            myPart.Set(cells=myPart.cells, name=self.fibersSetsName[-1]) # create set            
            # Create inner surface of fiber
            surface = myPart.faces.findAt(((a_in[0][0], a_in[0][1], depth/2.), )) #[myPart.edges.findAt(((point[0], point[1], 0.), )) for point in a_in]  # find inner face
            self.lumenSurfName.append('Lumen-'+self.fibersName[-1]) # append surface name
            myPart.Surface(side1Faces=surface, name=self.lumenSurfName[-1]) # create surface
            myPart.Set(faces=surface, name=self.lumenSurfName[-1]) # create set           
            # Create outer surface of fiber
            surface = myPart.faces.findAt(((a_out[0][0], a_out[0][1], depth/2.), )) #[myPart.edges.findAt(((point[0], point[1], 0.), )) for point in a_out] # find inner face
            self.fibersExtSurfName.append('ExtFiber-'+self.fibersName[-1]) # append surface name
            myPart.Surface(side1Faces=surface, name=self.fibersExtSurfName[-1]) # create surface
            myPart.Set(faces=surface, name=self.fibersExtSurfName[-1]) # create set
            # Create fiber instances
            myAssembly.Instance(name=self.fibersName[-1], part=myPart, dependent=ON)
        
        # Merge fibers in one single part
        myAssembly.InstanceFromBooleanMerge(name=self.fibersTotName , instances=([myAssembly.instances[self.fibersName[i]]
            for i in range(len(self.fibersName))] ), keepIntersections=ON, domain=GEOMETRY, originalInstances=DELETE)
        myAssembly.features.changeKey(fromName=self.fibersTotName +'-1', toName=self.fibersTotName)
        # Create fibers assembly set 
        myModel.parts[self.fibersTotName ].Set(cells=myModel.parts[self.fibersTotName].cells, name='Set-'+self.fibersTotName) # create set
        # Merge outer surfaces of fibers
        myPart = myModel.parts[self.fibersTotName]
        myPart.SurfaceByBoolean(name=self.fibersTotExtSurfName[0], surfaces=([myPart.surfaces[self.fibersExtSurfName[i]]
            for i in range(len(self.fibersExtSurfName))]))
        # Merge inner surfaces of fibers
        myPart.SurfaceByBoolean(name=self.fibersTotInnerSurfName, surfaces=([myPart.surfaces[self.lumenSurfName[i]]
            for i in range(len(self.lumenSurfName))]))
        
        ############### CREATE MATRIX ###############
        print("Creating matrix...\n")
        mask = self.mstruct.mask # get mask data
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
        s.setPrimaryObject(option=STANDALONE) # set sketch objetc in viewport
        mask = self.totuple(mask)
        s.Spline(points=mask) # draw external fiber
        for i in range(self.mstruct.fibers.size): # for each fiber 
            f2 = self.mstruct.fibers[i] # get fiber data
            a_in = self.totuple(f2.innerpts) # inner points
            s.Spline(points=a_in) # mask inner points
        myPart = myModel.Part(name=self.matrixName, dimensionality=THREE_D, type=DEFORMABLE_BODY) # create part
        myPart.BaseSolidExtrude(sketch=s, depth=depth) # extrude mask
        s.unsetPrimaryObject() # unset sketch objetc in viewport
        del myModel.sketches['__profile__'] # close sketch        
        
        # Create mask assembly
        myPart = myModel.parts[self.matrixName]
        myAssembly.Instance(name=self.matrixName, part=myPart, dependent=ON) # create instance
        
        # Create front and back surface of matrix
        min_x = min(point[0] for point in mask)
        max_x = max(point[0] for point in mask)
        min_y = min(point[1] for point in mask)
        max_y = max(point[1] for point in mask)
        print(min_x, max_x, min_y, max_y)
        delta = 10.          
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, -1e-5, max_x+delta, max_y+delta, 1e-5)
        print(surface)
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[0]) # create back surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[0]) # create set
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, depth-1e-5, max_x+delta, max_y+delta, depth+1e-5)      
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[1]) # create front surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[1]) # create set

        # Create surface of matrix
        surface = myPart.faces
        myPart.Surface(side1Faces=surface, name=self.matrixExtSurfName[0]) # create outer surface
        
        # Create front and back surface of fibers
        myPart = myModel.parts[self.fibersTotName]
        delta = 1e-5
        surface = myPart.faces.getByBoundingBox(min_x, min_y, 0.-delta, max_x, max_y, 0.+delta)
        myPart.Surface(side1Faces=surface, name=self.fibersFaceSurfName[0]) # create back surface
        myPart.Set(faces=surface, name=self.fibersFaceSurfName[0]) # create set
        surface = myPart.faces.getByBoundingBox(min_x, min_y, depth-delta, max_x, max_y, depth+delta)
        myPart.Surface(side1Faces=surface, name=self.fibersFaceSurfName[1]) # create front surface
        myPart.Set(faces=surface, name=self.fibersFaceSurfName[1]) # create set

        # Cut fibers from mask
        tempName = 'mask-temp'
        myAssembly.InstanceFromBooleanCut(cuttingInstances=(myAssembly.instances[self.fibersTotName], ), 
            instanceToBeCut=myAssembly.instances[self.matrixName], name=tempName, originalInstances=SUPPRESS)
        myAssembly.features[self.fibersTotName ].resume()
        del myModel.parts[self.matrixName]
        del myAssembly.features[self.matrixName]
        myModel.parts.changeKey(fromName=tempName, toName=self.matrixName)
        myAssembly.features.changeKey(fromName=tempName+'-1', toName=self.matrixName)     
        myAssembly.regenerate()
        
        # Create inner surface of matrix
        myPart = myModel.parts[self.matrixName]
        surface = myPart.faces
        myPart.Surface(side1Faces=surface, name=self.matrixInnerSurfName) # create surface
        myPart.SurfaceByBoolean(name=self.matrixInnerSurfName, operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixInnerSurfName], myPart.surfaces[self.matrixExtSurfName[0]], ))
        myPart.Set(faces=myPart.surfaces[self.matrixInnerSurfName].faces, name=self.matrixInnerSurfName) # create set
        
        # Create external surface of matrix
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[0], operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixExtSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[0]],
            myPart.surfaces[self.matrixFaceSurfName[1]],))
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[0]].faces, name=self.matrixExtSurfName[0]) # create set
        # Create matrix set
        myModel.parts[self.matrixName].Set(cells=myModel.parts[self.matrixName].cells, name='Set-'+self.matrixName) # create mask set        
        myAssembly.regenerate() # regenerate assembly

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
 
    def setup_single_material(self, matpath, moisture, matrix_dens, fiber_dens, alpha=[0., 0., 60., -60., 15., 75.], ligred=1., lthick=1., scale=1e-21):
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

        ### CALCULATE LOCAL STIFFNESS OF EACH CELL WALL LAYER
        # Micromaterial calc. = [E1, E2, E3, nu23, nu13, nu12, G23, G13, G12] where:
        #                                                                       - 1-direction is the longitudinal direction
        #                                                                       - 2-direction is tangential to the cell wall
        #                                                                       - 3-direction is normal to the cell wall
        # Cell wall coord. sys. in Abaqus = [E1, E2, E3, nu12, nu13, nu23, G12, G13, G23] where:
        #                                                                       - 1-direction is normal to the cell wall
        #                                                                       - 2-direction is tangential to the cell wall
        #                                                                       - 3-direction is the longitudinal direction
        #                                 = [E3, E2, E1, nu32, nu31, nu21, G23, G13, G12] from micromaterial calculation
        materialStiffMat = micromat.layer_mat_prop(moisture, matpath, 'mht_full',ligred)[0] # local stiffness matrices of cell wall layers
        
        ### SET MATRIX MATERIAL
        layer_thick = [0.175, 0.175, 0.125, 0.125, 0.645, 0.035]
        # Calculate local eng. const of cell wall (P+S1(+)+S1(-)+S2+S3)
        self.fiber_eng_const = micromat.laminate_3D_mat_prop(alpha, layer_thick, materialStiffMat) # [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12] 
        print(np.round(self.fiber_eng_const,3))
        print(np.round(micromat.convert_mat_prop(self.fiber_eng_const, 'o', 'e', 'c')[0],3))
               
        """eng_const = [self.fiber_eng_const[2], self.fiber_eng_const[1], self.fiber_eng_const[0],
                        self.fiber_eng_const[3]/self.fiber_eng_const[1]*self.fiber_eng_const[2],
                        self.fiber_eng_const[4]/self.fiber_eng_const[0]*self.fiber_eng_const[2],
                        self.fiber_eng_const[5]/self.fiber_eng_const[0]*self.fiber_eng_const[1],
                        self.fiber_eng_const[6], self.fiber_eng_const[7], self.fiber_eng_const[8]]
        print(eng_const)"""

        materialname = self.matrixName # material name
        sectionname = 'Section-'+materialname # section name
        self.create_material(materialname, self.fiber_eng_const, self.matrix_density) # create material
        self.assign_section(self.matrixName, self.matrixName, sectionname, materialname) # create and assign section
        myPart = myModel.parts[self.matrixName] # matrix part
        myRegion = myPart.sets['Set-'+self.matrixName] # matrix region
        # Set discrete material orientation = global coord. sys.
        myPart.MaterialOrientation(region=myRegion, 
            orientationType=GLOBAL, axis=AXIS_1, additionalRotationType=ROTATION_NONE, 
            localCsys=None, fieldName='', stackDirection=STACK_3)
        
        myAssembly.regenerate() # regenerate assembly

        return
    # end: setup_single_material

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
            primaryAxisRegion = myPart.surfaces[self.lumenSurfName[i]] # primary axis normal surface       
            myPart.MaterialOrientation(region=myRegion, 
                orientationType=DISCRETE, axis=AXIS_1, normalAxisDefinition=SURFACE, primaryAxisVector=(0.0, 0.0, 1.0), 
                normalAxisRegion=primaryAxisRegion, flipNormalDirection=False, normalAxisDirection=AXIS_1, 
                primaryAxisDefinition=VECTOR, primaryAxisDirection=AXIS_3, flipPrimaryDirection=False, 
                additionalRotationType=ROTATION_NONE, angle=0.0, additionalRotationField='', stackDirection=STACK_3)

        end = time.time() # stop timer
        # Print elapsed time to create and set materials
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))

        myAssembly.regenerate() # regenerate assembly
        
        return
    # end: setup_materials

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
        elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
        elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD, 
            secondOrderAccuracy=ON, distortionControl=DEFAULT)
        
        ## SET FIBERS MESH
        myPart = myModel.parts[self.fibersTotName] # merged fibers part
        # Set element type
        myRegion = myPart.sets['Set-'+self.fibersTotName ] # fibers region (set)
        myPart.setElementType(elemTypes=(elemType1, elemType2, elemType3), regions=myRegion)
        # Set mesh control
        myPart.setMeshControls(elemShape=HEX, technique=SWEEP, algorithm=ADVANCING_FRONT, regions=myPart.cells)        
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeFiber, size=self.meshSizeFiber)
        # Generate mesh
        myPart.generateMesh()

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set element type
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2, elemType3), regions=myRegion)
        # Set mesh control
        myPart.setMeshControls(elemShape=TET, technique=FREE, regions=myPart.cells)        
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

        start = time.time() # start timer
        # Element type
        elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
        elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD, 
            secondOrderAccuracy=ON, distortionControl=DEFAULT)

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set element type
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2, elemType3), regions=myRegion)
        # Set mesh control
        myPart.setMeshControls(elemShape=TET, technique=FREE, regions=myPart.cells)        
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSize, size=self.meshSize)
        # Generate mesh
        myPart.generateMesh()

        myAssembly.regenerate() # regenerate assembly

        end = time.time() # stop timer
        
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

    def create_surface_partition(self, meshSize, sides=False):
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

        print(xmin, xmax, ymin, ymax)
        print(xmin-delta, ymin-delta, 0., xmax+delta, ymin+delta)

        # Partition top surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymax-delta/2.)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[13], faces=myPart.faces)
        myPart.Surface(name=self.matrixExtSurfName[2], side1Faces=myPart.faces.getByBoundingBox(xmin-delta*5, ymax-delta, 0.-delta, xmax+delta*5, ymax+delta, self.lthick+delta))
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[2], operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixExtSurfName[2]], myPart.surfaces[self.matrixFaceSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[1]], ))
        # Partition bottom surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymin+delta/2.)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[17], faces=myPart.faces)
        myPart.Surface(name=self.matrixExtSurfName[1], side1Faces=myPart.faces.getByBoundingBox(xmin-delta*5, ymin-delta, 0.-delta, xmax+delta*5, ymin+delta, self.lthick+delta))
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[1], operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixExtSurfName[1]], myPart.surfaces[self.matrixFaceSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[1]], ))
        
        # Partition sides surfaces
        if sides:
            surfName = 'SidesSurf-Matrix'
            myPart.SurfaceByBoolean(name=surfName, operation=DIFFERENCE, surfaces=(myPart.surfaces[self.matrixExtSurfName[0]], 
                myPart.surfaces[self.matrixExtSurfName[1]],  myPart.surfaces[self.matrixExtSurfName[2]], ))        
            faces = myPart.surfaces[surfName].faces
            xmax = max(mask, key=lambda x: x[0])[0]
            xmin = min(mask, key=lambda x: x[0])[0]
            leftSide = ()
            rightSide = ()
            for i in range(len(faces)):
                if faces[i].pointOn[0][0] < (xmin+xmax)/2.:
                    leftSide += faces[i],
                else:
                    rightSide += faces[i],
            myPart.Surface(name=self.matrixExtSurfName[4], side1Faces=FaceArray(rightSide))
            myPart.Surface(name=self.matrixExtSurfName[3], side1Faces=FaceArray(leftSide))
     
        return
    # end:create_surface_partition

    def create_surfaces(self):
    #
        myModel = self.mdb.models[self.modelName] # model
        myPart = myModel.parts[self.matrixName] # matrix part
        mask = self.mstruct.mask # get mask data
        mask = self.totuple(mask)
        delta = 1e-3
        depth = 1.

        # Find points on surfaces
        ymax = max(row[1] for row in mask)
        ymin = min(row[1] for row in mask)
        xmin = min(row[0] for row in mask)
        xmax = max(row[0] for row in mask)

        # Top surface
        myPart.Surface(name=self.matrixExtSurfName[2], side1Faces=myPart.faces.getByBoundingBox(xmin-delta, ymax-delta, 0.-delta, xmax+delta, ymax+delta, depth+delta))
        # Bottom surface
        myPart.Surface(name=self.matrixExtSurfName[1], side1Faces=myPart.faces.getByBoundingBox(xmin-delta, ymin-delta, 0.-delta, xmax+delta, ymin+delta, depth+delta))
        # Left surface
        myPart.Surface(name=self.matrixExtSurfName[3], side1Faces=myPart.faces.getByBoundingBox(xmin-delta, ymin-delta, 0.-delta, xmin+delta, ymax+delta, depth+delta))
        # Right surface
        myPart.Surface(name=self.matrixExtSurfName[4], side1Faces=myPart.faces.getByBoundingBox(xmax-delta, ymin-delta, 0.-delta, xmax+delta, ymax+delta, depth+delta))

        """# Partition sides surfaces
        surfName = 'SidesSurf-Matrix'
        myPart.SurfaceByBoolean(name=surfName, operation=DIFFERENCE, surfaces=(myPart.surfaces[self.matrixExtSurfName[0]], 
            myPart.surfaces[self.matrixExtSurfName[1]],  myPart.surfaces[self.matrixExtSurfName[2]], ))        
        faces = myPart.surfaces[surfName].faces
        xmax = max(mask, key=lambda x: x[0])[0]
        xmin = min(mask, key=lambda x: x[0])[0]
        leftSide = ()
        rightSide = ()
        for i in range(len(faces)):
            if faces[i].pointOn[0][0] < (xmin+xmax)/2.:
                leftSide += faces[i],
            else:
                rightSide += faces[i],
        myPart.Surface(name=self.matrixExtSurfName[4], side1Faces=FaceArray(rightSide))
        myPart.Surface(name=self.matrixExtSurfName[3], side1Faces=FaceArray(leftSide))"""
     
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
        myJob.submit(consistencyChecking=OFF) # submit job
        start = time.time() # start timer
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
        # Create back RP
        self.refPointName.append('RP-Back') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., (self.tissueBox[3]+self.tissueBox[2])/2., -t))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create front RP
        self.refPointName.append('RP-Front') # append RP name
        depth = 1.     
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., (self.tissueBox[3]+self.tissueBox[2])/2., depth+t))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
               
        # Create sides node sets
        setName1 = ['Bottom-nodes', 'Top-nodes', 'Left-nodes', 'Right-nodes']
        for i in range(len(setName1)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName1[i], nodes=myPart.surfaces[self.matrixExtSurfName[1+i]].nodes)
        myPart.SetByBoolean(name=setName1[2], operation=DIFFERENCE, sets=(
            myPart.sets[setName1[2]], myPart.sets[setName1[0]], myPart.sets[setName1[1]],))
        myPart.SetByBoolean(name=setName1[3], operation=DIFFERENCE, sets=(
            myPart.sets[setName1[3]], myPart.sets[setName1[0]], myPart.sets[setName1[1]],))
        # Create back/front node sets
        setName2 = ['Back-nodes', 'Front-nodes']
        for i in range(len(setName2)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName2[i], nodes=myPart.surfaces[self.matrixFaceSurfName[i]].nodes)
        
        myPart.Set(name=self.matrixInnerSurfName+'-nodes', nodes=myPart.sets[self.matrixInnerSurfName].nodes)

        myPart.SetByBoolean(name=setName2[0], operation=DIFFERENCE, sets=(myPart.sets[setName2[0]], myPart.sets[setName1[0]],
                            myPart.sets[setName1[1]], myPart.sets[setName1[2]], myPart.sets[setName1[3]], myPart.sets[self.matrixInnerSurfName+'-nodes'],))    
        myPart.SetByBoolean(name=setName2[1], operation=DIFFERENCE, sets=(myPart.sets[setName2[1]], myPart.sets[setName1[0]],
                            myPart.sets[setName1[1]], myPart.sets[setName1[2]], myPart.sets[setName1[3]], myPart.sets[self.matrixInnerSurfName+'-nodes'],))        
        # Merge back/front matrix nodes with fibers nodes
        myPart = myModel.parts[self.fibersTotName]
        setName2 = ['Back-nodes', 'Front-nodes']
        for i in range(len(setName2)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName2[i], nodes=myPart.surfaces[self.fibersFaceSurfName[i]].nodes)
            myAssembly.SetByBoolean(name=setName2[i], sets=(
                myAssembly.allInstances[self.matrixName].sets[setName2[i]], 
                myAssembly.allInstances[self.fibersTotName].sets[setName2[i]], ))

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
        # Create back RP
        self.refPointName.append('RP-Back') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., (self.tissueBox[3]+self.tissueBox[2])/2., -t))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create front RP
        self.refPointName.append('RP-Front') # append RP name
        depth = 1.     
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., (self.tissueBox[3]+self.tissueBox[2])/2., depth+t))
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
        # Create back/front node sets
        setName2 = ['Back-nodes', 'Front-nodes']
        for i in range(len(setName2)):
            myPart.Set(name=setName2[i], nodes=myPart.surfaces[self.matrixFaceSurfName[i]].nodes)
        myPart.SetByBoolean(name=setName2[0], operation=DIFFERENCE, sets=(myPart.sets[setName2[0]], myPart.sets[setName1[0]],
                            myPart.sets[setName1[1]], myPart.sets[setName1[2]], myPart.sets[setName1[3]],))    
        myPart.SetByBoolean(name=setName2[1], operation=DIFFERENCE, sets=(myPart.sets[setName2[1]], myPart.sets[setName1[0]],
                            myPart.sets[setName1[1]], myPart.sets[setName1[2]], myPart.sets[setName1[3]],))        
        
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
            [0.0, 0.0, depth]
        ]
        
        # Calculate displacements
        # Array of strains
        disp_x = strain[0]*self.LatticeVector[0][0]/2.
        disp_y = strain[1]*self.LatticeVector[1][1]/2.
        disp_z = strain[2]*self.LatticeVector[2][2]/2.
        #################################################
        disp_yz_z = strain[3]*self.LatticeVector[1][1]/2.
        disp_yz_y = strain[3]*self.LatticeVector[2][2]/2.
        #################################################
        disp_xz_z = strain[4]*self.LatticeVector[0][0]/2.
        disp_xz_x = strain[4]*self.LatticeVector[2][2]/2.
        #################################################
        disp_xy_y = strain[5]*self.LatticeVector[0][0]/2.
        disp_xy_x = strain[5]*self.LatticeVector[1][1]/2.
        
        # Set displacements and coupling DOF
        if strain[0] != 0 and all(x == 0.0 for x in strain[1:]):  # Case 1: x-strain
            strains = [[UNSET, 0, 0, 0, 0, 0], [UNSET, 0, 0, 0, 0, 0],      # bottom, top
                       [-disp_x, 0, 0, UNSET, UNSET, UNSET], [disp_x, 0, 0, UNSET, UNSET, UNSET],   # left, right
                       [UNSET, 0, 0, UNSET, UNSET, UNSET], [UNSET, 0, 0, UNSET, UNSET, UNSET]]      # back, front
            
        elif strain[1] != 0 and all(x == 0.0 for x in [strain[0]] + strain[2:]):  # Case 2: y-strain
            strains = [[0, -disp_y, 0, UNSET, UNSET, UNSET], [0, disp_y, 0, UNSET, UNSET, UNSET],   # bottom, top
                       [0, UNSET, 0, UNSET, UNSET, UNSET], [0, UNSET, 0, UNSET, UNSET, UNSET],      # left, right
                       [0, UNSET, 0, UNSET, UNSET, UNSET], [0, UNSET, 0, UNSET, UNSET, UNSET]]      # back, front
            
        elif strain[2] != 0 and all(x == 0.0 for x in strain[:2] + strain[3:]):  # Case 3: z-strain
            strains = [[0, 0, UNSET, UNSET, UNSET, UNSET], [0, 0, UNSET, UNSET, UNSET, UNSET],      # bottom, top
                       [0, 0, UNSET, UNSET, UNSET, UNSET], [0, 0, UNSET, UNSET, UNSET, UNSET],      # left, right
                       [0, 0, -disp_z, UNSET, UNSET, UNSET], [0, 0, disp_z, UNSET, UNSET, UNSET]]   # back, front
            
        ##############################################################################################################

        elif strain[3] != 0 and all(x == 0.0 for x in strain[:3] + strain[4:]):  # Case 4: yz-strain
            strains = [[0, UNSET, disp_yz_z, UNSET, 0, 0], [0, UNSET, -disp_yz_z, UNSET, 0, 0],  # bottom, top
                       [0, UNSET, UNSET, UNSET, 0, 0], [0, UNSET, UNSET, UNSET, 0, 0],           # left, right
                       [0, disp_yz_y, UNSET, UNSET, 0, 0], [0, -disp_yz_y, UNSET, UNSET, 0, 0]]  # back, front
            
        elif strain[4] != 0 and all(x == 0.0 for x in strain[:4] + strain[5:]):  # Case 5: xz-strain
            strains = [[UNSET, 0, UNSET, 0, UNSET, 0], [UNSET, 0, UNSET, 0, UNSET, 0],          # bottom, top
                       [UNSET, 0, disp_xz_z, 0, UNSET, 0], [UNSET, 0, -disp_xz_z, 0, UNSET, 0], # left, right
                       [disp_xz_x, 0, UNSET, 0, UNSET, 0], [-disp_xz_x, 0, UNSET, 0, UNSET, 0]] # back, front         
            
        elif strain[5] != 0 and all(x == 0.0 for x in strain[:5]):  # Case 6: xy-strain
            strains = [[disp_xy_x, UNSET, 0, 0, 0, UNSET], [-disp_xy_x, UNSET, 0, 0, 0, UNSET], # bottom, top
                       [UNSET, disp_xy_y, 0, 0, 0, UNSET], [UNSET, -disp_xy_y, 0, 0, 0, UNSET], # left, right
                       [UNSET, UNSET, 0, 0, 0, UNSET], [UNSET, UNSET, 0, 0, 0, UNSET]]          # back, front 

        ##############################################################################################################

        elif strain[0] != 0 and strain[1] != 0 and all(x == 0.0 for x in strain[2:]):  # Case 7: x+y-strain
            strains = [[UNSET, -disp_y, 0, UNSET, UNSET, UNSET], [UNSET, disp_y, 0, UNSET, UNSET, UNSET],   # bottom, top
                       [-disp_x, UNSET, 0, UNSET, UNSET, UNSET], [disp_x, UNSET, 0, UNSET, UNSET, UNSET],   # left, right
                       [UNSET, UNSET, 0, UNSET, UNSET, UNSET], [UNSET, UNSET, 0, UNSET, UNSET, UNSET]]      # back, front 
            
        elif strain[0] != 0 and strain[2] != 0 and strain[1] == 0 and not any(strain[3:]):  # Case 8: x+z-strain
            strains = [[UNSET, 0, UNSET, UNSET, UNSET, UNSET], [UNSET, 0, UNSET, UNSET, UNSET, UNSET],      # bottom, top
                       [-disp_x, 0, UNSET, UNSET, UNSET, UNSET], [disp_x, 0, UNSET, UNSET, UNSET, UNSET],   # left, right
                       [UNSET, 0, -disp_z, UNSET, UNSET, UNSET], [UNSET, 0, disp_z, UNSET, UNSET, UNSET]]   # back, front 
            
        elif strain[1] != 0 and strain[2] != 0 and strain[0] == 0 and not any(strain[3:]):  # Case 9: y+z-strain
            strains = [[0, -disp_y, UNSET, UNSET, UNSET, UNSET], [0, disp_y, UNSET, UNSET, UNSET, UNSET],   # bottom, top
                       [0, UNSET, UNSET, UNSET, UNSET, UNSET], [0, UNSET, UNSET, UNSET, UNSET, UNSET],      # left, right
                       [0, UNSET, -disp_z, UNSET, UNSET, UNSET], [0, UNSET, disp_z, UNSET, UNSET, UNSET]]   # back, front 
            
        elif all(x!=0.0 for x in strain[:3]) and not any(strain[3:]):  # Case 10: x+y+z-strain
            strains = [[UNSET, -disp_y, UNSET, UNSET, UNSET, UNSET], [UNSET, disp_y, UNSET, UNSET, UNSET, UNSET],   # bottom, top
                       [-disp_x, UNSET, UNSET, UNSET, UNSET, UNSET], [disp_x, UNSET, UNSET, UNSET, UNSET, UNSET],      # left, right
                       [UNSET, UNSET, -disp_z, UNSET, UNSET, UNSET], [UNSET, UNSET, disp_z, UNSET, UNSET, UNSET]]   # back, front 
            
        else:
            raise ValueError("The strain case you provided is not implemented: ", strain)
        
        coup = [[ON if value != UNSET else OFF for value in row] for row in strains]
        strains = [[value if value != UNSET else 0 for value in row] for row in strains]

        # Define reference and set names
        constName = ['Couple-Bottom', 'Couple-Top', 'Couple-Left', 'Couple-Right', 'Couple-Back', 'Couple-Front']
        setName = ['Bottom-nodes', 'Top-nodes', 'Left-nodes', 'Right-nodes', 'Back-nodes', 'Front-nodes']

        # Apply Coupling constraints
        for i, (const, set_name) in enumerate(zip(constName, setName)):
            region1 = myAssembly.sets[self.refPointName[i]]
            region2 = (myInstance.sets[set_name] if i < 4 else myAssembly.sets[set_name])
            myModel.Coupling(name=const, controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
                localCsys=None, u1=coup[i][0], u2=coup[i][1], u3=coup[i][2], ur1=coup[i][3], ur2=coup[i][4], ur3=coup[i][5])

        # Create amplitude for displacements
        ampName = 'Amp-Disp'
        myModel.EquallySpacedAmplitude( name=ampName, timeSpan=STEP,
            smooth=SOLVER_DEFAULT, fixedInterval=1.0, begin=0.0, data=(0.0, 1.0))

        # Apply boundary conditions and history output requests
        for j, refName in enumerate(self.refPointName):
            myModel.DisplacementBC(name='BC-'+refName, createStepName=stepName, region=myAssembly.sets[refName],
                u1=strains[j][0], u2=strains[j][1], u3=strains[j][2], ur1=strains[j][3], ur2=strains[j][4], ur3=strains[j][5],
                amplitude=ampName, distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None)

            myModel.HistoryOutputRequest(name='HO-'+refName,
                createStepName=stepName, variables=('U1', 'U2', 'U3', 'RF1', 'RF2', 'RF3'),
                region=myAssembly.sets[refName], sectionPoints=DEFAULT, rebar=EXCLUDE, numIntervals=10)

        # Update the number of intervals for history output
        numInt = 10
        myModel.fieldOutputRequests['F-Output-1'].setValues(numIntervals=numInt, timeMarks=ON)
        myModel.historyOutputRequests['H-Output-1'].setValues(numIntervals=numInt)
        myAssembly.regenerate()

        return
    # end: def set_periodic_boundary_conditions

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

    def compute_3D_stiffness_matrix(self, eps, stepName, jobName, savepath):
    # Compute 2D stiffness matrix by running 4 elementary cases.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - StiffMat = stiffness matrix

        # Set elementary cases strains
        eps_x =     [eps,       eps,    eps,    0.5*eps,    eps,            0., 0., 0.]
        eps_y =     [0.5*eps,   eps,    eps,    eps,        0.5*eps,        0., 0., 0.]
        eps_z =     [0.,        0.,     eps,    0.5*eps,    0.5*eps,        0., 0., 0.]
        eps_yz =    [0.,        0.,     0.,     0.,         0.,             eps, 0., 0.]
        eps_xz =    [0.,        0.,     0.,     0.,         0.,             0., eps, 0.]
        eps_xy =    [0.,        0.,     0.,     0.,         0.,             0., 0., eps]

        # Copy model
        numElCases = 8
        modelName = self.create_model_copy(numElCases-1)
    
        # Run elementary cases
        for k in range(numElCases):
            strains = [eps_x[k], eps_y[k], eps_z[k], eps_yz[k], eps_xz[k], eps_xy[k]]
            # Set boundary conditions
            self.set_periodic_boundary_conditions(modelName[k], stepName, strains)
            # Run job
            self.submit_job(jobName+'-{}'.format(k+1), modelName[k])

        return #StiffMat
    # end: def compute_3D_stiffness_matrix

    def calculate_3D_stiffness_matrix(self, jobName, inc=-1):
    # Calculate 2D stiffness matrix from 2 elementary cases results.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - StiffMat = stiffness matrix

        numElCase = 8 # Number of elementary cases
        numComp = 6
        
        # Initialize arrays
        stress_mat = [[0.]*numComp for i in range(numElCase)]
        strain_mat = [[0.]*numComp for i in range(numElCase)]
        self.StiffMat = [[0.]*numComp for i in range(numComp)]
        
        # Get stresses and strains from odb files
        for i in range(numElCase):
            #print('Case {}:'.format(i+1))
            stress_mat[i], strain_mat[i] = self.calculate_stress_strain_from_RP(jobName=jobName+("-{}".format(i+1)), inc=inc)
            #print(stress_mat[i], strain_mat[i])
        
        # Calculate stiffness matrix coefficients
        # Shear diagonal elements
        k = 3
        for i in range(numElCase-3,numElCase): # last 3 elementary cases
            stress = stress_mat[i]
            strain = strain_mat[i]
            self.StiffMat[k][k] = stress[k]/strain[k]
            k += 1
        
        # RT plane elements
        i = 0
        j = 1 
        # Case 1
        stressi = [stress_mat[0][i],stress_mat[1][i]]
        straini = [strain_mat[0][i],strain_mat[1][i]]
        # Case 2
        stressj = [stress_mat[0][j],stress_mat[1][j]]
        strainj = [strain_mat[0][j],strain_mat[1][j]]            
        # Out-of-diagonal elements
        Cij = (stressi[0]*straini[1]-stressi[1]*straini[0])/(strainj[0]*straini[1]-straini[0]*strainj[1])           
        Cji = (stressj[0]*strainj[1]-stressj[1]*strainj[0])/(straini[0]*strainj[1]-strainj[0]*straini[1])    
        self.StiffMat[i][j] = (Cij+Cji)/2.
        self.StiffMat[j][i] = self.StiffMat[i][j]
        #print(Cij, Cji)           
        # Diagonal elements
        self.StiffMat[0][0] = stressi[0]/straini[0]-(Cij+Cji)/2.*strainj[0]/straini[0]
        self.StiffMat[1][1] = stressj[0]/strainj[0]-(Cij+Cji)/2.*straini[0]/strainj[0]
        
        # L direction elements
        stressi = []
        straini = []
        b = []
        C13 = []
        C23 = []
        for i in range(2, 5):
            stressi.append([stress_mat[i][0],stress_mat[i][1],stress_mat[i][2]])
            straini.append([strain_mat[i][0],strain_mat[i][1],strain_mat[i][2]])
            b.append(stress_mat[i][2])
            C13.append((stressi[-1][0]-straini[-1][0]*self.StiffMat[0][0]-straini[-1][1]*self.StiffMat[0][1])/straini[-1][2])
            C23.append((stressi[-1][1]-straini[-1][0]*self.StiffMat[1][0]-straini[-1][1]*self.StiffMat[1][1])/straini[-1][2])

        solution = np.linalg.solve(straini, b)
        
        self.StiffMat[2][0] = (solution[0]+np.mean(C13))/2.
        self.StiffMat[0][2] = self.StiffMat[2][0]
        self.StiffMat[2][1] = (solution[1]+np.mean(C23))/2.
        self.StiffMat[1][2] = self.StiffMat[2][0]
        self.StiffMat[2][2] = solution[2]

        """k = 2
        stress = stress_mat[k]
        strain = strain_mat[k]
        print(stress)
        print(strain)
        for i in range(2):
            self.StiffMat[i][2] = (stress[i]-strain[0]*self.StiffMat[i][0]-strain[1]*self.StiffMat[i][1])/strain[2]
            self.StiffMat[2][i] = self.StiffMat[i][2]
        self.StiffMat[2][2] = (stress[2]-strain[0]*self.StiffMat[2][0]-strain[1]*self.StiffMat[2][1])/strain[2]"""
        print('\nStiffness matrix:')
        print(np.round(np.matrix(self.StiffMat),3))

        print('\nCompliance matrix:')
        np.set_printoptions(suppress=True)
        self.CompMat = np.linalg.inv(self.StiffMat)
        print(np.round(np.matrix(self.CompMat),3))

        self.calculate_engineering_constants()
        

        return self.StiffMat
    # end: def calculate_3D_stiffness_matrix

    def calculate_stress_strain_from_RP(self, jobName, inc=-1):
    # Calculate stress and strain from the reaction forces and displacements of the reference points.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - stress = array of stress components
    # - strain = array of strain components

        numComp = 6

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
        for i in range(numComp):
            myRegion = 'Node ASSEMBLY.{}'.format(i+1)
            myHistory = myStep.historyRegions[myRegion]
            RF.append([myHistory.historyOutputs['RF1'].data[inc][1], myHistory.historyOutputs['RF2'].data[inc][1], myHistory.historyOutputs['RF3'].data[inc][1]])
            U.append([myHistory.historyOutputs['U1'].data[inc][1], myHistory.historyOutputs['U2'].data[inc][1], myHistory.historyOutputs['U3'].data[inc][1]])
        
        # Find bounding box
        nodeSet = myAssembly.instances['MATRIX'].nodeSets['SET-MATRIX']
        xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
        yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
        zCoord = [curNode.coordinates[2] for curNode in nodeSet.nodes]
        self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord), min(zCoord), max(zCoord)] # store matrix bounding box coordinates       
        self.LatticeVector = [[self.tissueBox[1]-self.tissueBox[0], 0., 0.], [0., self.tissueBox[3]-self.tissueBox[2], 0.], [0.0, 0.0, self.tissueBox[5]-self.tissueBox[4]]]
        
        # Calculate stresses from reaction forces
        stress[0] = (abs(RF[2][0])+abs(RF[3][0]))/(2*self.LatticeVector[1][1]*self.LatticeVector[2][2])
        stress[1] = (abs(RF[0][1])+abs(RF[1][1]))/(2*self.LatticeVector[0][0]*self.LatticeVector[2][2])
        stress[2] = (abs(RF[4][2])+abs(RF[5][2]))/(2*self.LatticeVector[0][0]*self.LatticeVector[1][1])       
        stress[3] = ((abs(RF[0][2])+abs(RF[1][2]))/(self.LatticeVector[0][0]*self.LatticeVector[2][2])+
                     (abs(RF[4][1])+abs(RF[5][1]))/(self.LatticeVector[1][1]*self.LatticeVector[0][0]))/4. # yz       
        stress[4] = ((abs(RF[2][2])+abs(RF[3][2]))/(self.LatticeVector[2][2]*self.LatticeVector[1][1])+
                     (abs(RF[4][0])+abs(RF[5][0]))/(self.LatticeVector[0][0]*self.LatticeVector[1][1]))/4. # xz
        
        stress[5] = ((abs(RF[0][0])+abs(RF[1][0]))/(self.LatticeVector[0][0]*self.LatticeVector[2][2])+
                     (abs(RF[2][1])+abs(RF[3][1]))/(self.LatticeVector[1][1]*self.LatticeVector[2][2]))/4. # xy

        # Calculate strains from displacements
        strain[0] = (abs(U[2][0])+abs(U[3][0]))/self.LatticeVector[0][0]
        strain[1] = (abs(U[0][1])+abs(U[1][1]))/self.LatticeVector[1][1]
        strain[2] = (abs(U[4][2])+abs(U[5][2]))/self.LatticeVector[2][2]
        strain[3] = ((abs(U[0][2])+abs(U[1][2]))/self.LatticeVector[1][1]+
                     (abs(U[4][1])+abs(U[5][1]))/self.LatticeVector[2][2]) # engineering yz            
        strain[4] = ((abs(U[2][2])+abs(U[3][2]))/self.LatticeVector[0][0]+
                     (abs(U[4][0])+abs(U[5][0]))/self.LatticeVector[2][2]) # engineering xz
        strain[5] = ((abs(U[0][0])+abs(U[1][0]))/self.LatticeVector[1][1]+
                     (abs(U[2][1])+abs(U[3][1]))/self.LatticeVector[0][0]) # engineering xy

        odb.close()

        #print(RF, stress)
        return stress, strain
    # end: def calculate_stress_strain_from_RP

    def calculate_engineering_constants(self):
        # Extract engineering constants from the compliance matrix
        E_T = 1 / self.CompMat[0, 0]
        E_R = 1 / self.CompMat[1, 1]
        E_L = 1 / self.CompMat[2, 2]

        nu_TR = -self.CompMat[0, 1] / self.CompMat[0, 0]
        nu_TL = -self.CompMat[0, 2] / self.CompMat[0, 0]
        nu_RL = -self.CompMat[1, 2] / self.CompMat[1, 1]

        nu_RT = -self.CompMat[1, 0] / self.CompMat[1, 1]
        nu_LT = -self.CompMat[2, 0] / self.CompMat[2, 2]
        nu_LR = -self.CompMat[2, 1] / self.CompMat[2, 2]

        G_TR = 1 / self.CompMat[5, 5]
        G_TL = 1 / self.CompMat[4, 4]
        G_RL = 1 / self.CompMat[3, 3]

        # Print the results using the % operator
        print("\nE_T = {:.3f}, E_R = {:.3f}, E_L = {:.3f}".format(E_T, E_R, E_L))
        print("nu_TR = {:.3f}, nu_TL = {:.3f}, nu_RL = {:.3f}".format(nu_TR, nu_TL, nu_RL))
        print("nu_RT = {:.3f}, nu_LT = {:.3f}, nu_LR = {:.3f}".format(nu_RT, nu_LT, nu_LR))
        print("G_TR = {:.3f}, G_TL = {:.3f}, G_RL = {:.3f}".format(G_TR, G_TL, G_RL))
    
        return
    # end: def calculate_engineering_constants
   
