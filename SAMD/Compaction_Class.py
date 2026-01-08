from abaqus import *
from abaqusConstants import *
from caeModules import *
from part import *
from mesh import *
from visualization import *
from copy import copy, deepcopy
from scipy.io import loadmat
import numpy as np
from visualization import *
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
        self.addPressure = False

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

    def set_general_constraints(self, modelName, stepName, stepTime=1, numIntervals=20, massScale=1e12):
    # Create step and set general interaction constraints:
    #   - tie interaction between matrix and fibers
    #   - self-contact to fiber inner surfaces
    #   - general contact to lateral surfaces
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - massScale = mass scaling factor (=1e12 by default)
    # - numIntervals = number of intervals for output request (=20 by default)
        
        myModel = self.mdb.models[modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        ## CREATE STEP
        self.stepTime = stepTime
        myModel.ExplicitDynamicsStep(name=stepName, timePeriod=self.stepTime, previous='Initial', improvedDtMethod=ON) # create step
        myStep = myModel.steps[stepName] # step
        myStep.setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, massScale, 0.0, None, 0, 0, 0.0, 0.0, 0, None), )) # set mass scaling
        
        ## CREATE FLUID CAVITY
        self.create_fluid_cavity(stepName, numIntervals=numIntervals)
        print("Set fluid cavity conditions.\n")

        ## CREATE INTERACTION CONSTRAINTS
        # Tie matrix (inner) and fibers (outer)
        contIntName = 'FiberMatrix-Tie' # interaction name
        region1 = myAssembly.instances[self.fibersTotName].surfaces[self.fibersTotExtSurfName[0]] # region 1 = fibers external surface
        region2 = myAssembly.instances[self.matrixName].surfaces[self.matrixInnerSurfName] # region 2 = matrix inner surface
        myModel.Tie(name=contIntName, master=region1, slave=region2, constraintEnforcement=SURFACE_TO_SURFACE, 
            positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, thickness=ON) # create tie interaction
        # Create contact interaction properties
        contPropName = 'Contact-Prop' # contact property name
        myModel.ContactProperty(contPropName) # create contact property
        myModel.interactionProperties[contPropName].TangentialBehavior(formulation=FRICTIONLESS) # set contact tangential behavior
        myModel.interactionProperties[contPropName].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT) # set contact normal behavior
        myModel.interactionProperties[contPropName].GeometricProperties(contactArea=self.lthick, padThickness=None) # set contact geometric properties (surfaces thickness)
        # Set self-contact for fibers (inner)
        for k in range(len(self.lumenSurfName)): # for each fiber
            myRegion = myAssembly.instances[self.fibersTotName ].surfaces[self.lumenSurfName[k]] # fiber inner surface
            contIntName = self.lumenSurfName[k]+'-Self-Contact' # interaction name
            myModel.SelfContactExp(createStepName='Initial', interactionProperty=contPropName,
                mechanicalConstraint=KINEMATIC, name=contIntName, surface=myRegion) # create self-contact interaction
        # Set  general contact for outer matrix and inner fiberes surfaces
        contIntName = 'ExtSurf-General-Contact' # interaction name
        myModel.ContactExp(name=contIntName, createStepName=stepName) # create general contact interaction
        region1 = myAssembly.instances[self.matrixName ].surfaces[self.matrixExtSurfName[0]] # region 1 = external surfaces
        myInteraction = myModel.interactions[contIntName] # interaction
        myInteraction.includedPairs.setValuesInStep(stepName=stepName, useAllstar=OFF, addPairs=((region1, SELF),  )) # set contact pairs
        myInteraction.contactPropertyAssignments.appendInStep(stepName=stepName, assignments=((GLOBAL, SELF, contPropName), )) # set contact property

        myAssembly.regenerate() # regenerate assembly

        ## CREATE OUTPUT REQUEST
        myModel.fieldOutputRequests['F-Output-1'].setValues(numIntervals=numIntervals)
        myModel.historyOutputRequests['H-Output-1'].setValues(numIntervals=numIntervals)
        # Energy
        hoName = 'Energy-Output'
        myModel.HistoryOutputRequest(name=hoName, createStepName=stepName,
            variables=('ALLIE', 'ALLKE'), numIntervals=numIntervals)
        # Volume
        foName = 'Volume-Output'
        myModel.FieldOutputRequest(name=foName, createStepName=stepName,
            variables=('EVOL',), numIntervals=numIntervals, timeMarks=ON)
        # Contact stress
        foName = 'Contact-Output'
        myModel.FieldOutputRequest(name=foName, createStepName=stepName,
            variables=('CSTRESS', ), timeMarks=ON, region=MODEL, exteriorOnly=OFF,
            sectionPoints=DEFAULT, rebar=EXCLUDE, numIntervals=numIntervals)

        return contIntName
    # end: set_general_constraints
    
    def create_fluid_cavity(self, stepName, numIntervals=20, density=1, expcoeff=-1000, bulkmodulus=1e-15, scale=1e-21):
    # Create fluid cavity interaction.
    # Args:
    # - stepName = name of the step
    # - density = density of the fluid [kg/m3] (=1e-3 by default)
    # - expcoeff = expansion coefficient of the fluid [1/K] (=-1000 by default)
    # - bulkmodulus = bulk modulus of the fluid [GPa] (=1e-18 by default)
    # - scale = density scale factor for unit consistency (=1e-21 by default: kg/m3 -> ton/microm3)
    # - numIntervals = number of intervals for output request (=20 by default)
        
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        # Create fluid cavity property
        contPropName = 'FluidCavity-Prop'
        myModel.FluidCavityProperty(name=contPropName, fluidDensity=density*scale, 
            useExpansion=True, expansionTable=((expcoeff, ), ), useBulkModulus=True, 
            bulkModulusTable=((bulkmodulus, ), ))
        # Create and couple RPs
        self.refPointNameCavity = [] # initialize array of reference point names
        fiber_in = [] # external fibers points
        RPKey = []
        coordRP = []
        for i in range(self.mstruct.fibers.size): # repeat for each fiber
            # Find fiber bounding box
            f1 = self.mstruct.fibers[i] # get fiber data
            a_in = f1.innerpts # external points
            fiber_in.append(a_in) # append external points
            fiber_temp = self.totuple(fiber_in[-1]) # get external coordinates
            fiber_x = [] # x-coordinates of fibers
            fiber_y = [] # y-coordinates of fibers
            for k in range(len(fiber_temp)): # for each point
                out_temp = fiber_temp[k] # point coordinates
                fiber_x.append(out_temp[0]) # append x
                fiber_y.append(out_temp[1]) # append y
            # Create reference point (RP)
            self.refPointNameCavity.append('RP-'+str(i+1)) # append RP name
            coordRP.append([min(fiber_x)+(max(fiber_x)-min(fiber_x))/2., min(fiber_y)+(max(fiber_y)-min(fiber_y))/2.])
            myAssembly.ReferencePoint(point=(coordRP[-1][0], coordRP[-1][1], 0.))
            if i==0:
                myAssembly.features.changeKey(fromName='RP-1', toName=self.refPointNameCavity[i])
            RPKey.append(myAssembly.referencePoints.keys()[0])
            region1=myAssembly.Set(referencePoints=(myAssembly.referencePoints[RPKey[-1]], ), name='Set-'+self.refPointNameCavity[-1])
            region2=myAssembly.instances[self.fibersTotName].surfaces[self.lumenSurfName[i]]
            intPropName = 'Fluid-Cavity-'+str(i)
            myModel.FluidCavity(name=intPropName, createStepName='Initial', cavityPoint=region1, 
                cavitySurface=region2, interactionProperty=contPropName, thickness=self.lthick)

        # Create set of RPs
        setName = 'Set-Cavity-RPs'
        refPoints=([myAssembly.referencePoints[RPKey[i]] for i in range(len(self.refPointNameCavity))])
        myAssembly.Set(referencePoints=(refPoints), name=setName)
        # Create cavity output
        myRegion = myModel.rootAssembly.sets[setName]
        hoName = 'Cavity-Output'
        myModel.HistoryOutputRequest(name=hoName, createStepName=stepName, variables=('PCAV', 'CVOL'),
            numIntervals=numIntervals, region=myRegion, sectionPoints=DEFAULT, rebar=EXCLUDE)

        return
    # end: create_fluid_cavity

    def create_rigid_foils(self, stepName, contIntName):
    #
                
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        # TOP FOIL
        # Create analytical rigid body
        foilName = 'Top-Foil'
        l = (self.tissueBox[1]-self.tissueBox[0])*1.5
        t = self.meshSizeMatrix
        tol = self.meshSizeMatrix/2.
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        s.rectangle(point1=(self.tissueBox[0]-(l-(self.tissueBox[1]-self.tissueBox[0]))/2.,self.tissueBox[3]+tol),
                    point2=(self.tissueBox[1]+(l-(self.tissueBox[1]-self.tissueBox[0]))/2.,self.tissueBox[3]+t+tol))
        myPart = myModel.Part(dimensionality=TWO_D_PLANAR, name=foilName, type=ANALYTIC_RIGID_SURFACE)
        myModel.parts[foilName].AnalyticRigidSurf2DPlanar(sketch=s)
        del myModel.sketches['__profile__']
        myAssembly.DatumCsysByDefault(CARTESIAN)
        myAssembly.Instance(dependent=ON, name=foilName, part=myPart)
        # Create surfaces
        foilSurfName = 'TopSurf-Foil'
        surface = myPart.edges.findAt(((self.tissueBox[0], self.tissueBox[3]+t+tol, 0.), ))
        myPart.Surface(side2Edges=surface, name=foilSurfName)       
        # Create rigid body constraint
        myRegion1 = myAssembly.sets[self.refPointName[0]]
        myRegion2 = myAssembly.instances[foilName].surfaces[foilSurfName] 
        myModel.RigidBody(name='Constraint-Top-foil', refPointRegion=myRegion1, surfaceRegion=myRegion2)         
        # Set general contact
        myRegion1 = myAssembly.instances[foilName].surfaces[foilSurfName]
        myRegion2 = myAssembly.instances[self.matrixName ].surfaces[self.matrixExtSurfName[0]]
        myRegion3 = myAssembly.instances[self.fibersTotName].surfaces[self.fibersTotExtSurfName[0]]
        myModel.interactions[contIntName].includedPairs.setValuesInStep(
            stepName=stepName, addPairs=((myRegion1, myRegion2),(myRegion1, myRegion3),))

        # BOTTOM FOIL
        foilName = 'Bottom-Foil'
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        s.rectangle(point1=(self.tissueBox[0]-(l-(self.tissueBox[1]-self.tissueBox[0]))/2.,self.tissueBox[2]-tol),
                    point2=(self.tissueBox[1]+(l-(self.tissueBox[1]-self.tissueBox[0]))/2.,self.tissueBox[2]-t-tol))
        myPart = myModel.Part(dimensionality=TWO_D_PLANAR, name=foilName, type=ANALYTIC_RIGID_SURFACE)
        myModel.parts[foilName].AnalyticRigidSurf2DPlanar(sketch=s)
        del myModel.sketches['__profile__']
        myAssembly.DatumCsysByDefault(CARTESIAN)
        myAssembly.Instance(dependent=ON, name=foilName, part=myPart)
        # Create surfaces
        foilSurfName = 'BottomSurf-Foil'
        surface = myPart.edges.findAt(((self.tissueBox[0], self.tissueBox[2]-t-tol, 0.), ))
        myPart.Surface(side1Edges=surface, name=foilSurfName)
        # Create rigid body constraint
        myRegion1 = myAssembly.sets[self.refPointName[1]]
        myRegion2 = myAssembly.instances[foilName].surfaces[foilSurfName] 
        myModel.RigidBody(name='Constraint-Bottom-foil', refPointRegion=myRegion1, surfaceRegion=myRegion2)
        # Set general contact
        myRegion1 = myAssembly.instances[foilName].surfaces[foilSurfName]
        myRegion2 = myAssembly.instances[self.matrixName ].surfaces[self.matrixExtSurfName[0]]
        myRegion3 = myAssembly.instances[self.fibersTotName].surfaces[self.fibersTotExtSurfName[0]]
        myModel.interactions[contIntName].includedPairs.setValuesInStep(
            stepName=stepName, addPairs=((myRegion1, myRegion2),))

        return
    # end: create_rigid_foils

    def mechanical_densification_bc(self, modelName, stepName, stepTime, bcValue, numIntervals=20, massScale=1e12):
    # Set boundary conditions for mechanical densification
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - bcValue = compaction level (from 0 to 1 = fraction of the tissue height)
    # - numIntervals = number of intervals for output request (=20 by default)

        self.addRPsRequest = True # request displacement and force outputs for RPs

        myModel = self.mdb.models[modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        myPart = myModel.parts[self.matrixName] # matrix part
        

        # Create step and set general contact interactions
        contIntName = self.set_general_constraints(modelName, stepName, stepTime=stepTime, numIntervals=numIntervals, massScale=massScale)
        print("Created general constraints.\n")

        # Find bounding box
        myPart.Set(nodes=myPart.nodes, name=self.matrixName+'-nodes') # create matrix nodes set
        nodeSet = myPart.sets[self.matrixName+'-nodes'] # matrix nodes set
        xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
        yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
        self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord)] # store matrix bounding box coordinates

        ## CREATE FOILS AND COUPLE REFERENCE POINTS
        self.refPointName = [] # initialize array of reference point names
        t = self.meshSizeMatrix
        # Create RP for TOP foil
        self.refPointName.append('RP-Top-foil') # append RP name
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[3]+t*2, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[0])
        myAssembly.Set(name=self.refPointName[0], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create RP for BOTTOM foil
        self.refPointName.append('RP-Bottom-foil') # append RP name
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[2]-t*2, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[1])
        myAssembly.Set(name=self.refPointName[1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create rigid foils
        self.create_rigid_foils(stepName, contIntName)
        print("Created rigid foils.\n")

        RPnode = True
        disp_nodes = [UNSET, 0, 0]
        coup_nodes = [ON, OFF, OFF]
        disp_foils = [0, 1, 0]
        ## COUPLE REFERENCE POINTS TO NODES 
        if RPnode:
            # Create RP for TOP nodes
            self.refPointName.append('RP-Top-nodes') # append RP name
            myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[3]*1.1, 0.))
            myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[2])
            myAssembly.Set(name=self.refPointName[2], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
            # Create RP for BOTTOM nodes
            self.refPointName.append('RP-Bottom-nodes') # append RP name
            myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[2]-self.tissueBox[3]*0.1, 0.))
            myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[3])
            myAssembly.Set(name=self.refPointName[3], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
                
            # Create set of top nodes and couple to RP
            setName = 'Top-nodes'
            delta = self.meshSizeMatrix/10.
            myPart.Set(name=setName, nodes=nodeSet.nodes.getByBoundingBox(min(xCoord)-delta, max(yCoord)-delta, 0, max(xCoord)+delta, max(yCoord)+delta, 0))
            myModel.Coupling(name='Constraint-Top-nodes', controlPoint=myAssembly.sets[self.refPointName[2]], 
                surface=myAssembly.instances[self.matrixName].sets[setName], influenceRadius=WHOLE_SURFACE,
                couplingType=KINEMATIC, localCsys=None, u1=coup_nodes[0], u2=coup_nodes[1], ur3=coup_nodes[2])
                
            # Create set of bottom nodes and couple to RP
            setName = 'Bottom-nodes'
            myPart.Set(name=setName, nodes=nodeSet.nodes.getByBoundingBox(min(xCoord)-delta, min(yCoord)-delta, 0, max(xCoord)+delta, min(yCoord)+delta, 0))
            myModel.Coupling(name='Constraint-Bottom-nodes', controlPoint=myAssembly.sets[self.refPointName[3]], 
                surface=myAssembly.instances[self.matrixName].sets[setName], influenceRadius=WHOLE_SURFACE,
                couplingType=KINEMATIC, localCsys=None, u1=coup_nodes[0], u2=coup_nodes[1], ur3=coup_nodes[2])
       
        # Set boundary conditions
        if bcValue>=0 and bcValue<=1: # check if the compaction level lies between 0 and 1 (included)
            
            # Apply displacement to top RP
            disp = (-(self.tissueBox[3]-self.tissueBox[2])*bcValue)/2. # displacement value
            ampName = 'Amp-Top-Disp' # amplitude name
            myModel.SmoothStepAmplitude(name=ampName, timeSpan=STEP, data=((
                0.0, 0.0), (self.stepTime, disp))) # create smooth step amplitude
            bcondName = 'Disp-Top-foil' # boundary condition name
            myRegion = myAssembly.sets[self.refPointName[0]] # region = top RP set
            myModel.DisplacementBC(amplitude=ampName, createStepName=stepName,
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                name=bcondName, region=myRegion, u1=disp_foils[0], u2=disp_foils[1], ur3=disp_foils[2]) # create displacement boundary condition
            
            # Apply displacement to top nodes
            """if RPnode:
                bcondName = 'Disp-Top-nodes' # boundary condition name
                myRegion = myAssembly.sets[self.refPointName[2]] # region = top RP set
                myModel.DisplacementBC(amplitude=ampName, createStepName=stepName, distributionType=UNIFORM,
                    fieldName='', fixed=OFF, localCsys=None, name=bcondName, region=myRegion, u1=disp_nodes[0], u2=disp_nodes[1])"""
            # Apply displacement to bottom RP
            disp = ((self.tissueBox[3]-self.tissueBox[2])*bcValue)/2. # displacement value
            ampName = 'Amp-Bottom-Disp' # amplitude name
            myModel.SmoothStepAmplitude(name=ampName, timeSpan=STEP, data=((0.0, 0.0), (self.stepTime, disp))) # create smooth step amplitude
            bcondName = 'Disp-Bottom-foil' # boundary condition name
            myRegion = myAssembly.sets[self.refPointName[1]] # region = bottom RP set
            myModel.DisplacementBC(amplitude=ampName, createStepName=stepName,
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                name=bcondName, region=myRegion, u1=disp_foils[0], u2=disp_foils[1], ur3=disp_foils[2]) # create displacement boundary condition
            
            """# Apply displacement to bottom nodes
            if RPnode:
                bcondName = 'Disp-Bottom-nodes' # boundary condition name
                myRegion = myAssembly.sets[self.refPointName[3]] # region = top RP set
                myModel.DisplacementBC(amplitude=ampName, createStepName=stepName, distributionType=UNIFORM,
                    fieldName='', fixed=OFF, localCsys=None, name=bcondName, region=myRegion, u1=disp_nodes[0], u2=disp_nodes[1])""" 
        
            # Request history output
            variables = ('U1', 'U2', 'RF1', 'RF2')
            for refPoint in self.refPointName:
                myModel.HistoryOutputRequest(name=refPoint+'-Output', createStepName=stepName, variables=variables,
                    numIntervals=numIntervals, region=myAssembly.sets[refPoint], sectionPoints=DEFAULT, rebar=EXCLUDE)

        else:
            raise ValueError("The target compaction level must be between 0 and 1. You passed {}.".format(bcValue))
        
        myAssembly.regenerate() # regenerate assembly

        return
    # end: mechanical_densification_bc

    def shear_mechanical_densification_bc(self, modelName, stepName, stepTime ,bcValue, shearAmp, shearFreq, numIntervals=20, massScale=1e12):
    # Set boundary conditions for shear-assisted mechanical densification
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - bcValue = compaction level (from 0 to 1 = fraction of the tissue height)
    # - shearAmp = amplitude of lateral periodic vibration
    # - shearFreq = frequency of lateral periodic vibration
    # - numIntervals = number of intervals for output request (=20 by default)

        myModel = self.mdb.models[modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        
        # Set mechanical densification boundary conditions
        self.mechanical_densification_bc(modelName, stepName, stepTime, bcValue, numIntervals=numIntervals, massScale=massScale)
        
        # Fix bottom
        #myModel.boundaryConditions['Disp-Bottom-RP'].setValues(u2=0) # fix vertical displacement of bottom foil
        #myModel.boundaryConditions['Disp-Bottom-nodes'].setValues(u2=0) # fix vertical displacement of bottom foil
        # Double displacement of top nodes
        #disp = (-(self.tissueBox[3]-self.tissueBox[2])*bcValue) # displacement value
        #myModel.amplitudes['Amp-Top-Disp'].setValues(data=((0.0, 0.0), (self.stepTime, disp)))

        # Create periodic amplitude
        ampName = 'Amp-Shear' # amplitude name
        myModel.PeriodicAmplitude(name=ampName, timeSpan=STEP, 
            frequency=shearFreq, start=0.0, a_0=0.0, data=((0.0, shearAmp), )) # create periodic amplitude
        # Apply shear movements
        bcondName = 'Shear-Top' # boundary condition name
        myRegion = myAssembly.sets[self.refPointName[2]]
        myModel.DisplacementBC(amplitude=ampName, createStepName=stepName,
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
            name=bcondName, region=myRegion, u1=1.0) # create displacement boundary condition
        
        bcondName = 'Shear-Bottom' # boundary condition name
        myRegion = myAssembly.sets[self.refPointName[3]]
        myModel.DisplacementBC(amplitude=ampName, createStepName=stepName,
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
            name=bcondName, region=myRegion, u1=-1.0) # create displacement boundary condition

        myAssembly.regenerate() # regenerate assembly

        return
    # end: shear_mechanical_densification_bc

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

    def create_surface_partition(self, meshSize):
    #
        myModel = self.mdb.models[self.modelName] # model
        myPart = myModel.parts[self.matrixName] # matrix part
        mask = self.mstruct.mask # get mask data
        mask = self.totuple(mask)
        delta = meshSize/3.

        # Find points on surfaces
        tolerance = 25.  # Define your tolerance
        close_points = self.find_points_within_tolerance(mask, tolerance)
        ymax, ymin = self.find_max_min_y(close_points)
        xmin = min(row[0] for row in mask)
        xmax = max(row[0] for row in mask)
        # Partition top surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymax-delta)
        myPart.PartitionEdgeByDatumPlane(datumPlane=myPart.datums[7], edges=myPart.edges)
        myPart.Surface(name=self.matrixExtSurfName[2], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymax-delta, 0., xmax+delta, ymax+delta))
        myPart.Set(name=self.matrixExtSurfName[2], edges=myPart.edges.getByBoundingBox(xmin-delta, ymax-delta, 0., xmax+delta, ymax+delta))
        # Partition bottom surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymin+delta)
        myPart.PartitionEdgeByDatumPlane(datumPlane=myPart.datums[11], edges=myPart.edges)
        myPart.Surface(name=self.matrixExtSurfName[1], side1Edges=myPart.edges.getByBoundingBox(xmin-delta, ymin-delta, 0., xmax+delta, ymin+delta))
        myPart.Set(name=self.matrixExtSurfName[1], edges=myPart.edges.getByBoundingBox(xmin-delta, ymin-delta, 0., xmax+delta, ymin+delta))
         
        return
    # end:create_surface_partition

    def vacuum_densification_bc(self, modelName, stepName, stepTime, pressValue, numIntervals=20, massScale=1e12):
    # Set boundary conditions for vacuum densification: 
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - bcValue = compaction level value (from 0 to 1, fraction of the tissue height) in 
    #   displacement control, or force value in force control
    # - dispControl = True for displacement control (default), False for force control

        myModel = self.mdb.models[modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        myInstance = myAssembly.instances[self.matrixName]

        # Create step and set general contact interactions
        self.set_general_constraints(modelName, stepName, stepTime=stepTime, numIntervals=numIntervals, massScale=massScale)
        print("Created general constraints.\n")

        # Apply pressure to top nodes   
        ampName = 'Amp-Press'
        myModel.SmoothStepAmplitude(name=ampName, timeSpan=STEP, data=((0.0, 0.0), (self.stepTime, pressValue))) # create smooth step amplitude
        loadName = 'Press-Top-nodes'
        myRegion = myInstance.surfaces[self.matrixExtSurfName[2]]
        myModel.Pressure(name=loadName, createStepName=stepName, region=myRegion, 
            distributionType=UNIFORM, field='', magnitude=1.0, amplitude=ampName)
        """myModel.SurfaceTraction(name=loadName, createStepName=stepName, region=myRegion, magnitude=1.0,
            amplitude=ampName, directionVector=((0.0, 0.0, 0.0), (0.0, -1.0, 0.0)), distributionType=UNIFORM,
            field='', localCsys=None, traction=GENERAL, follower=OFF)"""
        # Apply pressure to bottom nodes
        loadName = 'Press-Bottom-nodes'
        myRegion = myInstance.surfaces[self.matrixExtSurfName[1]]
        myModel.Pressure(name=loadName, createStepName=stepName, region=myRegion, 
            distributionType=UNIFORM, field='', magnitude=1.0, amplitude=ampName)
        """myModel.SurfaceTraction(name=loadName, createStepName=stepName, region=myRegion, magnitude=1.0,
            amplitude=ampName, directionVector=((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), distributionType=UNIFORM,
            field='', localCsys=None, traction=GENERAL, follower=OFF)"""
        
        # Create field output
        self.addPressure = True
        setName = 'Pressure-Set'
        myAssembly.SetByBoolean(name=setName, sets=(myInstance.sets[self.matrixExtSurfName[1]], 
            myInstance.sets[self.matrixExtSurfName[2]],))
        myModel.FieldOutputRequest(name='Pressure-Output', createStepName=stepName, variables=('P', ),
            numIntervals=numIntervals, region=myAssembly.sets[setName], sectionPoints=DEFAULT, rebar=EXCLUDE)
        
        myAssembly.regenerate() # regenerate assembly

        return
    # end: vacuum_densification_bc

    def self_densification_bc(self, modelName, stepName, stepTime, density, expcoeff, bulkmodulus, temperatureVar, numIntervals=20, scale=1e-21, massScale=1e12):
    # Set boundary conditions for self-densification: 
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - density = density of the fluid [kg/m3] (=1e-3 by default)
    # - expcoeff = expansion coefficient of the fluid [1/K] (=-1000 by default)
    # - bulkmodulus = bulk modulus of the fluid [GPa] (=1e-18 by default)
    # - temperatureVar = temperature variation
    # - scale = density scale factor for unit consistency (=1e-21 by default: kg/m3 -> ton/microm3)
    # - numIntervals = number of intervals for output request (=20 by default)

        myModel = self.mdb.models[modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        
        # Setup general boudary conditions
        self.set_general_constraints(modelName, stepName, stepTime=stepTime, massScale=massScale, numIntervals=numIntervals)
        
        # Setup fluid property
        contPropName = 'FluidCavity-Prop'
        myModel.interactionProperties[contPropName].setValues(fluidDensity=density*scale, 
            expansionTable=((expcoeff, ), ), bulkModulusTable=((bulkmodulus, ),))

        # Create tempearture field
        setName = 'Set-Cavity-RPs'
        ampName = 'Amp-Temperature'
        myModel.SmoothStepAmplitude(name=ampName, timeSpan=STEP, data=((0.0, 0.0), (stepTime, temperatureVar)))
        fieldName = 'Temperature-Variation'
        myRegion = myAssembly.sets[setName]
        myModel.Temperature(name=fieldName, createStepName=stepName, region=myRegion, distributionType=UNIFORM, 
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(1.0, ), amplitude=ampName)

        return
    #end: self_densification_bc

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

    def extract_info(self, jobName):
    # 
        #Split the string by '-' to get the parts
        parts = jobName.split('-')
        
        # Extract the structure name and test type
        structure = parts[0]
        test_type = parts[1]
        
        # Initialize default values
        moisture_content = None
        lignin_content = None
        rps = None
        amp = None

        # Check the test type and extract relevant info
        if test_type in ['MD', 'VD', 'SD']:
            # Extract moisture content (mc) and lignin content (lig)
            mc_match = search(r'mc(\d+)', jobName)
            lig_match = search(r'lig(\d+)', jobName)           
            if mc_match:
                moisture_content = int(mc_match.group(1))/100.
            if lig_match:
                lignin_content = int(lig_match.group(1))/100.
        
        elif test_type == 'SAMD':
            # Extract rounds per second (rps) and amplitude (amp)
            rps_match = search(r'rps(\d+)', jobName)
            amp_match = search(r'amp(\d+)', jobName)           
            if rps_match:
                rps = int(rps_match.group(1))
            if amp_match:
                amp = int(amp_match.group(1))            
            # Assign the known values for moisture content and lignin content
            moisture_content = 0.24
            lignin_content = 0.8

        if test_type in ['MD', 'SAMD']:
            self.addRPsRequest = True
        elif test_type == 'VD':
            self.addPressure = True

        # Return the extracted info as a dictionary
        info_dict = {
            'structure': structure,
            'test_type': test_type,
            'moisture_content': moisture_content,
            'lignin_stiffness': lignin_content}
        
        if test_type == 'SAMD':
            info_add = {
            'frequency': rps,
            'amplitude': amp}
            for key, value in info_add.items():
                info_dict[key] = value
                
        return info_dict
    # end: extract_info

    def save_results(self, jobName, savepath): 
    # Save results 
    
        """# Get struct data
        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab data structure
        self.mstruct = data['mstruct'] # get points data"""
        # Extract data from job name
        comp_info = self.extract_info(jobName)
        
        # Open odb file
        try:
            odb = openOdb(jobName + ".odb")
        except:
            #savepath = '/home/aferrara/Desktop/MD/LW1_MD'
            odb = openOdb(savepath+'/'+jobName + ".odb")
                    
        myAssembly = odb.rootAssembly
        stepName = odb.steps.keys()[0]
        myStep = odb.steps[stepName]
        t1 = myStep.timePeriod
        # Initialize arrays
        numInc = len(myStep.frames)-1
        tcurrent = np.zeros(numInc)
        volume1 = np.zeros(numInc)
        volume2 = np.zeros(numInc)
        cpress = np.zeros(numInc)
        strainen = np.zeros(numInc)
        kineticen = np.zeros(numInc)
        cvol = np.zeros(numInc)
        pcav = np.zeros(numInc)
        force1 = np.zeros(numInc)
        disp1 = np.zeros(numInc)
        edge1 = np.zeros(numInc)
        force2 = np.zeros(numInc)
        disp2 = np.zeros(numInc)
        edge2 = np.zeros(numInc)
        work = np.zeros(numInc)
        pressure = np.zeros(numInc)

        # Read results
        for inc in range(numInc):
            lastFrame = myStep.frames[inc]
            tinst = lastFrame.frameValue # time for the frame in the current step
            tcurrent[inc] = tinst
            
            # Tissue volume
            myField = lastFrame.fieldOutputs['EVOL']
            myRegion1 = myAssembly.instances[self.fibersTotName.upper()].elementSets['SET-'+self.fibersTotName.upper()]
            volumeValues1 = myField.getSubset(region=myRegion1).values
            myRegion2 = myAssembly.instances[self.matrixName.upper()].elementSets['SET-'+self.matrixName.upper()]
            volumeValues2 = myField.getSubset(region=myRegion2).values
            volume1[inc] = 0
            for vol in volumeValues1:
                volume1[inc] += vol.data
            volume2[inc] = 0
            for vol in volumeValues2:
                volume2[inc] += vol.data
            
            # Contact stress
            myField = lastFrame.fieldOutputs['CPRESS']
            cpressValues = myField.values
            cpress[inc] = 0
            count = 0
            for press in cpressValues:
                cpress[inc] += press.data
                if press.data != 0.:
                    count += 1.  
                if count == 0.:
                    count = 1
            cpress[inc] = cpress[inc]/count
            
            # Energy
            myHistory = myStep.historyRegions['Assembly ASSEMBLY']
            myOutput = myHistory.historyOutputs['ALLSE']
            strainen[inc] = myOutput.data[inc][1]
            myOutput = myHistory.historyOutputs['ALLKE']
            kineticen[inc] = myOutput.data[inc][1]
            
            # Cavity variables
            numFibers = len(myStep.historyRegions)-5 #self.mstruct.fibers.size
            for k in range(numFibers):
                myHistory = myStep.historyRegions['Node ASSEMBLY.{}'.format(k+1)]
                myOutput = myHistory.historyOutputs['CVOL']
                cvol[inc] += myOutput.data[inc][1]
                myOutput = myHistory.historyOutputs['PCAV']
                pcav[inc] += myOutput.data[inc][1]
            pcav[inc] = pcav[inc]/numFibers
            
            # RPs displacement and force for MD/SAMD
            work[0] = 0.
            if self.addRPsRequest:
                # Find bounding box
                nodeSet = myAssembly.instances['MATRIX'].nodeSets['SET-MATRIX']
                xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
                yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
                tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord)] # store matrix bounding box coordinates
                """
                # Find displacements and forces
                self.refPointName = ['Node ASSEMBLY.{}'.format(len(myStep.historyRegions)-2), 'Node ASSEMBLY.{}'.format(len(myStep.historyRegions)-1)]
                myHistory = myStep.historyRegions[self.refPointName[0]]
                myOutput = myHistory.historyOutputs['RF2']
                force1[inc] = myOutput.data[inc][1]
                myOutput = myHistory.historyOutputs['U2']
                disp1[inc] = myOutput.data[inc][1]
                edge1[inc] = disp1[inc]+tissueBox[3]

                myHistory = myStep.historyRegions[self.refPointName[1]]
                myOutput = myHistory.historyOutputs['RF2']
                force2[inc] = myOutput.data[inc][1]
                myOutput = myHistory.historyOutputs['U2']
                disp2[inc] = myOutput.data[inc][1]
                edge2[inc] = disp2[inc]+tissueBox[2]
                """
                # Find displacements and forces
                self.refPointName = ['Node ASSEMBLY.{}'.format(len(myStep.historyRegions)-4), 'Node ASSEMBLY.{}'.format(len(myStep.historyRegions)-3),
                                     'Node ASSEMBLY.{}'.format(len(myStep.historyRegions)-2), 'Node ASSEMBLY.{}'.format(len(myStep.historyRegions)-1)]
                myRF2 = []
                myU2 = []
                for refPoint in self.refPointName:
                    myHistory = myStep.historyRegions[refPoint]
                    myRF2.append(myHistory.historyOutputs['RF2'].data[inc][1])
                    myU2.append(myHistory.historyOutputs['U2'].data[inc][1])
                
                force1[inc] = myRF2[0]# + myRF2[2]
                force2[inc] = myRF2[1]# + myRF2[3]
                
                disp1[inc] = myU2[0]
                edge1[inc] = disp1[inc]+tissueBox[3]
                disp2[inc] = myU2[1]
                edge2[inc] = disp2[inc]+tissueBox[2]
                if inc > 0:
                    work[inc] = force1[inc]*(disp1[inc]-disp1[inc-1])+force2[inc]*(disp2[inc]-disp2[inc-1])

            if self.addPressure:
                myField = lastFrame.fieldOutputs['P']
                myRegion = myAssembly.elementSets['PRESSURE-SET']
                pressureValues = myField.getSubset(region=myRegion).values
                pressure[inc] = 0
                for press in pressureValues:
                    pressure[inc] += press.data
                pressure[inc] = pressure[inc]/len(pressureValues)

        # Store data in a dictionary
        tvol = volume1 + volume2
        porosity = np.divide(cvol,(cvol+tvol))
        data = {
            'step_time': t1, 
            'time': tcurrent,
            'tissue_volume': tvol,
            'cavity_volume': cvol,
            'porosity': porosity,
            'strain_energy': strainen,
            'kinetic_energy': kineticen,
            'contact_pressure': cpress,
            'cavity_pressure': pcav}
        for key, value in data.items():
            comp_info[key] = value
        
        if self.addRPsRequest: # if MD, SAMD
            data_add = {
                'force1': force1,
                'force2': force2,
                'displacement1': disp1,
                'displacement2': disp2,
                'edge1': edge1,
                'edge2': edge2,
                'work': work}
            for key, value in data_add.items():
                comp_info[key] = value

        if self.addPressure: # if VD
            data_add = {
                'pressure': pressure}
            for key, value in data_add.items():
                comp_info[key] = value
        
        # Save dictionary
        np.savez(savepath+jobName+'.npz',**comp_info)

        # Close odb file    
        odb.close()

        return comp_info
    # end: save_results


