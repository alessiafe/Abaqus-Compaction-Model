
# Micromaterial_Calculation.py
# Author Alessia Ferrara, May 2022
# Last modified October 2023
# Library of functions to compute the hygroelastic properties of the wood cell wall layers, as well as the equivalent properties as planar or 3D laminate composite.
#
# Original Matlab version on /nas22.ethz.ch/baug_ifb_woodmat/share/projects/SNF_Mech_Creep/Wood_cellular_model/01_microstructure
# D. F. Mora Mendez et al. 'Mechanical behavior of chemically modified Norway spruce: a generic hierarchical model for wood modifications' (2019)

import numpy as np
import math as math
from mpmath import *
mp.dps = 150 # increase precision to stabilize matrix operations (default = 15)

def laminate_3D_mat_prop(ply_angles, ply_thick, ply_materials):
# Calculate engineering constants of 3D laminate composite [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12]
# Args:
# - ply_angles = array of laminate angles
# - ply_thick = array of laminate thicknesses
# - ply_materials = array of laminate stiffness matrices
# Return:
# - lam_eng_const = total laminate engineering constants

    C = full_laminate_matrix(ply_angles, ply_thick, ply_materials) # laminate stiffness matrix
    lam_eng_const = convert_mat_prop(C, 'o', 'c', 'e')[0]

    return lam_eng_const
# End: def laminate_3D_mat_prop

def full_laminate_matrix(ply_angles, ply_thick, ply_materials):
# Calculate full stiffness matrix of 3D laminate composite
# Args:
# - ply_angles = array of laminate angles
# - ply_thick = array of laminate thicknesses
# - ply_materials = array of laminate stiffness matrices
# Return:
# - D = total laminate stiffness matrix

    Nplies = len(ply_angles) # number of plies
    tot_thick = sum(ply_thick) # total laminate thickness
    V = np.divide(ply_thick,tot_thick) # array of volume fractions
    C = np.empty([Nplies,6,6])
    # Rotate stiffness matrix from local to global sys (on 1-2 plane (laminate plane), around 3rd axis)
    for k in range(Nplies):
        C[k] = transform_stiffness_matrix(ply_materials[k],ply_angles[k],3) # rotate stiffness matrix
    D = np.empty([6,6]) # initialize total stiffness matrix
    # 3D laminate formulation by Chou et al. Elastic Constants of Layered Media (1972)
    for i in [0,1,2,5]:
        for j in [0,1,2,5]:
            D[i][j] = 0
            for k in range(Nplies):
                a = 0
                for l in range(Nplies):
                    a = a + V[l]*C[l][2][j]/C[l][2][2]
                b = 0
                for l in range(Nplies):
                    b = b + V[l]/C[l][2][2]
                D[i][j] = D[i][j] + V[k]*(C[k][i][j] - C[k][i][2]*C[k][2][j]/C[k][2][2] + C[k][i][2] * a/(C[k][2][2]*b))

    for i in [3, 4]:
        for j in [3, 4]:
            a = 0
            for k in range(Nplies):
                a = a + V[k]*C[k][i][j]/np.linalg.det([[C[k][3][3], C[k][3][4]], [C[k][4][3], C[k][4][4]]])
            b = 0
            for k in range(Nplies):
                for l in range(Nplies):
                    b = b + V[k]*V[l]*(C[k][3][3]*C[l][4][4]-C[k][3][4]*C[l][4][3])/(np.linalg.det([[C[k][3][3], C[k][3][4]], [C[k][4][3], C[k][4][4]]])*np.linalg.det([[C[l][3][3], C[l][3][4]], [C[l][4][3], C[l][4][4]]]))
            D[i][j] = a / b

    # Impose orthotropy
    D[[0,1,2],[5,5,5]] = 0.
    D[[5,5,5],[0,1,2]] = 0.
    D[[3,4],[4,3]] = 0.

    return D
# End: def full_laminate_matrix

def laminate_2D_mat_prop(ply_angles, ply_thick, ply_materials):
# Calculate engineering constants of plate laminate composite [E1,E2,nu12,nu21,G12]
# Args:
# - ply_angles = array of laminate angles
# - ply_thick = array of laminate thicknesses
# - ply_materials = array of laminate stiffness tensor 
# Return:
# - [E1, E2, nu12, nu21, G12] = plane laminate engineering constants

    ADBm = ADB_matrix(ply_angles, ply_thick, ply_materials) # calculate ADB matrix = [[A, B]; [B, D]]
    A = ADBm[0][0] # get A matrix
    B = ADBm[0][1] # get B matrix
    D = ADBm[1][1] # get D matrix
   
    # inverse of ADB = [[A, B]; [B, D]]^(-1) = [[A', B']; [B', D']]
    D_inv = np.linalg.inv(D-np.matmul(np.matmul(B, np.linalg.inv(A)), np.linalg.inv(B))) # D' matrix of inverse of ADB
    A_inv = np.linalg.inv(A)+(np.linalg.multi_dot([np.linalg.inv(A), B, D_inv, B, np.linalg.inv(A)])) # A' matrix of inverse of ADB

    tot_thick = sum(ply_thick) # laminate thickness
    # Calculate engineering constants
    E1 = 1./(tot_thick*A_inv[0,0])
    E2 = 1./(tot_thick*A_inv[1,1])
    G12 = 1./(tot_thick*A_inv[2,2])
    nu12 = -A_inv[1,0]/A_inv[0,0]
    nu21 = -A_inv[0,1]/A_inv[1,1]

    return [E1, E2, nu12, nu21, G12]
# End: def laminate_2D_mat_prop

def ADB_matrix(ply_angles, ply_thick, ply_materials):
# Calculate ABD matrix of a laminate.
# Args:
# - ply_angles = array of laminate angles
# - ply_thick = array of laminate thicknesses
# - ply_materials = array of laminate stiffness tensor
# Return:
# - ADBm = ADB matrix of plane laminate

    # Laminate definition
    Nplies = len(ply_angles) # number of plies
    tot_thick = sum(ply_thick) # total laminate thickness

    # Compute stiffness and compliance tensor for each layer
    layers_mat = np.empty(Nplies, dtype=object)
    layers_comp = np.empty(Nplies, dtype=object)
    for i in range(Nplies):
        [layers_comp[i], type] = convert_mat_prop(ply_materials[i],'a','c','ca') # compliance matrix
        [layers_mat[i], type] = convert_mat_prop(ply_materials[i],'o','c','e') # engineering constants
        layers_mat[i] = np.array(layers_mat[i])
    
    # Arrays initialization
    comp_mat = np.empty(Nplies, dtype=object)
    ps_comp_mat = np.empty(Nplies, dtype=object)
    ps_stiff_mat = np.empty(Nplies, dtype=object)
    gl_stiff_mat = np.empty(Nplies, dtype=object)
    zbar = np.empty(Nplies)
    A = np.empty([3,3])
    B = np.empty([3,3])
    D = np.empty([3,3])
    
    for i in range(Nplies):
        [comp_mat[i], type] = convert_mat_prop(layers_mat[i],'o','e','s') # local compliance matrix (fiber coord. sys)
        ps_comp_mat[i] = [[comp_mat[i][0][0], comp_mat[i][0][1], 0], [comp_mat[i][1][0], comp_mat[i][1][1], 0], [0, 0, comp_mat[i][5][5]]] # local reduced (plane stress) compliance matrix
        ps_stiff_mat[i] = np.linalg.inv(ps_comp_mat[i]) # local reduced (plane stress) stiffness matrix
        c = np.cos(ply_angles[i]*math.pi/180) # cos of ply angle [rad]
        s = np.sin(ply_angles[i]*math.pi/180) # sin of ply angle [rad]
        T1 = [[c*c, s*s, 2*c*s], [s*s, c*c, -2*c*s], [-c*s, c*s, c*c-s*s]] # Rotation matrix T1
        T2 = [[c*c, s*s, c*s], [s*s, c*c, -c*s], [-2*c*s, 2*c*s, c*c-s*s]] # Rotation matrix T2
        gl_stiff_mat[i] = np.matmul(np.matmul(np.linalg.inv(T1),ps_stiff_mat[i]),T2) # global reduced stiffness matrix (laminate coord. sys)
        zbar[i] = (tot_thick + ply_thick[i])/2 - sum(ply_thick[0:i+1]) # ply coordinate
        
        A = A + gl_stiff_mat[i] * ply_thick[i] # A matrix
        B = B + gl_stiff_mat[i] * ply_thick[i] * zbar[i] # B matrix
        D = D + gl_stiff_mat[i] * ply_thick[i] * (zbar[i]**2  + ply_thick[i]**2 / 12) # D matrix
        
    ADBm = [[A, B], [B, D]] # ADB matrix

    return ADBm
# End: def ABDmatrix

def layer_mat_prop(moisture,runpath,rule,scalelig=1.):
# Get material properties (stiffness matrix and array of hygroexpansion coefficients) to give to Abaqus as input
# Args:
# - moisture = moisture content (can be an vector)
# - runpath = path where to find the material properties data file
# - rule = choose a set of composite materials rule:
#   'mr2': mixing-rule for ML, P, S1, S2, S3 (HC+L -> +C)
#   'mht': Halpin-Tsai for CML, CML, S1, S2, S3 (C+HC -> +L)
#   'mht_full': full Halpin-Tsai for ML, P, S1, S2, S3 (C+HC -> +L)
#   'mr2_full': full mixing-rule for ML, P, S1, S2, S3 (HC+L -> +C)
# Return:
# - layers_stiff = array of stiffness matrices (one for each layer: ML, P, S1-, S1+, S2, S3)
# - layers_hygro = array of hygroexpansion coefficients arrays (one for each layer: ML, P, S1-, S1+, S2, S3)
 
   layers_mat = get_current_layer_materials(moisture,runpath,rule,scalelig) # get layers materials (9 eng. constants + 3 hygroexp. coeff.)
   
   # Add S1 layer with opposite MFA
   nl = len(layers_mat) # number of layers
   layers_mat = [layers_mat[0],layers_mat[1],layers_mat[2],layers_mat[3],layers_mat[4],[]] #np.zeros(len(layers_mat[0]))
   temp = np.empty(4, dtype=object)
   temp[0] = layers_mat[2]
   for i in range(3,nl+1):
       temp[i-2] = layers_mat[i]
       layers_mat[i] = temp[i-3]
   #end for
   
   # Compute stiffness tensor for each layer
   layers_stiff = [[]]*(nl+1) #np.empty(nl+1, dtype=object)
   for i in range(nl+1):
       layers_stiff[i] = convert_mat_prop(layers_mat[i][0:9],'o','e','c')[0]; # stiffness matrix
   #end for
   # Get hygroexpansion coefficient fot each layer
   layers_hygro = np.empty(nl+1, dtype=object)
   for i in range(nl+1):
       layers_hygro[i] = layers_mat[i][9:12]

   return [layers_stiff, layers_hygro]
# End: def ayer_mat_prop

def get_current_layer_materials(moisture,runpath,rule,scalelig=1.):
# Get the current materials properties for the different layers for different level of moisture content.
# Args:
# - moisture = moisture content (can be an vector)
# - runpath = path where to find the material properties data file
# - rule = choose a set of composite materials rule:
#   'mr2': mixing-rule for ML, P, S1, S2, S3 (HC+L -> +C)
#   'mht': Halpin-Tsai for CML, CML, S1, S2, S3 (C+HC -> +L)
#   'mht_full': full Halpin-Tsai for ML, P, S1, S2, S3 (HC+C -> +L)
#   'mr2_full': full mixing-rule for ML, P, S1, S2, S3 (HC+L -> +C)
# Return:
# - materials = array of m kx12 matrices:
#               - m = number of elements in moisture, each matrix kx12 corresponds to a value of moisture
#               - 12 engineering constants for each material in a row (9 elastic constants + 3 hygroexpansion coefficients)
#               - k = number of materials

    # Set rule to calculate composite materials
    if rule=='mr2':
            matidx = [8,9,10,11,12]
    if rule=='mht':
            matidx = [14,14,16,18,20]
    if rule=='mht_full':
            matidx = [36,37,38,39,40]    # Hering 17.01.2012
    if rule=='mr2_full':
            matidx = [26,27,28,29,30]    # Hering 17.01.2012
    if (rule!='mr2' and rule!='mht' and rule!='mht_full' and rule!='mr2_full'):
            raise ValueError ('Unknown mode')

    m = get_materials(runpath,moisture,scalelig)
    #print(m)
    materials = m[matidx,:]

    return materials
# End: def get_current_layer_materials

def get_materials(runpath,moisture,scalelig=1.):
# Calculate material properties for different layers corresponding to a given moisture content
# Args:
# - runpath = path where to find the material properties data file
# - moisture = moisture content at which the properties must be calculated
# Return:
# - mats = engineering constants and hygroexpansion coefficients of the materials (n x (9+3) matrix)
#   Structure of the mats cell vector:
#   1: cellulose (c)
#   2: hemicellulose (h)
#   3: lignin (l)
#   4: composite materials
#      ...

    # Load material properties (Stiffness matrices + moisture dependency + hygroexpansion coefficients)
    materials_data = np.load(runpath+'basic_components.npy', allow_pickle=True, encoding='latin1') # basic components
    elastic_prop = [materials_data[0], materials_data[1], materials_data[2]*scalelig]
    n =  len(elastic_prop) # num. of chemical consituents
    hygroexp_coeff = materials_data[3] # nx3 cell, each row contain the data about the dependency of the
    # elastic constants wrt the moisture (c(ub), ub(w)), and the hygroexpansion coefficients
    
    # Load composite material definitions
    materials = np.load(runpath+'materials_composite.npy', allow_pickle=True, encoding='latin1',fix_imports=True) # composite materials
    
    if(len(materials[0])<13): #if columns < 13
        sz=[len(materials), len(materials[0])]
        for i in range(sz[0]):
            materials[i]=np.pad(materials[i], (0, 13-sz[1]), 'constant') #np.concatenate(materials[i],np.zeros(14-sz[1])) # fill remaining columns with zero
    #end if

    ub = boundwater(moisture)
    c_deps = moistcoeff(ub)

    for i in range(n): # for each chemical constituent
        [eng_const,type] = convert_mat_prop(elastic_prop[i]*c_deps[i],'o','c','e')
        materials[i] = np.concatenate(([-9], eng_const, hygroexp_coeff[i][0]), axis=0)
                    # [-9 = orthotropic material,
                    #  stiff. matrix * moisture-dependent function -> E, G, nu,
                    #  hygroexpansion coefficients]
    #end for
    
    mats = compute_materials(materials) # compute 9 material constants + 3 hygroexpansion coefficients

    return mats
# End: def get_materials

def ub_hc(moisture):
# Calculate fraction of bound water content of hemicellulose
# Args:
# - moisture = moisture content
# Return:
# - ub_hc = fraction of bound water of hemicellulose
    ub = 2.6/3.6*0.2*((1-math.exp(-moisture*0.9/0.2))/(1-math.exp(-0.3*0.9/0.2)))
    ub = ub/(2.6/3.6*0.2*((1-math.exp(-0.28*0.9/0.2))/(1-math.exp(-0.3*0.9/0.2)))) #modified equation
    return ub
# End: ub_hc

def ub_l(moisture):
# Calculate fraction of bound water content of lignin
# Args:
# - moisture = moisture content
# Return:
# - ub_l = fraction of bound water of lignin
    ub = 1/3.6*0.2*((1-math.exp(-moisture*0.9/0.2))/(1-math.exp(-0.3*0.9/0.2)))
    ub = ub/(1/3.6*0.2*((1-math.exp(-0.28*0.9/0.2))/(1-math.exp(-0.3*0.9/0.2)))) #modified equation
    return ub
# End: ub_l

def c_hc(ub):
# Calculate stiffness matrix coefficient of hemicellulose
# Args:
# - ub = fraction of bound water
# Return:
# - cb_hc = stiffness matrix coefficient of hemicellulose
    #c = 6/2*(0.1+0.9*(1-ub))
    c = 3/2*(0.1+0.9*(1-ub)) #modified equation (h=3)
    return c
# End: c_hc

def c_l(ub):
# Calculate stiffness matrix coefficient of lignin
# Args:
# - ub = fraction of bound water
# Return:
# - cb_hc = stiffness matrix coefficient of lignin
    #c = (1+2*(1-ub))/2
    c = (448.33*pow(ub,8)-2291.3*pow(ub,7)+4756.1*pow(ub,6)-5120.1*pow(ub,5)+2997.5*pow(ub,4)-891.76*pow(ub,3)+100.93*pow(ub,2)-1.0292*ub+2.3744)*0.7 #modified equation
    return c
# End: c_l

def boundwater(moisture):
# Calculate fraction of bound water content partitioned between cellulose, hemicellulose, and lignin
# Args:
# - moisture = moisture content
# Return:
# - ub = fractions of bound water
    ub = [0, ub_hc(moisture), ub_l(moisture)]
    return ub
# End: boundwater

def moistcoeff(ub):
# Calculate stiffness matrix coefficient of cellulose, hemicellulose, lignin
# Args:
# - ub = fraction of bound water
# Return:
# - cb = stiffness matrix coefficients
    c_deps = [1, c_hc(ub[1]), c_l(ub[2])]
    return c_deps
# End: moistcoeff

def compute_materials(data):
# Compute materials engineering constants.
# Args:
# - data = n x (10+3) matrix where
#           - each row corresponds to a material
#           - the 1st element of each row determines how the rest is interpreted
#           - 9 engineering constants (E, G, nu)
#           - 3 hygroexpansion coefficients
#  The 1st element can assume one of the following values:
#   -2: the material is isotropic: the next 2 numbers are E and G. The rest of the row is ignored
#   -5: the material is orthotropic with transversal isotropy: the next 5 numbers are E1, E2, nu12, G12, G23. The rest of the row is ignored
#   -9: the material is orthotropic: the rest of the row is E1, E2, E3, nu23, nu13, nu12, G23, G13, G12
#   -1: the material is composed of two other materials: one for the fibers and the other for the matrix; it is assumed that the matrix and the
#       fiber material are both orthotropic, transversally isotropic wrt the same orientation: the next 4 numbers are the number of the material of
#       the fibers, the volume fraction of the fibers, the number of the material of the matrix and the volume fraction of the matrix
#   -10: the material is composed of two other materials: one for the fibers and the other for the matrix; it is assumed that the matrix and the
#        fiber material are both orthotropic, transversally isotropic wrt the same orientation: the next 6 numbers are the number of the material of
#        the fibers, the volume fraction of the fibers, the number of the material of the matrix and the volume fraction of the matrix and 2
#        geometry coefficients: for E and for G. The computation is done with Halpin-Tsai formulae.
#   a number between 0 and 1: the material constants are computed with the constants of other materials. Up to 5 materials can be combined this
#                             way. Each couple of numbers correspond to a material. The first number is the volume fraction and the second the number of
#                             the material (the row number in the material matrix, starting at 1). The materials are assumed to be isotropic.
# Return:
# - mm = n x 12 matrix containing:
#               - 9 engineering constants (E, G, nu)
#               - 3 hygroexpansion coeffients

    n = len(data) # num. rows = num. of materials
    mm = np.zeros((n,12)) # initialize output matrix

    for i in range(n): # for each material
        row = data[i] # row corresponding to a material
        # 1st element of data define how the rest is interpreted
        if row[0]==-2: # isotropic material
            E = row[1]
            G = row[2]
            nu = (E-2*G)/(2*G)
            values = np.concatenate(([E,E,E,nu,nu,nu,G,G,G],row[10:13]), axis=0)

        if row[0]==-5: # orthotropic with transversal isotropy material
            E1 = row[1]
            E2 = row[2]
            nu12 = row[3]
            G12 = row[4]
            G23 = row[5]
            G13 = G12
            E3 = E2
            nu13 = nu12
            nu23 = (E2-2*G23)/(2*G23)
            values = np.concatenate(([E1,E2,E3,nu23,nu13,nu12,G23,G13,G12],row[10:13]), axis=0)

        if row[0]==-9: # orthotropic material
            values = [row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12]]

        if row[0]==-1: # matrix+fibers orthotropic, transversally isotropic -> use get_E_nu_G
            matf = mm[int(row[1])-1]
            matm = mm[int(row[3])-1]
            vf = row[2]
            vm = row[4]
            vf = vf/(vf+vm)
            vm = 1-vf
            r = get_E_nu_G(matf[0],matf[1],matm[0],matm[1],matf[5],matm[5],matf[8],matf[6],matm[8],matm[6],vf,vm)
            epsl = (vf*matf[0]*matf[9]+vm*matm[0]*matm[9])/(vf*matf[0]+vm*matm[0])
            epst = matf[10]*math.sqrt(vf)+(1-math.sqrt(vf))*(1+vf*matm[10]*matf[0])/(vf*matf[0]+vm*matm[0])*matm[10] # changed Hering 17.01.2012 - C.C.Chamis Trans. thermal expansion
            values = [r[0],r[1],r[1],r[3],r[2],r[2],r[5],r[4],r[4],epsl,epst,epst]

        if row[0]==-10: # matrix+fibershygroexp_coeff[3][0] orthotropic, transversally isotropic -> use halpin_tsai2
            xe = row[5]
            xg = row[6]
            matf = mm[int(row[1])-1]
            matm = mm[int(row[3])-1]
            vf = row[2]
            vm = row[4]
            vf = vf/(vf+vm)
            vm = 1-vf
            r = halpin_tsai2(matf[0],matf[1],matm[0],matm[1],matf[5],matm[5],matf[8],matf[6],matm[8],matm[6],vf,vm,xe,xg)
            epsl = (vf*matf[0]*matf[9]+vm*matm[0]*matm[9])/(vf*matf[0]+vm*matm[0])
            epst = matf[10]*math.sqrt(vf)+(1-math.sqrt(vf))*(1+vf*matm[5]*(matf[0]/(vf*matf[0]+vm*matm[0])))*matm[10] # changed Hering 17.01.2012 - C.C.Chamis Transv. thermal expansion
            if(xg==0):
                epsltlt = 0
                for j in range(0,359,10):
                    angle = j*math.pi/180
                    epsltlt = epsltlt+(abs(epsl*r[0]*math.cos(angle))+abs(epst*r[2]*math.sin(angle)))/(r[0]+r[2]); # inserted by SH 19.01.2012 to average main swelling directions
                #end for
                epslt = epsltlt/36
                values = [r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],epslt,epslt,epst]
            else:
                values = [r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],epsl,epst,epst]
            #end if

        if row[0]==-11: # elastic matrix material + composite moisture dependency
            xe = row[5]
            xg = row[6]
            matf = mm[row[1]-1]
            matm = mm[row[3]-1]
            vf = row[2]
            vm = row[4]
            vf = vf/(vf+vm)
            vm = 1-vf
            r = mm[row[3],0:9]
            epsl = (vf*matf[0]*matf[9]+vm*matm[0]*matm[9])/(vf*matf[0]+vm*matm[0])
            epst = matf[10]*math.sqrt(vf)+(1-math.sqrt(vf))*(1+vf*matm[5]*(matf[0]/(vf*matf[0]+vm*matm[0])))*matm[10] # changed Hering 17.01.2012 - C.C.Chamis Transv. thermal expansion
            if(xg==0):
                epsltlt = 0
                for j in range(0,359,10):
                    angle = j*math.pi/180
                    epsltlt = epsltlt+(abs(epsl*r[0]*math.cos(angle))+abs(epst*r(3)*math.sin(angle)))/(r[0]+r[2]) # inserted by SH 19.01.2012 to average main swelling directions
                #end for
                epslt = epsltlt/36
                values = [r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],epslt,epslt,epst]
            else:
                values = [r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],epsl,epst,epst]
            #end if

        if (row[0]>0 and row[0]<1):
            c1 = mm[max(1,row[1])] # constants for the first material
            v = row[0]
            for j in range (3,10,2):
                c2 = mm[max(1,row[j]),:]
                v1 = v/(v+row[j-1])
                v2 = 1-v1
                c1[0] = c1[0]*v1+c2[0]*v2
                c1[6] = c1[6]*v1+c2[6]*v2
                c1[10:13] = c1[10:13]*v1+c2[10:13]*v2
                v = v/v1
            #end for
            nu = 0.5*c1[0]/c1[6]-1
            values = [c1[0],c1[0],c1[0],nu,nu,nu,c1[6],c1[6],c1[6],c1[10:13]]
        #end if
        mm[i] = values
    #end for

    return mm
# End: def compute_materials

def halpin_tsai2(Efl,Eft,Eml,Emt,nuf,num,Gfl,Gft,Gml,Gmt,Vf,Vm,xsi_e,xsi_g):
# Compute the stiffness matrix of the composite via Halpin-Tsai equations taking into account geometrical factors
# Args:
# - Efl = Young's modulus of fiber longitudinally to the fibers
# - Eft = Young's modulus of fiber transversally to the fibers
# - Eml = Young's modulus of matrix longitudinally to the fibers
# - Emt = Young's modulus of matrix transversally to the fibers
# - nuf = Poisson's modulus of fiber
# - num = Poisson's modulus of matrix
# - Gfl = shear modulus of fiber with longitudinal and transversal directions
# - Gft = shear modulus of fiber with both transversal directions
# - Gml = shear modulus of matrix longitudinal and transversal directions
# - Gmt = shear modulus of matrix with both transversal directions
# - Vf = Volume fraction of fiber
# - Vm = Volume fraction of matrix
# - xsi_e = geometrical factor for the Young modulus -> if xsi_e=2 and xsi_g=0: average in all directions (random orientation)
# - xsi_g = geometrical factor for the shear modulus
#        -> if xsi_g=0: the properties are averaged around direction 3
#        -> if xsi_g=1: the properties are averaged around direction 1 (in general not necessary)
# Return:
# - E = array containing longitudinal and transversal Young's moduli, Poisson's ratios and shear moduli of the material
# composed of the fiber and the matrix (6 elements)

    e_f = [Efl,Eft,nuf,np.NaN,Gfl,Gft] # fiber material constants
    e_m = [Eml,Emt,num,np.NaN,Gml,Gmt] # matrix material constants
    # Compute matrices
    [c_f,type_f] = convert_mat_prop(e_f,'t','e','c') # t = transv. isotropic, e = engineering constants, c = stiffness matrix
    [c_m,type_m] = convert_mat_prop(e_m,'t','e','c') # t = transv. isotropic, e = engineering constants, c = stiffness matrix

    c = halpin_tsai_hill(c_m,c_f,Vf) # composite stiffness matrix
    if (xsi_g==0): #average around the third axis
        cc = np.zeros((6,6))
        for i in range(0,359,10):
            cc = cc + transform_stiffness_matrix(c,i,3)
        # end for
        cc[[0,0,0,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5],[3,4,5,3,4,5,3,4,5,0,1,2,4,5,0,1,2,3,5,0,1,2,3,4]] = 0.0
        c = cc/36

        if (xsi_e==2): # average in all directions
            cct = [transform_stiffness_matrix(c,90,2),c,c]
            cc = np.zeros((6,6))
            for i in range(0,179,10):
                for j in range (3):
                    cc = cc + transform_stiffness_matrix(cct[j],i,j+1)
                # end for
            # end for
            cc[[0,0,0,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5],[3,4,5,3,4,5,3,4,5,0,1,2,4,5,0,1,2,3,5,0,1,2,3,4]] = 0.0
            c = cc/54
        # end if
    # end if
    
    if(xsi_g==1): #average around the first axis - in general not necessary - due to the nature aligned fiber
        cc = np.zeros((6,6))
        for i in range(0,359,10):
            cc = cc + transform_stiffness_matrix(c,i,1) # direction changed to 1 SH 19.01.2012 --> additional option for aligned fibers
        # end for
        cc[[0,0,0,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5],[3,4,5,3,4,5,3,4,5,0,1,2,4,5,0,1,2,3,5,0,1,2,3,4]] = 0.0;
        c = cc/36
    # end if
    
    if(xsi_g==2): #
        for rep in range(3):
            for k in range(3):
                cc=np.zeros((6,6))
                for i in range(0,359,10):
                    cc = cc + transform_stiffness_matrix(c,i,3-k)
                    # end for
                cc[[0,0,0,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5],[3,4,5,3,4,5,3,4,5,0,1,2,4,5,0,1,2,3,5,0,1,2,3,4]] = 0.0;
                c = cc/36
    # end if

    #INVC = np.linalg.inv(c) # compliance matrix
    [E,type_E]=convert_mat_prop(c,'o','c','e') # engineering constants

    return E
# End: def halpin_tsai2

def halpin_tsai_hill(c_matrix,c_fiber,v_f):
# Compute the stiffness matrix of the composite via Halpin-Tsai equations
# Args:
# - c_matrix = stiffness matrix (6x6) of the matrix
# - c_fiber = stiffness matrix (6x6) of the fiber
# - v_f = volume fraction of fibers
# Return:
# - c = stiffness matrix of the composite

    v_m = 1-v_f # Hill Notation
    kf = (c_fiber[1,1]+c_fiber[1,2])/2
    km = (c_matrix[1,1]+c_matrix[1,2])/2
    mf = (c_fiber[1,1]-c_fiber[1,2])/2
    mm = (c_matrix[1,1]-c_matrix[1,2])/2
    uf = c_fiber[4,4]
    um = c_matrix[4,4]
    nf = c_fiber[0,0]
    nm = c_matrix[0,0]
    lf = c_fiber[0,1]
    lm = c_matrix[0,1]
    zk = mm/km
    Mk = kf/km
    ek = (Mk-1)/(Mk+zk)
    k = km*(1+zk*ek*v_f)/(1-ek*v_f)
    zm = (km/mm)/(km/mm+2)
    Mm = mf/mm
    em = (Mm-1)/(Mm+zm)
    m = mm*(1+zm*em*v_f)/(1-em*v_f)
    zu = 1
    Mu = uf/um
    eu = (Mu-1)/(Mu+zu)
    u = um*(1+zu*eu*v_f)/(1-eu*v_f)
    nu12f = lf/2/kf
    nu12m = lm/2/km
    nu12 = nu12f*v_f+nu12m*v_m
    l = nu12*2*k
    e11f = nf-lf/2/kf
    e11m = nm-lm/2/km
    e11 = e11f*v_f+e11m*v_m
    n = e11+l/2/k
    c=np.zeros((6,6))
    c[[1,2,2],[0,0,1]] = [l,l,k-m] # left out-of-diagonal elements
    c = c + c.transpose() # left out-of-diagonal elements + right out-of-diagonal elements
    np.fill_diagonal(c, [n,k+m,k+m,m,u,u]) # diagonal elements

    return c
# End: def halpin_tsai_hill

def convert_mat_prop(data,symm,type1,type2):
# Convert material properties from one representation to another one
# Args:
# - data = matrix containing the data
# - symm = symmetry of the material:
#      - a (anisotropic)
#      - i (isotropic)
#      - t (transversally isotropic: direction 1 is different)
#      - o (orthotropic)
#      - g (orthotropic generalized plane strain)
# - type1, type2 = data types before and after conversion:
#      - ca (stiffness matrix for abaqus)
#      - c (stiffness matrix) -> 6x6 matrix, c44=c2323
#      - s (compliance matrix)
#      - e (engineering constants) -> [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12] if orthotropic, [E1,E2,nu12,nu23,G12,G23] if transversally isotropic
#      - ea (eng. const. for abaqus)
# Return:
# - data = converted data
# - types = cell array containing the types that have been used (the first element is type1 and the last type2, there may be other
# types in the middle.

    types = list(type1)
    c = [[], [], []]
    conv_names = np.array(['e2s','s2e','c2s','s2c','ca2c','c2ca','ea2e','e2ea'])
    while(not types[-1]==type2):
        c[0] = [types[-1],'2',type2];
        c[1] = [types[-1],'2','c'];
        c[2] = [types[-1],'2','s'];
        found = False
        idx = -1
        for i in range (1,len(c)+1):
            for j in range(len(conv_names)):
                if conv_names[j]==''.join(c[i-1]):
                    idx = j
            if (idx!=-1):
                [data,temp] = convert_to(data,symm,conv_names[idx])
                types.append(temp)
                found = True
                break
            #end if
        #end for
        if(not found):
            raise ValueError('Cannot find how to convert')
        #end if
    #end while
    return [data,types]
# End: def convert_mat_prop

def convert_to(data,symm,conv):
# Convert material properties from one representation to another one
# Args:
# - data = matrix containing the data
# - symm = symmetry of the material:
#      - a (anisotropic)
#      - i (isotropic)
#      - t (transversally isotropic: direction 1 is different)
#      - o (orthotropic)
#      - g (orthotropic generalized plane strain)
# - conv = conversion name
# Return:
# - data2 = converted data
# - t = cell array containing the types that have been used

    if conv=='ca2c':
        t = 'c'
        if(symm=='a'):
            c=np.zeros((6,6))
            c[0:3,0:3] = data[[0,1,3], [1,2,4], [3,4,5]]
            c[0:3,3:6] = data[[15,10,6],[16,11,76],[17,12,8]]
            c_temp = c[0:3,3:6]
            c[3:6,0:3] = c_temp.transpose()
            c[3:6,3:6] = data[[20,19,20],[19,14,13],[18,13,9]]
            data2 = c
        else:
            raise ValueError('Conversion not implemented for that symmetry')
    # end if 'ca2c'
    if conv=='c2ca':
        t = 'ca';
        if(symm=='a'):
            ca = np.zeros((21))
            ca[0:6] = data[[0,0,1,0,1,2],[0,1,1,2,2,2]]
            ca[6:10] = data[[0,1,2,5],[5,5,5,5]]
            ca[10:15] = data[[0,1,2,5,4],[4,4,4,4,4]]
            ca[15:21] = data[[0,1,2,5,4,3],[3,3,3,3,3,3]]
            data2 = ca
        else:
            raise ValueError('Conversion not implemented for that symmetry')
        # end if
    # end if 'c2ca'
    if conv=='c2s':
        t = 's'
        data2 = np.linalg.inv(data)
    # end if 'c2s'
    if conv=='s2c':
        t = 'c'
        data2 = np.linalg.inv(data)
    # end if 's2c'
    if conv=='s2e':
        t = 'e'
        if(symm=='t'):
            E1 = 1/data[0,0]
            E2 = 1/data[1,1]
            nu12 = -data[1,0]*E1
            nu23 = -data[2,1]*E2
            G12 = 1/data[5,5]
            G23 = 1/data[3,3]
            data2 = [E1,E2,nu12,nu23,G12,G23]
        elif(symm=='o'):
            E1 = 1/data[0,0]
            E2 = 1/data[1,1]
            E3 = 1/data[2,2]
            nu12 = -data[1,0]*E1
            nu23 = -data[2,1]*E2
            nu13 = -data[2,0]*E1
            G12 = 1/data[5,5]
            G13 = 1/data[4,4]
            G23 = 1/data[3,3]
            data2 = [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12]
        elif(symm=='g'):
            E1 = 1/data[0,0]
            E2 = 1/data[1,1]
            E3 = 1/data[2,2]
            nu12 = -data[1,0]*E1
            nu23 = -data[2,1]*E2
            nu13 = -data[2,0]*E1
            G12 = 1/data[5,5]
            data2 = [E1,E2,E3,nu23,nu13,nu12,G12]
        else:
            raise ValueError('Convertion not implemented for that symmetry')
    # end if 's2e'
    if conv=='e2s':
        t = 's'
        if(symm=='t'):
            if(np.sum(np.isnan(data))>1):
                raise ValueError('Not enough parameters given')
            data2 = np.zeros((6,6))
            if(np.isnan(data[1])):
                data[1] = 2*data[5]*(1+data[3])
            if(np.isnan(data[3])):
                data[3] = data[1]/data[4]/2-1
            data2[[0,0,1],[1,2,2]] = [-data[2]/data[0],-data[2]/data[0],-data[3]/data[1]]
            temp = data2[0:3,0:3]
            data2[0:3,0:3] = temp+temp.transpose()
            data2[[0,1,2],[0,1,2]] = [1/data[0],1/data[1],1/data[1]];
            data2[[3,4,5],[3,4,5]] = [2*(1+data[3])/data[1],1/data[4],1/data[4]];
        elif(symm=='o'):
            data2=np.zeros((6,6))
            data2[[0,0,1],[1,2,2]] = [-data[5]/data[0],-data[4]/data[0],-data[3]/data[1]]
            temp = data2[0:3,0:3]
            data2[0:3,0:3] = temp+temp.transpose()
            data2[[0,1,2],[0,1,2]] = [1/data[0],1/data[1],1/data[2]]
            data2[[3,4,5],[3,4,5]] = [1/data[6],1/data[7],1/data[8]]
        elif(symm=='g'):
            data2=np.zeros((4,4))
            data2[[0,0,1],[1,2,2]] = [-data[5]/data[0],-data[4]/data[0],-data[3]/data[1]]
            temp = data2[0:3,0:3]
            data2[0:3,0:3] = temp+temp.transpose()
            data2[[0,1,2],[0,1,2]] = [1/data[0],1/data[1],1/data[2]]
            data2[[3],[3]] = [1/data[6]]
        else:
            raise ValueError('Convertion not implemented for that symmetry')
    #end if 'e2s'
    if (conv!='e2s' and conv!='s2e' and conv!='c2s' and conv!='s2c' and conv!='ca2c' and conv!='c2ca'):
        raise ValueError('Unknown conversion')
    #end if
    return [data2,t]
# End: def convert_to

def get_E_nu_G(Efl,Eft,Eml,Emt,nuf,num,Gfl,Gft,Gml,Gmt,Vf,Vm):
# Calculate engineering constants of the composite material for given engineering constants
# of fiber and matrix
# Args:
# - Efl = Young's modulus of fiber longitudinally to the fibers
# - Eft = Young's modulus of fiber transversally to the fibers
# - Eml = Young's modulus of matrix longitudinally to the fibers
# - Emt = Young's modulus of matrix transversally to the fibers
# - nuf = Poisson's modulus of fiber
# - num = Poisson's modulus of matrix
# - Gfl = shear modulus of fiber with longitudinal and transversal directions
# - Gft = shear modulus of fiber with both transversal directions
# - Gml = shear modulus of matrix longitudinal and transversal directions
# - Gmt = shear modulus of matrix with both transversal directions
# - Vf = Volume fraction of fiber
# - Vm = Volume fraction of matrix
# Return:
# - E = array containing longitudinal and transversal Young's moduli, Poisson's ratios and
# shear moduli of the material composed of the fiber and the matrix
    
    nuflt = nuf
    nuftt = (Eft-2*Gft)/(2*Gft)
    numlt = num
    numtt = (Emt-2*Gmt)/(2*Gmt)
    El = Efl*Vf+Eml*Vm
    Et = Eft*Emt/(Eft*Vm+Emt*Vf)
    nult = nuflt*Vf+numlt*Vm
    Glt = Gfl*Gml/(Gfl*Vm+Gml*Vf)
    Gtt = Gft*Gmt/(Gft*Vm+Gmt*Vf)
    nutt = (Et-2*Gtt)/(2*Gtt)
    E = [El,Et,nutt,nult,Gtt,Glt]

    return E
# End: def get_E_nu_G

def transform_stiffness_matrix(c,angle_deg,direction):
# Transform the stiffness matrix (6x6) by rotating the coordinate system around the given
# axis by a given angle.
# Args:
# - c = stiffness matrix expressed in local coordinate system
# - angle_deg = angle of rotation (degrees)
# - direction = direction around which the system is rotated (1,2 or 3)
# Return:
# - c_rot = stiffness matrix expressed in th global coordinate system

    T = get_transformation_matrix(angle_deg,direction) # transformation matrix
    c_rot = np.dot(np.dot(T,c),T.transpose()) # stiffness matrix

    return c_rot
# End: def transform_stiffness_matrix

def get_transformation_matrix(angle_deg,direction):
# Calculate the transformation matrix for matrices 6x6 for a given angle and direction.
# Args:
# - angle_deg = angle of rotation (degrees)
# - direction = direction around which the system is rotated (1,2 or 3)
# Return:
# - T = transformation matrix

    angle = angle_deg*math.pi/180 # angle in radiant
    c = math.cos(angle) # cosine
    s = math.sin(angle) # sine
    cs = c*s
    cc = c*c
    ss = s*s
    T = np.zeros((6,6))
    T[0:2,0:2] = ([[cc, ss], [ss, cc]])
    T[2,2] = 1
    T[3:5,3:5] = ([[c,s], [-s,c]])
    T[5,0:2] = ([-cs, cs])
    T[0:2,5] = ([2*cs,-2*cs])
    T[5,5] = cc-ss
    if (direction<3):
        if (direction==1):
            t = T[[2,0,1,5,4,3],:]
            T = t[:,[2,0,1,5,4,3]]
        elif (direction==2):
            t = T[[0,2,1,4,5,3],:]
            T = t[:,[0,2,1,4,5,3]]
        # end if
    # end if

    return T
# End: def get_transformation_matrix
