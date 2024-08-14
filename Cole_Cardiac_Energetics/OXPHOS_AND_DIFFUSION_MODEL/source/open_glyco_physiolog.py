# Intialise OpenCMISS-Iron
from opencmiss.iron import iron
import numpy as np
import math
#import meshio
import time

starttime = time.time()

# Program Parameters

CoordinateSystemUserNumber=1
RegionUserNumber=2
BasisUserNumber=3
MeshUserNumber=4
DecompositionUserNumber=5
GeometricFieldUserNumber=6

ATPFieldUserNumber=7
ATPMaterialsFieldUserNumber=8
ATPEquationsSetUserNumber=9
ATPEquationsSetFieldUserNumber=10
ATPSourceFieldUserNumber=17

ADPFieldUserNumber=11
ADPMaterialsFieldUserNumber=12
ADPEquationsSetUserNumber=13
ADPEquationsSetFieldUserNumber=14
ADPSourceFieldUserNumber=23

AMPFieldUserNumber=24
AMPMaterialsFieldUserNumber=25
AMPEquationsSetUserNumber=26
AMPEquationsSetFieldUserNumber=27
AMPSourceFieldUserNumber=28

PCrFieldUserNumber=29
PCrMaterialsFieldUserNumber=30
PCrEquationsSetUserNumber=31
PCrEquationsSetFieldUserNumber=32
PCrSourceFieldUserNumber=33

CrFieldUserNumber=34
CrMaterialsFieldUserNumber=35
CrEquationsSetUserNumber=36
CrEquationsSetFieldUserNumber=37
CrSourceFieldUserNumber=38

PiFieldUserNumber=39
PiMaterialsFieldUserNumber=40
PiEquationsSetUserNumber=41
PiEquationsSetFieldUserNumber=42
PiSourceFieldUserNumber=43

OxyFieldUserNumber=69
OxyMaterialsFieldUserNumber=70
OxyEquationsSetUserNumber=71
OxyEquationsSetFieldUserNumber=72
OxySourceFieldUserNumber=73

DPsiFieldUserNumber=44
DPsiMaterialsFieldUserNumber=45
DPsiEquationsSetUserNumber=46
DPsiEquationsSetFieldUserNumber=47
DPsiSourceFieldUserNumber=48

ProblemUserNumber=15
ControlLoopNode=0
AnalyticFieldUserNumber=16
CellMLUserNumber=18
CellMLModelsFieldUserNumber=19
CellMLStateFieldUserNumber=20
CellMLIntermediateFieldUserNumber=21
CellMLParametersFieldUserNumber=22
BCFieldUserNumber=74

### Input Parameters ###
input_parameters = open('input.txt', 'r')

# Reading the input parameters
#node_file, elem_file = input_parameters.readline().split()
simplex_order = input_parameters.readline().split()

init_ATP, ATPDiffx, ATPDiffy, ATPDiffz = input_parameters.readline().split()
init_ADP, ADPDiffx, ADPDiffy, ADPDiffz = input_parameters.readline().split()
init_AMP, AMPDiffx, AMPDiffy, AMPDiffz = input_parameters.readline().split()
init_PCr, PCrDiffx, PCrDiffy, PCrDiffz = input_parameters.readline().split()
init_Cr, CrDiffx, CrDiffy, CrDiffz = input_parameters.readline().split()
init_Pi, PiDiffx, PiDiffy, PiDiffz = input_parameters.readline().split()
init_Oxy, OxyDiffx, OxyDiffy, OxyDiffz = input_parameters.readline().split()
init_DPsi ,DPsiDiffx, DPsiDiffy, DPsiDiffz = input_parameters.readline().split()
MITOVOLS= input_parameters.readline().split()

input_parameters.close()

# Simulation Runtime Inputs
startT          = float(0.0)
#This following two lines represent a short time to test, which can be commented out later and a longer time that will go until past steady state
#endT 		= float(20)
endT            = float(25001.0)
#endT           = float(10)
Tstep           = float(1.0)
ODE_TIME_STEP   = float(0.01)         #Simulation Time Parameters
outputFreq      = int(5000)           #Output Frequence of Results
#outputFreq      = int(1)

# total number of loops
n_loops = int((endT - startT)/Tstep)
n_output = int(n_loops/outputFreq) + 1 # +1 is for initial conditions


### Initialize ###
#------------------------------------------------------------------------------------------------
# COMPUTATIONAL NODE INFORMATION (PARALELLIZATION)
#------------------------------------------------------------------------------------------------
ErrorHandlingModes = iron.ErrorHandlingModeSet(2)
# Get the total number of processors (value given to the -np flag)
ComputationalNumberOfNodes = iron.ComputationalNumberOfNodesGet()
# Get the rank of the current process. Will range from 0 to numberOfComputationalNodes
ComputationalNodeNumber = iron.ComputationalNodeNumberGet()

print ('NumberofComputationalNodes:', ComputationalNumberOfNodes)
print ('ComputationalNodeNumber:', ComputationalNodeNumber)

#------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
#------------------------------------------------------------------------------------------------

# Two Dimensional Coordinate System
CoordinateSystem = iron.CoordinateSystem()
CoordinateSystem.CreateStart(CoordinateSystemUserNumber)
CoordinateSystem.DimensionSet(2) #The coordinate system is 3D by default;set it to be 2D.
CoordinateSystem.CreateFinish()

#------------------------------------------------------------------------------------------------
# REGION
#------------------------------------------------------------------------------------------------

# Start Region
Region = iron.Region()
Region.CreateStart(RegionUserNumber, iron.WorldRegion)
Region.CoordinateSystemSet(CoordinateSystem)
Region.LabelSet("Cell")
Region.CreateFinish()

#------------------------------------------------------------------------------------------------
# BASIS
#------------------------------------------------------------------------------------------------

# Simplex Basis Reaction Diffusion - Creation of a trilinear-simplex basis
Basis = iron.Basis()
Basis.CreateStart(BasisUserNumber)
Basis.TypeSet(iron.BasisTypes.SIMPLEX)
Basis.NumberOfXiSet(2)
#Basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*2
print ('Simplex_Order:', simplex_order)
Basis.CreateFinish()


#------------------------------------------------------------------------------------------------
# MESH
#------------------------------------------------------------------------------------------------
### Create the mesh ###

# Reading node file
node_file = open('CA001.1.node', 'r')

# Reading the mesh details from the first line of the node file
number_of_nodes, number_of_coords, number_of_attributes, boundary_marker = node_file.readline().split()
number_of_nodes = int(number_of_nodes)
number_of_coords=int(number_of_coords)
number_of_attributes=int(number_of_attributes)
boundary_marker=int(boundary_marker)

# Creating variables to store node number & boundary marker
NodeNums = [[0 for m in range(2)] for n in range(number_of_nodes)]
# Creating array to store x and y coordinates
NodeCoords = [[0 for m in range(number_of_coords+1)] for n in range(number_of_nodes)]

# Reading details from node file
for i in range(number_of_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeCoords[i][2], NodeNums[i][1] = node_file.readline().split()  #node number, x, y, structure type, boundary marker
    # Converting from string to appropriate type
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])
    NodeCoords[i][2] = float(NodeCoords[i][2])
    NodeNums[i][1] = int(NodeNums[i][1])
    #This above line is the attribute of the nodes (glycogen, mitochondria, myofibrils, etc)
node_file.close()

#Input element file
elem_file = open('CA001.1.ele', 'r')

#Reading the values of the first line
number_of_elements, nodes_per_ele, ele_attributes = elem_file.readline().split()
number_of_elements = int(number_of_elements)
nodes_per_ele = int(nodes_per_ele)
ele_attributes = int(ele_attributes)

# Creating variable to store the element map
ElemMap = [[0 for x in range(nodes_per_ele+ele_attributes)] for y in range(number_of_elements)]
Elemindex = [[0 for m in range(1)] for n in range(number_of_elements)]

#elements_list = [[0 for i in range(4)] for elem in range(number_of_elements)]
# Reading element data from elemfile
for i in range(number_of_elements):
    # Performing the mapping
    Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2] = elem_file.readline().split()
#    elements_list[i] = [Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2]]
elem_file.close()

##########

#Initialise Nodes
Nodes = iron.Nodes()
Nodes.CreateStart(Region, number_of_nodes)
print('Number of Nodes:', number_of_nodes)
Nodes.CreateFinish()

# Initialise Mesh
Mesh = iron.Mesh()
Mesh.CreateStart(MeshUserNumber, Region, number_of_coords)
Mesh.NumberOfElementsSet(number_of_elements)
print('Number of Elements:', number_of_elements)
Mesh.NumberOfComponentsSet(1)

# Initialise Elements
MeshElements = iron.MeshElements()
MeshElements.CreateStart(Mesh, 1, Basis)

for i in range(number_of_elements):
    element=Elemindex[i][0]
    NodesList = list(
      map(int,[ElemMap[i][0], ElemMap[i][1], ElemMap[i][2]]))
    MeshElements.NodesSet(int(element), NodesList)

MeshElements.CreateFinish()

Mesh.CreateFinish()


#------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
#------------------------------------------------------------------------------------------------

# Parallelization
Decomposition = iron.Decomposition()
Decomposition.CreateStart(DecompositionUserNumber, Mesh)
Decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
Decomposition.NumberOfDomainsSet(ComputationalNumberOfNodes)
Decomposition.CreateFinish()

#------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
#------------------------------------------------------------------------------------------------

# Geometric Field - Start to create a default geometric field on the region
GeometricField = iron.Field()
GeometricField.CreateStart(GeometricFieldUserNumber, Region)
GeometricField.LabelSet('Geometry')
GeometricField.MeshDecompositionSet(Decomposition)
GeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
GeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
GeometricField.CreateFinish()

###  Set geometric values from customized mesh  ###
for i in range(number_of_nodes):
    node = NodeNums[i][0]
    NodeDomain = Decomposition.NodeDomainGet(node,1)
    if NodeDomain == ComputationalNodeNumber :
       node_x = NodeCoords[i][0]*float(0.01)
       node_y = NodeCoords[i][1]*float(0.01)
       GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, node, 1, node_x)
       GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 2, node_y)
       GeometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)

# Update the geometric field to ensure all nodal values are updated accross all cpus
GeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
GeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Export mesh results
#Fields = iron.Fields()
#Fields.CreateRegion(Region)
#Fields.NodesExport("Mesh", "FORTRAN")
#Fields.ElementsExport("Mesh", "FORTRAN")
#Fields.Finalise()

### metabolites ###

############################## ATP ###########################################

# Create ATP reaction diffusion with constant source equations_set

ATPEquationsSet = iron.EquationsSet()
ATPEquationsSetField = iron.Field()

ATPEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

ATPEquationsSet.CreateStart(ATPEquationsSetUserNumber,Region,GeometricField, ATPEquationsSetSpecification,
                            ATPEquationsSetFieldUserNumber,ATPEquationsSetField)

#  Set the equations set to be a standard Diffusion no source problem
#  Finish creating the equations set
ATPEquationsSet.CreateFinish()

# Create the equations set dependent field variables
ATPField = iron.Field()
ATPEquationsSet.DependentCreateStart(ATPFieldUserNumber,ATPField)
ATPField.VariableLabelSet(iron.FieldVariableTypes.U,'ATP Field')
ATPField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'ATP DU_DN Field')
ATPEquationsSet.DependentCreateFinish()

ATPField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ATPField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the equations set material field variables
ATPMaterialsField = iron.Field()
ATPEquationsSet.MaterialsCreateStart(ATPMaterialsFieldUserNumber,ATPMaterialsField)
ATPMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'ATP Materials Field')

# Set to element based interpolation
ATPMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
ATPMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)

# Finish the equations set materials field variables
ATPEquationsSet.MaterialsCreateFinish()

for i in range(number_of_nodes):
    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0], 1)
    if NodeDomain == ComputationalNodeNumber:
        # If the node is mitochondria with a border that ISN'T 1, make diffusion 1/10'th of normal
        if NodeCoords[i][2] == 1 and NodeNums[i][1] != 1:
            ATPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(0.1)*float(ATPDiffx))
            ATPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(0.1)*float(ATPDiffy))
        # Else, make the diffusion the normal value for that metabolite
        else:
            ATPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(ATPDiffx))
            ATPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(ATPDiffy))

ATPMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ATPMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set up source field for reaction diffusion equation set.
# Note that for the split problem subtype, the source field is not useall.  

ATPSourceField = iron.Field()
ATPEquationsSet.SourceCreateStart(ATPSourceFieldUserNumber,ATPSourceField)
ATPSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'ATP Source Field')

# Finish the equations set source field variables
ATPEquationsSet.SourceCreateFinish()

ATPSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

ATPSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ATPSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)



# ############################# ADP ##########################################

ADPEquationsSet = iron.EquationsSet()
ADPEquationsSetField = iron.Field()


ADPEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

ADPEquationsSet.CreateStart(ADPEquationsSetUserNumber,Region,GeometricField, ADPEquationsSetSpecification,
                            ADPEquationsSetFieldUserNumber,ADPEquationsSetField)

ADPEquationsSet.CreateFinish()


# Create the equations set dependent field variables

ADPField = iron.Field()
ADPEquationsSet.DependentCreateStart(ADPFieldUserNumber,ADPField)
ADPField.VariableLabelSet(iron.FieldVariableTypes.U,'ADP Field')
ADPField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'ADP DU_DN Field')
ADPEquationsSet.DependentCreateFinish()

ADPField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ADPField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


# Create the equations set material field variables

ADPMaterialsField = iron.Field()
ADPEquationsSet.MaterialsCreateStart(ADPMaterialsFieldUserNumber, ADPMaterialsField)
ADPMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'ADP Materials Field')

# Set to element based interpolation

ADPMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
ADPMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
ADPEquationsSet.MaterialsCreateFinish()


for i in range(number_of_nodes):
    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0], 1)
    if NodeDomain == ComputationalNodeNumber:        
        # If the node is mitochondria with a border that ISN'T 1, make diffusion 1/10'th of normal
        if NodeCoords[i][2] == 1 and NodeNums[i][1] != 1:
            ADPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(0.1)*float(ADPDiffx))
            ADPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(0.1)*float(ADPDiffy))
        # Else, make the diffusion the normal value for that metabolite
        else:
            ADPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(ADPDiffx))
            ADPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(ADPDiffy))



ADPMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ADPMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
# Set up source field for reaction diffusion equation set.

ADPSourceField = iron.Field()
ADPEquationsSet.SourceCreateStart(ADPSourceFieldUserNumber, ADPSourceField)
ADPSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'ADP Source Field')
# Finish the equations set source field variables
ADPEquationsSet.SourceCreateFinish()

ADPSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))
# Set the equations set to be a standard Diffusion no source problem

# Finish creating the equations set

ADPSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ADPSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)



############################## AMP ###########################################


# Create AMP reaction diffusion with constant source equations_set

AMPEquationsSet = iron.EquationsSet()
AMPEquationsSetField = iron.Field()


AMPEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

AMPEquationsSet.CreateStart(AMPEquationsSetUserNumber,Region,GeometricField, AMPEquationsSetSpecification,
                            AMPEquationsSetFieldUserNumber, AMPEquationsSetField)

AMPEquationsSet.CreateFinish()

# Create the equations set dependent field variables

AMPField = iron.Field()
AMPEquationsSet.DependentCreateStart(AMPFieldUserNumber, AMPField)
AMPField.VariableLabelSet(iron.FieldVariableTypes.U,'AMP Field')
AMPField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'AMP DU_DN Field')
AMPEquationsSet.DependentCreateFinish()

AMPField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
AMPField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the equations set material field variables

AMPMaterialsField = iron.Field()
AMPEquationsSet.MaterialsCreateStart(AMPMaterialsFieldUserNumber, AMPMaterialsField)
AMPMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'AMP Materials Field')

# Set to element based interpolation

AMPMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
AMPMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
AMPEquationsSet.MaterialsCreateFinish()


for i in range(number_of_nodes):
    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0], 1)
    if NodeDomain == ComputationalNodeNumber:
        # If the node is mitochondria with a border that ISN'T 1, make diffusion 1/10'th of normal
        if NodeCoords[i][2] == 1 and NodeNums[i][1] != 1:
            AMPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(0.1)*float(AMPDiffx))
            AMPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(0.1)*float(AMPDiffy))
        # Else, make the diffusion the normal value for that metabolite
        else:
            AMPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(AMPDiffx))
            AMPMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(AMPDiffy))

AMPMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
AMPMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set up source field for reaction diffusion equation set.

AMPSourceField = iron.Field()
AMPEquationsSet.SourceCreateStart(AMPSourceFieldUserNumber, AMPSourceField)
AMPSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'AMP Source Field')
# Finish the equations set source field variables
AMPEquationsSet.SourceCreateFinish()

AMPSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

AMPSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
AMPSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set the equations set to be a standard Diffusion no source problem

# Finish creating the equations set


############################## PCr ###########################################


# Create PCr reaction diffusion with constant source equations_set

PCrEquationsSet = iron.EquationsSet()
PCrEquationsSetField = iron.Field()


PCrEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

PCrEquationsSet.CreateStart(PCrEquationsSetUserNumber,Region,GeometricField, PCrEquationsSetSpecification,
                            PCrEquationsSetFieldUserNumber, PCrEquationsSetField)

PCrEquationsSet.CreateFinish()

# Create the equations set dependent field variables

PCrField = iron.Field()
PCrEquationsSet.DependentCreateStart(PCrFieldUserNumber, PCrField)
PCrField.VariableLabelSet(iron.FieldVariableTypes.U,'PCr Field')
PCrField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'PCr DU_DN Field')
PCrEquationsSet.DependentCreateFinish()

PCrField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PCrField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


# Create the equations set material field variables

PCrMaterialsField = iron.Field()
PCrEquationsSet.MaterialsCreateStart(PCrMaterialsFieldUserNumber, PCrMaterialsField)
PCrMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'PCr Materials Field')

# Set to element based interpolation

PCrMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
PCrMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
PCrEquationsSet.MaterialsCreateFinish()

for i in range(number_of_nodes):
    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0], 1)
    if NodeDomain == ComputationalNodeNumber:
        # If the node is mitochondria with a border that ISN'T 1, make diffusion 1/10'th of normal
        if NodeCoords[i][2] == 1 and NodeNums[i][1] != 1:
            PCrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(0.1)*float(PCrDiffx))
            PCrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(0.1)*float(PCrDiffy))
        # Else, make the diffusion the normal value for that metabolite
        else:
            PCrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(PCrDiffx))
            PCrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(PCrDiffy))

PCrMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PCrMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set up source field for reaction diffusion equation set.

PCrSourceField = iron.Field()
PCrEquationsSet.SourceCreateStart(PCrSourceFieldUserNumber, PCrSourceField)
PCrSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'PCr Source Field')
# Finish the equations set source field variables
PCrEquationsSet.SourceCreateFinish()

PCrSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

PCrSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PCrSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
# Set the equations set to be a standard Diffusion no source problem

# Finish creating the equations set


############################## Cr ############################################


# Create Cr reaction diffusion with constant source equations_set

CrEquationsSet = iron.EquationsSet()
CrEquationsSetField = iron.Field()


CrEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

CrEquationsSet.CreateStart(CrEquationsSetUserNumber,Region,GeometricField, CrEquationsSetSpecification,
                            CrEquationsSetFieldUserNumber, CrEquationsSetField)

CrEquationsSet.CreateFinish()

# Create the equations set dependent field variables

CrField = iron.Field()
CrEquationsSet.DependentCreateStart(CrFieldUserNumber, CrField)
CrField.VariableLabelSet(iron.FieldVariableTypes.U,'Cr Field')
CrField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Cr DU_DN Field')
CrEquationsSet.DependentCreateFinish()

CrField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
CrField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the equations set material field variables

CrMaterialsField = iron.Field()
CrEquationsSet.MaterialsCreateStart(CrMaterialsFieldUserNumber, CrMaterialsField)
CrMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'Cr Materials Field')

# Set to element based interpolation

CrMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
CrMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
CrEquationsSet.MaterialsCreateFinish()

for i in range(number_of_nodes):
    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0], 1)
    if NodeDomain == ComputationalNodeNumber:
        # If the node is mitochondria with a border that ISN'T 1, make diffusion 1/10'th of normal
        if NodeCoords[i][2] == 1 and NodeNums[i][1] != 1:
            CrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(0.1)*float(CrDiffx))
            CrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(0.1)*float(CrDiffy))
        # Else, make the diffusion the normal value for that metabolite
        else:
            CrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(CrDiffx))
            CrMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(CrDiffy))

CrMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
CrMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


# Set up source field for reaction diffusion equation set.

CrSourceField = iron.Field()
CrEquationsSet.SourceCreateStart(CrSourceFieldUserNumber, CrSourceField)
CrSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'Cr Source Field')
# Finish the equations set source field variables
CrEquationsSet.SourceCreateFinish()

CrSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

CrSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
CrSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)



############################## Pi ############################################


# Create Pi reaction diffusion with constant source equations_set

PiEquationsSet = iron.EquationsSet()
PiEquationsSetField = iron.Field()


PiEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

PiEquationsSet.CreateStart(PiEquationsSetUserNumber,Region,GeometricField, PiEquationsSetSpecification,
                            PiEquationsSetFieldUserNumber, PiEquationsSetField)

PiEquationsSet.CreateFinish()

# Create the equations set dependent field variables

PiField = iron.Field()
PiEquationsSet.DependentCreateStart(PiFieldUserNumber, PiField)
PiField.VariableLabelSet(iron.FieldVariableTypes.U,'Pi Field')
PiField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Pi DU_DN Field')
PiEquationsSet.DependentCreateFinish()

PiField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PiField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the equations set material field variables

PiMaterialsField = iron.Field()
PiEquationsSet.MaterialsCreateStart(PiMaterialsFieldUserNumber, PiMaterialsField)
PiMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'Pi Materials Field')

# Set to element based interpolation

PiMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
PiMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
PiEquationsSet.MaterialsCreateFinish()

for i in range(number_of_nodes):
    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0], 1)
    if NodeDomain == ComputationalNodeNumber:
        # If the node is mitochondria with a border that ISN'T 1, make diffusion 1/10'th of normal
        if NodeCoords[i][2] == 1 and NodeNums[i][1] != 1:
            PiMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(0.1)*float(PiDiffx))
            PiMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(0.1)*float(PiDiffy))
        # Else, make the diffusion the normal value for that metabolite
        else:
            PiMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, float(PiDiffx))
            PiMaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 2, float(PiDiffy))

PiMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PiMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


# Set up source field for reaction diffusion equation set.

PiSourceField = iron.Field()
PiEquationsSet.SourceCreateStart(PiSourceFieldUserNumber, PiSourceField)
PiSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'Pi Source Field')
# Finish the equations set source field variables
PiEquationsSet.SourceCreateFinish()

PiSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

PiSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PiSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

############################## Oxy ###########################################


# Create Oxy reaction diffusion with constant source equations_set

OxyEquationsSet = iron.EquationsSet()
OxyEquationsSetField = iron.Field()


OxyEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

OxyEquationsSet.CreateStart(OxyEquationsSetUserNumber,Region,GeometricField, OxyEquationsSetSpecification,
                            OxyEquationsSetFieldUserNumber, OxyEquationsSetField)

OxyEquationsSet.CreateFinish()

OxyField = iron.Field()
OxyEquationsSet.DependentCreateStart(OxyFieldUserNumber, OxyField)
OxyField.VariableLabelSet(iron.FieldVariableTypes.U,'Oxy Field')
OxyField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Oxy DU_DN Field')
OxyEquationsSet.DependentCreateFinish()

OxyField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
OxyField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the equations set material field variables

OxyMaterialsField = iron.Field()
OxyEquationsSet.MaterialsCreateStart(OxyMaterialsFieldUserNumber, OxyMaterialsField)
OxyMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'Oxy Materials Field')

OxyMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
OxyMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
OxyEquationsSet.MaterialsCreateFinish()
#check > diffx, diffy??? in ALL
OxyMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, float(OxyDiffx)) 
OxyMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(OxyDiffx)) 

OxyMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
OxyMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set up source field for reaction diffusion equation set.

OxySourceField = iron.Field()
OxyEquationsSet.SourceCreateStart(OxySourceFieldUserNumber, OxySourceField)
OxySourceField.VariableLabelSet(iron.FieldVariableTypes.U,'Oxy Source Field')
# Finish the equations set source field variables
OxyEquationsSet.SourceCreateFinish()

OxySourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

OxySourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
OxySourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set the equations set to be a standard Diffusion no source problem   # Finish creating the equations set


############################## DPsi ##########################################

# Create DPsi reaction diffusion with constant source equations_set

DPsiEquationsSet = iron.EquationsSet()
DPsiEquationsSetField = iron.Field()


DPsiEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

DPsiEquationsSet.CreateStart(DPsiEquationsSetUserNumber,Region,GeometricField, DPsiEquationsSetSpecification,
                            DPsiEquationsSetFieldUserNumber, DPsiEquationsSetField)

DPsiEquationsSet.CreateFinish()

DPsiField = iron.Field()
DPsiEquationsSet.DependentCreateStart(DPsiFieldUserNumber, DPsiField)
DPsiField.VariableLabelSet(iron.FieldVariableTypes.U,'DPsi Field')
DPsiField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'DPsi DU_DN Field')
DPsiEquationsSet.DependentCreateFinish()

DPsiField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
DPsiField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the equations set material field variables

DPsiMaterialsField = iron.Field()
DPsiEquationsSet.MaterialsCreateStart(DPsiMaterialsFieldUserNumber, DPsiMaterialsField)
DPsiMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'DPsi Materials Field')

DPsiMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
DPsiMaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
DPsiEquationsSet.MaterialsCreateFinish()
#check > diifx, diffy??? in ALL
DPsiMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0)) 
DPsiMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, float(0.0))
#specify parameter, zeroval = (0.0_dp) above

DPsiMaterialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
DPsiMaterialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set up source field for reaction diffusion equation set.

DPsiSourceField = iron.Field()
DPsiEquationsSet.SourceCreateStart(DPsiSourceFieldUserNumber, DPsiSourceField)
DPsiSourceField.VariableLabelSet(iron.FieldVariableTypes.U,'DPsi Source Field')
# Finish the equations set source field variables
DPsiEquationsSet.SourceCreateFinish()

DPsiSourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(0.0))

DPsiSourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
DPsiSourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Set the equations set to be a standard Diffusion no source problem   # Finish creating the equations set

print('Now start to set up CellML fields')
##############################################################################

# Start to set up CellML Fields
# Create the CellML environment
CellML = iron.CellML()
CellML.CreateStart(CellMLUserNumber,Region)
#import cellml model
MitochondriaIndex = CellML.ModelImport('mitochondria_new.cellml')
MyofibrilIndex = CellML.ModelImport('myofibril_new.cellml')

CellML.VariableSetAsKnown(MitochondriaIndex, 'general_constants/param')
CellML.VariableSetAsKnown(MyofibrilIndex, 'ATP/param')
CellML.VariableSetAsWanted(MitochondriaIndex, 'ANT_flux/V_ANT')
CellML.VariableSetAsWanted(MyofibrilIndex, 'H_ATP/H_ATP')
CellML.VariableSetAsWanted(MyofibrilIndex, 'v_CK/v_CK')
CellML.VariableSetAsWanted(MitochondriaIndex, 'v_MiCK/v_MiCK')
CellML.VariableSetAsWanted(MitochondriaIndex, 'Electron_flux_complex_I/V_C1')
CellML.VariableSetAsWanted(MitochondriaIndex, 'Electron_flux_complex_III/V_C3')
CellML.VariableSetAsWanted(MitochondriaIndex, 'Electron_flux_complex_IV/V_C4')
CellML.VariableSetAsWanted(MitochondriaIndex, 'ATP_synthesis_flux/V_F1')
CellML.VariableSetAsWanted(MitochondriaIndex, 'dATP_x_dt/ATP_xx')
CellML.VariableSetAsWanted(MitochondriaIndex, 'dADP_x_dt/ADP_xx')
CellML.VariableSetAsWanted(MitochondriaIndex, 'dPi_x_dt/Pi_xx')
CellML.VariableSetAsWanted(MitochondriaIndex, 'Proton_motive_force/dG_H')

# Finish the CellML environment
CellML.CreateFinish()

### Start the creation of CellML OpenCMISS field maps ###

CellML.FieldMapsCreateStart()
# dependent field, solve the dae, and then put the result of the dae into the source field.
CellML.CreateFieldToCellMLMap(ATPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'ATPi/ATPi', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'ATPi/ATPi', iron.FieldParameterSetTypes.VALUES,
                              ATPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(ADPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'ADPi/ADPi', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'ADPi/ADPi', iron.FieldParameterSetTypes.VALUES,
                              ADPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(AMPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'AMPi/AMPi', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'AMPi/AMPi', iron.FieldParameterSetTypes.VALUES,
                              AMPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(PCrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'PCri/PCri', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'PCri/PCri', iron.FieldParameterSetTypes.VALUES,
                              PCrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(CrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'Cri/Cri', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'Cri/Cri', iron.FieldParameterSetTypes.VALUES,
                              CrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(PiField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'Pii/Pii', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'Pii/Pii', iron.FieldParameterSetTypes.VALUES,
                              PiField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(OxyField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'dO2_dt/O2', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'dO2_dt/O2', iron.FieldParameterSetTypes.VALUES,
                              OxyField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(DPsiField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'dPsi_dt/dPsi', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MitochondriaIndex, 'dPsi_dt/dPsi', iron.FieldParameterSetTypes.VALUES,
                              DPsiField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(ATPSourceField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MitochondriaIndex, 'general_constants/param', iron.FieldParameterSetTypes.VALUES)
CellML.CreateFieldToCellMLMap(ATPSourceField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'ATP/param', iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(ATPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'ATP/ATP', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MyofibrilIndex, 'ATP/ATP', iron.FieldParameterSetTypes.VALUES,
                              ATPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(ADPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'ADP/ADP', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MyofibrilIndex, 'ADP/ADP', iron.FieldParameterSetTypes.VALUES,
                              ADPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(AMPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'AMP/AMP', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MyofibrilIndex, 'AMP/AMP', iron.FieldParameterSetTypes.VALUES,
                              AMPField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(PCrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'PCr/PCr', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MyofibrilIndex, 'PCr/PCr', iron.FieldParameterSetTypes.VALUES,
                              PCrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(CrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'Cr/Cr', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MyofibrilIndex, 'Cr/Cr', iron.FieldParameterSetTypes.VALUES,
                              CrField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)

CellML.CreateFieldToCellMLMap(PiField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                              MyofibrilIndex, 'Pi/Pi', iron.FieldParameterSetTypes.VALUES)
CellML.CreateCellMLToFieldMap(MyofibrilIndex, 'Pi/Pi', iron.FieldParameterSetTypes.VALUES,
                              PiField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)


# Finish the creation of CellML OpenCMISS field maps
CellML.FieldMapsCreateFinish()


### Start the creation of the CellML models field ###
#This field is an integer field that stores which nodes have which cellml model

CellMLModelsField = iron.Field()
CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber, CellMLModelsField)
CellML.ModelsFieldCreateFinish()

# The CellMLModelsField is an integer field that stores which model is being used by which node.
# By default all field parameters have default model value of 1, i.e. the first model. But, this command below is for example purposes

# Assigning the bufferNryr cellml model (model 1) for all nodes.

CellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, int(2))

for i in range(number_of_nodes):
    if NodeCoords[i][2] == 1:
       NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
       if NodeDomain == ComputationalNodeNumber :
           print(" Setting first model for node ",NodeNums[i][0],"\n")
           CellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, int(1))
    elif NodeCoords[i][2] < 10:
        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
        if NodeDomain == ComputationalNodeNumber :
            print(" Setting no model for node ",NodeNums[i][0],"\n")
            CellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, int(0))

CellMLModelsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
CellMLModelsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

### Start the creation of the CellML state field ###

CellMLStateField = iron.Field()
CellML.StateFieldCreateStart(CellMLStateFieldUserNumber,CellMLStateField)
CellML.StateFieldCreateFinish()

### Start the creation of CellML parameters field ###

CellMLParametersField = iron.Field()
CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber, CellMLParametersField)
CellML.ParametersFieldCreateFinish()

### Start the creation of CellML intermediate field ###

CellMLIntermediateField = iron.Field()
CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber, CellMLIntermediateField)
CellML.IntermediateFieldCreateFinish()

# Set initial value of the dependent field/state variable,
# Set initial value of the dependent field/state variable,

# Cole D changes start:
# This code is purely to test and see if mito produce, myo consume, and diffusion is happening by changing initial metabolite values
# 1 = Mitochondria, 10 = Myofibrils, 2 = Glycogen, 3 = 'Other'


#for i in range(number_of_nodes):
#    if NodeCoords[i][2] == 1:
#       NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#       if NodeDomain == ComputationalNodeNumber :
#           ATPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ATP))
#    elif NodeCoords[i][2] >= 10:
#        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#        if NodeDomain == ComputationalNodeNumber :
#           ATPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, 1*float(init_ATP))
#    elif NodeCoords[i][2] == 2:
#        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#        if NodeDomain == ComputationalNodeNumber :
#           ATPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ATP))
#    elif NodeCoords[i][2] == 3:
#        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#        if NodeDomain == ComputationalNodeNumber :
#           ATPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ATP))

#for i in range(number_of_nodes):
#    if NodeCoords[i][2] == 1:
#       NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#       if NodeDomain == ComputationalNodeNumber :
#           ADPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ADP))
#    elif NodeCoords[i][2] >= 10:
#        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#        if NodeDomain == ComputationalNodeNumber :
#           ADPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ADP))
#    elif NodeCoords[i][2] == 2:
#        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#        if NodeDomain == ComputationalNodeNumber :
#           ATPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ADP))
#    elif NodeCoords[i][2] == 3:
#        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#        if NodeDomain == ComputationalNodeNumber :
#           ATPField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, NodeNums[i][0], 1, .000001*float(init_ADP))

# Cole D changes:
# don't forget to comment these next two lines out.

#ATPField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, .000001*float(init_ATP))
#ADPField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, .000001*float(init_ADP))

# Cole D changes:
# don't forget to comment these next two lines back in after testing of initial values is done.
ATPField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_ATP))
ADPField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_ADP))
AMPField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_AMP))
PCrField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_PCr))
CrField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_Cr))
PiField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_Pi))
OxyField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_Oxy))
DPsiField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, float(init_DPsi))


DPsiField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
DPsiField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

ATPField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ATPField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

ADPField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
ADPField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

AMPField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
AMPField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

PCrField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PCrField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

CrField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
CrField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

PiField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
PiField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

### Create the equations set equations for ATP ###
ATPEquations = iron.Equations()
ATPEquationsSet.EquationsCreateStart(ATPEquations)

# Set the equations matrices sparsity type
ATPEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
ATPEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)

# Set the equations set output
#ATPEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
#ATPEquations.OutputTypeSet(iron.EquationsOutputTypes.TIMING)
#ATPEquations.OutputTypeSet(iron.EquationsOutputTypes.MATRIX)
#ATPEquations.OutputTypeSet(iron.EquationsOutputTypes.ELEMENT_MATRIX)
# Finish the equations set equations
ATPEquationsSet.EquationsCreateFinish()

# Create the equations set equations for ADP
ADPEquations = iron.Equations()
ADPEquationsSet.EquationsCreateStart(ADPEquations)
ADPEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
ADPEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
ADPEquationsSet.EquationsCreateFinish()

# Create the equations set equations for AMP
AMPEquations = iron.Equations()
AMPEquationsSet.EquationsCreateStart(AMPEquations)
AMPEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
AMPEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
AMPEquationsSet.EquationsCreateFinish()

# Create the equations set equations for PCr
PCrEquations = iron.Equations()
PCrEquationsSet.EquationsCreateStart(PCrEquations)
PCrEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
PCrEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
PCrEquationsSet.EquationsCreateFinish()

# Create the equations set equations for Cr
CrEquations = iron.Equations()
CrEquationsSet.EquationsCreateStart(CrEquations)
CrEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
CrEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
CrEquationsSet.EquationsCreateFinish()

# Create the equations set equations for Pi
PiEquations = iron.Equations()
PiEquationsSet.EquationsCreateStart(PiEquations)
PiEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
PiEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
PiEquationsSet.EquationsCreateFinish()

# Create the equations set equations for Oxy
OxyEquations = iron.Equations()
OxyEquationsSet.EquationsCreateStart(OxyEquations)
OxyEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
OxyEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
OxyEquationsSet.EquationsCreateFinish()

# Create the equations set equations for DPsi
DPsiEquations = iron.Equations()
DPsiEquationsSet.EquationsCreateStart(DPsiEquations)
DPsiEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
DPsiEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
DPsiEquationsSet.EquationsCreateFinish()

#________________________________________________________________________________________________
# REACTION DIFFUSION PROBLEM AND SOLVER
#________________________________________________________________________________________________

#------------------------------------------------------------------------------------------------
# PROBLEM - REACTION DIFFUSION
#------------------------------------------------------------------------------------------------

### Create the problem ###
Problem = iron.Problem()
ProblemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                             iron.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                             iron.ProblemSubtypes.CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT]
# Set the problem to be a strang split reaction diffusion problem
Problem.CreateStart(ProblemUserNumber, ProblemSpecification)
Problem.CreateFinish()

# Create the problem control loop
Problem.ControlLoopCreateStart()
# Get the control loop
ControlLoop = iron.ControlLoop()
Problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], ControlLoop)
# start time, stop time, time increment
ControlLoop.TimesSet(startT,endT,Tstep)

# Control Loop Outputs
ControlLoop.TimeOutputSet(outputFreq)
#controlLoop.LoadOutputSet(1)
ControlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
ControlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)
# Set the times
#ControlLoop.TimesSet(float(0.0), float(2.0), float(1.0))    #UPDATE TIME STEP SOLVEEEE
ControlLoop.TimeOutputSet(outputFreq)
Problem.ControlLoopCreateFinish()


Problem.SolversCreateStart()
#------------------------------------------------------------------------------------------------
# SOLVER - REACTION DIFFUSION
#------------------------------------------------------------------------------------------------
#    1st Solver --> DAE for cellml ODE
#         |
#         v
#    2nd Solver --> Dynamic for PDE
#         |
#         v
#    3rd Solver --> DAE for cellml ODE

# First solver is a DAE solver
Solver = iron.Solver()
Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, Solver)
Solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
Solver.DAETimeStepSet(float(0.01))
Solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

# Second solver is the dynamic solver for solving the parabolic equation
Solver = iron.Solver()
LinearSolver = iron.Solver()
Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 2, Solver)
# Set theta - backward vs forward time step parameter
Solver.DynamicThetaSet([float(1.0)])
# Solver.OutputTypeSet(iron.SolverOutputTypes.NONE)
Solver.OutputTypeSet(iron.SolverOutputTypes.TIMING)
# Get the dynamic linear solver from the solver
Solver.DynamicLinearSolverGet(LinearSolver)
LinearSolver.LibraryTypeSet(iron.SolverLibraries.LAPACK)
LinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
LinearSolver.LinearDirectTypeSet(iron.DirectLinearSolverTypes.LU)

#LinearSolver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
#LinearSOlver.LinearIterativeMaximumIterationsSet(10000)
#LinearSolver.OutputTypeSet(iron.SolverOutputTypes.NONE)
#LinearSolver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
#LinearSolver.OutputTypeSet(iron.SolverOutputTypes.TIMING)
#LinearSolver.OutputTypeSet(iron.SolverOutputTypes.SOLVER)

# Third solver is another DAE solver
Solver = iron.Solver()
Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 3, Solver)
Solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
Solver.DAETimeStepSet(float(0.01))
#Solver.DAETimeStepSet(ODE_TIME_STEP)
Solver.OutputTypeSet(iron.SolverOutputTypes.NONE)
#Solver.OutputTypeSet(iron.SolverOutputTypes.TIMING)
#Solver.OutputTypeSet(iron.SolverOutputTypes.SOLVER)
#Solver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)

# Finish the creation of the problem solver
Problem.SolversCreateFinish()

### Start the creation of the problem solver CellML equations ###
Problem.CellMLEquationsCreateStart()

# Get the first solver.      Get the CellML Equations
Solver = iron.Solver()
Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, Solver)
CellMLEquations = iron.CellMLEquations()
Solver.CellMLEquationsGet(CellMLEquations)
# Add in the CellML environment
CellMLIndex = CellMLEquations.CellMLAdd(CellML)

# Get the third solver.      Get the CellML equations
Solver = iron.Solver()
Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 3, Solver)
CellMLEquations = iron.CellMLEquations()
Solver.CellMLEquationsGet(CellMLEquations)
# Add in the CellML environment
CellMLIndex = CellMLEquations.CellMLAdd(CellML)

# Finish the creation of the problem solver CellML equations
Problem.CellMLEquationsCreateFinish()

# Start the creation of the problem solver equations
Problem.SolverEquationsCreateStart()

# Get the second solver.      Get the solver equations
Solver = iron.Solver()
Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 2, Solver)
SolverEquations = iron.SolverEquations()
Solver.SolverEquationsGet(SolverEquations)
# Set the solver equations sparsity
SolverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)

# Add in the equation set
ATPEquationSetIndex = SolverEquations.EquationsSetAdd(ATPEquationsSet)
ADPEquationSetIndex = SolverEquations.EquationsSetAdd(ADPEquationsSet)
AMPEquationSetIndex = SolverEquations.EquationsSetAdd(AMPEquationsSet)
PCREquationSetIndex = SolverEquations.EquationsSetAdd(PCrEquationsSet)
CrEquationSetIndex = SolverEquations.EquationsSetAdd(CrEquationsSet)
PiEquationSetIndex = SolverEquations.EquationsSetAdd(PiEquationsSet)
OxyEquationSetIndex = SolverEquations.EquationsSetAdd(OxyEquationsSet)
#DPsiEquationSetIndex = SolverEquations.EquationsSetAdd(DPsiEquationsSet)

# Finish the creation of the problem solver equations
Problem.SolverEquationsCreateFinish()


##############################################################################

print('Set up boundary conditions')
BoundaryConditions = iron.BoundaryConditions()
SolverEquations.BoundaryConditionsCreateStart(BoundaryConditions)

CONDITION = iron.BoundaryConditionsTypes.FIXED
BCVALUE = float(0.0)

for i in range(number_of_nodes):
    if NodeNums[i][1] == 1:
        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
        if NodeDomain == ComputationalNodeNumber :
           BoundaryConditions.SetNode(ATPField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, BCVALUE)
           BoundaryConditions.SetNode(ADPField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, BCVALUE)
           BoundaryConditions.SetNode(AMPField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, BCVALUE)
           BoundaryConditions.SetNode(PCrField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, BCVALUE)
           BoundaryConditions.SetNode(CrField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, BCVALUE)
           BoundaryConditions.SetNode(PiField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, BCVALUE)

    if NodeNums[i][1] == 1:
        NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
        if NodeDomain == ComputationalNodeNumber :
           BoundaryConditions.SetNode(OxyField, iron.FieldVariableTypes.U, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, NodeNums[i][0], 1,
                                   CONDITION, float(init_Oxy))

SolverEquations.BoundaryConditionsCreateFinish()
##############################################################################
Fields = iron.Fields()
Fields.CreateRegion(Region)
Fields.NodesExport("Diffusion", "FORTRAN")
Fields.ElementsExport("Diffusion", "FORTRAN")
Fields.Finalise()


print('sooolllllvvveeeeeeeee')

# Solve the problem

# FOR loop
#UPDATE TIME loop
# SOLUTION FIELDS



Problem.Solve()

Fields = iron.Fields()
Fields.CreateRegion(Region)
Fields.NodesExport("Diffusion1", "FORTRAN")
Fields.ElementsExport("Diffusion1", "FORTRAN")
Fields.Finalise()
#
# Obtain the nodeal values and their coordinates
#points = np.array(NodeCoords)
# Obtain the connectivity matrix between the nodes at each element
# Note that meshio uses a 0-index node labeling system, while Opencmiss uses a 1-index node labeling system. So all node numbers need to be reduced by 1.
# elements_list[0] represents the element number
# elements_list[1:5] represents the 4 nodes associated with that element
#cells_array = np.array(elements_list)[:][1:4]
#cells = [("tetra", cells_array)]

#ATPFieldNodes = np.zeros(number_of_nodes)
#for i in range(number_of_nodes):
#    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#    if NodeDomain == ComputationalNodeNumber :
#       node = i + 1
#       ATPFieldNodes[i] = ATPField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)
#meshio.write_points_cells("output_ATP.vtk",points,cells,point_data={"solution":ATPFieldNodes})

#ADPFieldNodes = np.zeros(number_of_nodes)
#for i in range(number_of_nodes):
#    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#    if NodeDomain == ComputationalNodeNumber :
#       node = i + 1
#       ADPFieldNodes[i] = ADPField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)
#meshio.write_points_cells("output_ADP.vtk",points,cells,point_data={"solution":ADPFieldNodes})

#AMPFieldNodes = np.zeros(number_of_nodes)
#for i in range(number_of_nodes):
#    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#    if NodeDomain == ComputationalNodeNumber :
#       node = i + 1
#       AMPFieldNodes[i] = AMPField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)
#meshio.write_points_cells("output_AMP.vtk",points,cells,point_data={"solution":AMPFieldNodes})

#PCrFieldNodes = np.zeros(number_of_nodes)
#for i in range(number_of_nodes):
#    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#    if NodeDomain == ComputationalNodeNumber :
#       node = i + 1
#       PCrFieldNodes[i] = PCrField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)
#meshio.write_points_cells("output_PCr.vtk",points,cells,point_data={"solution":PCrFieldNodes})

#PiFieldNodes = np.zeros(number_of_nodes)
#for i in range(number_of_nodes):
#    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#    if NodeDomain == ComputationalNodeNumber :
#       node = i + 1
#       PiFieldNodes[i] = PiField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)
#meshio.write_points_cells("output_Pi.vtk",points,cells,point_data={"solution":PiFieldNodes})

#CrFieldNodes = np.zeros(number_of_nodes)
#for i in range(number_of_nodes):
#    NodeDomain = Decomposition.NodeDomainGet(NodeNums[i][0],1)
#    if NodeDomain == ComputationalNodeNumber :
#       node = i + 1
#       CrFieldNodes[i] = CrField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, node, 1)
#meshio.write_points_cells("output_Cr.vtk",points,cells,point_data={"solution":CrFieldNodes})

# Export results

#Fields = iron.Fields()
#Fields.CreateRegion(Region)
#Fields.NodesExport("Diffusion", "FORTRAN")
#Fields.ElementsExport("Diffusion", "FORTRAN")
#Fields.Finalise()

print('Successfully completed')

iron.Finalise()

endtime=time.time()

# Calculate the total execution time
execution_time = endtime - starttime

# Output execution time to a text file
with open("execution_time.txt", "w") as file:
        file.write(f"Total execution time: {execution_time} seconds") 
#  END

