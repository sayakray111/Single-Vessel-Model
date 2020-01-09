#!/usr/bin/env python
import os
import numpy as np
import time
import placentagen as pg
from opencmiss.iron import iron
el_type = 2

nel_x = 10
nel_y = 10
nel_z = 10

export_mesh = False
filename_mesh = 'expected-results/placenta_mesh'

export_results = True
export_directory = 'output'

porosity = 1.0
perm_over_vis = 1.0 # permiability over vicosity

initial_velocity = 1
initial_pressure = 1

# Set problem parameters
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    materialFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,13)

#Decide number of computational nodes and mesh properties
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "DarcyRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-quadratic Lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.numberOfXi = 3
if(el_type == 2):
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
    numberOfNodesXi = 3
    numberOfGaussXi = 3
elif(el_type ==1):
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
    numberOfNodesXi = 2
    numberOfGaussXi = 2
basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*3)
basis.quadratureLocalFaceGaussEvaluate = True
basis.CreateFinish()

# Create a generated mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [1,1,1]
generatedMesh.numberOfElements = [nel_x,nel_y,nel_z]
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()
#geometricField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
#geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
#Update the geometric field parameters
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
generatedMesh.GeometricParametersCalculate(geometricField)

# Create standard Darcy equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
        iron.EquationsSetTypes.DARCY_EQUATION,
        iron.EquationsSetSubtypes.STANDARD_DARCY]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 1, 1)

dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()


materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,porosity)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,perm_over_vis)

dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,1,initial_velocity)
#dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,2,initial_velocity)
#dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,3,initial_velocity)
#dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,initial_pressure)

#Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Darcy equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
        iron.ProblemTypes.DARCY_EQUATION,
        iron.ProblemSubtypes.STANDARD_DARCY]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
#solver.outputType = iron.SolverOutputTypes.SOLVER
solver.outputType = iron.SolverOutputTypes.NONE
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-13
solver.linearIterativeRelativeTolerance = 1.0E-13
problem.SolversCreateFinish()

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create Boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
nodes = iron.Nodes()
region.NodesGet(nodes)
Total_Nodes = nodes.NumberOfNodesGet()

##Putting Boundary conditions
f_val = []
Node_Number = []
f_val1 = []
f_val2 = []
f_val3 = []
z_val = 0
x_val = 0
y_val = 0
k = 0
for i in range(1,Total_Nodes):
    z_val = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,i,3)
    x_val = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,i,1)
    y_val = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,i,2)
    if(z_val == 0 or z_val == 1 or y_val == 0 or y_val == 1 or x_val == 0 or x_val == 1):
        Node_Number.append(i)
        f_val.append(x_val*y_val*z_val)
        f_val1.append(y_val*z_val)
        f_val2.append(x_val*z_val)
        f_val3.append(y_val*x_val)
        k = k+1
for j in range(0,k):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Node_Number[j],1,iron.BoundaryConditionsTypes.FIXED,-f_val1[j])
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Node_Number[j],2,iron.BoundaryConditionsTypes.FIXED,-f_val2[j])
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Node_Number[j],3,iron.BoundaryConditionsTypes.FIXED,-f_val3[j])
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Node_Number[j],4,iron.BoundaryConditionsTypes.FIXED,f_val[j])


boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Total_Nodes,1,iron.BoundaryConditionsTypes.FIXED,1)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Total_Nodes,2,iron.BoundaryConditionsTypes.FIXED,1)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Total_Nodes,3,iron.BoundaryConditionsTypes.FIXED,1)
boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,Total_Nodes,4,iron.BoundaryConditionsTypes.FIXED,1)
    
solverEquations.BoundaryConditionsCreateFinish()
start_time = time.time()
#print(dependentField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,27,1))
# Solve the problem
problem.Solve()
end_time = time.time()

print ('Total time for solve = '+ str((end_time-start_time)/60.0) + ' mins')


#if(export_results):
    #if not os.path.exists(export_directory):
        #os.makedirs(export_directory)
        ## Export results
    #export_file = export_directory + '/StaticDarcy'
    #fields = iron.Fields()
    #fields.CreateRegion(region)
    #fields.NodesExport(export_file,"FORTRAN")
    #fields.ElementsExport(export_file,"FORTRAN")
    #fields.Finalise()


iron.Finalise()
raise SystemExit

