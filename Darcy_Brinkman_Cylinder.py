#> This example program solves a FSI problem in a multi-block tube using OpenCMISS.
#>
#> By Chris Bradley
#>
#>

#================================================================================================================================
#  Start Program
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import math, numpy, csv, time, sys, os, pdb
from opencmiss.iron import iron
from mpl_toolkits import mplot3d
#import matplotlib.pyplot as plt

LINEAR = 1
QUADRATIC = 2

numberOfSquareElements =2
numberOfArmElements = 2
numberOfLengthElements = 5

pipeRadius = 5.#12
lengthSize = 40.#24
squareSizeRatio = 0.500

pressure_bc = False
Pdrop = 100.
Pbase =0
v_z = 21.
beta = 0.4
porosity = .4
viscosity = 3e-3
vol_frac = 1-porosity
d = .5
perm = d**2. * (1.-vol_frac)**3./(180.*vol_frac**2.)
perm_visc = perm/viscosity

print(porosity,perm_visc)
uInterpolation = QUADRATIC
pInterpolation = LINEAR

progressDiagnostics = True
debug = False
#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()############################
computationalNodeNumber = iron.ComputationalNodeNumberGet()#############################
#================================================================================================================================
#  Should not need to change anything below here.
#================================================================================================================================

if numberOfLengthElements == 0:
    numberOfDimensions = 2
else:
    numberOfDimensions = 3

fluidCoordinateSystemUserNumber     = 1
  
fluidRegionUserNumber = 2

uBasisUserNumber = 1
pBasisUserNumber = 2

fluidMeshUserNumber     = 1
  
fluidDecompositionUserNumber     = 1
  
fluidGeometricFieldUserNumber = 11
 
#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

#iron.OutputSetOn("Testing")

#================================================================================================================================
#  Coordinate Systems
#================================================================================================================================

if (progressDiagnostics):
    print(' ')
    print('Coordinate systems ...')

# Create a RC coordinate system for the fluid region
fluidCoordinateSystem = iron.CoordinateSystem()
fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber)
fluidCoordinateSystem.DimensionSet(3)
fluidCoordinateSystem.CreateFinish()

if (progressDiagnostics):
    print('Coordinate systems ... Done')
  
#================================================================================================================================
#  Regions
#================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

# Create a fluid region
fluidRegion = iron.Region()
fluidRegion.CreateStart(fluidRegionUserNumber,iron.WorldRegion)
fluidRegion.label = 'FluidRegion'
fluidRegion.coordinateSystem = fluidCoordinateSystem
fluidRegion.CreateFinish()

if (progressDiagnostics):
    print('Regions ... Done')

#================================================================================================================================
#  Bases
#================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')
    
numberOfNodesXi = uInterpolation+1
numberOfGaussXi = uInterpolation+1

uBasis = iron.Basis()
uBasis.CreateStart(uBasisUserNumber)
uBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
uBasis.numberOfXi = 3
if (uInterpolation == LINEAR):
    uBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
elif (uInterpolation == QUADRATIC):
    uBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
else:
    print('Invalid u interpolation')
    exit()
uBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
uBasis.quadratureLocalFaceGaussEvaluate = True
uBasis.CreateFinish()

pBasis = iron.Basis()
pBasis.CreateStart(pBasisUserNumber)
pBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
pBasis.numberOfXi = 3
pBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
pBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
pBasis.quadratureLocalFaceGaussEvaluate = True
pBasis.CreateFinish()

numberOfLocalNodes = numberOfNodesXi*numberOfNodesXi*numberOfNodesXi
numberOfLocalInterfaceNodes = numberOfNodesXi*numberOfNodesXi
localNodeIdx000 = 0
localNodeIdx100 = numberOfNodesXi-1
localNodeIdx010 = numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx110 = numberOfNodesXi*numberOfNodesXi-1
localNodeIdx001 = numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx101 = numberOfNodesXi-1+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx011 = numberOfNodesXi*(numberOfNodesXi-1)+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx111 = numberOfLocalNodes-1

if (progressDiagnostics):
    print('Basis functions ... Done')
  
#================================================================================================================================
#  Mesh
#================================================================================================================================

numberOfFluidNodesPerBlock = numberOfSquareElements*(numberOfNodesXi-1)*(numberOfArmElements*(numberOfNodesXi-1)+1)
numberOfFluidElementsPerBlock = numberOfSquareElements*numberOfArmElements
numberOfFluidNodesPerLength = 4*numberOfFluidNodesPerBlock+ \
                     (numberOfSquareElements*(numberOfNodesXi-1)-1)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
numberOfFluidElementsPerLength = 4*numberOfFluidElementsPerBlock+numberOfSquareElements*numberOfSquareElements
numberOfFluidNodes = numberOfFluidNodesPerLength*(numberOfLengthElements*(numberOfNodesXi-1)+1)
numberOfFluidElements = numberOfFluidElementsPerLength*numberOfLengthElements

if (debug):
    print('  Mesh Parameters:')
    print('    numberOfSquareElements: %d' % (numberOfSquareElements))
    print('    numberOfArmElements: %d' % (numberOfArmElements))
    print('    numberOfLengthElements: %d' % (numberOfLengthElements))
    print('    numberOfWallElements: %d' % (numberOfWallElements))
    print('    numberOfNodesXi: %d' % (numberOfNodesXi))
    print('    numberOfFluidNodesPerBlock: %d' % (numberOfFluidNodesPerBlock))
    print('    numberOfElementPerBlock: %d' % (numberOfFluidElementsPerBlock))
    print('    numberOfFluidNodesPerLength: %d' % (numberOfFluidNodesPerLength))
    print('    numberOfFluidElementsPerLength: %d' % (numberOfFluidElementsPerLength))
    print('    numberOfFluidNodes: %d' % (numberOfFluidNodes))
    print('    numberOfFluidElements: %d' % (numberOfFluidElements))
        
fluidNodes = iron.Nodes()
fluidNodes.CreateStart(fluidRegion,numberOfFluidNodes)
fluidNodes.CreateFinish()

fluidMesh = iron.Mesh()
fluidMesh.CreateStart(fluidMeshUserNumber,fluidRegion,3)
fluidMesh.NumberOfElementsSet(numberOfFluidElements)
fluidMesh.NumberOfComponentsSet(2)

if (debug):
    print('  Fluid Elements:')

fluidUElements = iron.MeshElements()
fluidUElements.CreateStart(fluidMesh,1,uBasis)
fluidPElements = iron.MeshElements()
fluidPElements.CreateStart(fluidMesh,2,pBasis)

for zElementIdx in range(1,max(numberOfLengthElements+1,2)):
    #Handle the arm blocks first
    previousBlock = 4
    for blockIdx in range(1,5):
        for yElementIdx in range(1,numberOfArmElements+1):
            for xElementIdx in range(1,numberOfSquareElements+1):
                localNodes = [0]*numberOfLocalNodes
                elementNumber = xElementIdx+(yElementIdx-1)*numberOfSquareElements+(blockIdx-1)*numberOfSquareElements*numberOfArmElements+\
                                (zElementIdx-1)*numberOfFluidElementsPerLength
                if (xElementIdx == 1):
                    localNodes[localNodeIdx000] = (previousBlock-1)*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)+ \
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                 (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx100] = (blockIdx-1)*numberOfFluidNodesPerBlock+numberOfNodesXi-1+\
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                 (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                else:
                    localNodes[localNodeIdx000] = (blockIdx-1)*numberOfFluidNodesPerBlock+(xElementIdx-2)*(numberOfNodesXi-1)+(numberOfNodesXi-2)+1+\
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSquareElements*(numberOfNodesXi-1))+\
                                                 (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+numberOfNodesXi-1
                localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1)*(numberOfNodesXi-1)
                localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + numberOfSquareElements*(numberOfNodesXi-1)*(numberOfNodesXi-1)
                localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                if(uInterpolation == QUADRATIC):
                    localNodes[1] = localNodes[localNodeIdx100] - 1
                    localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1)
                    localNodes[4] = localNodes[1] + numberOfSquareElements*(numberOfNodesXi-1)
                    localNodes[5] = localNodes[4] + 1
                    localNodes[7] = localNodes[localNodeIdx110] - 1
                    localNodes[9] = localNodes[0]+numberOfFluidNodesPerLength
                    localNodes[10] = localNodes[1]+numberOfFluidNodesPerLength
                    localNodes[11] = localNodes[2]+numberOfFluidNodesPerLength
                    localNodes[12] = localNodes[3]+numberOfFluidNodesPerLength
                    localNodes[13] = localNodes[4]+numberOfFluidNodesPerLength
                    localNodes[14] = localNodes[5]+numberOfFluidNodesPerLength
                    localNodes[15] = localNodes[6]+numberOfFluidNodesPerLength
                    localNodes[16] = localNodes[7]+numberOfFluidNodesPerLength
                    localNodes[17] = localNodes[8]+numberOfFluidNodesPerLength
                    localNodes[19] = localNodes[10]+numberOfFluidNodesPerLength
                    localNodes[21] = localNodes[12]+numberOfFluidNodesPerLength
                    localNodes[22] = localNodes[13]+numberOfFluidNodesPerLength
                    localNodes[23] = localNodes[14]+numberOfFluidNodesPerLength
                    localNodes[25] = localNodes[16]+numberOfFluidNodesPerLength
                linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                               localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
                if (debug):
                    print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                    if (uInterpolation == QUADRATIC):
                        print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                fluidPElements.NodesSet(elementNumber,linearNodes)
                fluidUElements.NodesSet(elementNumber,localNodes)
        previousBlock = blockIdx
    #Handle the square block
    if (numberOfSquareElements==1):
        elementNumber = elementNumber + 1
        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx010] = 2*numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx110] = numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        if(uInterpolation == QUADRATIC):
            localNodes[1] = localNodes[localNodeIdx100] - 1
            localNodes[3] = localNodes[localNodeIdx000] - 1
            localNodes[4] = localNodes[localNodeIdx100] + 1
            localNodes[5] = localNodes[localNodeIdx110] - 1
            localNodes[7] = localNodes[localNodeIdx010] - 1
        localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                       localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
        if (uInterpolation == QUADRATIC):
            localNodes[9] = localNodes[0]+numberOfFluidNodesPerLength
            localNodes[10] = localNodes[1]+numberOfFluidNodesPerLength
            localNodes[11] = localNodes[2]+numberOfFluidNodesPerLength
            localNodes[12] = localNodes[3]+numberOfFluidNodesPerLength
            localNodes[13] = localNodes[4]+numberOfFluidNodesPerLength
            localNodes[14] = localNodes[5]+numberOfFluidNodesPerLength
            localNodes[15] = localNodes[6]+numberOfFluidNodesPerLength
            localNodes[16] = localNodes[7]+numberOfFluidNodesPerLength
            localNodes[17] = localNodes[8]+numberOfFluidNodesPerLength
            localNodes[19] = localNodes[10]+numberOfFluidNodesPerLength
            localNodes[21] = localNodes[12]+numberOfFluidNodesPerLength
            localNodes[22] = localNodes[13]+numberOfFluidNodesPerLength
            localNodes[23] = localNodes[14]+numberOfFluidNodesPerLength
            localNodes[25] = localNodes[16]+numberOfFluidNodesPerLength
        if (debug):
            print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                  (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
            if (uInterpolation == QUADRATIC):
                print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                                
        fluidPElements.NodesSet(elementNumber,linearNodes)
        fluidUElements.NodesSet(elementNumber,localNodes)
    else:
        for yElementIdx in range(1,numberOfSquareElements+1):
            for xElementIdx in range(1,numberOfSquareElements+1):
                localNodes = [0]*numberOfLocalNodes
                elementNumber = 4*numberOfFluidElementsPerBlock+xElementIdx+(yElementIdx-1)*numberOfSquareElements+\
                                (zElementIdx-1)*numberOfFluidElementsPerLength
                if (yElementIdx == 1):
                    if (xElementIdx == 1):
                        #Bottom-left
                        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 3*numberOfFluidNodesPerBlock+numberOfArmElements*(numberOfNodesXi-1)*\
                                                      numberOfSquareElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 3*numberOfFluidNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-2)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements*(numberOfNodesXi-1)
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Bottom-right
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)-(numberOfNodesXi-2)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx110] - 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                        elif(uInterpolation == cubic):
                            print("Not implemented.")
                            exit()
                    else:
                        #Bottom
                        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-2)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]+(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                        elif(uInterpolation == cubic):
                            print("Not implemented.")
                            exit()
                elif (yElementIdx == numberOfSquareElements):
                    if (xElementIdx == 1):
                        #Top-left
                        localNodes[localNodeIdx000] = 2*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 2*numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = 2*numberOfFluidNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[1] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] + 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Top-right
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+\
                                                      (numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+\
                                                      (numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx110] - 1
                            localNodes[7] = localNodes[localNodeIdx010] - 1
                    else:
                        #Top
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+\
                                                      (xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 2*numberOfFluidNodesPerBlock-(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]-(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx010] - 1
                else:
                    if (xElementIdx == 1):
                        #Left
                        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock-(yElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]-(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx100]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements*(numberOfNodesXi-1) 
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Right
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = numberOfSquareElements*(numberOfNodesXi-1)*numberOfArmElements*(numberOfNodesXi-1)+\
                                                      (yElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx100]+(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx100] + 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                    else:
                        #Middle
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]+(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                               localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
                if (uInterpolation == QUADRATIC):
                    localNodes[9] = localNodes[0]+numberOfFluidNodesPerLength
                    localNodes[10] = localNodes[1]+numberOfFluidNodesPerLength
                    localNodes[11] = localNodes[2]+numberOfFluidNodesPerLength
                    localNodes[12] = localNodes[3]+numberOfFluidNodesPerLength
                    localNodes[13] = localNodes[4]+numberOfFluidNodesPerLength
                    localNodes[14] = localNodes[5]+numberOfFluidNodesPerLength
                    localNodes[15] = localNodes[6]+numberOfFluidNodesPerLength
                    localNodes[16] = localNodes[7]+numberOfFluidNodesPerLength
                    localNodes[17] = localNodes[8]+numberOfFluidNodesPerLength
                    localNodes[19] = localNodes[10]+numberOfFluidNodesPerLength
                    localNodes[21] = localNodes[12]+numberOfFluidNodesPerLength
                    localNodes[22] = localNodes[13]+numberOfFluidNodesPerLength
                    localNodes[23] = localNodes[14]+numberOfFluidNodesPerLength
                    localNodes[25] = localNodes[16]+numberOfFluidNodesPerLength
                if (debug):
                    print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                    if (uInterpolation == QUADRATIC):
                        print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                
                fluidPElements.NodesSet(elementNumber,linearNodes)
                fluidUElements.NodesSet(elementNumber,localNodes)

fluidUElements.CreateFinish()
fluidPElements.CreateFinish()

fluidMesh.CreateFinish()

if(progressDiagnostics):
    print('Meshes ... Done')    
debug=False
#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')
    
# Create a decomposition for the fluid mesh
fluidDecomposition = iron.Decomposition()
fluidDecomposition.CreateStart(fluidDecompositionUserNumber,fluidMesh)
fluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
fluidDecomposition.CalculateFacesSet(True)
fluidDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')
    
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')
    
# Start to create a default (geometric) field on the fluid region
fluidGeometricField = iron.Field()
fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber,fluidRegion)
# Set the decomposition to use
fluidGeometricField.MeshDecompositionSet(fluidDecomposition)
# Set the scaling to use
fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
fluidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidGeometry')
# Set the domain to be used by the field components.
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
# Finish creating the second field
fluidGeometricField.CreateFinish()


if (progressDiagnostics):
    print('Geometric Field ... Done')
    
if (progressDiagnostics):
    print('Geometric Parameters ...')

armSize = (1.0-squareSizeRatio)*pipeRadius
squareSize = pipeRadius-armSize

if (debug):
    print('  Fluid Nodes:')
maxz = numberOfLengthElements*(numberOfNodesXi-1)+1
top_nodes = []
bottom_nodes = []
side_nodes = []
side_nodes_top = []
side_nodes_bottom = []
for zNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):#1 to 4
    
    #Handle the arm blocks first
    previousBlock = 4
    for blockIdx in range(1,5):
        for yNodeIdx in range(1,numberOfArmElements*(numberOfNodesXi-1)+2):# 1 to 4
            for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):#1 to 3
                nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                             (zNodeIdx-1)*numberOfFluidNodesPerLength#nodeNumber go from 1 upto all nodenumber
               
                nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)#all nodeDomain is zero
                
                if (nodeDomain == computationalNodeNumber):
                    if((zNodeIdx == 1)and(yNodeIdx == 1)):
                        side_nodes_top.append(nodeNumber)
                    elif(zNodeIdx == 1):
                        top_nodes.append(nodeNumber)
                    elif((zNodeIdx == maxz)and(yNodeIdx == 1)):
                        side_nodes_bottom.append(nodeNumber)
                    elif(zNodeIdx == maxz):
                        bottom_nodes.append(nodeNumber)
                    elif(yNodeIdx == 1):
                        side_nodes.append(nodeNumber)

                    zPosition = float(zNodeIdx-1)*lengthSize/float(numberOfLengthElements*(numberOfNodesXi-1))
                    if (yNodeIdx == numberOfArmElements*(numberOfNodesXi-1)+1):#if it is the last loop
                        #On the square
                        if (blockIdx == 1):
                            xPosition = squareSize - xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            
                        elif (blockIdx == 2):
                            xPosition = -xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = (numberOfSquareElements*(numberOfNodesXi-1)-xNodeIdx)*squareSize/ \
                                        (numberOfSquareElements*(numberOfNodesXi-1))
                            
                        elif (blockIdx == 3):
                            xPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize
                            yPosition = -xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            
                        elif (blockIdx == 4):
                            xPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = -(numberOfSquareElements*(numberOfNodesXi-1)-xNodeIdx)*squareSize/ \
                                        (numberOfSquareElements*(numberOfNodesXi-1))
                            
                    else:
                        #In the arm
                        #Work out the r, theta position
                        theta = xNodeIdx*math.pi/(2*numberOfSquareElements*(numberOfNodesXi-1))+(blockIdx-1)*math.pi/2.0
                        fraction = 1.0/(abs(math.sin(theta))+abs(math.cos(theta)))
                        armRadius = armSize+squareSize*(1.0-fraction)
                        radius = (numberOfArmElements*(numberOfNodesXi-1)-yNodeIdx+1)*armRadius/ \
                                 (numberOfArmElements*(numberOfNodesXi-1))+squareSize*fraction
                        xPosition = radius*math.cos(theta)
                        yPosition = radius*math.sin(theta)
                        
                    zPosition = float(zNodeIdx-1)*lengthSize/float(numberOfLengthElements*(numberOfNodesXi-1))
                    
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Position         = [ %.2f, %.2f, %.2f ]' % (zPosition,xPosition,yPosition))

    #Now handle square
    for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = 4*numberOfFluidNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)+\
                         (zNodeIdx-1)*numberOfFluidNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = squareSize*((xNodeIdx-1)-(yNodeIdx-1))/(numberOfSquareElements*(numberOfNodesXi-1))
                yPosition = squareSize*((xNodeIdx-1)+(yNodeIdx-1))/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize
                zPosition = float(zNodeIdx-1)*lengthSize/float(numberOfLengthElements*(numberOfNodesXi-1))
                if(zNodeIdx == 1):
                    top_nodes.append(nodeNumber)
                elif(zNodeIdx == maxz):
                    bottom_nodes.append(nodeNumber)
                #print ('Z2',zPosition)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Position         = [ %.2f, %.2f, %.2f ]' % (zPosition,xPosition,yPosition))
                        
# Update fields            
fluidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


if (progressDiagnostics):
    print('Geometric Parameters ... Done')

# Export results
fluidFields = iron.Fields()
fluidFields.CreateRegion(fluidRegion)
fluidFields.NodesExport("output/TubeFluid","FORTRAN")
fluidFields.ElementsExport("output/TubeFluid","FORTRAN")
fluidFields.Finalise()
#raise SystemExit
#########################################################################################################
#-----------------------------------------------------------------------------------------------------------
#EQUATION SETS
#-----------------------------------------------------------------------------------------------------------
(fluidcoordinateSystemUserNumber,
    fluidRegionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    fluidMeshUserNumber,
    fluidDecompositionUserNumber,
    fluidGeometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    materialFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,13)
'''
# Create a tri-linear lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [3]*3
basis.quadratureLocalFaceGaussEvaluate = True
basis.CreateFinish()
'''
# Create standard Diffusion equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
        iron.EquationsSetTypes.DARCY_EQUATION,
        iron.EquationsSetSubtypes.DARCY_BRINKMAN]
equationsSet.CreateStart(equationsSetUserNumber,fluidRegion,fluidGeometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#DEPENDENT FIELD
#-----------------------------------------------------------------------------------------------------------

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,1)
equationsSet.DependentCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#MATERIAL FIELD
#-----------------------------------------------------------------------------------------------------------

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

## I believe this will change the diffusion coeff
# k and mu?
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,porosity)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,perm_visc)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,beta)

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.01)
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,0.01)
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,0.01)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#PROBLEM
#-----------------------------------------------------------------------------------------------------------

# Create Diffusion equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
        iron.ProblemTypes.DARCY_EQUATION,
        iron.ProblemSubtypes.DARCY_BRINKMAN]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVER
#-----------------------------------------------------------------------------------------------------------

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

## Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = iron.Nodes()
fluidRegion.NodesGet(nodes)


for i in range (0,len(top_nodes)):
    if(pressure_bc):
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(top_nodes[i]),4,iron.BoundaryConditionsTypes.FIXED,Pdrop)
    else:
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(top_nodes[i]),3,iron.BoundaryConditionsTypes.FIXED,v_z)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(top_nodes[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(top_nodes[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)

for i in range (0,len(side_nodes_top)):
    if(pressure_bc):
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_top[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_top[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_top[i]),4,iron.BoundaryConditionsTypes.FIXED,Pdrop)
    else:
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_top[i]),3,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_top[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_top[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)
for i in range (0,len(bottom_nodes)):
    if(pressure_bc):
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(bottom_nodes[i]),4,iron.BoundaryConditionsTypes.FIXED,Pbase)
    else:
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(bottom_nodes[i]),3,iron.BoundaryConditionsTypes.FIXED,v_z)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(bottom_nodes[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(bottom_nodes[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(bottom_nodes[i]),4,iron.BoundaryConditionsTypes.FIXED,Pbase)
for i in range (0,len(side_nodes_bottom)):
    if(pressure_bc):
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),4,iron.BoundaryConditionsTypes.FIXED,Pbase)
    else:
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),3,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes_bottom[i]),4,iron.BoundaryConditionsTypes.FIXED,Pbase)
for i in range (0,len(side_nodes)):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes[i]),3,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes[i]),2,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(side_nodes[i]),1,iron.BoundaryConditionsTypes.FIXED,0.0)



solverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVE
#-----------------------------------------------------------------------------------------------------------

problem.Solve()

#-----------------------------------------------------------------------------------------------------------
#OUTPUT
#-----------------------------------------------------------------------------------------------------------


# Ensure output directories exist
if not os.path.exists('./output'):
    os.makedirs('./output')

# Export results
fields = iron.Fields()
fields.CreateRegion(fluidRegion)
fields.NodesExport("output/StaticDarcy","FORTRAN")
fields.ElementsExport("output/StaticDarcy","FORTRAN")
fields.Finalise()

#dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 1) 
#fig = plt.figure()
a#x = fig.add_subplot(111, projection='3d')




if(pressure_bc):
  velocity = (Pdrop-Pbase)*(perm_visc)/(lengthSize*porosity)
else:
  pressure_drop = lengthSize*v_z*porosity/(perm_visc)

Total_nodes = nodes.NumberOfNodesGet()
x=numpy.zeros(Total_nodes)
y=numpy.zeros(Total_nodes)
z=numpy.zeros(Total_nodes)
p=numpy.zeros(Total_nodes)
vx=numpy.zeros(Total_nodes)
vy=numpy.zeros(Total_nodes)
vz=numpy.zeros(Total_nodes)

for i in range(0,Total_nodes):
    x[i] = fluidGeometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 1)
    y[i] = fluidGeometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 2)
    z[i] = fluidGeometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 3)
    p[i] =dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 4)
    vz[i]=dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 3)
    vy[i]=dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 2)
    vx[i]=dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(i+1), 1)
   
if(pressure_bc):
  print(numpy.mean(vz),velocity)
else:
  pressure_error = (numpy.max(p)-pressure_drop)
  print('Pressure error: ' + str(pressure_error) + ' (Pa)')
  print('Predicted driving pressure ' + str(numpy.max(p)))

if(pressure_bc):
    pic = ax.scatter(x,y,z,c=vz,s=50)
else:
    pic = ax.scatter(x,y,z,c=p,s=50)

#plt.colorbar(pic)
#plt.show()

iron.Finalise()