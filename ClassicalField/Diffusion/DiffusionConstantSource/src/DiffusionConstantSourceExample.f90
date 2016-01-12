!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example ClassicalField/Diffusion/DiffusionConstantSource/src/DiffusionConstantSourceExample.f90
!! Example program to solve a diffusion equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Diffusion/DiffusionConstantSource/history.html
!<

!> Main program
PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=14

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR
  
  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField, EquationsSetField, DependentField, MaterialsField, SourceField, AnalyticField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver, LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
 
#define IS_ANALYTIC 1

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: Err
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  NUMBER_GLOBAL_X_ELEMENTS=10
  NUMBER_GLOBAL_Y_ELEMENTS=10
  NUMBER_GLOBAL_Z_ELEMENTS=10
  NUMBER_OF_DOMAINS=1


  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the coordinate system to be 2D
      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
    ELSE
      !Set the coordinate system to be 3D
      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
    ENDIF
    !Finish the creation of the coordinate system
    CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)


    !Start the creation of the region
    CALL cmfe_Region_Initialise(Region,Err)
    CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL cmfe_Region_CreateFinish(Region,Err)

    !Start the creation of a basis (default is trilinear lagrange)
    CALL cmfe_Basis_Initialise(Basis,Err)
    CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear Lagrange basis
      CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
    ELSE
      !Set the basis to be a trilinear Lagrange basis
      CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BASIS,Err)

    !Start the creation of a generated mesh in the region
    CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular x*y*z mesh
    CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
    !Set the default basis
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
        & NUMBER_GLOBAL_Z_ELEMENTS],Err)
    ENDIF    
    !Finish the creation of a generated mesh in the region
    CALL cmfe_Mesh_Initialise(Mesh,Err)
    CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

    !Create a decomposition
    CALL cmfe_Decomposition_Initialise(Decomposition,Err)
    CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    !Finish the decomposition
    CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

       
    !Update the geometric field parameters
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
!  ENDIF

  !IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR
  
  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
  & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE], &
  & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem
  
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set source field variables
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)

#if IS_ANALYTIC
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(AnalyticField,Err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1, &
    & AnalyticFieldUserNumber, &
    & AnalyticField,Err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)
#endif

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

!  write(*,*) 'cmfe_BoundaryConditions_Initialise'
!  !Create the equations set boundary conditions
!   CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
!   CALL cmfe_EquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
!   !Set the first node to 0.0 and the last node to 1.0
! !   FirstNodeNumber=1
! !   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
! !     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
! !   ELSE
! !     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
! !   ENDIF
! !   CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,FirstNodeNumber,1, &
! !     & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
! !   CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,LastNodeNumber,1, &
! !     & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
!   !Finish the creation of the equations set boundary conditions
!   CALL cmfe_EquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

  !Create the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_DIFFUSION_EQUATION_TYPE, &
    & CMFE_PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)


  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,ControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,1.0001_CMISSRP,0.01_CMISSRP,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)


  !Start the creation of the problem solvers
!  
! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
! 

  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,10000,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

#if IS_ANALYTIC
  write(*,*) 'cmfe_BoundaryConditions_Initialise analytic'
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
!  CALL cmfe_BoundaryConditions_AddConstant( &
!     BoundaryConditions, & !/*<The boundary conditions to add the constant to. */,                                                               const cmfe_BoundaryConditionsType
!     field, &              !/*<The dependent field to set the boundary condition on. */,                                                         const cmfe_FieldType
!     variableType, &       !/*<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes */,  const int 
!     componentNumber, &    !/*<The component number of the dependent field to set the boundary condition at. */,                                 const int
!     condition, &          !/*<The boundary condition type to set \see OPENCMISS_BoundaryConditionsTypes,OPENCMISS */,                           const int
!     value)                !/*<The value of the boundary condition to add. */);                                                                  const double
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
#else
  write(*,*) 'cmfe_BoundaryConditions_Initialise'
  !Create the equations set boundary conditions
   CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
   CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
   !Set the first node to 0.0 and the last node to 1.0
   FirstNodeNumber=1
   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
   ELSE
     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
   ENDIF
!    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,FirstNodeNumber,1, &
!      CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
!    CALL cmfe_BoundaryConditions_AddConstant(BoundaryConditions, DependentField, CMFE_FIELD_U_VARIABLE_TYPE, &
!      1, CMFE_BOUNDARY_CONDITION_FIXED, 1.0, Err)
 !   CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,LastNodeNumber,1, &
 !     & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

#endif

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

#if IS_ANALYTIC
  !Output Analytic analysis
  CALL cmfe_AnalyticAnalysis_Output(DependentField,"DiffusionAnalytics_x10_y10_z10_L_T1",Err)
#endif

  EXPORT_FIELD=.FALSE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
#if IS_ANALYTIC
    CALL cmfe_Fields_NodesExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#else
    CALL cmfe_Fields_NodesExport(Fields,"DiffusionConstantSource","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"DiffusionConstantSource","FORTRAN",Err)
#endif
    CALL cmfe_Fields_Finalise(Fields,Err)

  ENDIF

  !CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
 
END PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE
