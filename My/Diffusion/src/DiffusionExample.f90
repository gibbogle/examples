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

!> \example ClassicalField/Diffusion/Diffusion/src/DiffusionExample.f90
!! Example program to solve a diffusion equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Diffusion/Diffusion/history.html
!<

MODULE IOSTUFF
  
  USE OPENCMISS
  USE OpenCMISS_Iron

  IMPLICIT NONE

!GB
    INTEGER(CMISSIntg),parameter :: BASIS_TYPE_SIMPLEX = 1
    INTEGER(CMISSIntg),parameter :: BASIS_TYPE_LAGRANGE = 2

  CONTAINS

  ! --------------------------------------------------------------------------------------------------------------
  !>Reads in a mesh and returns a mesh object. Basically handles all Basis and Mesh creation calls (comment them out entirely)
  ! --------------------------------------------------------------------------------------------------------------
  SUBROUTINE READ_MESH(Filename,basis_type, MeshUserNumber,Region, Mesh,Bases,Nodes,Elements)
    character(len=*), intent(in) :: Filename
    INTEGER(CMISSIntg), intent(in) :: MeshUserNumber, basis_type
    type(cmfe_RegionType), intent(in) :: Region
    type(cmfe_MeshType), intent(inout) :: Mesh
    type(cmfe_BasisType), allocatable, intent(out) :: Bases(:)
    TYPE(cmfe_MeshElementsType), allocatable, intent(out) :: Elements(:)
    TYPE(cmfe_NodesType), intent(out) :: Nodes
    !Local variables
    INTEGER(CMISSIntg),parameter :: fid=77
    character*6, parameter :: MAXFMT='(A255)'
    character(len=255) :: word
    INTEGER(CMISSIntg) :: NumberOfMeshDimensions,NumberOfMeshComponents
    INTEGER(CMISSIntg) :: NumberOfNodes,NumberOfElements,NumberOfBases
    INTEGER(CMISSIntg) :: MeshComponentNumber,MyComputationalNode,i,Err
    INTEGER(CMISSIntg) :: InterpolationType
    INTEGER(CMISSIntg) :: compn,basisn,gaussn, basis_order,el(64),lnn
    
    CALL cmfe_ComputationalNodeNumberGet(MyComputationalNode,Err)

    ! only if root process
    if (MyComputationalNode==0) then
      ! open file
      open(unit=fid, file=Filename, status='old', action='read', err=998)
      CALL cmfe_Mesh_Initialise(Mesh,Err)

      compn=0; basisn=0
      ! read header info
      do
        read(fid,FMT=MAXFMT,end=776) word
        ! skip blanks/comments
        if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
        ! proper keywords start with a #
        if (trim(adjustl(word))=="#number_of_mesh_dimensions") then
          read(fid,*,end=777,err=999) NumberOfMeshDimensions
          CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
          write(*,*) 'NumberOfMeshDimensions: ',NumberOfMeshDimensions
        elseif (trim(adjustl(word))=="#number_of_mesh_components") then
          read(fid,*,end=777,err=999) NumberOfMeshComponents
          CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
          allocate(Elements(NumberOfMeshComponents))
          write(*,*) 'NumberOfMeshComponents: ',NumberOfMeshComponents
        elseif (trim(adjustl(word))=="#number_of_nodes") then
          read(fid,*,end=777,err=999) NumberOfNodes
          !Define nodes for the mesh
          CALL cmfe_Nodes_Initialise(Nodes,Err)
          CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
          CALL cmfe_Nodes_CreateFinish(Nodes,Err)  
        elseif (trim(adjustl(word))=="#number_of_elements") then
          read(fid,*,end=777,err=999) NumberOfElements
          CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)  
        elseif (trim(adjustl(word))=="#number_of_bases") then
          read(fid,*,end=777,err=999) NumberOfBases
          allocate(Bases(NumberOfBases))
        elseif (trim(adjustl(word))=="#mesh_component") then
          ! start of a mesh component block
          compn=compn+1
          if (compn>NumberOfMeshComponents) then
            write(*,*) "READ_MESH: incorrect number of mesh components are defined"
            close(fid)
            return
          endif
          read(fid,*,end=777,err=999) MeshComponentNumber
          write(*,*) 'MeshComponentNumber: ',MeshComponentNumber
          do
            read(fid,FMT=MAXFMT,end=777,err=999) word
            if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
            if (trim(adjustl(word))/="#basis_order") then
              write(*,*) "READ_MESH: #basis_order must follow the mesh_component definition"
              close(fid)
              return
            endif
            exit
          enddo
          read(fid,*,end=777,err=999) basis_order
          write(*,*) 'basis_order: ',basis_order
          basisn=basisn+1
          InterpolationType=0; gaussn=0; lnn=0
          if (basis_type == BASIS_TYPE_LAGRANGE) then
              ! set up basis: hardcoded for lagrange basis type
              CALL cmfe_Basis_Initialise(Bases(basisn),Err)
              CALL cmfe_Basis_CreateStart(basisn,Bases(basisn),Err) 
              CALL cmfe_Basis_TypeSet(Bases(basisn),cmfe_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
              CALL cmfe_Basis_NumberOfXiSet(Bases(basisn),NumberOfMeshDimensions,Err)
              select case (basis_order)
              case (1)
                InterpolationType=cmfe_BASIS_LINEAR_LAGRANGE_INTERPOLATION
                gaussn=3; lnn=8
              case (2)
                InterpolationType=cmfe_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
                gaussn=3; lnn=27
              case (3)
                InterpolationType=cmfe_BASIS_CUBIC_LAGRANGE_INTERPOLATION
                gaussn=3; lnn=64
              end select
              CALL cmfe_Basis_InterpolationXiSet(Bases(basisn),[InterpolationType,InterpolationType,InterpolationType],Err)
              CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Bases(basisn),[gaussn,gaussn,gaussn],Err)  
          elseif (basis_type == BASIS_TYPE_SIMPLEX) then
              CALL cmfe_Basis_Initialise(Bases(basisn),Err)
              CALL cmfe_Basis_CreateStart(basisn,Bases(basisn),Err) 
              CALL cmfe_Basis_TypeSet(Bases(basisn),cmfe_BASIS_SIMPLEX_TYPE,Err)
              CALL cmfe_Basis_NumberOfXiSet(Bases(basisn),NumberOfMeshDimensions,Err)
              InterpolationType=cmfe_BASIS_LINEAR_SIMPLEX_INTERPOLATION
              gaussn = 3
              write(*,*) 'cmfe_Basis_InterpolationXiSet'
              if (NumberOfMeshDimensions == 1) then
                  lnn = 2
                  CALL cmfe_Basis_InterpolationXiSet(Bases(basisn),[InterpolationType],Err)
              elseif (NumberOfMeshDimensions == 2) then
                  lnn = 3
                  CALL cmfe_Basis_InterpolationXiSet(Bases(basisn),[InterpolationType,InterpolationType],Err)
                  write(*,*) 'did cmfe_Basis_InterpolationXiSet'
              elseif (NumberOfMeshDimensions == 3) then
                  lnn = 4
                  CALL cmfe_Basis_InterpolationXiSet(Bases(basisn),[InterpolationType,InterpolationType,InterpolationType],Err)
                  write(*,*) 'cmfe_Basis_QuadratureOrderSet'
                  CALL cmfe_Basis_QuadratureOrderSet(Bases(basisn),3,Err)
              endif
          endif
          if (NumberOfMeshDimensions == 3) then
              CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(Bases(basisn),.true.,Err)
          endif
          CALL cmfe_Basis_CreateFinish(Bases(basisn),Err)
          ! element definition now
          do
            read(fid,FMT=MAXFMT,end=777,err=999) word
            if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
            if (trim(adjustl(word))/="#elements") then
              write(*,*) "READ_MESH: #elements must follow the basis definition"
              close(fid)
              return
            endif
            exit
          enddo
          write(*,*) 'Creating elements'
          CALL cmfe_MeshElements_Initialise(Elements(compn),Err)
          CALL cmfe_MeshElements_CreateStart(Mesh,compn,Bases(basisn),Elements(compn),Err)
          do i=1,NumberOfElements
            read(fid,*,end=777,err=999) el(1:lnn)
            write(*,*) 'cmfe_MeshElements_NodesSet: elem: ',i,el(1:lnn)
            CALL cmfe_MeshElements_NodesSet(Elements(compn),i,el(1:lnn),Err)
          enddo
          CALL cmfe_MeshElements_CreateFinish(Elements(compn),Err)
        elseif (trim(adjustl(word))=="#number_of_surfaces") then
          write(*,*) "READ_MESH: not implemented yet"
          close(fid)
          return
        endif
      enddo
    endif

776 CALL cmfe_Mesh_CreateFinish(Mesh,Err)
    close(fid)
    return ! happy return
777 write(*,*) "READ_MESH: unexpected end of file encountered"
    close(fid)
    return
998 write(*,*) "READ_MESH: could not open file "//Filename
    close(fid)
    return
999 write(*,*) "READ_MESH: error reading line"
    close(fid)
    return
  END SUBROUTINE READ_MESH

! --------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------
  SUBROUTINE READ_NODES(Filename,GeometricField)
    character(len=*), intent(in) :: Filename
    type(cmfe_FieldType), intent(inout) :: GeometricField
    !Local variables
    INTEGER(CMISSIntg),parameter :: fid=79
    character*6, parameter :: MAXFMT='(A255)'
    character(len=255) :: word
    INTEGER(CMISSIntg) :: NumberOfNodes,NumberOfCoordinateDimensions
    INTEGER(CMISSIntg) :: MyComputationalNode,i,j,Err
    REAL(CMISSRP) :: coord(3)

    ! skipping all dimension, size or otherwise error checks

    CALL cmfe_ComputationalNodeNumberGet(MyComputationalNode,Err)    

    if (MyComputationalNode==0) then
      ! open the file
      open(unit=fid,file=Filename,status='old',action='read',err=998)

      ! read some header info
      do
        read(fid,FMT=MAXFMT,end=776,err=999) word
        ! skip blanks/comments
        if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
        if (trim(adjustl(word))=="#number_of_nodes") then
          read(fid,*,end=777,err=999) NumberOfNodes
        elseif (trim(adjustl(word))=="#number_of_coordinate_dimensions") then
          read(fid,*,end=777,err=999) NumberOfCoordinateDimensions
        elseif (trim(adjustl(word))=="#nodes") then
          do i=1,NumberOfNodes
            read(fid,*,err=999) coord(1:NumberOfCoordinateDimensions)
            do j=1,NumberOfCoordinateDimensions
              CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
                  1,1,i,j,coord(j),Err)
            enddo
          enddo
          exit
        endif
      enddo
    endif

776 close(fid)
    return ! happy return
777 write(*,*) "READ_NODES: unexpected end of file encountered"
    close(fid)
    return
998 write(*,*) "READ_NODES: could not open file "//Filename
    close(fid)
    return
999 write(*,*) "READ_NODES: error reading line"
    close(fid)
    return
  END SUBROUTINE READ_NODES

END MODULE IOSTUFF


! --------------------------------------------------------------------------------------------------------------
!> Main program
! --------------------------------------------------------------------------------------------------------------
PROGRAM DIFFUSIONEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  use iostuff    ! GB

#ifdef WIN32
  USE IFQWIN
#endif

#define GENERATE_MESH 0
#define USE_SOURCE 1

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP
  
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
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=14

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  INTEGER(CMISSIntg) :: MPI_IERROR
  LOGICAL :: CASE_2D
  INTEGER(CMISSIntg) :: node  !GB
  character*(255) :: meshname
 
    !CMISS variables

! GB
#if GENERATE_MESH
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh 
  TYPE(cmfe_FieldType) :: AnalyticField 
#else
  INTEGER(CMISSIntg) :: FirstNodeNumber,FirstNodeDomain,LastNodeNumber,LastNodeDomain,NumberOfNodes
  INTEGER(CMISSIntg) :: INTERPOLATION_TYPE
  TYPE(cmfe_BasisType),allocatable :: Bases(:)
  TYPE(cmfe_MeshElementsType),allocatable :: Elements(:) 
  TYPE(cmfe_NodesType) :: Nodes
  REAL(CMISSRP) :: source_value
#endif

!  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField, SourceField  !GB added sourcefield
  TYPE(cmfe_FieldsType) :: Fields
!  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver, LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions

  LOGICAL :: EXPORT_FIELD

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
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

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

#if GENERATE_MESH
  NUMBER_GLOBAL_X_ELEMENTS=10
  NUMBER_GLOBAL_Y_ELEMENTS=10
  NUMBER_GLOBAL_Z_ELEMENTS=0
  CASE_2D = (NUMBER_GLOBAL_Z_ELEMENTS == 0)
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
#else
!   meshname = 'example1'
   write(*,*) 'Enter the base name of the mesh files (e.g. zzzz if the files are zzzz-mesh, zzzz-nodes)'
   read(*,'(a)') meshname
   CASE_2D = .false.    ! GB
   INTERPOLATION_TYPE = 0
#endif
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diagnostic levels on for testing

  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
!  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
  IF (CASE_2D) then
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
  
#if GENERATE_MESH
  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
!  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
  IF (CASE_2D) then
    !Set the basis to be a bilinear Lagrange basis
    !CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
!    CALL cmfe_Basis_InterpolationXiSet(Basis,[3,3],Err)
!    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[4,4],Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

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
#else
  write(*,*) 'READ_MESH'
  call READ_MESH(trim(meshname)//'-mesh',BASIS_TYPE_SIMPLEX, MeshUserNumber,Region, Mesh, Bases, Nodes, Elements)
  write(*,*) 'did READ_MESH'
#endif

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
!  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
  IF (.not.CASE_2D) then
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  
#if GENERATE_MESH
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
#else
  write(*,*) 'READ_NODES'
  CALL READ_NODES(trim(meshname)//'-nodes',GeometricField)
  write(*,*) 'did READ_NODES'
#endif
 
#if USE_SOURCE
  !Create the diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
!    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE], &
    & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Set the equations set to be a standard Diffusion constant source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

#else

!  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
  & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE],EquationsSetFieldUserNumber, &
  & EquationsSetField,EquationsSet,Err)
!  !Set the equations set to be a standard Diffusion problem
  
!  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)
#endif

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

#if USE_SOURCE
  ! Copied from ReactionDiffusionConstantSource1D.f90
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)
  !Set up source field for reaction diffusion equation set. Setting source field to be 0.5 for this test across
  !across all nodes
  !CALL cmfe_FieldComponentValuesInitialise(SourceField,cmfe_FieldUVariableType,cmfe_FieldValuesSetType, &
  ! & 1,0.5_CMISSRP,Err)
  !Set up source field for reaction diffusion as a nodally varying field. Zeros at the two ends, and 0.5 at middle node
  !node=2
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,NumberOfNodes,Err)
  write(*,*) 'NumberOfNodes: ',NumberOfNodes
  do node = 1, NumberOfNodes
    source_value = 0.0_CMISSRP
    if (node == 2) source_value = 5.0_CMISSRP
      CALL cmfe_Field_ParameterSetUpdateNode(RegionUserNumber,SourceFieldUserNumber, &
   &   CMFE_FIELD_U_VARIABLE_TYPE, &
   &   CMFE_FIELD_VALUES_SET_TYPE, &
   &   1,1,node,1,source_value,Err)
  enddo
#endif

#if GENERATE_MESH && !USE_SOURCE
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(AnalyticField,Err)
!  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN  
  IF (CASE_2D) then
    CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1, &
      & AnalyticFieldUserNumber, &
      & AnalyticField,Err)
  ELSE
    WRITE(*,'(A)') "Three dimensions is not implemented."
    STOP
  ENDIF
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)
#endif
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Create the equations set boundary conditions
  !Find the first and last dof numbers and ranks
!   NULLIFY(FIELD_VARIABLE)
!   CALL FIELD_VARIABLE_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
!   DEPENDENT_DOF_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
!   first_global_dof=1
!   first_local_dof=0
!   first_local_rank=0
!   DO rank_idx=1,DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%NUMBER_OF_DOMAINS
!     IF(DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
!       first_local_dof=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%LOCAL_NUMBER(rank_idx)
!       first_local_rank=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%DOMAIN_NUMBER(rank_idx)
!       EXIT
!     ENDIF
!   ENDDO !rank_idx  
!   NULLIFY(FIELD_VARIABLE)
!   CALL FIELD_VARIABLE_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
!   DEPENDENT_DOF_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
!   last_global_dof=DEPENDENT_DOF_MAPPING%NUMBER_OF_GLOBAL
!   last_local_dof=0
!   last_local_rank=0
!   DO rank_idx=1,DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%NUMBER_OF_DOMAINS
!     IF(DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
!       last_local_dof=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%LOCAL_NUMBER(rank_idx)
!       last_local_rank=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%DOMAIN_NUMBER(rank_idx)
!       EXIT
!     ENDIF
!   ENDDO !rank_idx
!   NULLIFY(BOUNDARY_CONDITIONS)
!   CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
!   IF(MY_COMPUTATIONAL_NODE_NUMBER==first_local_rank) &
!     & CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,first_local_dof, &
!     & BOUNDARY_CONDITION_FIXED,1.0_DP,ERR,ERROR,*999)
!   IF(MY_COMPUTATIONAL_NODE_NUMBER==last_local_rank) &
!     & CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD_DELUDELN_VARIABLE_TYPE,last_local_dof, &
!     & BOUNDARY_CONDITION_FIXED,1.0_DP,ERR,ERROR,*999)
!   CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Create the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_DIFFUSION_EQUATION_TYPE, &
    & CMFE_PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
#if GENERATE_MESH
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,1.001_CMISSRP,0.001_CMISSRP,Err)
#else
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,1.01_CMISSRP,0.01_CMISSRP,Err)
#endif
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)


  !Start the creation of the problem solvers

! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
! 
!   CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)
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
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Create the solver equations boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
#if GENERATE_MESH
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
#else
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  ENDIF
#endif
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

#if GENERATE_MESH && !USE_SOURCE
  !Output Analytic analysis
  Call cmfe_AnalyticAnalysis_Output(DependentField,"DiffusionAnalytics_x4_y4_q_T1",Err)
  meshname = "Diffusion_x4_y4_q_T1"
#endif

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,trim(meshname),"FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,trim(meshname),"FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

  ENDIF
  
  !Output timing summary
  !CALL TIMING_SUMMARY_OUTPUT(ERR,ERROR,*999)

  !Calculate the stop times and write out the elapsed user and system times
!   CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
!   CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)
! 
!   CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
!     & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)
!   
  !CALL CMFE_FINALISE(ERR,ERROR,*999)
  !CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM DIFFUSIONEXAMPLE
