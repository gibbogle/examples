!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/Laplace/Fortran/src/LaplaceExample.f90
!! Example program to solve a Laplace equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/Laplace/history.html
!!
!<

MODULE IOSTUFF
  
  USE OPENCMISS
  IMPLICIT NONE

  CONTAINS

  ! --------------------------------------------------------------------------------------------------------------
  !>Reads in a mesh and returns a mesh object. Basically handles all Basis and Mesh creation calls (comment them out entirely)
  SUBROUTINE READ_MESH(Filename,MeshUserNumber,Region, Mesh,Bases,Nodes,Elements)
    character(len=*), intent(in) :: Filename
    INTEGER(CMISSIntg), intent(in) :: MeshUserNumber
    type(CMISSRegionType), intent(in) :: Region
    type(CMISSMeshType), intent(inout) :: Mesh
    type(CMISSBasisType), allocatable, intent(out) :: Bases(:)
    TYPE(CMISSMeshElementsType), allocatable, intent(out) :: Elements(:)
    TYPE(CMISSNodesType), intent(out) :: Nodes
    !Local variables
    INTEGER(CMISSIntg),parameter :: fid=77
    character*6, parameter :: MAXFMT='(A255)'
    character(len=255) :: word
    INTEGER(CMISSIntg) :: NumberOfMeshDimensions,NumberOfMeshComponents
    INTEGER(CMISSIntg) :: NumberOfNodes,NumberOfElements,NumberOfBases
    INTEGER(CMISSIntg) :: MeshComponentNumber,MyComputationalNode,i,Err
    INTEGER(CMISSIntg) :: InterpolationType
    INTEGER(CMISSIntg) :: compn,basisn,gaussn, basis_order,el(64),lnn
    
    CALL CMISSComputationalNodeNumberGet(MyComputationalNode,Err)

    ! only if root process
    if (MyComputationalNode==0) then
      ! open file
      open(unit=fid, file=Filename, status='old', action='read', err=998)
      CALL CMISSMesh_Initialise(Mesh,Err)

      compn=0; basisn=0
      ! read header info
      do
        read(fid,FMT=MAXFMT,end=776) word
        ! skip blanks/comments
        if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
        ! proper keywords start with a #
        if (trim(adjustl(word))=="#number_of_mesh_dimensions") then
          read(fid,*,end=777,err=999) NumberOfMeshDimensions
          CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
        elseif (trim(adjustl(word))=="#number_of_mesh_components") then
          read(fid,*,end=777,err=999) NumberOfMeshComponents
          CALL CMISSMesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
          allocate(Elements(NumberOfMeshComponents))
        elseif (trim(adjustl(word))=="#number_of_nodes") then
          read(fid,*,end=777,err=999) NumberOfNodes
          !Define nodes for the mesh
          CALL CMISSNodes_Initialise(Nodes,Err)
          CALL CMISSNodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
          CALL CMISSNodes_CreateFinish(Nodes,Err)  
        elseif (trim(adjustl(word))=="#number_of_elements") then
          read(fid,*,end=777,err=999) NumberOfElements
          CALL CMISSMesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)  
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
          ! set up basis: hardcoded for lagrange basis type
          basisn=basisn+1
          CALL CMISSBasis_Initialise(Bases(basisn),Err)
          CALL CMISSBasis_CreateStart(basisn,Bases(basisn),Err) 
          CALL CMISSBasis_TypeSet(Bases(basisn),CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
          CALL CMISSBasis_NumberOfXiSet(Bases(basisn),NumberOfMeshDimensions,Err)
          select case (basis_order)
          case (1)
            InterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
            gaussn=3; lnn=8
          case (2)
            InterpolationType=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
            gaussn=3; lnn=27
          case (3)
            InterpolationType=CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION
            gaussn=3; lnn=64
          end select
          CALL CMISSBasis_InterpolationXiSet(Bases(basisn),[InterpolationType,InterpolationType,InterpolationType],Err)
          CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Bases(basisn),[gaussn,gaussn,gaussn],Err)  
          CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(Bases(basisn),.true.,Err)
          CALL CMISSBasis_CreateFinish(Bases(basisn),Err)
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
          CALL CMISSMeshElements_Initialise(Elements(compn),Err)
          CALL CMISSMeshElements_CreateStart(Mesh,compn,Bases(basisn),Elements(compn),Err)
          do i=1,NumberOfElements
            read(fid,*,end=777,err=999) el(1:lnn)
            CALL CMISSMeshElements_NodesSet(Elements(compn),i,el(1:lnn),Err)
          enddo
          CALL CMISSMeshElements_CreateFinish(Elements(compn),Err)
        elseif (trim(adjustl(word))=="#number_of_surfaces") then
          write(*,*) "READ_MESH: not implemented yet"
          close(fid)
          return
        endif
      enddo
    endif

776 CALL CMISSMesh_CreateFinish(Mesh,Err)
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

  SUBROUTINE READ_NODES(Filename,GeometricField)
    character(len=*), intent(in) :: Filename
    type(CMISSFieldType), intent(inout) :: GeometricField
    !Local variables
    INTEGER(CMISSIntg),parameter :: fid=79
    character*6, parameter :: MAXFMT='(A255)'
    character(len=255) :: word
    INTEGER(CMISSIntg) :: NumberOfNodes,NumberOfCoordinateDimensions
    INTEGER(CMISSIntg) :: MyComputationalNode,i,j,Err
    REAL(CMISSDP) :: coord(3)

    ! skipping all dimension, size or otherwise error checks

    CALL CMISSComputationalNodeNumberGet(MyComputationalNode,Err)    

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
            CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,i,j, &
              & coord(j),Err)
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

  ! --------------------------------------------------------------------------------------------------------------

  !>Reads in a field. Only works for a material field currently, all components using same interpolation. handles full field definition code block.
  SUBROUTINE READ_FIELD(Filename,FieldUserNumber,Region,GeometricField, Field)
    character(len=*), intent(in) :: Filename
    INTEGER(CMISSIntg), intent(in) :: FieldUserNumber
    type(CMISSRegionType), intent(in) :: Region
    type(CMISSFieldType), intent(in) :: GeometricField
    type(CMISSFieldType), intent(out) :: Field
    !Local variables
    TYPE(CMISSDecompositionType) :: Decomposition
    INTEGER(CMISSIntg), parameter :: fid=78
    character*6, parameter :: MAXFMT='(A255)'
    character(len=255) :: word,field_type,interpolation_type,data_type
    INTEGER(CMISSIntg) :: NumberOfVariables,NumberOfComponents
    INTEGER(CMISSIntg) :: varn,mcompn,VariableType,InterpolationType,DataType
    INTEGER(CMISSIntg) :: MyComputationalNode,i,ind,Err
    INTEGER(CMISSIntg) :: data_int(100)
    REAL(CMISSDP) :: data_dp(100)

    CALL CMISSComputationalNodeNumberGet(MyComputationalNode,Err)    

    if (MyComputationalNode==0) then
      ! Open the file see if everything is fine
      open(unit=fid,file=Filename,status='old',action='read',err=998)

      ! initial field setup
!       CALL CMISSField_Initialise(Field,Err)
!       CALL CMISSField_CreateStart(FieldUserNumber,Region,Field,Err)
!       ! --- set material type here then do the below ---
!       CALL CMISSDecomposition_Initialise(Decomposition,Err)
!       CALL CMISSField_MeshDecompositionGet(GeometricField,Decomposition,Err)
!       CALL CMISSField_MeshDecompositionSet(Field,Decomposition,Err)
!       CALL CMISSField_GeometricFieldSet(Field,GeometricField,Err)

      do
        read(fid,FMT=MAXFMT,end=776) word
        ! skip blanks/comments
        if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
        ! proper keywords start with a #
        if (trim(adjustl(word))=="#field_type") then
          read(fid,FMT=MAXFMT,end=777) field_type
          select case (trim(field_type))
          case ("material")
            CALL CMISSField_Initialise(Field,Err)
            CALL CMISSField_CreateStart(FieldUserNumber,Region,Field,Err)
            CALL CMISSField_TypeSet(Field,CMISS_FIELD_MATERIAL_TYPE,Err)
            CALL CMISSDecomposition_Initialise(Decomposition,Err)
            CALL CMISSField_MeshDecompositionGet(GeometricField,Decomposition,Err)
            CALL CMISSField_MeshDecompositionSet(Field,Decomposition,Err)
            CALL CMISSField_GeometricFieldSet(Field,GeometricField,Err)
          case default
            write(*,*) "READ_FIELD: field types other than material has not been implemented"
            close(fid)
            return
          end select
        elseif (trim(adjustl(word))=="#number_of_variables") then
          read(fid,*,end=777) NumberOfVariables
          CALL CMISSField_NumberOfVariablesSet(Field,NumberOfVariables,Err)
        elseif (trim(adjustl(word))=="#variable") then
          read(fid,*,end=777) varn
          select case (varn)
          case (1)
            VariableType=CMISS_FIELD_U_VARIABLE_TYPE
          case (2)
            VariableType=CMISS_FIELD_V_VARIABLE_TYPE
          case (3)
            VariableType=CMISS_FIELD_U1_VARIABLE_TYPE
          case (4)
            VariableType=CMISS_FIELD_U2_VARIABLE_TYPE
          case default
            write(*,*) "READ_FIELD: more than 4 variables, take care of it"
          end select
          ! start a variable block
          do
            read(fid,FMT=MAXFMT,end=777) word
            ! skip blanks/comments
            if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
            if (trim(adjustl(word))=="#number_of_components") then
              read(fid,*,end=777,err=999) NumberOfComponents
              CALL CMISSField_NumberOfComponentsSet(Field,VariableType,NumberOfComponents,Err) 
            elseif (trim(adjustl(word))=="#interpolation_type") then
              read(fid,FMT=MAXFMT,end=777,err=999) interpolation_type
              select case (trim(adjustl(interpolation_type)))
              case ("constant")
                InterpolationType=CMISS_FIELD_CONSTANT_INTERPOLATION
              case ("elemental")
                InterpolationType=CMISS_FIELD_ELEMENT_BASED_INTERPOLATION
              case ("nodal")
                InterpolationType=CMISS_FIELD_NODE_BASED_INTERPOLATION
              case default
                write(*,*) "READ_FIELD: unsupported interpolation type encountered"
                close(fid)
                return
              end select
              ! hardcoded: assume all components are interpolated the same way
              do i=1,NumberOfComponents
                CALL CMISSField_ComponentInterpolationSet(Field,VariableType,i,InterpolationType,Err)
              enddo
            elseif (trim(adjustl(word))=="#mesh_component") then
              read(fid,*,end=777,err=999) mcompn
!               if (InterpolationType==CMISS_FIELD_NODE_BASED_INTERPOLATION) then
                do i=1,NumberOfComponents
                  CALL CMISSField_ComponentMeshComponentSet(Field,VariableType,i,mcompn,Err)
                enddo
!               endif
            elseif (trim(adjustl(word))=="#data_type") then
              read(fid,*,end=777,err=999) data_type
              select case (trim(data_type))
              case ("dp")
                DataType=CMISS_FIELD_DP_TYPE
              case ("intg")
                DataType=CMISS_FIELD_INTG_TYPE
              case default
                write(*,*) "READ_FIELD: unsupported data type encountered"
                close(fid)
                return
              end select
              CALL CMISSField_DataTypeSet(Field,VariableType,DataType,Err)
              CALL CMISSField_CreateFinish(Field,Err)
            elseif (trim(adjustl(word))=="#data") then
              ! don't know how many lines there will be - just go until blank line
              do
                read(fid,FMT=MAXFMT,end=777,err=999) word
                if (len_trim(word)==0) exit ! exit if blank line
                select case (DataType)
                case (CMISS_FIELD_DP_TYPE)
                  read(word,*,end=777,err=999) ind,data_dp(1:NumberOfComponents)
                  do i=1,NumberOfComponents
                    select case (InterpolationType)
                    case (CMISS_FIELD_CONSTANT_INTERPOLATION)
                      CALL CMISSField_ParameterSetUpdateConstant(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,i,data_dp(i),Err)
                    case (CMISS_FIELD_ELEMENT_BASED_INTERPOLATION)
                      CALL CMISSField_ParameterSetUpdateElement(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,ind,i,data_dp(i),Err)
                    case (CMISS_FIELD_NODE_BASED_INTERPOLATION)
                      CALL CMISSField_ParameterSetUpdateNode(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,1,1,ind,i,data_dp(i), &
                        & Err)
                    end select
                  enddo
                case (CMISS_FIELD_INTG_TYPE)
                  read(word,*,end=777,err=999) ind,data_int(1:NumberOfComponents)
                  do i=1,NumberOfComponents
                    select case (InterpolationType)
                    case (CMISS_FIELD_CONSTANT_INTERPOLATION)
                      CALL CMISSField_ParameterSetUpdateConstant(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,i,data_int(i),Err)
                    case (CMISS_FIELD_ELEMENT_BASED_INTERPOLATION)
                      CALL CMISSField_ParameterSetUpdateElement(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,ind,i,data_int(i), &
                        & Err)
                    case (CMISS_FIELD_NODE_BASED_INTERPOLATION)
                      CALL CMISSField_ParameterSetUpdateNode(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,1,1,ind,i,data_int(i), &
                        & Err)
                    end select
                  enddo
                end select
              enddo
            elseif (trim(adjustl(word))=="#data_all") then  ! debug-option to initialise all field values in one go, NO INDEX
              read(fid,FMT=MAXFMT,end=777,err=999) word
              select case (DataType)
              case (CMISS_FIELD_DP_TYPE)
                read(word,*,end=777,err=999) data_dp(1:NumberOfComponents)
                do i=1,NumberOfComponents
                  CALL CMISSField_ComponentValuesInitialise(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE, &
                    & i,data_dp(i),Err)
                enddo
              case (CMISS_FIELD_INTG_TYPE)
                read(word,*,end=777,err=999) data_int(1:NumberOfComponents)
                do i=1,NumberOfComponents
                  CALL CMISSField_ComponentValuesInitialise(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE, &
                    & i,data_int(i),Err)
                enddo
              end select
            endif
          enddo
        endif
      enddo
    endif

776 close(fid)
    return ! happy return
777 write(*,*) "READ_FIELD: unexpected end of file encountered"
    close(fid)
    return
998 write(*,*) "READ_FIELD: could not open file "//Filename
    close(fid)
    return
999 write(*,*) "READ_FIELD: error reading line"
    close(fid)
    return
  END SUBROUTINE READ_FIELD

  ! --------------------------------------------------------------------------------------------------------------
  !>read in a bunch of node indices (which form a surface)
  SUBROUTINE READ_SURFACE(Filename,surface)
    character(len=*), intent(in) :: Filename
    INTEGER(CMISSIntg), allocatable, intent(inout) :: surface(:)
    !Local variables
    INTEGER(CMISSIntg),parameter :: fid=79
    character*6,parameter :: MAXFMT='(A255)'
    character(len=255) :: word
    INTEGER(CMISSIntg) :: nn,i,MyComputationalNode,Err

    CALL CMISSComputationalNodeNumberGet(MyComputationalNode,Err)

    if (MyComputationalNode==0) then
      ! open file
      open(unit=fid,file=Filename,status='old',action='read',err=998)
      do
        read(fid,FMT=MAXFMT,end=776) word
        ! skip blanks/comments
        if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
        ! proper keywords
        if (trim(adjustl(word))=="#number_of_nodes") then
          read(fid,*,end=777,err=999) nn
          if (allocated(surface)) deallocate(surface)
          allocate(surface(nn))
        elseif (trim(adjustl(word))=="#nodes") then
          ! read the numbers
          do i=1,nn
            read(fid,*,end=777,err=999) surface(i)
          enddo
        endif
      enddo
    endif

776 close(fid)
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

  END SUBROUTINE READ_SURFACE

END MODULE IOSTUFF

!--------------------------------------------------------------------------------------------------------------
!> Main program
!--------------------------------------------------------------------------------------------------------------
PROGRAM LAPLACEEXAMPLE

  USE OPENCMISS
  USE MPI
  use iostuff    ! GB


#ifdef WIN32
  USE IFQWIN
#endif

#define GENERATE_MESH 0

  IMPLICIT NONE

! This has been modified by Gib to read the mesh from -mesh and -nodes files

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename

  !CMISS variables

#if GENERATE_MESH
  TYPE(CMISSBasisType) :: Basis    ! GB
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh   ! GB
#else
  TYPE(CMISSBasisType),allocatable :: Bases(:)
  TYPE(CMISSMeshElementsType),allocatable :: Elements(:)    ! GB
#endif
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,EquationsSetField,DependentField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain
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

#if GENERATE_MESH
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 4) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Y_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Z_ELEMENTS<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification 
    NUMBER_GLOBAL_X_ELEMENTS=1
    NUMBER_GLOBAL_Y_ELEMENTS=3
    NUMBER_GLOBAL_Z_ELEMENTS=1
!    INTERPOLATION_TYPE=1
    
    INTERPOLATION_TYPE=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION    
    
  ENDIF
#else
   NUMBER_GLOBAL_X_ELEMENTS=11
   NUMBER_GLOBAL_Y_ELEMENTS=22
   NUMBER_GLOBAL_Z_ELEMENTS = 33    ! GB
#endif

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  CALL CMISSRandomSeedsSet(9999,Err)
  
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Laplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE
  
  CALL CMISSOutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
    
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)

! Restrict to 3D    ! GB
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF

  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"LaplaceRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

#if GENERATE_MESH
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    NUMBER_OF_GAUSS_XI=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)
   
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
#else

! e.g. Filename = 'LV242-mesh'  reads only the mesh, nodes are read later
  write(*,*) 'READ_MESH'
  call READ_MESH('LV242-mesh',MeshUserNumber,Region, Mesh, Bases, Nodes, Elements)
  write(*,*) 'did READ_MESH'
#endif

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
 
  !Destory the mesh now that we have decomposed it
  !CALL CMISSMesh_Destroy(Mesh,Err)
 
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

#if GENERATE_MESH
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
#else

  write(*,*) 'READ_NODES'
  CALL READ_NODES('LV242-nodes',GeometricField)
  write(*,*) 'did READ_NODES'
#endif

  !Create the Standard Laplace Equations set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.5_CMISSDP, &
    & Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMISS_PROBLEM_STANDARD_LAPLACE_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  
!  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
!  CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSDP,Err)
!  CALL CMISSSolver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSDP,Err)

  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  
  !CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_LAPACK_LIBRARY,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_SUPERLU_LIBRARY,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_PASTIX_LIBRARY,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSRegion_NodesGet(Region,Nodes,Err)
  CALL CMISSNodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Export results
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"Laplace","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"Laplace","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM LAPLACEEXAMPLE



