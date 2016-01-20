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
!GB            CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,i,j, &
            CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1,1,i,j, &
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
#if 0

  !>Reads in a field. Only works for a material field currently, all components using same interpolation. handles full field definition code block.
  SUBROUTINE READ_FIELD(Filename,FieldUserNumber,Region,GeometricField, Field)
    character(len=*), intent(in) :: Filename
    INTEGER(CMISSIntg), intent(in) :: FieldUserNumber
    type(cmfe_RegionType), intent(in) :: Region
    type(cmfe_FieldType), intent(in) :: GeometricField
    type(cmfe_FieldType), intent(out) :: Field
    !Local variables
    TYPE(cmfe_DecompositionType) :: Decomposition
    INTEGER(CMISSIntg), parameter :: fid=78
    character*6, parameter :: MAXFMT='(A255)'
    character(len=255) :: word,field_type,interpolation_type,data_type
    INTEGER(CMISSIntg) :: NumberOfVariables,NumberOfComponents
    INTEGER(CMISSIntg) :: varn,mcompn,VariableType,InterpolationType,DataType
    INTEGER(CMISSIntg) :: MyComputationalNode,i,ind,Err
    INTEGER(CMISSIntg) :: data_int(100)
    REAL(CMISSSP) :: data_dp(100)

    CALL cmfe_ComputationalNodeNumberGet(MyComputationalNode,Err)    

    if (MyComputationalNode==0) then
      ! Open the file see if everything is fine
      open(unit=fid,file=Filename,status='old',action='read',err=998)

      ! initial field setup
!       CALL cmfe_Field_Initialise(Field,Err)
!       CALL cmfe_Field_CreateStart(FieldUserNumber,Region,Field,Err)
!       ! --- set material type here then do the below ---
!       CALL cmfe_Decomposition_Initialise(Decomposition,Err)
!       CALL cmfe_Field_MeshDecompositionGet(GeometricField,Decomposition,Err)
!       CALL cmfe_Field_MeshDecompositionSet(Field,Decomposition,Err)
!       CALL cmfe_Field_GeometricFieldSet(Field,GeometricField,Err)

      do
        read(fid,FMT=MAXFMT,end=776) word
        ! skip blanks/comments
        if (len_trim(word)==0 .or. index(adjustl(word),"!")==1) cycle
        ! proper keywords start with a #
        if (trim(adjustl(word))=="#field_type") then
          read(fid,FMT=MAXFMT,end=777) field_type
          select case (trim(field_type))
          case ("material")
            CALL cmfe_Field_Initialise(Field,Err)
            CALL cmfe_Field_CreateStart(FieldUserNumber,Region,Field,Err)
            CALL cmfe_Field_TypeSet(Field,CMISS_FIELD_MATERIAL_TYPE,Err)
            CALL cmfe_Decomposition_Initialise(Decomposition,Err)
            CALL cmfe_Field_MeshDecompositionGet(GeometricField,Decomposition,Err)
            CALL cmfe_Field_MeshDecompositionSet(Field,Decomposition,Err)
            CALL cmfe_Field_GeometricFieldSet(Field,GeometricField,Err)
          case default
            write(*,*) "READ_FIELD: field types other than material has not been implemented"
            close(fid)
            return
          end select
        elseif (trim(adjustl(word))=="#number_of_variables") then
          read(fid,*,end=777) NumberOfVariables
          CALL cmfe_Field_NumberOfVariablesSet(Field,NumberOfVariables,Err)
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
              CALL cmfe_Field_NumberOfComponentsSet(Field,VariableType,NumberOfComponents,Err) 
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
                CALL cmfe_Field_ComponentInterpolationSet(Field,VariableType,i,InterpolationType,Err)
              enddo
            elseif (trim(adjustl(word))=="#mesh_component") then
              read(fid,*,end=777,err=999) mcompn
!               if (InterpolationType==CMISS_FIELD_NODE_BASED_INTERPOLATION) then
                do i=1,NumberOfComponents
                  CALL cmfe_Field_ComponentMeshComponentSet(Field,VariableType,i,mcompn,Err)
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
              CALL cmfe_Field_DataTypeSet(Field,VariableType,DataType,Err)
              CALL cmfe_Field_CreateFinish(Field,Err)
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
                      CALL cmfe_Field_ParameterSetUpdateConstant(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,i,data_dp(i),Err)
                    case (CMISS_FIELD_ELEMENT_BASED_INTERPOLATION)
                      CALL cmfe_Field_ParameterSetUpdateElement(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,ind,i,data_dp(i),Err)
                    case (CMISS_FIELD_NODE_BASED_INTERPOLATION)
                      CALL cmfe_Field_ParameterSetUpdateNode(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,1,1,ind,i,data_dp(i), &
                        & Err)
                    end select
                  enddo
                case (CMISS_FIELD_INTG_TYPE)
                  read(word,*,end=777,err=999) ind,data_int(1:NumberOfComponents)
                  do i=1,NumberOfComponents
                    select case (InterpolationType)
                    case (CMISS_FIELD_CONSTANT_INTERPOLATION)
                      CALL cmfe_Field_ParameterSetUpdateConstant(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,i,data_int(i),Err)
                    case (CMISS_FIELD_ELEMENT_BASED_INTERPOLATION)
                      CALL cmfe_Field_ParameterSetUpdateElement(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,ind,i,data_int(i), &
                        & Err)
                    case (CMISS_FIELD_NODE_BASED_INTERPOLATION)
                      CALL cmfe_Field_ParameterSetUpdateNode(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE,1,1,ind,i,data_int(i), &
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
                  CALL cmfe_Field_ComponentValuesInitialise(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE, &
                    & i,data_dp(i),Err)
                enddo
              case (CMISS_FIELD_INTG_TYPE)
                read(word,*,end=777,err=999) data_int(1:NumberOfComponents)
                do i=1,NumberOfComponents
                  CALL cmfe_Field_ComponentValuesInitialise(Field,VariableType,CMISS_FIELD_VALUES_SET_TYPE, &
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

    CALL cmfe_ComputationalNodeNumberGet(MyComputationalNode,Err)

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

#endif

END MODULE IOSTUFF