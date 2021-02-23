module ExternalModelAlquimiaMod

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use decompMod                    , only : bounds_type
  use EMI_DataMod, only : emi_data_list, emi_data
  use clm_varctl                   , only : iulog

  use ExternalModelBaseType        , only : em_base_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use EMI_CNCarbonStateType_Constants
  use EMI_CNNitrogenStateType_Constants
  use EMI_CNCarbonFluxType_Constants

   use AlquimiaContainers_module, only : AlquimiaSizes,AlquimiaProblemMetaData,AlquimiaProperties,&
            AlquimiaState,AlquimiaAuxiliaryData,AlquimiaAuxiliaryOutputData, AlquimiaEngineStatus, &
            AlquimiaEngineFunctionality,AlquimiaGeochemicalCondition
            
   use alquimia_fortran_interface_mod, only : AlquimiaFortranInterface
   use iso_c_binding, only : c_ptr

  implicit none

  type, public, extends(em_base_type) :: em_alquimia_type
    ! Initialization data needed
    integer :: index_l2e_init_col_type
    integer :: index_l2e_init_lunit_index
    integer :: index_l2e_init_col_active
    integer :: index_l2e_init_col_grid_ind
    integer :: index_l2e_init_lunit_type
    integer :: index_l2e_init_col_wgt
    integer :: index_l2e_init_state_watsatc ! Porosity
    integer :: index_l2e_init_state_temperature_soil
    integer :: index_l2e_init_state_h2osoi_liq
    integer :: index_l2e_init_state_h2osoi_ice
    
    ! Solve data needed
    integer :: index_l2e_col_active
    integer :: index_l2e_lunit_type
    integer :: index_l2e_lunit_index
    integer :: index_l2e_state_h2osoi_liq
    integer :: index_l2e_state_h2osoi_ice
    integer :: index_l2e_state_decomp_cpools
    integer :: index_l2e_state_decomp_npools
    integer :: index_l2e_state_temperature_soil
    integer :: index_l2e_soil_pool_decomp_k
    integer :: index_l2e_state_nh4
    integer :: index_l2e_state_no3
    
    ! Solve data returned to land model
    integer :: index_e2l_state_decomp_cpools
    integer :: index_e2l_state_decomp_npools
    integer :: index_e2l_flux_hr
    integer :: index_e2l_state_nh4
    integer :: index_e2l_state_no3
    
    ! Chemistry engine: Should be one per thread
    type(AlquimiaFortranInterface)       :: chem
    type(AlquimiaEngineStatus)    :: chem_status
    type(c_ptr)                   :: chem_engine
    
    ! Chemistry metadata
    type(AlquimiaSizes)           :: chem_sizes
    type(AlquimiaProblemMetaData) :: chem_metadata
    
    ! Persistent chemical properties and state
    type(AlquimiaProperties), pointer :: chem_properties(:,:)  ! Dimensions of column, depth
    type(AlquimiaState)     , pointer :: chem_state(:,:)       ! Column,depth. Contains a list of species in the structure
    type(AlquimiaAuxiliaryData), pointer   :: chem_aux_data(:,:)
    type(AlquimiaAuxiliaryOutputData), pointer   :: chem_aux_output(:,:)
    
    ! Initial condition
    type(AlquimiaGeochemicalCondition) :: chem_ic
    
    ! State variables
    real(r8), pointer, dimension(:,:)    :: water_density, porosity, temperature, aqueous_pressure
    
    ! Mapping between ELM and alquimia decomp pools
    integer, pointer, dimension(:)       :: carbon_pool_mapping
    integer, pointer, dimension(:)       :: nitrogen_pool_mapping
    integer, pointer, dimension(:)       :: pool_reaction_mapping
    integer                              :: CO2_pool_number
    integer                              :: NH4_pool_number,NO3_pool_number
    
   contains
     procedure, public :: Populate_L2E_Init_List  => EMAlquimia_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EMAlquimia_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EMAlquimia_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EMAlquimia_Populate_E2L_List
     procedure, public :: Init                    => EMAlquimia_Init
     procedure, public :: Solve                   => EMAlquimia_Solve
  end type em_alquimia_type


  real(r8),parameter :: min_dt = 0.1

contains

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during initialization stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index
    
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    ! id                                             = L2E_COLUMN_TYPE
    ! call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    ! this%index_l2e_init_col_type              = index
    
    id                                             = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_lunit_index              = index
    
    id                                             = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active              = index
    
    ! id                                             = L2E_COLUMN_GRIDCELL_INDEX
    ! call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    ! this%index_l2e_init_col_grid_ind              = index
    
    id                                             = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_lunit_type              = index
    
    id                                             = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_watsatc              = index
    
    id                                             = L2E_STATE_TSOIL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_temperature_soil              = index
    
    ! id                                             = L2E_COLUMN_DZ
    ! call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    ! this%index_l2e_init_col_dz              = index
    
    deallocate(em_stages)
    
    write(iulog,*)'L2EInit List:'
    call l2e_init_list%PrintInfo()

  end subroutine EMAlquimia_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from external model to land
    ! model during initialization stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_init_list

    ! write(iulog,*)'EMAlquimia_Populate_E2L_Init_List must be extended by a child class.'
    ! call endrun(msg=errMsg(__FILE__, __LINE__))
    ! write(iulog,*)'EMAlquimia_Populate_E2L_Init_List is empty.'
    write(iulog,*)'E2LInit List:'
    call e2l_init_list%PrintInfo()

  end subroutine EMAlquimia_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_Alquimia_SOLVE_STAGE


    id                                             = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_lunit_index              = index
    
    id                                             = L2E_COLUMN_ACTIVE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_active              = index
    
    id                                             = L2E_LANDUNIT_TYPE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_lunit_type              = index
    
    ! Liquid water
    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    ! Solid  water
    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index
    
    ! Carbon pools
    id                                             = L2E_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_decomp_cpools              = index
    
    ! Nitrogen pools
    id                                             = L2E_STATE_NITROGEN_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_decomp_npools              = index

    ! Soil temperature
    id                                             = L2E_STATE_TSOIL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_temperature_soil              = index
    
    ! Decomposition rate constants
    id                                             = L2E_FLUX_SOIL_POOL_DECOMP_K
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_soil_pool_decomp_k              = index
    
    id                                             = L2E_STATE_NH4_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_nh4              = index
    
    id                                             = L2E_STATE_NO3_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_no3              = index


    deallocate(em_stages)
    
    write(iulog,*)'L2E List:'
    call l2e_list%PrintInfo()

  end subroutine EMAlquimia_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from external model to land
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_list
    
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    ! Updated Carbon pools
    ! May want to change this to rates of change instead?
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_ALQUIMIA_SOLVE_STAGE
    
    id                                             = E2L_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_decomp_cpools              = index
    
    ! Nitrogen pools
    id                                             = E2L_STATE_NITROGEN_POOLS_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_decomp_npools              = index

    ! Heterotrophic respiration flux
    id                                             = E2L_FLUX_HETEROTROPHIC_RESP_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_hr              = index
    
    id                                             = E2L_STATE_NH4_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_nh4              = index
    
    id                                             = E2L_STATE_NO3_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_no3              = index

    deallocate(em_stages)

    write(iulog,*)'E2L List:'
    call e2l_list%PrintInfo()

  end subroutine EMAlquimia_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialize an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    use alquimia_fortran_interface_mod, only : AllocateAlquimiaEngineStatus, &
                                            AllocateAlquimiaProblemMetaData,&
                                            AllocateAlquimiaState,&
                                            AllocateAlquimiaProperties,&
                                            AllocateAlquimiaAuxiliaryData,&
                                            AllocateAlquimiaAuxiliaryOutputData, &
                                            AllocateAlquimiaGeochemicalCondition
                                            
                                                
    use iso_c_binding, only :  C_BOOL, C_CHAR, C_INT, c_f_pointer
    use c_f_interface_module, only : c_f_string_ptr, f_c_string_ptr
    use AlquimiaContainers_module, only : kAlquimiaMaxStringLength

    use clm_varpar, only : nlevdecomp, ndecomp_pools
    use clm_varcon, only : istcrop,istsoil
    use clm_varctl, only : alquimia_inputfile,alquimia_engine_name,alquimia_IC_name,alquimia_CO2_name,&
                           alquimia_NO3_name,alquimia_NH4_name,alquimia_handsoff
    use CNDecompCascadeConType, only : decomp_cascade_con
    
    use PFloTranAlquimiaInterface_module, only : PrintSizes,PrintProblemMetaData, ProcessCondition,PrintState

    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent (in)   :: bounds_clump
    
    
    ! Local variables
    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: lun_type(:)
    integer                              :: c,j,l,ii,jj
    real(r8) , pointer                   :: porosity(:,:),temperature(:,:)
    
    character(len=kAlquimiaMaxStringLength) :: alq_poolname
    type (c_ptr), pointer :: name_list(:)
    logical :: found_pool
    integer :: pool_num
    
    
    ! Should read this from a namelist
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: inputfile
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: engine_name
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: IC_name
    
    
    logical(C_BOOL) :: hands_off
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    type(AlquimiaEngineFunctionality) :: chem_engine_functionality

    
    print *, 'Entering Alquimia setup'
    
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active          , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_lunit_index  , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_lunit_type       , lun_type     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_lunit_type       , lun_type     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_state_watsatc       , porosity     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_state_temperature_soil , temperature     )
    
    

    ! These are parameters of Alquimia engine. Need to read them in somehow
    inputfile   = alquimia_inputfile
    engine_name = alquimia_engine_name
    IC_name     = alquimia_IC_name  ! Name of initial condition
    hands_off = alquimia_handsoff  ! hands_off = .false. allows/requires rate constants, mineral rate const, CEC, complexation site density, and isotherms to be passed through alquimia 
    
    ! Allocate memory for status container
    call AllocateAlquimiaEngineStatus(this%chem_status)
    ! Point Alquimia interface to correct subroutines (based on engine that was specified in engine_name)
    call this%chem%CreateInterface(engine_name, this%chem_status)
    
    ! Print out the result of the interface creation call
    call c_f_string_ptr(this%chem_status%message,status_message)
    if(this%chem_status%error /= 0) then
      call endrun(msg='Alquimia error: '//status_message)
    endif
    
    ! Set up the engine and get the storage requirements
    call this%chem%Setup(inputfile, hands_off, this%chem_engine, this%chem_sizes, chem_engine_functionality, this%chem_status)
    ! Print out the result of the interface creation call
    call c_f_string_ptr(this%chem_status%message,status_message)
    if(this%chem_status%error /= 0) then
      call endrun(msg='Alquimia error: '//status_message)
    endif
    
    write(iulog,'(a,L1)') 'Alquimia hands off mode should be ',hands_off
    call PrintSizes(this%chem_sizes)
    
    ! Allocate memory for chemistry data
    call AllocateAlquimiaProblemMetaData(this%chem_sizes, this%chem_metadata)
    
    call this%chem%GetProblemMetaData(this%chem_engine, this%chem_metadata, this%chem_status)
    if(this%chem_status%error /= 0) then
      call c_f_string_ptr(this%chem_status%message,status_message)
      call endrun(msg='Alquimia error: '//status_message)
    endif
    
    call printproblemmetadata(this%chem_metadata)
    
    ! Map out the location of pertinent pools in Alquimia data structure
    ! Assumes that organic matter pools in PFLOTRAN are named the same as decomp_pool_name_history
    ! Currently we are not mapping any non-CTC pools. Need to add another data structure for other chemicals if we want to save them in ELM restart/history or use in other ways
    print *,'Alquimia carbon pool mapping:'
    allocate(this%carbon_pool_mapping(ndecomp_pools))
    call c_f_pointer(this%chem_metadata%primary_names%data, name_list, (/this%chem_sizes%num_primary/))
    do ii=1, ndecomp_pools
      if(decomp_cascade_con%floating_cn_ratio_decomp_pools(ii)) then
        alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))//'C'
      else
        alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))
      endif
      pool_num = find_alquimia_pool(alq_poolname,name_list,this%chem_sizes%num_primary)
      if(pool_num>0) then 
        print '(a, i3, 1X,a7, a, i3, 1X, a)','ELM pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii)),' <-> Alquimia pool',pool_num,trim(alq_poolname)
      else
        print *,'WARNING: No match for pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii))
      endif
      this%carbon_pool_mapping(ii)=pool_num
    enddo
    
    pool_num = find_alquimia_pool(alquimia_CO2_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      print '(a,6x,a,i3,1x,a)','CO2 production', '<-> Alquimia pool',pool_num,trim(alquimia_CO2_name)
    else
      print '(a,i3,1X,a)','WARNING: No match for pool',ii,trim(alquimia_CO2_name)
    endif
    this%CO2_pool_number = pool_num
    
    
    print *,'Alquimia nitrogen pool mapping:'
    allocate(this%nitrogen_pool_mapping(ndecomp_pools))
    do ii=1, ndecomp_pools
      alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))//'N'
      pool_num = find_alquimia_pool(alq_poolname,name_list,this%chem_sizes%num_primary)
      if(pool_num>0) then 
        print '(a, i3, 1X,a7, a, i3, 1X, a)','ELM pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii)),' <-> Alquimia pool',pool_num,trim(alq_poolname)
      elseif  (decomp_cascade_con%floating_cn_ratio_decomp_pools(ii)) then
        print '(a,i3,1X,a)','WARNING: No match for pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii))
      endif
      this%nitrogen_pool_mapping(ii)=pool_num
    enddo
    
    pool_num = find_alquimia_pool(alquimia_NH4_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      print '(a,6x,a,i3,1x,a)','NH4', '<-> Alquimia pool',pool_num,trim(alquimia_NH4_name)
    else
      print '(a,i3,1X,a)','WARNING: No match for pool',ii,trim(alquimia_NH4_name)
    endif
    this%NH4_pool_number = pool_num
    
    pool_num = find_alquimia_pool(alquimia_NO3_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      print '(a,6x,a,i3,1x,a)','NO3', '<-> Alquimia pool',pool_num,trim(alquimia_NO3_name)
    else
      print '(a,i3,1X,a)','WARNING: No match for pool',ii,trim(alquimia_NO3_name)
    endif
    this%NO3_pool_number = pool_num
    ! print *,this%carbon_pool_mapping
    ! print *,this%nitrogen_pool_mapping
    

    ! Need to map out reactions as well
    allocate(this%pool_reaction_mapping(ndecomp_pools))
    call c_f_pointer(this%chem_metadata%aqueous_kinetic_names%data, name_list, (/this%chem_metadata%aqueous_kinetic_names%size/))
    print *,'Alquimia reactions:'
    do ii=1,this%chem_metadata%aqueous_kinetic_names%size
      call c_f_string_ptr(name_list(ii),alq_poolname)
      print *,trim(alq_poolname)
    enddo
    do ii=1, ndecomp_pools
      if(decomp_cascade_con%cascade_receiver_pool(ii)>0) then
        alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))//' decay to '// trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(ii)))//' (SOMDEC sandbox)'
      else
        alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))//' decay to CO2 (SOMDEC sandbox)'
      endif
      pool_num = find_alquimia_pool(alq_poolname,name_list,this%chem_metadata%aqueous_kinetic_names%size)
      if(pool_num>0) then 
        print '(a, i3, 1X,a7,a, i3, 1X, a)','ELM pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii)),' decay <-> Alquimia reaction',pool_num,trim(alq_poolname)
      else
        print '(a,i3,1x,a,1x,a)','WARNING: No match for reaction:',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii)),': ',trim(alq_poolname)
      endif
      this%pool_reaction_mapping(ii)=pool_num
    enddo
    ! Todo: Keep track of other (non-SOMDEC) reactions too somehow
    
    ! Initial condition. The zero length for constraints suggest that it will be read in from input file?
    call AllocateAlquimiaGeochemicalCondition(len_trim(ic_name,C_INT),0,0,this%chem_ic)
    call f_c_string_ptr(ic_name,this%chem_ic%name,len_trim(ic_name)+1)
    
    ! These are duplicated for each cell
    allocate(this%chem_state(nlevdecomp,bounds_clump%begc:bounds_clump%endc))
    allocate(this%chem_properties(nlevdecomp,bounds_clump%begc:bounds_clump%endc))
    allocate(this%chem_aux_data(nlevdecomp,bounds_clump%begc:bounds_clump%endc))
    allocate(this%chem_aux_output(nlevdecomp,bounds_clump%begc:bounds_clump%endc))
    
    allocate(this%porosity(nlevdecomp,bounds_clump%begc:bounds_clump%endc))
    
    do c = bounds_clump%begc, bounds_clump%endc
      l = col_landunit(c)
      if ((col_active(c) == 1).and. &
           (lun_type(l) == istsoil  .or. lun_type(l) == istcrop)) then
        do j = 1, nlevdecomp    ! Should be changed when we start doing a whole column at once
            call AllocateAlquimiaState(this%chem_sizes, this%chem_state(j,c))
            call AllocateAlquimiaProperties(this%chem_sizes, this%chem_properties(j,c))
            call AllocateAlquimiaAuxiliaryData(this%chem_sizes, this%chem_aux_data(j,c))
            call AllocateAlquimiaAuxiliaryOutputData(this%chem_sizes, this%chem_aux_output(j,c))
            
            ! Initialize the state for the cell
            this%chem_properties(j,c)%volume = 1.0 ! l2e_volume(j,c)
            this%chem_properties(j,c)%saturation = 1.0 ! l2e_saturation(j,c)
            this%chem_state(j,c)%water_density = 999.9720
            this%porosity(j,c) = porosity(j,c)
            this%chem_state(j,c)%porosity = this%porosity(j,c) ! l2e_porosity(j,c)
            this%chem_state(j,c)%aqueous_pressure = 101325.0
            this%chem_state(j,c)%temperature = temperature(j,c) - 273.15

            call this%chem%ProcessCondition(this%chem_engine, this%chem_ic, this%chem_properties(j,c), this%chem_state(j,c), &
                                           this%chem_aux_data(j,c), this%chem_status)
            if(this%chem_status%error /= 0) then
              call c_f_string_ptr(this%chem_status%message,status_message)
              call endrun(msg='Alquimia error in ProcessCondition: '//status_message)
            endif
            
            
        enddo
      endif
    enddo
    ! call PrintState(this%chem_state(1,1))
    
    ! At this point, Alquimia transport driver also allocates boundary conditions
    ! but not sure how we handle that in this case
    
    ! Here, any chemistry info (CEC, ion exchange, etc to chem_state and isotherm kd, langmuir, etc would be copied to chem_properties)

  end subroutine EMAlquimia_Init
  
  integer function find_alquimia_pool(pool_name,name_list,n_names) result(pool_number)
    use AlquimiaContainers_module, only : kAlquimiaMaxStringLength
    use c_f_interface_module, only : c_f_string_ptr

    implicit none
    
    character(*),intent(in) :: pool_name
    type (c_ptr), pointer,intent(in) :: name_list(:)
    integer, intent(in) :: n_names
    
    integer :: jj
    character(len=kAlquimiaMaxStringLength) :: alq_poolname
    
    
    pool_number=-1
    
    do jj=1, n_names
      call c_f_string_ptr(name_list(jj),alq_poolname)
      if(trim(alq_poolname) == trim(pool_name)) then
        pool_number=jj
        exit
      endif
    enddo
    
  end function find_alquimia_pool

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    use clm_varpar, only : nlevdecomp,ndecomp_pools
    use clm_varcon, only : istcrop,istsoil,catomw,natomw
    use AlquimiaContainers_module, only : AlquimiaEngineStatus, kAlquimiaMaxStringLength
    use alquimia_fortran_interface_mod, only :  ReactionStepOperatorSplit, GetAuxiliaryOutput
    use, intrinsic :: iso_c_binding, only : C_CHAR, c_double, c_f_pointer
    use c_f_interface_module, only : c_f_string_ptr
    use PFloTranAlquimiaInterface_module, only : printState
    
    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none
    !
    ! !ARGUMENTS
    class(em_alquimia_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt ! s
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump
    
    
    ! Local variables
    integer                              :: c,j,l,poolnum
    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: lun_type(:)
    integer                              :: max_cuts
    real(r8) , pointer, dimension(:,:,:)    :: soilcarbon_l2e,soilcarbon_e2l 
    real(r8) , pointer, dimension(:,:,:)    :: soilnitrogen_l2e,soilnitrogen_e2l 
    real(r8) , pointer, dimension(:,:,:)    :: decomp_k
    real(r8) , pointer, dimension(:,:)    :: hr_e2l , temperature, h2o_liq, h2o_ice
    real(r8) , pointer, dimension(:,:)  :: no3_e2l,no3_l2e,nh4_e2l,nh4_l2e
    real(r8) :: CO2_before
    
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    procedure(ReactionStepOperatorSplit), pointer :: engine_ReactionStepOperatorSplit
    procedure(GetAuxiliaryOutput), pointer   :: engine_getAuxiliaryOutput
    real (c_double), pointer :: alquimia_data(:)

    ! write(iulog,*) 'Alquimia solving step!'
    
    ! Pass data from ELM
    
    ! Column filters
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_active          , col_active   )
    call l2e_list%GetPointerToInt1D(this%index_l2e_lunit_index  , col_landunit )
    call l2e_list%GetPointerToInt1D(this%index_l2e_lunit_type       , lun_type     )
    
    ! C and N pools
    call l2e_list%GetPointerToReal3D(this%index_l2e_state_decomp_cpools , soilcarbon_l2e)
    call l2e_list%GetPointerToReal3D(this%index_l2e_state_decomp_npools , soilnitrogen_l2e)
    
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_no3 , no3_l2e)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_nh4 , nh4_l2e)

    ! Abiotic factors
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_temperature_soil , temperature  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq, h2o_liq)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice, h2o_ice)
    
    ! Pool turnover rate constants calculated in ELM, incorporating T and moisture effects (1/s)
    call l2e_list%GetPointerToReal3D(this%index_l2e_soil_pool_decomp_k, decomp_k)

    ! C and N pools
    call e2l_list%GetPointerToReal3D(this%index_e2l_state_decomp_cpools , soilcarbon_e2l)
    call e2l_list%GetPointerToReal3D(this%index_e2l_state_decomp_npools , soilnitrogen_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_hr , hr_e2l)
    
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_no3 , no3_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_nh4 , nh4_e2l)

    
     ! Run the reactions engine for a step. Alquimia works on one cell at a time
     ! TODO: Transport needs to be integrated somehow. 
     do c = bounds_clump%begc, bounds_clump%endc
       l = col_landunit(c)
       if ((col_active(c) == 1).and. &
            (lun_type(l) == istsoil  .or. lun_type(l) == istcrop)) then
         do j = 1, nlevdecomp  
             ! Run the reactions one time step
             this%chem_state(j,c)%water_density = 999.9720
             this%chem_state(j,c)%porosity = this%porosity(j,c) ! l2e_porosity(j,c)
             this%chem_state(j,c)%aqueous_pressure = 101325.0
             this%chem_state(j,c)%temperature = temperature(j,c) - 273.15
             ! Set volume and saturation?
             ! Saturation determines volume of water and concentration of aqueous species
             ! Not sure if volume actually matters
             
             ! Set rate constants in properties.aqueous_kinetic_rate_cnst if in hands_on mode
             ! Should we check if we are in hands_on mode?
             call c_f_pointer(this%chem_properties(j,c)%aqueous_kinetic_rate_cnst%data, alquimia_data, (/this%chem_properties(j,c)%aqueous_kinetic_rate_cnst%size/))
             
             do poolnum=1,ndecomp_pools
               if(this%pool_reaction_mapping(poolnum)>0) then
                 ! write(iulog,'(a,i3,i3,a,e12.4)') 'Setting pool',poolnum,this%pool_reaction_mapping(poolnum),'kinetic rate const to',decomp_k(c,j,poolnum)
                 alquimia_data(this%pool_reaction_mapping(poolnum)) = decomp_k(c,j,poolnum)
               endif
             enddo
             ! write(iulog,'(a,8e12.4)') 'Rate constants: ',decomp_k(c,j,1:8) ;
             
             ! Set soil carbon and nitrogen from land model
             ! Have not done unit coversions yet!!!
             call c_f_pointer(this%chem_state(j,c)%total_immobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
             do poolnum=1,ndecomp_pools
               if(this%carbon_pool_mapping(poolnum)>0) &  ! Convert soil C from gC/m3 to mol C/m3. Assumes pool is defined as immobile, not aqueous
                 alquimia_data(this%carbon_pool_mapping(poolnum)) = soilcarbon_l2e(c,j,poolnum)/catomw
               ! Separate N pool only exists if floating CN ratio
               if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) &
                  alquimia_data(this%nitrogen_pool_mapping(poolnum)) = soilnitrogen_l2e(c,j,poolnum)/natomw
             enddo
             
             CO2_before  = alquimia_data(this%CO2_pool_number)
            
            ! Copy dissolved nitrogen species. Units need to be converted from gN/m3 to M/L. Currently assuming saturated porosity
            call c_f_pointer(this%chem_state(j,c)%total_mobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
             if(this%NO3_pool_number>0) alquimia_data(this%NO3_pool_number) = no3_l2e(c,j)/natomw*1000.0/this%chem_state(j,c)%porosity
             if(this%NH4_pool_number>0) alquimia_data(this%NH4_pool_number) = nh4_l2e(c,j)/natomw*1000.0/this%chem_state(j,c)%porosity

              call run_onestep(this, j,c, dt,0,max_cuts)
              ! print *,"Converged after",max_cuts,"cuts"
                                 
                                 
              ! Set updated land model values
              ! Have not done unit coversions yet!!!
              call c_f_pointer(this%chem_state(j,c)%total_immobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
              do poolnum=1,ndecomp_pools
                if(this%carbon_pool_mapping(poolnum)>0) &
                  soilcarbon_e2l(c,j,poolnum) = alquimia_data(this%carbon_pool_mapping(poolnum))*catomw
                ! Separate N pool only exists if floating CN ratio
                if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) then
                   soilnitrogen_e2l(c,j,poolnum) = alquimia_data(this%nitrogen_pool_mapping(poolnum))*natomw
                 elseif (this%carbon_pool_mapping(poolnum)>0) then
                   ! Calculate from CN ratio and C pool
                   soilnitrogen_e2l(c,j,poolnum) = soilcarbon_e2l(c,j,poolnum)/decomp_cascade_con%initial_cn_ratio(poolnum)
                endif
              enddo
              ! This needs to account for whether CO2 is stored in a mobile or immobile PFLOTRAN pool because units are different
              if(this%CO2_pool_number>0) hr_e2l(c,j) = (alquimia_data(this%CO2_pool_number)-CO2_before)*catomw/dt
              
              call c_f_pointer(this%chem_state(j,c)%total_mobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
              if(this%NO3_pool_number>0) no3_e2l(c,j) = alquimia_data(this%NO3_pool_number)*natomw/(1000.0/this%chem_state(j,c)%porosity)
              if(this%NH4_pool_number>0) nh4_e2l(c,j) = alquimia_data(this%NH4_pool_number)*natomw/(1000.0/this%chem_state(j,c)%porosity)
              
         enddo
       endif
     enddo
     

     ! Alquimia here calls GetAuxiliaryOutput which copies data back to interface arrays. We should do that here for EMI arrays
     ! Again, need to convert units back to ELM style, keeping track of what kind of species we are using so units are correct

  end subroutine EMAlquimia_Solve
  

  
  recursive subroutine run_onestep(this,j,c,dt,num_cuts,max_cuts)
    
    use, intrinsic :: iso_c_binding, only : C_CHAR
    use c_f_interface_module, only : c_f_string_ptr
    use AlquimiaContainers_module, only : kAlquimiaMaxStringLength
    
    implicit none
    
    class(em_alquimia_type)              :: this
    integer,intent(out)                  :: max_cuts
    integer,intent(in)                   :: num_cuts,j,c
    real(r8),intent(in)                  :: dt
    
    real(r8) :: actual_dt,porosity
    character(512) :: msg
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    integer :: ncuts2,ncuts,ii
    
    max_cuts = num_cuts
    actual_dt = dt/(2**num_cuts)
    
    ncuts=0
    ncuts2=0
    
    porosity=this%chem_state(j,c)%porosity
    call this%chem%ReactionStepOperatorSplit(this%chem_engine, actual_dt, this%chem_properties(j,c), this%chem_state(j,c), &
                                           this%chem_aux_data(j,c), this%chem_status)
    ! Reset porosity because Pflotran tends to mess it up
    this%chem_state(j,c)%porosity=porosity
    ! print *,'Converged =',this%chem_status%converged,"ncuts =",num_cuts
    if (this%chem_status%converged) then
      ! Success. Can get aux output and finish execution of the subroutine
      ! Get auxiliary output
      call this%chem%getAuxiliaryOutput(this%chem_engine, this%chem_properties(j,c), this%chem_state(j,c), &
                                  this%chem_aux_data(j,c), this%chem_aux_output(j,c), this%chem_status)
      if(this%chem_status%error /= 0) then
        call c_f_string_ptr(this%chem_status%message,status_message)
        call endrun(msg='Alquimia error in ReactionStepOperatorSplit: '//status_message)
      endif
      
    else ! Solve did not converge. Cut timestep, and bail out if too short
      if(actual_dt/2 < min_dt) then
        call c_f_string_ptr(this%chem_status%message,status_message)
        write(msg,*) "Error: Alquimia ReactionStepOperatorSplit failed to converge after",num_cuts,"cuts to dt =",actual_dt,'s'
        call endrun(msg=msg)
      else
        ! If we are not at minimum timestep yet, cut and keep going
        ! Need to run the step two times because we have cut the timestep in half
        call run_onestep(this, j,c, dt,num_cuts+1,ncuts)
        if(ncuts>max_cuts) max_cuts=ncuts
        ! print *,'Converged =',this%chem_status%converged,"ncuts =",ncuts,'(Substep 1)'
        
        ! The second one starts from the maximum number of cuts from the first one so it doesn't waste time retrying a bunch of failed timestep lengths
         do ii=1,2**(max_cuts-(num_cuts+1))
           call run_onestep(this, j,c, dt,ncuts,ncuts2)
           if(ncuts2>max_cuts) max_cuts=ncuts2
        !   print *,'Converged =',this%chem_status%converged,"ncuts =",ncuts2,'. Substep 2 +',ii
         enddo
        ! call run_onestep(this, j,c, dt,num_cuts+1,ncuts)
        ! if(ncuts>max_cuts) max_cuts=ncuts
        ! print *,'Converged =',this%chem_status%converged,"ncuts =",ncuts,'(Substep 2)'
      endif
    endif
      
      

  end subroutine run_onestep
  

end module ExternalModelAlquimiaMod
