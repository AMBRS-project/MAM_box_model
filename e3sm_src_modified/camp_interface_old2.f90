!> \file
!> The mam4_camp_interface module (MAM4 ↔ CAMP).
!>
!> Minimal, runtime-configurable interface:
!>  - CAMP config path comes from env var MAM4_CAMP_CONFIG (fallback: "mam4_config.json")
!>  - No external IC/emission/photolysis files are required for this test path
!>  - GMD/GSD are pushed into CAMP each step (modal-binned-mass update)
!>  - Gas and aerosol species are mapped in/out each step
!>  - Persistent state carries un-mapped species across timesteps

module mam4_camp_interface
  use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef MAM4_USE_CAMP
  use mam4_state
  use camp_camp_core
  use camp_camp_state
  use camp_chem_spec_data
  use camp_aero_rep_data
  use camp_aero_rep_modal_binned_mass
  use camp_constants
  use camp_util                      ! for assert/assert_msg, string_t
  use json_module
  use camp_solver_stats
  use camp_mechanism_data
  use camp_rxn_data
  use camp_rxn_emission
  use camp_rxn_photolysis
#endif
  implicit none

#ifdef MAM4_USE_CAMP
  logical :: first_step = .true.
  real(kind=r8), allocatable :: persistent_state(:)
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MAM4_USE_CAMP
  !> Run the CAMP module for the current MAM4 state (one time step).
  subroutine mam4_camp_interface_solve(env_state, aero_state, gas_state, del_t)
    ! Inputs
    type(env_state_t),  intent(inout) :: env_state
    type(aero_state_t), intent(inout) :: aero_state
    type(gas_state_t),  intent(inout) :: gas_state
    real(kind=r8),      intent(in)    :: del_t

    ! CAMP core objects
    type(camp_core_t),  pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    type(mechanism_data_t), pointer :: mechanism
    type(solver_stats_t), target :: solver_stats

    ! Aero-rep and modal geometry
    character(len=4),  parameter :: aero_rep_key = "MAM4"
    character(len=16), parameter :: mode_names(4) = (/ "accumulation    ", &
                                                       "aitken          ", &
                                                       "coarse          ", &
                                                       "primary_carbon  " /)
    integer(kind=i_kind) :: mode(4)
    integer :: n

    class(aero_rep_data_t), pointer :: aero_rep_ptr
    class(aero_rep_modal_binned_mass_t), pointer :: aero_rep_data
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

    ! mechanism list + config path
    type(string_t), allocatable :: mechs(:)
    character(len=512) :: cfg
    integer :: nlen, istat

    ! --- Initialize CAMP core from env var (fallback local file) ---
    cfg = 'mam4_config.json'
    call get_environment_variable('MAM4_CAMP_CONFIG', cfg, length=nlen, status=istat)
    if (istat == 0 .and. nlen > 0) cfg = cfg(1:nlen)
    camp_core => camp_core_t(trim(cfg))
    call camp_core%initialize()

    ! --- Pick the first mechanism defined in the config ---
    mechs = camp_core%mechanism_names()
    call assert(260845179, camp_core%get_mechanism(trim(mechs(1)%string), mechanism))

    ! --- Get aerosol representation and push modal geometry (GMD/GSD) ---
    call assert(209301925, camp_core%get_aero_rep(aero_rep_key, aero_rep_ptr))
    select type (aero_rep_ptr)
    type is (aero_rep_modal_binned_mass_t)
      call camp_core%initialize_update_object(aero_rep_ptr, update_data_GMD)
      call camp_core%initialize_update_object(aero_rep_ptr, update_data_GSD)
      call camp_core%solver_initialize()
      camp_state => camp_core%new_state()

      ! Carry previous step's full state (species not mapped from MAM are preserved)
      if (first_step) then
        first_step = .false.
        allocate(persistent_state(size(camp_state%state_var)))
        camp_state%state_var = 0.0_r8
        persistent_state     = 0.0_r8
      else
        camp_state%state_var = persistent_state
      end if

      ! Update modal geometry each step from MAM state
      do n = 1, 4
        call assert_msg(431201141, aero_rep_ptr%get_section_id(mode_names(n), mode(n)), "Could not get mode ID")
        call update_data_GMD%set_GMD(mode(n), aero_state%GMD(n))
        call update_data_GSD%set_GSD(mode(n), aero_state%GSD(n))
        call camp_core%update_data(update_data_GMD)
        call camp_core%update_data(update_data_GSD)
      end do

    class default
      write(*,*) "mam4_camp_interface: unsupported aerosol representation"
      stop 3
    end select

    ! --- Map environment + gas + aerosol from MAM into CAMP ---
    call env_state_set_camp_env_state(env_state, camp_state)
    call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)
    call mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)

    ! --- Solve chemistry/partitioning ---
    call camp_core%solve(camp_state, del_t, solver_stats=solver_stats)
    if (solver_stats%solver_flag .ne. 0) then
      write(*,*) 'CAMP solver failed with code ', solver_stats%solver_flag
      stop
    end if

    ! --- Map results back into MAM state ---
    call mam4_camp_interface_get_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)
    call gas_state_get_camp_conc(gas_state, camp_state, camp_core)

    ! --- Persist full state for the next step ---
    persistent_state = camp_state%state_var

    ! --- Cleanup ---
    deallocate(camp_core, camp_state)
  end subroutine mam4_camp_interface_solve
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Push aerosol-phase species from MAM → CAMP state
  subroutine mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_data)
#ifdef MAM4_USE_CAMP
    type(aero_state_t),               intent(inout) :: aero_state
    type(camp_core_t),    pointer,    intent(inout) :: camp_core
    type(camp_state_t),   pointer,    intent(inout) :: camp_state
    class(aero_rep_data_t), pointer,  intent(inout) :: aero_rep_data

    type(string_t), allocatable :: names(:)
    integer :: i, id, mode_id
    character(len=16) :: spec, mode_name
    logical :: mode_flag

    allocate(names(size(aero_rep_data%unique_names())))
    names = aero_rep_data%unique_names()

    do i = 1, size(names)
      select type (aero_rep_data)
      type is (aero_rep_modal_binned_mass_t)
        id        = aero_rep_data%spec_state_id(trim(names(i)%string))
        spec      = aero_rep_data%spec_name(trim(names(i)%string))
        mode_name = mode_extract(trim(names(i)%string), '.')
        mode_flag = aero_rep_data%get_section_id(trim(mode_name), mode_id)
      end select
      if (.not. mode_flag) then
        write(*,*) 'mam4_camp_interface_set_camp_aerosol: mode not found'
        stop
      end if

      select case (trim(spec))
      case('ASO4')
        camp_state%state_var(id) = aero_state%qso4(mode_id)
      case('APOC')
        camp_state%state_var(id) = aero_state%qpom(mode_id)
      case('SOA')
        camp_state%state_var(id) = aero_state%qsoa(mode_id)
      case('AEC')
        camp_state%state_var(id) = aero_state%qbc(mode_id)
      case('ASOIL')
        camp_state%state_var(id) = aero_state%qdst(mode_id)
      case('ANA')
        camp_state%state_var(id) = aero_state%qncl(mode_id)
      !case('AH2O')
      !  camp_state%state_var(id) = aero_state%qaerwat(mode_id)
      end select
    end do

    deallocate(names)
#endif
  end subroutine mam4_camp_interface_set_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pull aerosol-phase species from CAMP → MAM state
  subroutine mam4_camp_interface_get_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_data)
#ifdef MAM4_USE_CAMP
    type(aero_state_t),               intent(inout) :: aero_state
    type(camp_core_t),    pointer,    intent(inout) :: camp_core
    type(camp_state_t),   pointer,    intent(inout) :: camp_state
    class(aero_rep_data_t), pointer,  intent(inout) :: aero_rep_data

    type(string_t), allocatable :: names(:)
    integer :: i, id, mode_id
    character(len=16) :: spec, mode_name
    logical :: mode_flag

    allocate(names(size(aero_rep_data%unique_names())))
    names = aero_rep_data%unique_names()

    do i = 1, size(names)
      select type (aero_rep_data)
      type is (aero_rep_modal_binned_mass_t)
        id        = aero_rep_data%spec_state_id(trim(names(i)%string))
        spec      = aero_rep_data%spec_name(trim(names(i)%string))
        mode_name = mode_extract(trim(names(i)%string), '.')
        mode_flag = aero_rep_data%get_section_id(trim(mode_name), mode_id)
      end select
      if (.not. mode_flag) then
        write(*,*) 'mam4_camp_interface_get_camp_aerosol: mode not found'
        stop
      end if

      select case (trim(spec))
      case('ASO4')
        aero_state%qso4(mode_id) = camp_state%state_var(id)
      case('APOC')
        aero_state%qpom(mode_id) = camp_state%state_var(id)
      case('SOA')
        aero_state%qsoa(mode_id) = camp_state%state_var(id)
      case('AEC')
        aero_state%qbc(mode_id)  = camp_state%state_var(id)
      case('ASOIL')
        aero_state%qdst(mode_id) = camp_state%state_var(id)
      case('ANA')
        aero_state%qncl(mode_id) = camp_state%state_var(id)
      !case('AH2O')
      !  aero_state%qaerwat(mode_id) = camp_state%state_var(id)
      end select
    end do

    deallocate(names)
#endif
  end subroutine mam4_camp_interface_get_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Utility: extract mode token before separator (e.g., "accumulation" from "accumulation.mixed.ASO4")
  character(len=16) function mode_extract(str, sep)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: sep
    integer :: i, char_count

    char_count = 0
    do i = 1, len(str)
      if (str(i:i) .eq. sep) exit
      char_count = char_count + 1
    end do
    mode_extract = str(1:char_count)
  end function mode_extract

end module mam4_camp_interface
