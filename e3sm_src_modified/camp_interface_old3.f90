!> Interface between MAM4 box model and CAMP
module mam4_camp_interface

  use shr_kind_mod,                 only: r8 => shr_kind_r8
#ifdef MAM4_USE_CAMP
  use mam4_state
  use camp_camp_core
  use camp_camp_state
  use camp_chem_spec_data
  use camp_aero_rep_data
  use camp_aero_rep_modal_binned_mass
  use camp_constants
  use camp_util
  use json_module
  use camp_solver_stats
  use camp_mechanism_data
  use camp_rxn_data
  use camp_rxn_emission
  use camp_rxn_photolysis
#endif

  implicit none

#ifdef MAM4_USE_CAMP
  logical, save :: first_step = .true.
  real(kind=r8), allocatable, save :: persistent_state(:)
#endif

contains

#ifdef MAM4_USE_CAMP
  !---------------------------------------------------------------------------
  !> Run CAMP for the current MAM4 state over a single time step
  subroutine mam4_camp_interface_solve(env_state, aero_state, gas_state, del_t)

    ! CAMP objects
    type(camp_core_t),  pointer :: camp_core
    type(camp_state_t), pointer :: camp_state

    ! MAM4 states
    type(env_state_t),  intent(inout) :: env_state
    type(aero_state_t), intent(inout) :: aero_state
    type(gas_state_t),  intent(inout) :: gas_state

    ! Timestep
    real(kind=r8), intent(in) :: del_t

    ! Locals
    integer :: n, i, i_ic
    character(len=4),  parameter :: aero_rep_key = "MAM4"
    character(len=16), parameter :: mode_names(4) = (/ "accumulation    ", &
                                                       "aitken          ", &
                                                       "coarse          ", &
                                                       "primary_carbon  " /)
    real(kind=r8),     parameter :: pi = 3.14159265358979323846_r8
    character(len=255), allocatable :: ic_spec(:)
    character(len=255), allocatable :: aero_ic_names(:)
    character(len=255) :: cwd
    character(len=255) :: aero_mode, aero_spec, mode_name, aero_name
    real(kind=r8) :: aero_mass_frac, aero_dens, aero_vol(4), tmpdens(4), rtmpdens(4)
    type(string_t), allocatable :: names(:), aero_names(:), mechs(:)
    integer(kind=i_kind) :: mode(4)
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    class(aero_rep_modal_binned_mass_t), pointer :: aero_rep_data
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

    type(mechanism_data_t), pointer :: mechanism
    type(rxn_update_data_emission_t), allocatable :: q_update(:)
    type(rxn_update_data_photolysis_t), allocatable :: j_update(:)
    real(kind=r8), allocatable :: j(:), q(:), ic(:)
    integer, allocatable :: i_j(:), i_q(:), aero_ids(:)
    integer :: n_emis, n_phot, i_emis, i_phot, aero_id
    class(rxn_data_t), pointer :: rxn
    type(chem_spec_data_t), pointer :: chem_spec_data
    logical :: flag

    type(solver_stats_t), target :: solver_stats

    !--- Config path (env override MAM4_CAMP_CONFIG)
    character(len=512) :: cfg
    integer :: nlen, istat
    cfg = 'mam4_config.json'
    call get_environment_variable('MAM4_CAMP_CONFIG', cfg, length=nlen, status=istat)
    if (istat == 0 .and. nlen > 0) cfg = cfg(1:nlen)
    camp_core => camp_core_t(trim(cfg))
    call camp_core%initialize()

    !--- Count emission and photolysis reactions in the (first) mechanism
    mechs = camp_core%mechanism_names()
    call assert(260845179, camp_core%get_mechanism(trim(mechs(1)%string), mechanism))

    n_emis = 0; n_phot = 0
    do i = 1, mechanism%size()
      rxn => mechanism%get_rxn(i)
      select type (rxn)
      class is (rxn_photolysis_t); n_phot = n_phot + 1
      class is (rxn_emission_t);   n_emis = n_emis + 1
      end select
    end do

    if (n_emis > 0) allocate(q_update(n_emis), i_q(n_emis))
    if (n_phot > 0) allocate(j_update(n_phot), i_j(n_phot))

    !--- Build update objects for emis & photolysis
    i_emis = 0; i_phot = 0
    do i = 1, mechanism%size()
      rxn => mechanism%get_rxn(i)
      select type (rxn)
      class is (rxn_photolysis_t)
        i_phot = i_phot + 1
        i_j(i_phot) = i
        call camp_core%initialize_update_object(rxn, j_update(i_phot))
      class is (rxn_emission_t)
        i_emis = i_emis + 1
        i_q(i_emis) = i
        call camp_core%initialize_update_object(rxn, q_update(i_emis))
      end select
    end do

    !--- Read emission rates (flat test file, like Duncan’s driver)
    if (n_emis > 0) then
      allocate(q(n_emis))
      open(1, file='emis_q.dat', status='old')
      do i_emis = 1, n_emis
        read(1, *) i, q(i_emis)
      end do
      close(1)
      do i_emis = 1, n_emis
        rxn => mechanism%get_rxn(i_q(i_emis))
        select type (rxn)
        class is (rxn_emission_t)
          call q_update(i_emis)%set_rate(q(i_emis))
        end select
      end do
    end if

    !--- Optional photolysis (disabled unless you fill phot_j.dat)
    if (n_phot > 0) then
      allocate(j(n_phot))
      ! open(2, file='phot_j.dat', status='old')
      ! do i_phot = 1, n_phot
      !   read(2, *) i, j(i_phot)
      ! end do
      ! close(2)
      ! do i_phot = 1, n_phot
      !   rxn => mechanism%get_rxn(i_j(i_phot))
      !   select type (rxn)
      !   class is (rxn_photolysis_t)
      !     flag = rxn%property_set%get_real('base rate', j(i_phot))
      !     if (.not. flag) stop 'Photolysis rate not found'
      !     call j_update(i_phot)%set_rate(j(i_phot))
      !   end select
      ! end do
    end if

    !--- Initialize aerosol-rep update objects, solver, and state
    call assert(209301925, camp_core%get_aero_rep(aero_rep_key, aero_rep_ptr))
    select type (aero_rep_ptr)
    type is (aero_rep_modal_binned_mass_t)
      call camp_core%initialize_update_object(aero_rep_ptr, update_data_GMD)
      call camp_core%initialize_update_object(aero_rep_ptr, update_data_GSD)
      call camp_core%solver_initialize()
      camp_state => camp_core%new_state()
      do n = 1, 4
        call assert_msg(431201141, aero_rep_ptr%get_section_id(mode_names(n), mode(n)), 'Could not get mode ID')
        call update_data_GMD%set_GMD(mode(n), aero_state%GMD(n))
        call update_data_GSD%set_GSD(mode(n), aero_state%GSD(n))
        call camp_core%update_data(update_data_GMD)
        call camp_core%update_data(update_data_GSD)
      end do
      if (n_emis > 0) then
        do i_emis = 1, n_emis
          call camp_core%update_data(q_update(i_emis))
        end do
      end if
      if (n_phot > 0) then
        do i_phot = 1, n_phot
          call camp_core%update_data(j_update(i_phot))
        end do
      end if
    class default
      stop 'MAM4 CAMP: unsupported aerosol representation'
    end select

    !--- Pull chem spec data and set environment
    if (.not. camp_core%get_chem_spec_data(chem_spec_data)) stop 'CAMP chem_spec_data missing'
    call env_state_set_camp_env_state(env_state, camp_state)

    allocate(names(chem_spec_data%size()))
    allocate(aero_names(size(aero_rep_ptr%unique_names())))
    names = chem_spec_data%get_spec_names()
    aero_names = aero_rep_ptr%unique_names()

    if (first_step) then
      first_step = .false.

      camp_state%state_var = 0.0_r8

      allocate(persistent_state(size(camp_state%state_var)))
      persistent_state = 0.0_r8

      !--- Initial gas state from file (name, value) per line
      allocate(ic_spec(chem_spec_data%size()), ic(chem_spec_data%size()))
      open(3, file='ic_sulfate_condensation.dat', status='old')
      do i_ic = 1, size(names)
        read(3, *) ic_spec(i_ic), ic(i_ic)
        if (chem_spec_data%gas_state_id(trim(ic_spec(i_ic))) > 0) then
          camp_state%state_var(chem_spec_data%gas_state_id(trim(ic_spec(i_ic)))) = ic(i_ic)
        end if
      end do
      close(3)

      !--- Initialize aerosol-phase state from mass fractions
      !    Convert modal lognormal to total volume, then split by composition mf
      aero_vol   = (pi/6.0_r8) * aero_state%numc * aero_state%GMD**3 &
                 * exp(4.5_r8 * log(aero_state%GSD) * log(aero_state%GSD))
      rtmpdens = 0.0_r8

      call getcwd(cwd)
      allocate(aero_ids(size(aero_names)))
      allocate(aero_ic_names(size(aero_names)))
      allocate(aero_state%mf_aer(size(aero_names)))

      open(4, file=trim(cwd)//'/aero_mass_fracs.dat', status='old')
      do i_ic = 1, size(aero_names)
        read(4, *) aero_mode, aero_spec, aero_mass_frac, aero_dens
        do i = 1, size(aero_names)
          aero_id   = aero_rep_ptr%spec_state_id(trim(aero_names(i)%string))
          mode_name = mode_extract(trim(aero_names(i)%string), '.')
          aero_name = aero_rep_ptr%spec_name(trim(aero_names(i)%string))
          if ( trim(mode_name) == trim(aero_mode) .and. trim(aero_name) == trim(aero_spec) ) then
            aero_ids(i_ic)      = aero_id
            aero_ic_names(i_ic) = aero_names(i)%string
            do n = 1, 4
              if (trim(aero_mode) == trim(mode_names(n))) then
                rtmpdens(n) = rtmpdens(n) + aero_mass_frac / aero_dens
                exit
              end if
            end do
            exit
          end if
        end do
        aero_state%mf_aer(i_ic) = aero_mass_frac
      end do
      close(4)
      tmpdens = 1.0_r8 / rtmpdens

      do i = 1, size(aero_names)
        mode_name = mode_extract(trim(aero_ic_names(i)), '.')
        do n = 1, 4
          if (trim(mode_name) == trim(mode_names(n))) then
            camp_state%state_var(aero_ids(i)) = aero_vol(n) * tmpdens(n) * aero_state%mf_aer(i)
            exit
          end if
        end do
      end do

      call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)
      call mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)

      deallocate(aero_ids, aero_ic_names, ic, ic_spec)

    else
      camp_state%state_var = persistent_state
      call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)
      call mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)
    end if

    !--- Solve chemistry
    call camp_core%solve(camp_state, del_t, solver_stats = solver_stats)
    if (solver_stats%solver_flag /= 0) then
      write(*,*) 'CAMP solver failed with code ', solver_stats%solver_flag
      stop
    end if

    !--- Pull updated aerosol and gas back to MAM4
    call mam4_camp_interface_get_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)
    call gas_state_get_camp_conc(gas_state, camp_state, camp_core)

    !--- Persist full state_var between time steps
    persistent_state = camp_state%state_var

    ! Cleanup
    deallocate(camp_core, camp_state)
    if (allocated(j_update)) deallocate(j_update, j, i_j)
    if (allocated(q_update)) deallocate(q_update, q, i_q)
    deallocate(names, aero_names)

  end subroutine mam4_camp_interface_solve
  !---------------------------------------------------------------------------

  !> Push MAM4 aerosol modal masses into CAMP state
  subroutine mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_data)
    type(aero_state_t), intent(inout) :: aero_state
    type(camp_core_t),  pointer, intent(inout) :: camp_core
    type(camp_state_t), pointer, intent(inout) :: camp_state
    class(aero_rep_data_t), pointer :: aero_rep_data

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
      if (.not. mode_flag) stop 'Mode not found (set_camp_aerosol)'

      select case (trim(spec))
      case('ASO4'); camp_state%state_var(id) = aero_state%qso4(mode_id)
      case('APOC'); camp_state%state_var(id) = aero_state%qpom(mode_id)
      case('SOA');  camp_state%state_var(id) = aero_state%qsoa(mode_id)
      case('AEC');  camp_state%state_var(id) = aero_state%qbc(mode_id)
      case('ASOIL');camp_state%state_var(id) = aero_state%qdst(mode_id)
      case('ANA');  camp_state%state_var(id) = aero_state%qncl(mode_id)
      end select
    end do

    deallocate(names)
  end subroutine mam4_camp_interface_set_camp_aerosol

  !> Pull CAMP aerosol modal masses back into MAM4
  subroutine mam4_camp_interface_get_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_data)
    type(aero_state_t), intent(inout) :: aero_state
    type(camp_core_t),  pointer, intent(inout) :: camp_core
    type(camp_state_t), pointer, intent(inout) :: camp_state
    class(aero_rep_data_t), pointer :: aero_rep_data

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
      if (.not. mode_flag) stop 'Mode not found (get_camp_aerosol)'

      select case (trim(spec))
      case('ASO4'); aero_state%qso4(mode_id) = camp_state%state_var(id)
      case('APOC'); aero_state%qpom(mode_id) = camp_state%state_var(id)
      case('SOA');  aero_state%qsoa(mode_id) = camp_state%state_var(id)
      case('AEC');  aero_state%qbc(mode_id)  = camp_state%state_var(id)
      case('ASOIL');aero_state%qdst(mode_id) = camp_state%state_var(id)
      case('ANA');  aero_state%qncl(mode_id) = camp_state%state_var(id)
      end select
    end do

    deallocate(names)
  end subroutine mam4_camp_interface_get_camp_aerosol

  !> Extract the section/mode name before the first '.' in a label
  character(len=16) function mode_extract(str, sep)
    character(len=*), intent(in) :: str, sep
    integer :: i, char_count

    char_count = 0
    do i = 1, len(str)
      if (str(i:i) == sep) exit
      char_count = char_count + 1
    end do
    mode_extract = str(1:char_count)
  end function mode_extract

#endif  ! MAM4_USE_CAMP

end module mam4_camp_interface
