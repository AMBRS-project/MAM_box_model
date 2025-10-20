subroutine mam4_camp_interface_solve(env_state, aero_state, gas_state, del_t)
#ifdef MAM4_USE_CAMP
use camp_rxn_data               ! for types
implicit none
!> Environment / states
type(env_state_t),  intent(inout) :: env_state
type(aero_state_t), intent(inout) :: aero_state
type(gas_state_t),  intent(inout) :: gas_state
real(kind=r8),      intent(in)    :: del_t

! CAMP objects
type(camp_core_t),  pointer :: camp_core
type(camp_state_t), pointer :: camp_state
type(mechanism_data_t), pointer :: mechanism
type(solver_stats_t), target :: solver_stats

! Aero-rep + modal geometry updates
character(len=4),  parameter :: aero_rep_key = "MAM4"
character(len=16), parameter :: mode_names(4) = (/ "accumulation    ", &
                                                  "aitken          ", &
                                                  "coarse          ", &
                                                  "primary_carbon  " /)
integer(kind=i_kind) :: mode(4)
class(aero_rep_data_t), pointer :: aero_rep_ptr
class(aero_rep_modal_binned_mass_t), pointer :: aero_rep_data
type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

! mechanism pick + config path
type(string_t), allocatable :: mechs(:)
character(len=512) :: cfg
integer :: n, nlen, istat

! --- init CAMP core from env var (fallback local file) ---
cfg = 'mam4_config.json'
call get_environment_variable('MAM4_CAMP_CONFIG', cfg, length=nlen, status=istat)
if (istat == 0 .and. nlen > 0) cfg = cfg(1:nlen)
camp_core => camp_core_t(trim(cfg))
call camp_core%initialize()

! --- mechanism selection (first listed) ---
mechs = camp_core%mechanism_names()
call assert(260845179, camp_core%get_mechanism(trim(mechs(1)%string), mechanism))

! --- aerosol representation + modal geometry from MAM state ---
call assert(209301925, camp_core%get_aero_rep(aero_rep_key, aero_rep_ptr))
select type (aero_rep_ptr)
 type is (aero_rep_modal_binned_mass_t)
   call camp_core%initialize_update_object(aero_rep_ptr, update_data_GMD)
   call camp_core%initialize_update_object(aero_rep_ptr, update_data_GSD)
   call camp_core%solver_initialize()
   camp_state => camp_core%new_state()

   do n = 1, 4
     call assert_msg(431201141, aero_rep_ptr%get_section_id(mode_names(n), mode(n)), "Could not get mode ID")
     call update_data_GMD%set_GMD(mode(n), aero_state%GMD(n))
     call update_data_GSD%set_GSD(mode(n), aero_state%GSD(n))
     call camp_core%update_data(update_data_GMD)
     call camp_core%update_data(update_data_GSD)
   end do
class default
   write(*,*) "wrong aerosol rep"
   stop 3
end select

! --- map env + gas + aerosol from MAM into CAMP ---
call env_state_set_camp_env_state(env_state, camp_state)
call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)
call mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)

! --- persistent carryover across steps ---
if (first_step) then
 first_step = .false.
 allocate(persistent_state(size(camp_state%state_var)))
 persistent_state = camp_state%state_var
else
 camp_state%state_var = persistent_state
end if

! --- solve chemistry/partitioning ---
call camp_core%solve(camp_state, del_t, solver_stats=solver_stats)
if (solver_stats%solver_flag .ne. 0) then
 write(*,*) 'Solver failed with code ', solver_stats%solver_flag
 stop
end if

! --- map back to MAM state ---
call mam4_camp_interface_get_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)
call gas_state_get_camp_conc(gas_state, camp_state, camp_core)

persistent_state = camp_state%state_var

deallocate(camp_core, camp_state)
#endif
end subroutine mam4_camp_interface_solve