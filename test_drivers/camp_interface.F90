!> \file
!> The mam4_camp_interface module.

!> An interface between MAM4 and the CAMP
module mam4_camp_interface
#ifdef MAM4_USE_CAMP
    use shr_kind_mod, only: r8 => shr_kind_r8
    use mam4_state
    use camp_camp_core
    use camp_camp_state
    use camp_chem_spec_data
    use camp_aero_rep_data
    use camp_aero_rep_modal_binned_mass
    use camp_constants
    use camp_util
    ! use json_module
    use camp_solver_stats
    use camp_mechanism_data
    use camp_rxn_data
    use camp_rxn_emission
    use camp_rxn_photolysis
    use camp_rxn_first_order_loss

    implicit none

    public

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Run the CAMP module for the current MAM4 state
    subroutine mam4_camp_interface_solve( env_state, aero_state, gas_state, del_t)

    !> CAMP core
    type(camp_core_t), pointer :: camp_core
    !> CAMP state
    type(camp_state_t), pointer :: camp_state
    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol state
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state
    !> Time step (s)
    real(kind=r8), intent(in) :: del_t
    integer n, i, i_ic
    character(len=4), parameter :: aero_rep_key = "MAM4"
    character(len=16), parameter :: mode_names(4) = (/ "accumulation    ", &
                                                        "aitken          ", &
                                                        "coarse          ", &
                                                        "primary_carbon  " /)
    real(kind=r8), parameter :: pi = 3.14159265358979323846_r8
    character(len=255), allocatable :: ic_spec(:), aero_ic_names(:)
    character(len=255) :: aero_mode, aero_spec, mode_name, aero_name
    real(kind=r8) :: aero_mass_frac, aero_dens, aero_vol(4), tmpdens(4), rtmpdens(4)
    type(string_t), allocatable :: mech_names(:), aero_names(:)
    integer :: mode(4)
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    type(aero_rep_modal_binned_mass_t), pointer :: aero_rep_data
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD
    type(mechanism_data_t), pointer :: mechanism
    type(rxn_update_data_emission_t), allocatable :: q_update(:) !> Emissions
    type(rxn_update_data_photolysis_t), allocatable :: j_update(:) !> Photolysis
    type(rxn_update_data_first_order_loss_t), allocatable :: loss_update(:) !> Loss
    real(kind=r8), allocatable :: j(:), q(:), loss(:), ic(:), aero_mass_fracs(:)
    integer, allocatable :: i_j(:), i_q(:), i_l(:), aero_ids(:)
    integer n_emis, n_phot, i_emis, i_phot, persistent_id, aero_id, gas_id, dummy
    integer n_loss, i_loss
    class(rxn_data_t), pointer :: rxn
    type(chem_spec_data_t), pointer :: chem_spec_data
    logical :: flag
    type(solver_stats_t), target :: solver_stats
    integer, parameter :: gas_kind = 1, &
                            variable_kind = 1, &
                            constant_kind = 2, &
                            tracer_kind = 1
    character(len=255) :: config_key, mech_key

    namelist /camp_config/ config_key
    namelist /camp_mech/ mech_key

    open(0, file='camp_nml', status='old')
        read(0, camp_config)
        read(0, camp_mech)
    close(0)

    camp_core => camp_core_t( trim(config_key) )
    call camp_core%initialize()

    !> Count emission, photolysis, and loss reactions
    n_emis = 0
    n_phot = 0
    n_loss = 0

    call assert(260845179, camp_core%get_mechanism( trim(mech_key), mechanism ))
    do i = 1, mechanism%size()
        rxn => mechanism%get_rxn(i)
        select type (rxn)
            type is (rxn_photolysis_t)
                n_phot = n_phot + 1
            type is (rxn_emission_t)
                n_emis = n_emis + 1
            type is (rxn_first_order_loss_t)
                n_loss = n_loss + 1
        end select
    end do

    if (n_emis > 0) allocate(q_update(n_emis), i_q(n_emis))
    if (n_phot > 0) allocate(j_update(n_phot), i_j(n_phot))
    if (n_loss > 0) allocate(loss_update(n_loss), i_l(n_loss))

    !> Initialize emission and photolysis update objects
    i_emis = 0
    i_phot = 0
    i_loss = 0
    do i = 1, mechanism%size()
        rxn => mechanism%get_rxn(i)
        select type (rxn)
            type is (rxn_photolysis_t)
                i_phot = i_phot + 1
                i_j(i_phot) = i
                call camp_core%initialize_update_object(rxn, j_update(i_phot))
            type is (rxn_emission_t)
                i_emis = i_emis + 1
                i_q(i_emis) = i
                call camp_core%initialize_update_object(rxn, q_update(i_emis))
            type is (rxn_first_order_loss_t)
                i_loss = i_loss + 1
                i_l(i_loss) = i
                call camp_core%initialize_update_object(rxn, loss_update(i_loss))
        end select
    end do

    !> Set emission rates
    if (n_emis > 0) then
        allocate(q(n_emis))
        do i_emis = 1, n_emis
            rxn => mechanism%get_rxn(i_q(i_emis))
            select type (rxn)
                type is (rxn_emission_t)
                    flag = rxn%property_set%get_real('base rate', q(i_emis))
                    if ( .not.flag ) then
                        write(*,*) 'Emission rate not found'
                        stop
                    end if
                    call q_update(i_emis)%set_rate(q(i_emis))
            end select
        end do
    end if

    !> Set photolysis rates
    if (n_phot > 0) then
        allocate(j(n_phot))
        do i_phot = 1, n_phot
            rxn => mechanism%get_rxn(i_j(i_phot))
            select type (rxn)
                type is (rxn_photolysis_t)
                    flag = rxn%property_set%get_real('base rate', j(i_phot))
                    if ( .not.flag ) then
                        write(*,*) 'Photolysis rate not found'
                        stop
                    end if
                    call j_update(i_phot)%set_rate(j(i_phot))
            end select
        end do
    end if

    !> Set loss rates
    if (n_loss > 0) then
        allocate(loss(n_loss))
        do i_loss = 1, n_loss
            rxn => mechanism%get_rxn(i_l(i_loss))
            select type (rxn)
                type is (rxn_first_order_loss_t)
                    flag = rxn%property_set%get_real('base rate', loss(i_loss))
                    if ( .not.flag ) then
                        write(*,*) 'Loss rate not found'
                        stop
                    end if
                    call loss_update(i_loss)%set_rate(loss(i_loss))
            end select
        end do
    end if

    !> Initialize modal aero update objects and solver; update rates
    call assert(209301925, camp_core%get_aero_rep(aero_rep_key, aero_rep_ptr))
    select type (aero_rep_ptr)
        type is (aero_rep_modal_binned_mass_t)
        call camp_core%initialize_update_object(aero_rep_ptr, &
                                                    update_data_GMD)
        call camp_core%initialize_update_object(aero_rep_ptr, &
                                                    update_data_GSD)
        call camp_core%solver_initialize()
        camp_state => camp_core%new_state()
        do n = 1, 4
            call assert_msg(431201141, &
                    aero_rep_ptr%get_section_id(mode_names(n), mode(n)), &
                    "Could not get mode ID")
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
            write(*,*) "wrong aerosol rep"
            stop 3
    end select

    if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
        write(*,*) "Something's gone wrong!"
        stop 3
    end if

    ! Set the CAMP environmental state.
    call env_state_set_camp_env_state(env_state, camp_state)

    allocate( mech_names(chem_spec_data%size(variable_kind,gas_kind) + chem_spec_data%size(constant_kind,gas_kind)), &
                aero_names( size( aero_rep_ptr%unique_names( tracer_type=tracer_kind ) ) ) )

    mech_names = [chem_spec_data%get_spec_names(variable_kind,gas_kind), chem_spec_data%get_spec_names(constant_kind,gas_kind)]
    aero_names = aero_rep_ptr%unique_names( tracer_type=tracer_kind )

    if ( first_step ) then
    
        first_step = .false.

        camp_state%state_var = 0.0_r8
        
        allocate( persistent_state(size(camp_state%state_var)), &
                    ic_spec(chem_spec_data%size(variable_kind,gas_kind) + chem_spec_data%size(constant_kind,gas_kind)), &
                    ic(chem_spec_data%size(variable_kind,gas_kind) + chem_spec_data%size(constant_kind,gas_kind)) )

        persistent_state = 0.0_r8

        !> Load initial gas concentrations
        open(3, file = 'ic.dat', status='old')
        do i_ic = 1, size(mech_names)
            read(3, *) ic_spec(i_ic), ic(i_ic)
            ! if (chem_spec_data%gas_state_id( trim( ic_spec(i_ic)) ) > 0) then
            camp_state%state_var( chem_spec_data%gas_state_id( trim( ic_spec(i_ic)) ) ) = ic(i_ic)
            ! end if
        end do
        close(3)

        !> Calculate mode volume concentrations for converting mass fractions into mass concentrations
        aero_vol = (pi/6.0_r8) * aero_state%numc * aero_state%GMD**3 * exp(4.5_r8 * log(aero_state%GSD)**2)
        rtmpdens = 0.0_r8

        !> Read initial aerosol mass fractions, map to mechanism species, and convert to mass concentrations

        allocate( aero_ids( size( aero_names ) ), aero_state%mf_aer( size( aero_names ) ), aero_ic_names( size( aero_names ) ) )
        open(4, file = 'aero_mass_fracs.dat', status='old')
        do i_ic = 1, size(aero_names)
            read(4, *) aero_mode, aero_spec, aero_mass_frac, aero_dens
            do i = 1, size(aero_names)
                aero_id = aero_rep_ptr%spec_state_id( trim(aero_names(i)%string) )
                mode_name = mode_extract( trim(aero_names(i)%string), '.' )
                aero_name = aero_rep_ptr%spec_name( trim(aero_names(i)%string) )
                if ( trim(mode_name) == trim(aero_mode) .and. trim(aero_name) == trim(aero_spec) ) then
                    aero_ids(i_ic) = aero_id
                    aero_ic_names(i_ic) = aero_names(i)%string
                    do n = 1, 4
                        if ( trim(aero_mode) == trim(mode_names(n)) ) then
                            !> Harmonic mean particulate density, reciprocal
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

        !> Map aerosol mass concentrations to camp state object
        do i = 1, size(aero_names)
            mode_name = mode_extract( trim(aero_ic_names(i)), '.' )
            do n = 1, 4
                if ( trim(mode_name) == trim(mode_names(n)) ) then
                    camp_state%state_var( aero_ids(i) ) = aero_vol(n) * tmpdens(n) * aero_state%mf_aer(i)
                    ! write(*,*) trim(aero_ic_names(i)), camp_state%state_var( aero_ids(i) )
                    exit
                end if
            end do
        end do

        call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)
        call mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)

        deallocate( aero_ids, aero_ic_names, ic, ic_spec, aero_state%mf_aer )
        
    else

        !> Load persistent CAMP state
        camp_state%state_var = persistent_state

        call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)
        call mam4_camp_interface_set_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)
            
    end if

    ! do i = 1, size(aero_names)
    !     write(*,*) aero_names(i)%string
    ! end do
    
    ! Solve the multi-phase chemical system
    call camp_core%solve(camp_state, del_t, solver_stats = solver_stats)
    !call solver_stats%print()
    if (solver_stats%solver_flag /= 0) then
        write(*,*) 'Solver failed with code ', solver_stats%solver_flag
        stop
    end if

    call mam4_camp_interface_get_camp_aerosol(aero_state, camp_core, camp_state, aero_rep_ptr)
    ! Update the MAM4 gas-phase state
    call gas_state_get_camp_conc(gas_state, camp_state, camp_core)

    !> Set persistent CAMP state
    persistent_state = camp_state%state_var

    deallocate( camp_core, camp_state )
    if ( n_phot > 0 ) deallocate( j_update, j, i_j )
    if ( n_emis > 0 ) deallocate( q_update, q, i_q )
    deallocate( mech_names, aero_names )

    end subroutine mam4_camp_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Set the CAMP aerosol-phase species and mass concentrations
    subroutine mam4_camp_interface_set_camp_aerosol(aero_state, &
        camp_core, camp_state, aero_rep_data)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> CAMP core.
    type(camp_core_t), pointer, intent(inout) :: camp_core
    !> CAMP state.
    type(camp_state_t), pointer, intent(inout) :: camp_state
    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_data
    type(string_t), allocatable :: names(:)
    integer i, id, mode_id
    character(len=255) :: spec, mode_name
    logical :: mode_flag
    integer, parameter :: tracer_kind = 1

    allocate( names(size(aero_rep_data%unique_names( tracer_type=tracer_kind ))) )
    names = aero_rep_data%unique_names( tracer_type=tracer_kind )

    do i = 1, size(names)
        select type (aero_rep_data)
            type is (aero_rep_modal_binned_mass_t)
                id = aero_rep_data%spec_state_id( trim(names(i)%string) )
                spec = aero_rep_data%spec_name( trim(names(i)%string) )
                mode_name = mode_extract( trim(names(i)%string), '.' )
                mode_flag = aero_rep_data%get_section_id( trim(mode_name), mode_id )             
        if (.not.mode_flag) then
            write(*,*) 'Mode not found'
            stop
        end if        
            select case( trim(spec) )
                case('SO4')
                    if (lso4(mode_id) > 0) camp_state%state_var( id ) = aero_state%qso4(mode_id)
                    ! write(*,*) trim(mode_name)//'.'//trim(spec), camp_state%state_var(id)
                case('POM')
                    if (lpom(mode_id) > 0) camp_state%state_var( id ) = aero_state%qpom(mode_id)
                case('SOA')
                    if (lsoa(mode_id) > 0) camp_state%state_var( id ) = aero_state%qsoa(mode_id)
                case('BC')
                    if (lbc(mode_id) > 0) camp_state%state_var( id ) = aero_state%qbc(mode_id)
                case('DST')
                    if (ldst(mode_id) > 0) camp_state%state_var( id ) = aero_state%qdst(mode_id)
                case('NCL')
                    if (lncl(mode_id) > 0) camp_state%state_var( id ) = aero_state%qncl(mode_id)
                case('MOM')
                    if (lmom(mode_id) > 0) camp_state%state_var( id ) = aero_state%qmom(mode_id)
                case('H2O_aq')
                    camp_state%state_var( id ) = aero_state%qaerwat(mode_id)
            end select
        end select
    end do

    deallocate(names)

    end subroutine mam4_camp_interface_set_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Get the CAMP aerosol-phase species and mass concentrations
    subroutine mam4_camp_interface_get_camp_aerosol(aero_state, &
        camp_core, camp_state, aero_rep_data)
        
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> CAMP core.
    type(camp_core_t), pointer, intent(inout) :: camp_core
    !> CAMP state.
    type(camp_state_t), pointer, intent(inout) :: camp_state
    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_data
    type(string_t), allocatable :: names(:)
    integer i, id, mode_id
    character(len=255) :: spec, mode_name
    logical :: mode_flag
    integer, parameter :: tracer_kind = 1

    allocate(names(size(aero_rep_data%unique_names(tracer_type=tracer_kind))))
    names = aero_rep_data%unique_names(tracer_type=tracer_kind)

    do i = 1, size(names)
        ! write(*,*) names(i)%string
        select type (aero_rep_data)
            type is (aero_rep_modal_binned_mass_t)
                id = aero_rep_data%spec_state_id( trim(names(i)%string) )
                spec = aero_rep_data%spec_name( trim(names(i)%string) )
                mode_name = mode_extract( trim(names(i)%string), '.' )
                mode_flag = aero_rep_data%get_section_id( trim(mode_name), mode_id )
        if (.not.mode_flag) then
            write(*,*) 'Mode not found'
            stop
        end if
            select case( trim(spec) )
                case('SO4')
                    if (lso4(mode_id) > 0) aero_state%qso4(mode_id) = camp_state%state_var( id )
                    ! write(*,*) trim(mode_name)//'.'//trim(spec), aero_state%qso4(mode_id)
                case('POM')
                    if (lpom(mode_id) > 0) aero_state%qpom(mode_id) = camp_state%state_var( id )
                case('SOA')
                    if (lsoa(mode_id) > 0) aero_state%qsoa(mode_id) = camp_state%state_var( id )
                case('BC')
                    if (lbc(mode_id) > 0) aero_state%qbc(mode_id) = camp_state%state_var( id )
                case('DST')
                    if (ldst(mode_id) > 0) aero_state%qdst(mode_id) = camp_state%state_var( id )
                case('NCL')
                    if (lncl(mode_id) > 0) aero_state%qncl(mode_id) = camp_state%state_var( id )
                case('MOM')
                    if (lmom(mode_id) > 0) aero_state%qmom(mode_id) = camp_state%state_var( id )
                case('H2O_aq')
                    aero_state%qaerwat(mode_id) = camp_state%state_var( id )
            end select
        end select 
    end do
    
    deallocate(names)

    end subroutine mam4_camp_interface_get_camp_aerosol

    character(len=16) function mode_extract(str, sep)
        !> Read mode name from CAMP unique aerosol
        !> name pattern.

        character(len=*) :: str, sep
        integer i, char_count

        char_count = 0
        do i = 1, len(str)
            if (str(i:i) == sep) exit
            char_count = char_count + 1
        end do

        mode_extract = str(1:char_count)

    end function mode_extract

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
end module mam4_camp_interface
    