module mam4_camp

    use shr_kind_mod, only: r8 => shr_kind_r8
    use constituents, only: cnst_get_ind
    use chem_mods, only: imozart, adv_mass
    use mo_tracname, only: solsym
    use physconst, only: mwdry
    ! use rad_constituents, only: sigmag_amode_rc
    use modal_aero_data, only: sigmag_amode, specmw_amode
    use modal_aero_initialize_data, only: search_list_of_names

    use camp_camp_core
    use camp_camp_state
    use camp_chem_spec_data

    use camp_aero_rep_data
    use camp_aero_rep_modal_binned_mass

    use camp_mechanism_data
    use camp_rxn_data
    use camp_rxn_emission

    use camp_solver_stats
    use camp_util

    implicit none

    type mam4_env_state_t
        ! Temperature [K].
        real( kind=r8 ) :: temp
        ! Relative humidity [1].
        real( kind=r8 ) :: RH_CLEA
        ! Ambient pressure [Pa].
        real( kind=r8 ) :: press
    end type mam4_env_state_t

    type( mam4_env_state_t ) :: mam4_env_state
    real( r8), allocatable, save :: persistent_camp_state_var( : )
    logical, save :: persistence_initialized = .false.

    real( r8 ) :: dgn( 4 )
    real( r8 ) :: aircon
    integer :: mdo_gasaerexch

    character( len=14 ), parameter :: modename_aliases( 4 ) = (/ 'accumulation  ', &
                                                                 'aitken        ', &
                                                                 'coarse        ', &
                                                                 'primary_carbon' /)
    character( len=1 ), parameter :: mode_ints_as_chars( 4 ) = (/ '1', '2', '3', '4' /)
    character( len=3 ), parameter :: specname_aliases( 9 ) = (/ 'SO4', &
                                                                'NH4', &
                                                                'NO3', &
                                                                'POM', &
                                                                'SOA', &
                                                                'BC ', &
                                                                'NCL', &
                                                                'DST', &
                                                                'MOM' /)
    character( len=3 ), parameter :: specname_lowcaps( 9 ) = (/ 'so4', &
                                                                'nh4', &
                                                                'no3', &
                                                                'pom', &
                                                                'soa', &
                                                                'bc ', &
                                                                'ncl', &
                                                                'dst', &
                                                                'mom' /)

    contains

    character( len=16 ) function read_mode( str, sep )

        character( len=* ) :: str, sep
        integer i, char_count

        char_count = 0
        do i = 1, len( str )
            if ( str( i:i ) == sep ) exit
            char_count = char_count + 1
        end do

        read_mode = str( 1:char_count )

    end function read_mode

    subroutine set_camp_env_state( camp_state )

        type( camp_state_t ), pointer, intent( out ) :: camp_state

        real( kind=r8 ) temp, RH_CLEA, press

        namelist /met_input/ temp, press, RH_CLEA
        open( 0, file='namelist', status='old' )
            read( 0,met_input )
        close( 0 )

        mam4_env_state%temp = temp
        mam4_env_state%press = press
        mam4_env_state%RH_CLEA = RH_CLEA

        camp_state%env_states(1)%val%temp = mam4_env_state%temp
        camp_state%env_states(1)%val%rel_humid = mam4_env_state%RH_CLEA
        camp_state%env_states(1)%val%pressure = mam4_env_state%press

        return

    end subroutine set_camp_env_state

    subroutine set_camp_gas_state( mam4,camp_state,chem_spec_data )

        real( kind=r8 ), intent( in ) :: mam4( : ) ! vmr
        type( camp_state_t ), pointer, intent( out ) :: camp_state

        type( chem_spec_data_t ), pointer :: chem_spec_data

        integer( kind=i_kind ) :: i
        integer :: j, offset

        real( kind=r8 ), parameter :: t_steam = 373.15_r8
        real( kind=r8 ) :: a, water_vp
        real( kind=r8 ), parameter :: confac = 1.0e6_r8 ! [mol mol-1] -> [ppm]

        a = 1.0_r8 - t_steam / mam4_env_state%temp
        a = ( ( ( -0.1299_r8 * a - 0.6445_r8 ) * a - 1.976_r8 ) * a + 13.3185_r8 ) * a
        water_vp = 101325.0_r8 * exp( a )
        camp_state%state_var( chem_spec_data%gas_state_id( 'H2O' ) ) = &
            mam4_env_state%RH_CLEA * water_vp * confac / mam4_env_state%press

        offset = imozart - 1

        do i = 1_i_kind, size( camp_state%state_var )
            call cnst_get_ind( chem_spec_data%gas_state_name( i ),j,.false. )
            if ( j > 0 ) then
                camp_state%state_var( i ) = confac * mam4( j-offset )
            end if
        end do

        return

    end subroutine set_camp_gas_state

    subroutine get_camp_gas_state( mam4,camp_state,chem_spec_data )

        real( kind=r8 ), intent( out ) :: mam4( : ) ! vmr
        type( camp_state_t ), pointer, intent( in ) :: camp_state

        type( chem_spec_data_t ), pointer :: chem_spec_data

        integer( kind=i_kind ) :: i
        integer :: j, offset, count = 0

        real( kind=r8 ), parameter :: confac = 1.0e-6_r8 ! [ppm] -> [mol mol-1]

        offset = imozart - 1

        do i = 1_i_kind, size( camp_state%state_var )
            call cnst_get_ind( chem_spec_data%gas_state_name( i ),j,.false. )
            if ( j > 0 ) then
                mam4( j-offset ) = confac * camp_state%state_var( i )
            end if
        end do

        if ( mdo_gasaerexch == 0 ) then
            call cnst_get_ind( 'SOAG',j,.false. )
            mam4( j-offset ) = 0.0_r8
            do i = 1_i_kind, size( camp_state%state_var )
                if ( index( chem_spec_data%gas_state_name( i ), 'SOAG' ) > 0 ) then
                    if ( chem_spec_data%gas_state_name( i ) /= 'SOAG' ) then
                        mam4( j-offset ) = mam4( j-offset ) + confac * camp_state%state_var( i )
                    end if
                end if
            end do
        end if

        return

    end subroutine get_camp_gas_state

    subroutine set_camp_aero_state( mam4,camp_state,aero_rep_data )

        real( kind=r8 ), intent( in ) :: mam4( : )
        type( camp_state_t ), pointer, intent( out ) :: camp_state

        class( aero_rep_data_t ), pointer :: aero_rep_data

        type( string_t ), allocatable :: names( : )
        character( len=1 ) :: mode_id

        integer( kind=i_kind ) :: i
        integer :: j, l, n, offset, id

        ! offset = imozart - 1
        offset = 0

        select type ( aero_rep_data )
            type is ( aero_rep_modal_binned_mass_t )
                allocate( names( size( aero_rep_data%unique_names() ) ) )
                names = aero_rep_data%unique_names()
                do i = 1_i_kind, size( names )
                    call search_list_of_names( &
                        read_mode( trim( names( i )%string ), '.' ), &
                        n, &
                        modename_aliases, &
                        size( modename_aliases ) &
                    )
                    if ( n < 0 ) then
                        write(*,*) 'Cannot identify mode'
                        stop 3
                    end if
                    mode_id = mode_ints_as_chars( n )
                    call search_list_of_names( &
                        aero_rep_data%spec_name( trim( names( i )%string ) ), &
                        l, &
                        specname_aliases, &
                        size( specname_aliases ) &
                        )
                    if ( l < 0 ) then
                        write(*,*) 'Cannot identify species'
                        stop 3
                    end if
                    call search_list_of_names( &
                        trim( specname_lowcaps( l ) ) // '_a' // mode_id, &
                        j, &
                        solsym, &
                        size( solsym ) &
                    )
                    if ( j < 0 ) then
                        write(*,*) 'Cannot map species'
                        stop 3
                    end if
                    id = aero_rep_data%spec_state_id( trim( names( i )%string ) )
                    ! camp_state%state_var( id ) = mam4( j ) * aircon
                    ! camp_state%state_var( id ) = mam4( j ) * aircon * 0.012011_r8
                    ! camp_state%state_var( id ) = mam4( j ) * aircon * ( 1e-3_r8 * specmw_amode( l ) )
                    camp_state%state_var( id ) = mam4( j-offset ) * aircon * ( 1e-3_r8 * adv_mass( j-offset ) )
                end do
        end select
        deallocate( names )

        return

    end subroutine set_camp_aero_state

    subroutine get_camp_aero_state( mam4,camp_state,aero_rep_data )

        real( kind=r8 ), intent( out ) :: mam4( : )
        type( camp_state_t ), pointer, intent( in ) :: camp_state

        class( aero_rep_data_t ), pointer :: aero_rep_data

        type( string_t ), allocatable :: names( : )
        character( len=1 ) :: mode_id

        integer( kind=i_kind ) :: i
        integer :: j, l, n, idens, offset, id

        ! offset = imozart - 1
        offset = 0
        
        select type ( aero_rep_data )
            type is ( aero_rep_modal_binned_mass_t )
                allocate( names( size( aero_rep_data%unique_names() ) ) )
                names = aero_rep_data%unique_names()
                do i = 1_i_kind, size( names )
                    call search_list_of_names( &
                        read_mode( trim( names( i )%string ), '.' ), &
                        n, &
                        modename_aliases, &
                        size( modename_aliases ) &
                    )
                    if ( n < 0 ) then
                        write(*,*) 'Cannot identify mode'
                        stop 3
                    end if
                    mode_id = mode_ints_as_chars( n )
                    call search_list_of_names( &
                        aero_rep_data%spec_name( trim( names( i )%string ) ), &
                        l, &
                        specname_aliases, &
                        size( specname_aliases ) &
                        )
                    if ( l < 0 ) then
                        write(*,*) 'Cannot identify species'
                        stop 3
                    end if
                    call search_list_of_names( &
                        trim( specname_lowcaps( l ) ) // '_a' // mode_id, &
                        j, &
                        solsym, &
                        size( solsym ) &
                    )
                    if ( j < 0 ) then
                        write(*,*) 'Cannot map species'
                        stop 3
                    end if
                    id = aero_rep_data%spec_state_id( trim( names( i )%string ) )
                    ! mam4( j ) = camp_state%state_var( id ) / aircon
                    ! mam4( j ) = camp_state%state_var( id ) / 0.012011_r8 / aircon
                    ! mam4( j ) = camp_state%state_var( id ) / ( 1e-3_r8 * specmw_amode( l ) ) / aircon
                    mam4( j-offset ) = camp_state%state_var( id ) / ( 1e-3_r8 * adv_mass( j-offset ) ) / aircon
                end do
        end select
        deallocate( names )

        return

    end subroutine get_camp_aero_state

    subroutine solve_camp_chemistry( vmr,dt )

        real( kind=r8 ), intent( inout ) :: vmr( : )
        real( kind=r8 ), intent( in ) :: dt

        type( camp_core_t ), pointer :: camp_core
        type( camp_state_t ), pointer :: camp_state
        
        type( chem_spec_data_t ), pointer :: chem_spec_data
        class( aero_rep_data_t ), pointer :: aero_rep_ptr

        character( len=255 ) :: config_key, mech_key
        type( mechanism_data_t ), pointer :: mechanism
        class( rxn_data_t ), pointer :: rxn

        type( rxn_update_data_emission_t ), allocatable :: emission_update( : )
        type( aero_rep_update_data_modal_binned_mass_GMD_t ) :: update_data_GMD
        type( aero_rep_update_data_modal_binned_mass_GSD_t ) :: update_data_GSD
        
        type( solver_stats_t ), target :: solver_stats

        integer( kind=i_kind ) :: i
        integer :: j, n
        real( kind=r8 ) :: rate
        logical :: flag

        namelist /camp_config/ config_key, mech_key
        open( 2, file='namelist', status='old' )
            read( 2,camp_config )
        close( 2 )

        camp_core => camp_core_t( trim( config_key ) )
        call camp_core%initialize()

        if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
            write(*,*) "Something's gone wrong!"
            stop 3
        end if

        if ( .not. camp_core%get_mechanism( trim( mech_key ),mechanism ) ) then
            write(*,*) "Missing mechanism!"
            stop 3
        end if

        if ( .not. camp_core%get_aero_rep( 'Modal/binned', aero_rep_ptr ) ) then
            write(*,*) "Missing aerosol representation!"
            stop 3
        end if

        select type ( aero_rep_ptr )
            type is ( aero_rep_modal_binned_mass_t )
                call camp_core%initialize_update_object( aero_rep_ptr, &
                                                         update_data_GMD )
                call camp_core%initialize_update_object( aero_rep_ptr, &
                                                         update_data_GSD )
        end select

        allocate( emission_update( mechanism%size() ) )
        do i = 1_i_kind, mechanism%size()
            rxn => mechanism%get_rxn( i )
            select type ( rxn )
                type is ( rxn_emission_t )
                    call camp_core%initialize_update_object( rxn,emission_update( i ) )
            end select
        end do
        
        call camp_core%solver_initialize()
        camp_state => camp_core%new_state()

        select type ( aero_rep_ptr )
            type is ( aero_rep_modal_binned_mass_t )
                do j = 1, 4
                    flag = aero_rep_ptr%get_section_id( trim( modename_aliases( j ) ), n )
                    if ( .not. flag ) continue
                    call update_data_GMD%set_GMD( n, dgn( j ) )
                    call update_data_GSD%set_GSD( n, sigmag_amode( j ) )
                    call camp_core%update_data( update_data_GMD )
                    call camp_core%update_data( update_data_GSD )
                end do
        end select

        do i = 1_i_kind, mechanism%size()
            rxn => mechanism%get_rxn( i )
            select type ( rxn )
                type is ( rxn_emission_t )
                    flag = rxn%property_set%get_real( 'rate',rate )
                    call emission_update( i )%set_rate( rate )
                    call camp_core%update_data( emission_update( i ) )
            end select
        end do

        if ( persistence_initialized ) camp_state%state_var = persistent_camp_state_var

        call set_camp_env_state( camp_state )

        call set_camp_gas_state( vmr,camp_state,chem_spec_data )

        call set_camp_aero_state( vmr,camp_state,aero_rep_ptr )

        call camp_core%solve( camp_state,dt,solver_stats=solver_stats )
        if ( solver_stats%solver_flag /= 0 ) then
            call solver_stats%print()
            write(*,*) 'Solver failed with code ', solver_stats%solver_flag
            stop
        end if

        call get_camp_aero_state( vmr,camp_state,aero_rep_ptr )

        call get_camp_gas_state( vmr,camp_state,chem_spec_data )

        if ( .not. persistence_initialized ) then
            allocate( persistent_camp_state_var( size( camp_state%state_var ) ) )
            persistence_initialized = .true.
        end if

        if ( persistence_initialized ) persistent_camp_state_var = camp_state%state_var

        deallocate( camp_core,camp_state )

        return

    end subroutine solve_camp_chemistry

end module mam4_camp
