module mam4_camp

    use shr_kind_mod, only: r8 => shr_kind_r8
    ! use physconst, only: mwdry, r_universal
    ! use constituents, only: pcnst, cnst_get_ind
    use constituents, only: cnst_get_ind
    ! use chem_mods, only: gas_pcnst, adv_mass, imozart
    use chem_mods, only: imozart
    ! use mo_tracname, only: solsym

    use camp_camp_core
    use camp_camp_state
    use camp_chem_spec_data

    use camp_mechanism_data
    use camp_rxn_data
    use camp_rxn_emission
    ! use camp_rxn_factory

    use camp_solver_stats
    ! use camp_constants
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

    ! integer :: mdo_camp_chem

    contains

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
        ! type( camp_core_t ), pointer :: camp_core

        type( chem_spec_data_t ), pointer :: chem_spec_data

        integer( kind=i_kind ) :: i
        integer :: j, offset

        real( kind=r8 ), parameter :: t_steam = 373.15_r8
        real( kind=r8 ) :: a, water_vp
        real( kind=r8 ), parameter :: confac = 1.0e6_r8 ! [mol mol-1] -> [ppm]

        ! if ( .not. camp_core%get_chem_spec_data( chem_spec_data ) ) then
        !     write(*,*) 'Failed to retrieve chem_spec_data!'
        !     stop 3
        ! end if

        a = 1.0_r8 - t_steam / mam4_env_state%temp
        a = ( ( ( -0.1299_r8 * a - 0.6445_r8 ) * a - 1.976_r8 ) * a + 13.3185_r8 ) * a
        water_vp = 101325.0_r8 * exp( a )
        camp_state%state_var( chem_spec_data%gas_state_id( 'H2O' ) ) = &
            mam4_env_state%RH_CLEA * water_vp * 1.0e6_r8 / mam4_env_state%press

        offset = imozart - 1

        j = -1
        ! do i = 1,chem_spec_data%size( spec_phase=chem_spec_gas_phase )
        do i = 1_i_kind, size( camp_state%state_var )
            call cnst_get_ind( chem_spec_data%gas_state_name( i ),j,.false. )
            if ( j /= -1 ) then
                ! camp_state%state_var( &
                !     chem_spec_data%gas_state_id( &
                !         chem_spec_data%gas_state_name( i ) &
                !     ) &
                ! ) = confac * mam4( j-offset )
                camp_state%state_var( i ) = confac * mam4( j-offset )
            end if
        end do

        print *, 'Before CAMP, H2SO4 =', camp_state%state_var( chem_spec_data%gas_state_id( 'H2SO4' ) )

        return

    end subroutine set_camp_gas_state

    subroutine get_camp_gas_state( mam4,camp_state,chem_spec_data )

        real( kind=r8 ), intent( out ) :: mam4( : ) ! vmr
        type( camp_state_t ), pointer, intent( in ) :: camp_state
        ! type( camp_core_t ), pointer :: camp_core

        type( chem_spec_data_t ), pointer :: chem_spec_data

        integer( kind=i_kind ) :: i
        integer :: j, offset

        real( kind=r8 ), parameter :: confac = 1.0e-6_r8 ! [ppm] -> [mol mol-1]

        ! if ( .not. camp_core%get_chem_spec_data( chem_spec_data ) ) then
        !     write(*,*) 'Failed to retrieve chem_spec_data!'
        !     stop 3
        ! end if

        offset = imozart - 1

        j = -1
        ! do i = 1,chem_spec_data%size( spec_phase=chem_spec_gas_phase )
        do i = 1_i_kind, size( camp_state%state_var )
            call cnst_get_ind( chem_spec_data%gas_state_name( i ),j,.false. )
            if ( j /= -1 ) then
                ! mam4( j-offset ) = confac * camp_state%state_var( &
                !     chem_spec_data%gas_state_id( &
                !         chem_spec_data%gas_state_name( i ) &
                !     ) &
                ! )
                mam4( j-offset ) = confac * camp_state%state_var( i )
            end if
        end do

        print *, 'After CAMP, H2SO4 =', camp_state%state_var( chem_spec_data%gas_state_id( 'H2SO4' ) )

        return

    end subroutine get_camp_gas_state

    subroutine solve_camp_chemistry( vmr,dt )

        real( kind=r8 ), intent( inout ) :: vmr( : ) ! vmr
        real( kind=r8 ), intent( in ) :: dt

        type( camp_core_t ), pointer :: camp_core
        type( camp_state_t ), pointer :: camp_state
        
        type( chem_spec_data_t ), pointer :: chem_spec_data

        character( len=255 ) :: config_key, mech_key
        type( mechanism_data_t ), pointer :: mechanism
        ! type( rxn_factory_t ) :: rxn_factory
        class( rxn_data_t ), pointer :: rxn

        type( rxn_update_data_emission_t ) :: emission_update
        
        type( solver_stats_t ), target :: solver_stats

        integer( kind=i_kind ) :: i
        real( kind=r8 ) :: rate
        logical :: flag

        namelist /camp_config/ config_key, mech_key
        open( 1, file='namelist', status='old' )
            read( 1,camp_config )
        close( 1 )

        camp_core => camp_core_t( trim( config_key ) )
        call camp_core%initialize()

        if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
            write(*,*) "Something's gone wrong!"
            stop 3
        end if

        ! call assert( 260845179,camp_core%get_mechanism( trim( mech_key ),mechanism ) )
        if ( .not. camp_core%get_mechanism( trim( mech_key ),mechanism ) ) then
            write(*,*) "Missing mechanism!"
            stop 3
        end if

        do i = 1_i_kind, mechanism%size()
            rxn => mechanism%get_rxn( i )
            select type ( rxn )
                class is ( rxn_emission_t )
                    call camp_core%initialize_update_object( rxn,emission_update )
            end select
        end do
        
        call camp_core%solver_initialize()
        camp_state => camp_core%new_state()

        call set_camp_env_state( camp_state )

        call set_camp_gas_state( vmr,camp_state,chem_spec_data )

        do i = 1_i_kind, mechanism%size()
            rxn => mechanism%get_rxn( i )
            ! call emit_gas( camp_core,emission_update,rxn )
            select type ( rxn )
                class is ( rxn_emission_t )
                    flag = rxn%property_set%get_real( 'rate',rate )
                    call emission_update%set_rate( rate )
                    call camp_core%update_data( emission_update )
            end select
        end do

        call camp_core%solve( camp_state,dt,solver_stats=solver_stats )
        if ( solver_stats%solver_flag /= 0 ) then
            write(*,*) 'Solver failed with code ', solver_stats%solver_flag
            stop
        end if

        call get_camp_gas_state( vmr,camp_state,chem_spec_data )

        deallocate( camp_core,camp_state )

        return

    end subroutine solve_camp_chemistry

end module mam4_camp
