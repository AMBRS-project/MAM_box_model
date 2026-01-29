module mam4_camp
#ifdef ENABLE_CAMP

    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst, only: mwdry, r_universal
    use constituents, only: pcnst, cnst_get_ind
    use chem_mods, only: gas_pcnst, adv_mass
    use mo_tracname, only: solsym

    use camp_camp_core
    use camp_camp_state
    use camp_constants
    use camp_chem_spec_data

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

    character(len=*), parameter :: mam4_modes( 4 ) = (/ 'accumulation', &
                                                        'aitken', &
                                                        'coarse', &
                                                        'primary_carbon' /)

    contains

    subroutine set_camp_env_state( camp_env_state )

        type( camp_state_t ), pointer, intent( out ) :: camp_state

        camp_state%env_states(1)%val%temp = mam4_env_state%temp
        camp_state%env_states(1)%val%rel_humid = mam4_env_state%RH_CLEA
        camp_state%env_states(1)%val%pressure = mam4_env_state%press

        return

    end subroutine set_camp_env_state

    subroutine set_camp_gas_state( mam4,camp_state,camp_core )

        real( kind=r8 ), intent( in ) :: mam4( : ) ! vmr
        type( camp_state_t ), pointer, intent( out ) :: camp_state
        type( camp_core_t ), pointer :: camp_core

        type( chem_spec_data_t ), pointer :: chem_spec_data

        integer :: i,  j

        real( kind=r8 ), parameter :: confac = 1.0e6_r8 ! [mol mol-1] -> [ppm]

        if ( .not. camp_core%get_chem_spec_data( chem_spec_data ) ) then
            write(*,*) 'Failed to retrieve chem_spec_data!'
            stop 3
        end if

        do i = 1,chem_spec_data%size()
            call cnst_get_ind( chem_spec_data%gas_state_name( i ),j,.false. )
            if ( j /= -1 ) then
                camp_state%state_var( &
                    chem_spec_data%gas_state_id( &
                        chem_spec_data%gas_state_name( i ) &
                    ) &
                ) = confac * mam4( j )
            end if
        end do

        return

    end subroutine set_camp_gas_state

    subroutine get_camp_gas_state( mam4,camp_state,camp_core )

        real( kind=r8 ), intent( out ) :: mam4( : ) ! vmr
        type( camp_state_t ), pointer, intent( in ) :: camp_state
        type( camp_core_t ), pointer :: camp_core

        type( chem_spec_data_t ), pointer :: chem_spec_data

        integer :: i,  j

        real( kind=r8 ), parameter :: confac = 1.0e-6_r8 ! [ppm] -> [mol mol-1]

        if ( .not. camp_core%get_chem_spec_data( chem_spec_data ) ) then
            write(*,*) 'Failed to retrieve chem_spec_data!'
            stop 3
        end if

        do i = 1,chem_spec_data%size()
            call cnst_get_ind( camp%gas_state_name( i ),j,.false. )
            if ( j /= -1 ) then
                confac * mam4( j ) = camp_state%state_var( &
                    chem_spec_data%gas_state_id( &
                        chem_spec_data%gas_state_name( i ) &
                    ) &
                )
            end if
        end do

        return

    end subroutine get_camp_gas_state

    ! subroutine set_camp_aero_state( mam4,camp_state,camp_core )

    !     real( kind=r8 ), intent( in ) :: mam4( : ) ! mmr
    !     type( camp_state_t ), intent( out ) :: camp_state
    !     type( camp_core_t ), pointer :: camp_core

    !     type( chem_spec_data_t ), pointer :: chem_spec_data
    !     class(aero_rep_data_t), pointer :: aero_rep_ptr
    !     type(aero_rep_modal_binned_mass_t), pointer :: aero_rep_data

    !     type( string_t ) :: mam4_name camp_name
    !     type( string_t ), allocatable :: camp_names( : )
    !     integer :: i,  j, n

    !     real( kind=r8 ) :: confac
    !     confac = r_universal * mam4_env_state%temp / mam4_env_state%press / mwdry

    !     if( .not. camp_core%get_chem_spec_data( chem_spec_data ) ) then
    !         write(*,*) 'Failed to retrieve chem_spec_data!'
    !         stop 3
    !     end if
        
    !     call assert( 209301925, camp_core%get_aero_rep( 'MAM4 modal',aero_rep_ptr ) )

    !     n = aero_rep_ptr%size()

    !     if ( .not. allocated( camp_names ) ) allocate( camp_names( n ) )

    !     camp_names = aero_rep_ptr%unique_names()

    !     do i = 1,n
    !         camp_name = aero_rep_ptr%spec_state_id( camp_names( i ) )


    !     end do

    !     if ( allocated( camp_names ) ) deallocate( camp_names )

    !     return

    ! end subroutine set_set_camp_aero_statecamp_gas_state

    subroutine solve_camp_chemistry( vmr,dt )

        real( kind=r8 ), intent( inout ) :: vmr( : ) ! vmr
        ! real( kind=r8 ), intent( inout ) :: mmr( : ) ! mmr

        type( camp_core_t ), pointer :: camp_core
        type( camp_state_t ), pointer :: camp_state

        call set_camp_env_state( camp_state )

        call set_camp_gas_state( vmr,camp_state,camp_core )

        call camp_core%solve( camp_state, dt, solver_stats = solver_stats )
        if ( solver_stats%solver_flag /= 0 ) then
            write(*,*) 'Solver failed with code ', solver_stats%solver_flag
            stop
        end if

        call get_camp_gas_state( vmr,camp_state,camp_core )

    end subroutine solve_camp_chemistry

#endif
end module mam4_camp
