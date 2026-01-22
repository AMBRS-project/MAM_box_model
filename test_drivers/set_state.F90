!> \file
!> The mam4_state module.

!> The state_t structures and associated subroutines.

module mam4_state
#ifdef MAM4_USE_CAMP
      use camp_camp_core
      use camp_camp_state
      use camp_constants
      use camp_chem_spec_data
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

      type env_state_t
        !> Temperature (K).
        real(kind=r8) :: temp
        !> Ambient pressure (Pa).
        real(kind=r8) :: press
        !> Relative humidity (0-1).
        real(kind=r8) :: RH_CLEA
      end type env_state_t
      type aero_state_t
        !> Mass fractions (0-1).
        real(kind=r8) :: qso4(4)
        real(kind=r8) :: qpom(4)
        real(kind=r8) :: qsoa(4)
        real(kind=r8) :: qbc(4)
        real(kind=r8) :: qdst(4)
        real(kind=r8) :: qncl(4)
        real(kind=r8) :: qaerwat(4)
        !> Geometric mean diameter and standard deviation
        real(kind=r8) :: GMD(4), GSD(4)
        !> Aerosol densities and mass fractions
        real(kind=r8), allocatable :: dens_aer(:), mf_aer(:)
        !> Number concentration by mode
        real(kind=r8) :: numc(4)
      end type aero_state_t
      type gas_state_t
        !> Gas mixing ratios.
        real(kind=r8), allocatable :: vmr(:)
      end type gas_state_t
      integer :: mam_idx_so2g, mam_idx_h2so4g, mam_idx_soag
      real(kind=r8) :: to_kgperm3
      logical :: first_step
      !character(len=16), allocatable :: persistent_spec(:)
      real(kind=r8), allocatable :: persistent_state(:)

      contains

              subroutine env_state_set_camp_env_state(env_state, &
                                                      camp_state)

                type(env_state_t), intent(inout) :: env_state
                type(camp_state_t), intent(inout) :: camp_state

                camp_state%env_states(1)%val%rel_humid = &
                        env_state%RH_CLEA
                call camp_state%env_states(1)%set_temperature_K(env_state%temp)
                call camp_state%env_states(1)%set_pressure_Pa(env_state%press)

              end subroutine env_state_set_camp_env_state

              subroutine gas_state_set_camp_conc(camp_core, &
                                                 gas_state, &
                                                 env_state, &
                                                 camp_state)
                
                type(gas_state_t), intent(inout) :: gas_state
                type(env_state_t), intent(inout) :: env_state
                type(camp_state_t), intent(inout) :: camp_state
                type(camp_core_t), pointer :: camp_core

                integer(kind=i_kind) :: idx_SO2, &
                                        idx_SULRXN, &
                                        idx_SULF, &
                                        idx_SOAG
                integer(kind=i_kind) :: idx_M, &
                                        idx_H2, &
                                        idx_N2, &
                                        idx_O2, &
                                        idx_CH4, &
                                        idx_H2O
                real(kind=r8) :: weight(2)

                type(chem_spec_data_t), pointer :: chem_spec_data
            
                if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
                    write(*,*) "Something's gone wrong!"
                    stop 3
                end if

                idx_SO2 = chem_spec_data%gas_state_id("SO2")
                idx_SULRXN = chem_spec_data%gas_state_id("SULRXN")
                idx_SULF = chem_spec_data%gas_state_id("SULF")
                idx_SOAG = chem_spec_data%gas_state_id("SOAG")

                idx_M = chem_spec_data%gas_state_id("M")
                idx_H2 = chem_spec_data%gas_state_id("H2")
                idx_N2 = chem_spec_data%gas_state_id("N2")
                idx_O2 = chem_spec_data%gas_state_id("O2")
                idx_CH4 = chem_spec_data%gas_state_id("CH4")
                idx_H2O = chem_spec_data%gas_state_id("H2O")

                camp_state%state_var(idx_H2O) = &
                        env_state_rh_to_y(env_state)

                weight = (/ camp_state%state_var(idx_SULF), camp_state%state_var(idx_SULRXN) /) &
                       / ( camp_state%state_var(idx_SULF) + camp_state%state_var(idx_SULRXN) )
                camp_state%state_var(idx_SO2) =  gas_state%vmr(mam_idx_so2g)
                camp_state%state_var(idx_SULF) = weight(1) * gas_state%vmr(mam_idx_h2so4g)
                camp_state%state_var(idx_SULRXN) = weight(2) * gas_state%vmr(mam_idx_h2so4g)
                camp_state%state_var(idx_SOAG) = gas_state%vmr(mam_idx_soag)

                !camp_state%state_var(idx_M) =  1.0e+6_r8
                !camp_state%state_var(idx_H2) = 0.56_r8
                !camp_state%state_var(idx_N2) =  0.7808e+6_r8
                !camp_state%state_var(idx_O2) = 0.2095e+6_r8
                !camp_state%state_var(idx_CH4) = 1.85_r8
              
              end subroutine gas_state_set_camp_conc

              real(kind=r8) function env_state_rh_to_y(env_state)

                type(env_state_t), intent(inout) :: env_state

                real(kind=r8), parameter :: t_steam = 373.15_r8
                real(kind=r8) :: a, water_vp

                a = 1.0_r8 - t_steam / env_state%temp
                a = (((-0.1299_r8 * a - 0.6445_r8) * a - 1.976_r8) * a + 13.3185_r8) * a
                water_vp = 101325.0_r8 * exp(a)  ! (Pa)
                env_state_rh_to_y = env_state%RH_CLEA * water_vp * 1.0e6_r8 &
                        / env_state%press ! (ppm)

              end function env_state_rh_to_y

              subroutine gas_state_get_camp_conc(gas_state, &
                                                 camp_state, camp_core)
                
                type(gas_state_t), intent(inout) :: gas_state
                type(camp_state_t), intent(inout) :: camp_state
                type(camp_core_t), pointer :: camp_core

                integer(kind=i_kind) :: idx_SO2, &
                                        idx_SULRXN, &
                                        idx_SULF, &
                                        idx_SOAG

                type(chem_spec_data_t), pointer :: chem_spec_data
            
                if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
                    write(*,*) "Something's gone wrong!"
                    stop 3
                end if

                idx_SO2 = chem_spec_data%gas_state_id("SO2")
                idx_SULRXN = chem_spec_data%gas_state_id("SULRXN")
                idx_SULF = chem_spec_data%gas_state_id("SULF")
                idx_SOAG = chem_spec_data%gas_state_id("SOAG")

                gas_state%vmr(mam_idx_so2g) = camp_state%state_var(idx_SO2)
                gas_state%vmr(mam_idx_h2so4g) = camp_state%state_var(idx_SULRXN) &
                                              + camp_state%state_var(idx_SULF)
                gas_state%vmr(mam_idx_soag) = camp_state%state_var(idx_SOAG)
              
              end subroutine gas_state_get_camp_conc
#endif
end module mam4_state
