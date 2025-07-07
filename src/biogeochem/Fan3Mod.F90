module Fan3Mod

!-------------------------------------------------------------------------------
! !DESCRIPTION:
!
! This module implements the physical parameterizations of the FANv3 (Flow of
! Agricultureal Nitrogen version 3) process model, it includes the numerical 
! solutions of Nitrification and Denitrification in the first soil layer in CLM.
! The final target of FANv3 is to evaluate the NOx and N2O emission from the 
! agricultural soil and connect the nitrogen pool of FAN with CLM. 
!
!                                         Jinmu Luo, Nov 11 2022, Ithaca
!-------------------------------------------------------------------------------

use decompMod                        , only : bounds_type
use shr_kind_mod                     , only : r8 => shr_kind_r8
use shr_infnan_mod                   , only : isnan => shr_infnan_isnan
use shr_infnan_mod                   , only : isinf => shr_infnan_isinf
use clm_varctl                       , only : iulog
use spmdMod                          , only : masterproc
use abortutils                       , only : endrun
use clm_time_manager                 , only : get_step_size_real
use clm_time_manager                 , only : get_curr_date
use clm_time_manager                 , only : get_nstep_since_startup_or_lastDA_restart_or_pause
use clm_varcon                       , only : nitrif_n2o_loss_frac

!
implicit none
private

public :: eval_nitrification_n2o_nox
public :: eval_denitrification_n2o_nox
public :: eval_nitrate_loss
public :: eval_nh4_upward
public :: update_nitrate_pool
public :: eval_cr

! indexes for array of different flux (used in Fan2CTSM.F90)
integer, parameter, public :: n2o_n_id    = 1
integer, parameter, public :: nox_n_id    = 2
integer, parameter, public :: n2o_dn_id   = 3
integer, parameter, public :: nox_dn_id   = 4
integer, parameter, public :: n2_dn_id    = 5
integer, parameter, public :: no3flux_id  = 6
integer, parameter, public :: runoff_id   = 7
integer, parameter, public :: atm_id      = 8
integer, parameter, public :: deepsoil_id = 9

integer, parameter, public :: num_fan3fluxes  = 9

! Check button used in this module
logical, parameter, private :: v2log = .FALSE.

! Canopy reduction switch
logical, parameter, public :: use_canopy_reduction = .TRUE.
logical, parameter, public :: use_upward_diffusion = .TRUE.
!
real(r8), private, parameter :: skip_size = 3600.0_r8   ! Time steps to skip the balance check at startup (sec)
integer,  private            :: skip_steps = -999       ! Number of time steps to skip the balance check for

contains


!-------------------------------------------------------------------------------
! NO3 SOIL CONCENTRATION, SURFACE CONCENTRATION, RESISTENCE: SCSRFR
!-------------------------------------------------------------------------------
function eval_Raq(theta, thetasat, soil_depth) result(Raq)
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil)   
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil)
  real(r8), intent(in) :: soil_depth                       ! soil layer depth in CLM, 0.02(upward), 0.03(downward), meter
  real(r8), parameter :: Daq = 1.7e-9_r8                   ! Molecular diffusivity of NO3- in water, m2/s
  real(r8) :: tortuosity_water                             ! Tortuosity for aqueous-phase diffusion, unitless
  real(r8) :: Raq                                          ! resistances between soil layer, s/m 

  ! notice that julius used power of 7/3 here. 
  tortuosity_water = theta**(10_r8/3_r8) / (thetasat**2_r8)
  Raq = soil_depth / (2_r8 * tortuosity_water * Daq)
end function eval_Raq

function eval_no3_aqsoil(N_no3, theta) result(no3_aqsoil)
  real(r8), intent(in) :: N_no3(:)                         ! nitrate pools(age) of fanv3, gN/m2
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil) 
  real(r8), parameter :: soil_depth = 0.02_r8              ! first soil layer depth in CLM, 0.02 meter
  real(r8) :: no3_aqsoil(size(N_no3))                      ! no3 aqueous concentration in soil, gN/m3 in water

  no3_aqsoil = N_no3 / (theta * soil_depth)

  ! If column in extremely dry conditions
  where(isinf(no3_aqsoil)) no3_aqsoil = 0.0_r8
end function eval_no3_aqsoil

function eval_no3_aqsrf(N_no3, qr, theta, thetasat) result(no3_aqsrf)
  real(r8), intent(in) :: N_no3(:)                         ! nitrate pools(age) of fanv3, gN/m2
  real(r8), intent(in) :: qr                               ! surface runoff from CLM, m/s
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil)   
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil)
  real(r8), parameter :: soil_depth = 0.02_r8              ! first layer depth of soil in CLM, meter  
  real(r8) :: Raq                                          ! resistances between soil layer, s/m 
  real(r8) :: no3_aqsoil(size(N_no3))                      ! no3 aqueous concentration in soil, gN/m3 in water
  real(r8) :: no3_aqsrf(size(N_no3))                       ! no3 aqueous concentration at surface, gN/m3 in water

  Raq = eval_Raq(theta, thetasat, soil_depth)
  no3_aqsoil = eval_no3_aqsoil(N_no3, theta)
  no3_aqsrf = no3_aqsoil / (1_r8 + qr*Raq)
end function eval_no3_aqsrf

function eval_sgd(theta, thetasat) result(sg_diff)
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil)
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil)
  real(r8) :: sg_diff                                      ! soil gas diffusivity, unitless
 
  ! Huang et al.,(2015)
  ! sg_diff = 0.209*(1 - theta/thetasat)**(4/3)

  ! Julius et al.,(2020)
  sg_diff = (abs(theta-thetasat))**(10_r8/3_r8)/(thetasat**2)  
  
end function eval_sgd

! ---------------------- SCSRFR FUNCION END HERE -------------------------------


! ------------------------------------------------------------------------------
! NITRIFICATION PRODUCTION OF NOx and N2O: NP
! ------------------------------------------------------------------------------
subroutine eval_nitrification_n2o_nox(nitrify_flux, theta, thetasat, sg_diff, & 
           no3_flux, n2o_n, nox_n, RNOx)

  real(r8), intent(in) :: nitrify_flux(:)                  ! nitrification fluxes (age) from fanv2, gN/m2/s
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil) 
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil)
  real(r8), intent(in) :: sg_diff                          ! soil gas diffusivity from CLM, unitless
 
  real(r8), intent(out) :: no3_flux(:)                     ! production of NO3, gN/m2/s
  real(r8), intent(out) :: n2o_n(:)                        ! n2o production rate, gN/m2/s
  real(r8), intent(out) :: nox_n(:)                        ! nox production rate, gN/m2/s
  real(r8), intent(out) :: RNOx                            ! ratio of no/n2o (Parton et al.,1996, 2001), unitless

  real(r8) :: sg_diff_fan, wfps, nitrif_lost_as_n2o 
  real(r8), parameter :: pi = 3.1415926_r8      


  sg_diff_fan = eval_sgd(theta, thetasat)
 
  RNOx = 15.2_r8 + 35.5_r8*atan(0.68_r8*pi*(10.0_r8*sg_diff - 1.86_r8)) / pi
  wfps = max(min(theta/thetasat, 1._r8), 0._r8)
  nitrif_lost_as_n2o = 0.0006 + (0.01 - 0.0006)/( 1 + exp( -( -6.27 + 24.71*wfps)) ) 
  nitrif_lost_as_n2o = min(max(0.0006, nitrif_lost_as_n2o), 0.01) 
  no3_flux = nitrify_flux / (1_r8 + nitrif_lost_as_n2o + RNOx*nitrif_lost_as_n2o)
  n2o_n = nitrif_lost_as_n2o * no3_flux
  nox_n = RNOx * n2o_n

end subroutine eval_nitrification_n2o_nox

! ----------------------- NP FUNCION END HERE ----------------------------------


! ------------------------------------------------------------------------------
! DENITRIFICATION PRODUCTION OF NOx AND N2O : DNP
! Question remain: FNO3 and ratio_NO3_CO2 are evaluated by fertilizer type, not 
! total NO3 pool, it could make N2O and NOx emission smaller?
! ------------------------------------------------------------------------------
subroutine eval_denitrification_n2o_nox(N_no3, theta, thetasat, RNOx, FCO2, & 
           anaerobic_frac, sg_diff, soil_co2, n2o_dn, nox_dn, n2_dn)

  real(r8), intent(in) :: N_no3(:)                         ! nitrate pools(age) of fanv3, gN/m2
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil) 
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil)
  real(r8), intent(in) :: RNOx                             ! ratio of no/n2o (Parton et al.,1996, 2001), unitless
 
  real(r8), intent(in) :: FCO2                             ! maximum potential dn rates based on CO2, ugC/gsoil/day
  real(r8), intent(in) :: anaerobic_frac                   ! anoxic fraction of soils, unit?
  real(r8), intent(in) :: sg_diff                          ! soil gas diffusivity from CLM, unitless
  real(r8), intent(in) :: soil_co2                         ! soil co2 concentration, ugC/g(soil)/day
  
  real(r8) :: FNO3                                         ! maximum potential dn rates based on NO3-, ugN/gsoil/day
  real(r8) :: ratio_NO3_CO2                                ! ratio of no3 and co2 concentration in soil, unitless
  real(r8) :: fr_NC                                        ! function of the ration of no3 and co2, unitless
  real(r8) :: Rn2n2o                                       ! n2/n2o, unitless
  real(r8) :: k1                                           
  real(r8) :: soil_dens                                    ! bulk density of soil, kg/m3
  real(r8) :: no3_aqsoil(size(N_no3))                      ! no3- aqueous concentration in soil, gN/m3 in water
  real(r8) :: no3_total                                    ! summation of age pool, gN/m3 in water
 
  real(r8), intent(out) :: n2o_dn(:)                       ! n2o production rate, gN/m2/s
  real(r8), intent(out) :: nox_dn(:)                       ! nox production rate, gN/m2/s
  real(r8), intent(out) :: n2_dn(:)                        ! n2 production rate, gN/m2/s

  ! !LOCAL VARIABLES:
  real(r8) :: sg_diff_fan                                  ! soil gas diffusivity from FAN, unitless 
  real(r8) :: wfps                                         ! water-filled pore space, unitless
  real(r8) :: fr_WFPS                                      ! total water limitation function, unit?
  real(r8), parameter :: soil_part_dens = 2650.0_r8        ! soil particle density, kg/m3

  no3_aqsoil = eval_no3_aqsoil(N_no3, theta)
  no3_total = sum(no3_aqsoil)
  soil_dens = soil_part_dens * (1.0_r8 - thetasat)
 
  ! [gN/m3] to [ugN/gSoil], evaluate the total emission
  FNO3 = 1.15_r8 * (no3_total * 1.e3_r8 / soil_dens)**0.57_r8
 
  ! ratio of concentration of NO3 and CO2, unit in [ugN/gSoil / ugC/gSoil/day] 
  if ( soil_co2 > 1.0e9_r8 ) then 
     ratio_NO3_CO2 = (no3_total * 1.e3_r8 / soil_dens)/soil_co2
  else
     ratio_NO3_CO2 = 100._r8
  endif

  sg_diff_fan = eval_sgd(theta, thetasat)
  
  k1 = max(1.7_r8, 38.4_r8 - 350_r8*sg_diff)  
  fr_NC = max(k1*0.16_r8, k1*exp(-0.8_r8*ratio_NO3_CO2))

  wfps= max(min(theta/thetasat, 1._r8), 0._r8) * 100._r8
  fr_WFPS = max(0.1_r8, 0.015_r8 * wfps - 0.32_r8)

  Rn2n2o = fr_NC * fr_WFPS

  ! divide the total emission to age pool proportionally 
  if (no3_total /= 0.0_r8) then 
     n2o_dn = min(FNO3, FCO2) * max(min(anaerobic_frac, 1.0_r8), 0.0_r8) / & 
              (1_r8+Rn2n2o+RNOx) * no3_aqsoil/no3_total
  else
     n2o_dn = min(FNO3, FCO2) * max(min(anaerobic_frac, 1.0_r8), 0.0_r8) / & 
              (1_r8+Rn2n2o+RNOx) * 1/size(N_no3) 
  end if

  ! ugN/gsoil/day to gN/m2/sec
  n2o_dn = n2o_dn/1000/3600/24 * soil_dens * 0.02_r8  
  nox_dn = RNOx * n2o_dn
  n2_dn = Rn2n2o * n2o_dn

  ! Make sure Nan and infinite values are avoided here
  where(isnan(n2o_dn) .or. isinf(n2o_dn)) n2o_dn = 0.0_r8
  where(isnan(nox_dn) .or. isinf(nox_dn)) nox_dn = 0.0_r8
  where(isnan(n2_dn) .or. isinf(n2_dn)) n2_dn = 0.0_r8

end subroutine eval_denitrification_n2o_nox

! ---------------------- DNP FUNCION END HERE ----------------------------------



! ------------------------------------------------------------------------------
! NITRITE LOASS IN FANv3 : NLOSS
! ------------------------------------------------------------------------------
! Surface runoff, Diffusion and Flux to Atmosphere
! Runoff is a liner funcion, evaluate it in total or different has same results
! N loss should not bigger than the pool

subroutine eval_nitrate_loss(N_no3, qr, theta, thetasat, dt, n2o_dn, nox_dn, & 
                             n2_dn, no3_aqdeep, runoff, atm, deepsoil)

  real(r8), intent(in) :: N_no3(:)                         ! nitrate pools(age) of fanv3, gN/m2
  real(r8), intent(in) :: qr                               ! surface runoff from CLM, m/s
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil)   
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil
  real(r8), intent(in) :: dt                               ! time step of the model, second
  real(r8), intent(inout) :: n2o_dn(:)                     ! n2o production rate, gN/m2/s
  real(r8), intent(inout) :: nox_dn(:)                     ! nox production rate, gN/m2/s 
  real(r8), intent(inout) :: n2_dn(:)                      ! n2 production rate, gN/m2/s
  real(r8), intent(in) :: no3_aqdeep                       ! no3 concentration from deep layer of CLM, gN/m3 in water

  real(r8), parameter :: soil_depth = 0.03                 ! soil layer depth in CLM, 0.02(upward), 0.03(downward), meter

  real(r8) :: no3_aqsoil(size(N_no3))                      ! no3- aqueous concentration in soil, gN/m3 in water
  real(r8) :: no3_total                                    ! summation of age pool, gN/m3 in water
  real(r8) :: no3_aqsrf(size(N_no3))                       ! no3 aqueous concentration at surface, gN/m3 in water
  real(r8) :: Raq                                          ! resistances between soil layer, s/m
  real(r8) :: df(size(N_no3))

  real(r8), intent(out) :: runoff(:)                       ! surface runoff, gN/m2/s in water
  real(r8), intent(out) :: atm(:)                          ! N flux to the atmosphere from denitrification, gN/m2/s
  real(r8), intent(out) :: deepsoil(:)                     ! downward(positve)/upward(negative) diffusion, gN/m2/s
 
  ! temporary variables
  real(r8) :: all_flux(size(N_no3)), gas_flux(size(N_no3)) 
 
  no3_aqsoil = eval_no3_aqsoil(N_no3, theta)
  no3_total  = sum(no3_aqsoil) 
  no3_aqsrf = eval_no3_aqsrf(N_no3, qr, theta, thetasat) 
  Raq = eval_Raq(theta, thetasat, soil_depth)

  runoff = max(0.0_r8, qr * no3_aqsrf)
  atm = max(0.0_r8, n2_dn + n2o_dn + nox_dn)
 
  ! diffusion is evaluated in total then divide into age pool
  ! For now, we only evaluate the downward diffusion of nitrate
  if (no3_total /= 0.0_r8) then
     deepsoil = (no3_total - 0.0_r8)/Raq * no3_aqsoil/no3_total
  else
     deepsoil = (no3_total - 0.0_r8)/Raq * 1/size(N_no3)
  end if 

  ! Blance check, if necessary, rescalling the flux
  df = N_no3 - (runoff + atm + deepsoil)*dt 
  
  if (any(df < 0.0_r8)) then
     df = df/dt
     !write(iulog,*) "Output fluxes are larger the pool reservations! rescalling!"
     if (any(deepsoil < 0.0_r8)) then
        ! pump the nitrate from clm, rescale the emission and runoff
        where (df > 0.0_r8) df = 0.0_r8
        all_flux = runoff + atm
        where (all_flux == 0.0_r8) all_flux = 1e-25_r8
        gas_flux = n2_dn + n2o_dn + nox_dn
        runoff = max(0.0_r8, runoff + df*runoff/all_flux)
        where (gas_flux == 0.0_r8) gas_flux = 1e-25_r8
        nox_dn = max(0.0_r8, nox_dn + df*atm/all_flux * nox_dn/gas_flux)
        n2o_dn = max(0.0_r8, n2o_dn + df*atm/all_flux * n2o_dn/gas_flux)
        n2_dn = max(0.0_r8, n2_dn + df*atm/all_flux * n2_dn/gas_flux)
     else
        ! rescale all positive flux
        where (df > 0.0_r8) df = 0.0_r8
        all_flux = runoff + atm + deepsoil
        gas_flux = n2_dn + n2o_dn + nox_dn
        where (all_flux == 0.0_r8) all_flux = 1e-25_r8
        runoff = max(0.0_r8, runoff + df * runoff/all_flux)
        deepsoil = max(0.0_r8, deepsoil + df * deepsoil/all_flux)
        where (gas_flux == 0.0_r8) gas_flux = 1e-25_r8
        nox_dn = max(0.0_r8, nox_dn + df * atm/all_flux * nox_dn/gas_flux)
        n2o_dn = max(0.0_r8, n2o_dn + df * atm/all_flux * n2o_dn/gas_flux)
        n2_dn = max(0.0_r8, n2_dn + df * atm/all_flux * n2_dn/gas_flux)
     end if

  end if

  ! Make sure Nan and infinite values are avoided here
  where(isnan(deepsoil) .or. isinf(deepsoil)) deepsoil = 0.0_r8
  where(isnan(runoff) .or. isinf(runoff)) runoff = 0.0_r8
  where(isnan(n2o_dn) .or. isinf(n2o_dn)) n2o_dn = 0.0_r8
  where(isnan(nox_dn) .or. isinf(nox_dn)) nox_dn = 0.0_r8
  where(isnan(n2_dn) .or. isinf(n2_dn)) n2_dn = 0.0_r8
  atm = n2_dn + n2o_dn + nox_dn

end subroutine eval_nitrate_loss


! ---------------------- NLOSS FUNCION END HERE --------------------------------



! ------------------------------------------------------------------------------
! MASS BALANCE OF TOTAL NITRITE IN FANv3 : UPDATE
! ------------------------------------------------------------------------------
! FANv2    --> FANv3    : differnt age pool of ammonium
! FANv3    --> CLM, CAM : downward diffusion, flux to atmosphere
! CLM, CAM --> FANv3    : upward diffusion, (depostion in the future)

subroutine update_nitrate_pool(N_no3, I_no3, dt, runoff, atm, deepsoil, no3_flux)

  real(r8), intent(inout) :: N_no3(:)                      ! nitrate pools(age) of fanv3, gN/m2
  real(r8), intent(in) :: I_no3                            ! fertilization rate of no3, gN/m2/s
  real(r8), intent(in) :: dt                               ! time step of the model, second

  real(r8), intent(in) :: runoff(:)                        ! surface runoff, gN/m2/s in water
  real(r8), intent(in) :: atm(:)                           ! N flux to the atmosphere, gN/m2/s
  real(r8), intent(in) :: deepsoil(:)                      ! diffusion to/from deep soil, gN/m2/s
  
  real(r8), intent(in) :: no3_flux(:)                      ! recalculated no3 flux(fanv3), gN/m2/s
  
  ! Private variables
  real(r8) :: N_no3_old(size(N_no3))


  N_no3_old = N_no3
  ! N loss   
  N_no3 = N_no3 - (runoff + atm + deepsoil)*dt

  ! N gaining, add the fertilizers to the last age pool
  N_no3(size(N_no3)) = N_no3(size(N_no3)) + I_no3*dt 
  N_no3 = N_no3 + no3_flux * dt 
 
  ! Make sure nitrate concentration larger than zero
  if ( any( N_no3 < -1.00E-15_r8) ) then
     write(iulog, *), 'Negative pool!'
     write(iulog, *), N_no3_old
     write(iulog, *), N_no3
  end if 
  N_no3 = max(0.0_r8, N_no3)

end subroutine update_nitrate_pool

! ---------------------- UPDATE FUNCION END HERE -------------------------------


! ------------------------------------------------------------------------------
! CANOPY REDUCTION FUNCTION FOR NO IN FANv3 
! ------------------------------------------------------------------------------
subroutine eval_cr(bounds, filter_soilc, num_soilc, lai_patch, CR_col)
  !
  ! This subroutine calculate the above canopy reduction coefficient
  ! Is is called in CNPhenologyMod
  !
  use GridcellType         , only : grc
  use LandunitType         , only : lun
  use ColumnType           , only : col
  use PatchType            , only : patch
  use landunit_varcon      , only : istcrop, istsoil
  use abortutils           , only : endrun 
 
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: filter_soilc(:)
  integer, intent(in) :: num_soilc
  real(r8), intent(in) :: lai_patch(bounds%begp:bounds%endp)     ! Unburying one side Leaf area index, m2/m2
  real(r8), intent(inout) :: CR_col(bounds%begc:bounds%endc)     ! Canopy reduction coefficient, unitless
 
  ! local variables
  real(r8) :: SAI                                          ! Stomatal area index, m2/m2
  real(r8) :: SL_ratio, CR_patch
  integer :: lat, g, l, c, fc, p 

  ! Empirical model of NOx emission reduction by canopy  (Yienger and Levy et al.,1995) 
  do fc = 1, num_soilc
     c = filter_soilc(fc)

     CR_col(c) = 1.0_r8
     l = col%landunit(c)
     if (.not. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) cycle

     ! Agriculture plant has the same SAI/LAI ratio across the global
     CR_col(c) = 0.0_r8 
     if ( lun%itype(l) == istcrop ) then
        SL_ratio = 0.008_r8
        do p = col%patchi(c), col%patchf(c)
           SAI = SL_ratio * lai_patch(p)
           CR_patch = ( exp(-8.75*SAI) + exp(-0.24*lai_patch(p)) ) / 2
           CR_col(c) = CR_col(c) + CR_patch*patch%wtcol(p) 
        end do ! patch cycle
     end if

     ! Natural vegetation has different ratio, see Table 6 in (Yienger and Levy et al.,1995)
     if ( lun%itype(l) == istsoil ) then
        g = col%gridcell(c)
        lat = grc%latdeg(g)
        do p = col%patchi(c), col%patchf(c)
           if (abs(lat)>30) then
              select case(p)
                case(1, 2);    SL_ratio = 0.003_r8     ! Coniferous forest
                case default;  SL_ratio = 0.005_r8     ! Others
              end select
           else
              select case(p)
                case(9:11);    SL_ratio = 0.01_r8      ! Woodland
                case(3, 7, 8); SL_ratio = 0.015_r8     ! Drought-deciduous forest
                case(4, 6);    SL_ratio = 0.015_r8     ! Rain forest
                case default;  SL_ratio = 0.005_r8     ! Others
              end select
           end if

           SAI = SL_ratio * lai_patch(p)
           CR_patch = ( exp(-8.75*SAI) + exp(-0.24*lai_patch(p)) ) / 2
           CR_col(c) = CR_col(c) + CR_patch*patch%wtcol(p)
        end do ! patch cycle
     end if

     CR_col(c) = max(CR_col(c), 0.0_r8)
     CR_col(c) = min(CR_col(c), 1.0_r8)

  end do ! column cycle

end subroutine eval_cr


! ------------------ CANOPY REDUCTION FUNCTION ENDS HERE -----------------------


! ------------------------------------------------------------------------------
! AMMONIUM UPWARD DIFFUSION
! ------------------------------------------------------------------------------
subroutine eval_nh4_upward(nh4_totn, nh4_vr, dt, theta, thetasat, flux_upward)
  real(r8), intent(in) :: nh4_totn                         ! totoal nh4 concentration in fan, gN/m2
  real(r8), intent(in) :: nh4_vr                           ! inorganic nh4 concentration in CLM, gN/m2
  real(r8), intent(in) :: dt                               ! time step in clm, seconds
  real(r8), intent(in) :: theta                            ! soil water content, m3(water)/m3(soil)
  real(r8), intent(in) :: thetasat                         ! theta at saturation, m3(water)/m3(soil

  real(r8), intent(out) :: flux_upward                     ! upward diffusion flux, gN/m2/sec
  real(r8), parameter :: soil_depth = 0.03                 ! soil depth, m
  real(r8) :: Raq                                          ! resistances between soil layer, s/m
  ! some temporary variables 
  real(r8) :: nh4_up, nh4_down

  Raq = eval_Raq(theta, thetasat, soil_depth)

  ! in defalut condition, nh4_vr > nh4_totn
  flux_upward = (nh4_vr - nh4_totn)/Raq

  ! nh4 concentration at first layer should always small than in deep
  nh4_up = nh4_totn + flux_upward*dt
  nh4_down = nh4_vr - flux_upward*dt

  if (nh4_down < nh4_up) then
     flux_upward = flux_upward - (nh4_up - nh4_down)/dt
     flux_upward = max(0.0_r8, flux_upward)
  end if 
 
end subroutine eval_nh4_upward


! ------------------  AMMONIUM UPWARD DIFFUSION ENDS HERE ----------------------

 
end module Fan3Mod
