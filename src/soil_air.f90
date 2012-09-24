! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module soil_air

  !! > this module needs a summary < !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: saxton_parameters, soil_conductivity, soil_porosity, soil_resistance, &
            soil_water_potential, water_retention, water_uptake_layer
  ! Variables..
  public :: drainlayer, liquid, slayer, soilpor, unsat


  integer :: slayer,    &  ! Soil layer variable used when passing from a subroutine to a function which is acting
                           !  on a specific soil layer.
     water_retention_pass  ! Count of number of iterations gone through for soil layer in the Saxton equations

  real    :: drainlayer, & ! The field capacity of the specific layer being worked on in the s/r soil_balance
                 liquid, & ! Liquid fraction of the water present within a given soil layer, used in s/r soil_balance
                soilpor, & ! Soil porosity of a given soil layer, used in s/r soil_balance
                  unsat    ! Unsaturated volume of the soil layer below the current one being modelled (m3 m-2).
                           !  Calculated each timestep.

  save

contains
  !
  !----------------------------------------------------------------------
  !
  subroutine saxton_parameters()

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: cond1, cond2, cond3, potA, potB

    implicit none

    integer :: i
    real    :: A, B, CC, D, E, F, G, H, J, K, mult1, mult2, mult3, P, Q, R, T, U, V

    mult1 = 100.
    mult2 = 2.778e-6
    mult3 = 1000.0
    A = -4.396    ;  B = -0.0715   ; CC = -4.880e-4 ; D = -4.285e-5
    E = -3.140    ;  F = -2.22e-3  ;  G = -3.484e-5 ; H =  0.332
    J = -7.251e-4 ;  K =  0.1276   ;  P = 12.012    ; Q = -7.551e-2
    R = -3.895    ;  T =  3.671e-2 ;  U = -0.1103   ; V =  8.7546e-4

    do i = 1 , grid%core
       potA(i)  = exp( A + B * soil_frac_clay(i) + CC * soil_frac_sand(i) * soil_frac_sand(i) + &
                    D * soil_frac_sand(i) * soil_frac_sand(i) * soil_frac_clay(i) ) * 100.
       potB(i)  = E + F * soil_frac_clay(i) * soil_frac_clay(i) + G * soil_frac_sand(i) * soil_frac_sand(i) * soil_frac_clay(i)
       cond1(i) = mult2
       cond2(i) = P + Q * soil_frac_sand(i)
       cond3(i) = R + T * soil_frac_sand(i) + U * soil_frac_clay(i) + V * soil_frac_clay(i) * soil_frac_clay(i)
    enddo

  end subroutine saxton_parameters
  !
  !----------------------------------------------------------------------
  !
  real function soil_conductivity( wf )

    ! Used in the soil drainage integrator. !
    ! Returns a single-point value.         !
    ! 'slayer' is a module variable that    !
    !  provides the soil-layer number.      !

    use gv_soil_structure, only: cond1, cond2, cond3

    implicit none

    ! arguments..
    real, intent(in) :: wf ! fraction of water in soils

    if ( wf .lt. 0.05 ) then    ! Avoid floating-underflow fortran error
      soil_conductivity = 1e-30
    else
      soil_conductivity = cond1(slayer) * exp( cond2(slayer) + cond3(slayer) / wf )
      ! Soil conductivity (m s-1 )
    endif

  end function soil_conductivity
  !
  ! ---------------------------------------------------------------------
  !
  subroutine soil_porosity()

   ! Porosity is estimated from Saxton equations. !

    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: porosity

    implicit none

    ! local variables..
    integer :: i
    real    :: H, J, K

    ! saxton params relevant to porosity..
    H = 0.332  ;  J = -7.251e-4  ;  K = 0.1276

    ! loop over soil layers..
    do i = 1 , grid%core
      porosity(i) = H + J * soil_frac_sand(i) + K * log10( soil_frac_clay(i) )
    end do

  end subroutine soil_porosity
  !
  ! ---------------------------------------------------------------------
  !
  subroutine soil_resistance()

    !> subroutine description <!

    use gv_meteo,              only: head
    use gv_scale_declarations, only: pi, grid
    use gv_soil_structure,     only: abovebelow, conduc, rooted_layers, root_length, root_mass, &
                                     root_radius, rootresist, soilR, soilR1, soilR2, thickness

    implicit none

    ! local variables..
    integer :: i
    real    :: Lsoil, rs, rs2

    soilR = 0.

    ! Calculate soil-root hydraulic resistance
    do i = 1 , rooted_layers
      Lsoil = conduc(i) / head    !converts from ms-1 to m2 s-1 MPa-1
      if ( Lsoil .lt. 1e-35 ) then    !prevent floating point error
        soilR(i) = 1e35
      else 
        rs  = sqrt( 1. / (root_length(i) * pi ) )
        rs2 = log( rs / root_radius ) / ( 2.0 * pi * root_length(i) * thickness(i) * Lsoil )
        soilR1(i) = rs2 * 1E-6 * 18 * 0.001    ! convert from MPa s m2 m-3 to MPa s m2 mmol-1
        !second component of below ground resistance related to root hydraulics
        !rrcheck=rootresist/(root_mass(i)*thickness(i)/abovebelow)
        soilR2(i) = rootresist / ( root_mass(i) * thickness(i) / abovebelow )
        soilR(i)  = soilR1(i) + soilR2(i)
      endif
    enddo

  end subroutine soil_resistance
  !
  ! ---------------------------------------------------------------------
  !
  subroutine soil_water_potential()

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: potA, potB, rooted_layers, SWP, waterfrac

    implicit none

    ! local variables..
    integer :: i
    real    :: soil_WP

    do i = 1 , rooted_layers
      if ( waterfrac(i) .gt. 0. ) then
        soil_WP = -0.001 * potA(i) * waterfrac(i)**potB(i)   !  Soil water potential (MPa)
      else
        soil_WP = -9999.
      endif
      SWP(i) = soil_WP
    end do

  end subroutine soil_water_potential
  !
  !----------------------------------------------------------------------
  !
  subroutine water_retention()

    ! field capacity calculations for saxton eqns !

    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: field_capacity
    use math_tools,            only: zbrent

    implicit none

    ! local variables..
    integer        :: i
    real           :: x1, x2
    real,parameter :: xacc = 0.0001

    do i = 1 , grid%core
      x1   = 0.1    ! low guess
      x2   = 0.7    ! high guess
      water_retention_pass = i
      ! field capacity is water content at which SWP = -10 kPa
      field_capacity( i ) = zbrent( 'water_retention:water_retention_saxton_eqns' , water_retention_saxton_eqns , x1 , x2 , xacc )
    enddo

  end subroutine water_retention
  !
  !----------------------------------------------------------------------
  !
  subroutine water_uptake_layer( total_est_evap )

    ! Find from which layer water is withdrawn !

    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: fraction_uptake, iceprop, rooted_layers, soilR, SWP, weighted_SWP
    use gv_veg,                only: canopy_soil_resistance, lai, minlwp
    use log_tools

    implicit none

    ! arguments..
    real,intent(out) :: total_est_evap

    ! local variables..
    integer              :: i, j, p
    real                 :: est_evap(grid%core), frac
    integer, allocatable :: AR2(:) ! put values in array

    ! -- calculations begin below --

    j = size( shape( est_evap ) ) ! Get number of dimensions in array
    allocate( AR2(j) )            ! Allocate AR1 to number of dimensions in array

    total_est_evap = 0.  ;  weighted_SWP    = 0.
    est_evap       = 0.  ;  fraction_uptake = 0.

    do i = 1 , rooted_layers
      ! estimate max transpiration from gradient-gravity / soil resistance.
      est_evap(i) = ( SWP(i) - minlwp ) / soilR(i)
      est_evap(i) = max( 0. , est_evap(i) )       ! no negative 
      if ( iceprop(i) .gt. 0. ) est_evap(i) = 0. ! no uptake from frozen soils
    enddo
    total_est_evap = sum( est_evap )

    ! weighted soil water potential
    if ( total_est_evap .gt. 0. ) then
      ! Water was evaporated from some layers..
      do i = 1 , rooted_layers
        weighted_SWP = weighted_SWP + swp(i) * est_evap(i)
        ! fraction of total et taken from layer i...
        fraction_uptake(i) = est_evap(i) / total_est_evap
      enddo
      weighted_SWP = weighted_SWP / total_est_evap
    else
      ! No water was evaporated..
      fraction_uptake(i) = 1. / rooted_layers
    endif

    if ( nint ( sum( fraction_uptake ) ) .ne. 1. ) then
       call write_log( 'The sum of uptake fraction is not (nearly) equal to 1 '&
                     //' in water_uptake_layer' , msg_warning , __FILE__ , __LINE__ )
    endif
    if ( ( fraction_uptake(1) .gt. 1 ) .or. ( fraction_uptake(1) .lt. 0. ) ) then
       call write_log( 'Problem with the uptake fraction (either >1 or 0<)' , &
                       msg_warning , __FILE__ , __LINE__ )
    endif

    canopy_soil_resistance = 0.    ! reset
    do p = 1 , grid%canopy
      frac = lai(p) / sum(lai)
      do i = 1 , rooted_layers
        if ( frac .gt. 0. ) then 
          ! soil resistance for each canopy layer is related to leaf area
          canopy_soil_resistance(p) = canopy_soil_resistance(p) + 1. / ( soilR(i) / frac )
        else
          canopy_soil_resistance(p) = 0.001
        endif
      enddo
      canopy_soil_resistance(p) = 1. / canopy_soil_resistance(p)
    enddo

  end subroutine water_uptake_layer
  !
  !----------------------------------------------------------------------
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !----------------------------------------------------------------------
  !
  real function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    use gv_soil_structure, only: potA, potB

    implicit none

    ! arguments..
    real, intent(in) :: xin

    ! local variables..
    real ::soil_wp

    ! calculate the soil water potential (MPa)..
    soil_WP = -0.001 * potA( water_retention_pass ) * xin**potB( water_retention_pass )

    water_retention_saxton_eqns = 1000. * soil_wp + 10.    ! 10 kPa represents air-entry swp

  end function water_retention_saxton_eqns
  !
  !----------------------------------------------------------------------
  !
end module soil_air
!
!------------------------------------------------------------------------
