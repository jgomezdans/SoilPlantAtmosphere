! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_initialise

  !! This module gives initial values to variables at the !!
  !!  start of a 'from-scratch' SPA simulation.           !!
  !! The alternative is to start SPA from some previous   !!
  !!  state, which does not require these steps.          !!

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: initialise_soils, initialise_veg

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first...
  !
  !----------------------------------------------------------------------
  !
  subroutine initialise_soils( initialwater )

    ! Go through the steps necessary to initialise the soils from scratch !

    use gv_hourscale,          only: hourts
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: soil_temp, thickness, waterfrac, watericemm, wettingbot, wettingtop
    use soil_air,              only: saxton_parameters, soil_porosity, water_retention

    implicit none

    ! arguments..
    real,intent(out) :: initialwater

    ! local variables..
    integer          :: i

    call saxton_parameters ! calculate the conductivity parameters
    call water_retention

    watericemm    = 0.
    wettingbot    = 0.
    wettingtop    = 0.
    hourts        = soil_temp(1) ! initial soil temp estimate for light routine
    wettingbot(1) = thickness(1) ! top layer is saturated

    call soil_porosity      ! calculate the soil porosity

    ! v large mesophyll conductance - ACi curves have not been Epron adjusted
    initialwater  = 0.
    do i = 1 , grid%soil
      initialwater = initialwater + 1e3 * ( waterfrac(i) * thickness(i) )
    enddo

  end subroutine initialise_soils
  !
  !----------------------------------------------------------------------
  !
  subroutine initialise_veg( )

    ! Calculate the initial leaf-water potential, which requires !
    ! knowing the soil water potential, which in turn requires   !
    ! knowing about the conductance and resistance of the soil.  !

    use gv_meteo,              only: head
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: conduc, rooted_layers, root_length, root_mass, waterfrac, weighted_SWP
    use gv_veg,                only: canopy_height, lafrac, lai, layer_height, LMA, LWPstore, stock_foliage
    use soil_air,              only: slayer, soil_conductivity, soil_resistance, &
                                     soil_water_potential, water_uptake_layer

    implicit none

    ! local variables..
    real    :: dummy_output
    integer :: i

    canopy_height = layer_height(1)

    ! Calculate soil-conductance. This requires the saxton parameters
    ! to have been calculated already (which happens in init_soils)
    do slayer = 1 , grid%core   ! loop over layers..
      conduc(slayer) = soil_conductivity( waterfrac(slayer) )
    enddo

    ! Need to provide some (made-up!) initial values..
    root_length   = 0.1
    root_mass     = 0.1
    rooted_layers = grid%core
    if ( stock_foliage .eq. 0. ) stock_foliage = 0.1
    lai           = lafrac * stock_foliage / LMA

    ! Calculate resistance to flow (we need soilR)
    call soil_resistance()
    ! and the soil water potential..
    call soil_water_potential()
    ! ..use to calculate the weighted soil-water-potential..
    call water_uptake_layer( dummy_output )
    ! And finally used to get the initial leaf-water potential..
    do i = 1 , grid%canopy     
       LWPstore(i) = weighted_SWP - head * layer_height(i)
    enddo
 
  end subroutine initialise_veg
  !
  !----------------------------------------------------------------------
  !
end module spa_initialise
!
!------------------------------------------------------------------------
!
