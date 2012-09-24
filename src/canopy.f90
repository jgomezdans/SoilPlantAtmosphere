! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module canopy

  !! >this module requires a summary< !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: timestep_calcs

contains
  !
  !----------------------------------------------------------------------
  !
  subroutine timestep_calcs( time )

    ! CANOPY models (previously called day_sunshade_soilwater) !
    ! runs first by timestep, and within each step then runs   !
    ! through each layer.                                      !
    ! (calculate daily carbon gain for one layer)              !

    use allocate_carbon,       only: alloc_carbon_to_pools
    use gv_clim,               only: coa, gbw
    use gv_hourscale,          only: totet
    use gv_meteo,              only: la, psil, psis, temp
    use gv_scale_declarations, only: grid, time_holder
    use gv_soil_structure,     only: weighted_SWP
    use gv_veg,                only: avN, co2amb, GPP, gppt, lafrac, lai, lma, LWPstore, nfrac, nla, respt, &
                                     stock_foliage, totevap, totn, transt
    use light,                 only: leaf_sun, solar
    use soil_functions,        only: soil_processes

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    integer :: i
    real    :: frac_shade,  & ! part of leaf not hit by sunlight
               frac_sun,    & ! part of leaf hit by sunlight
               psil_store,  & ! temporary storage of psil value
               la_sun,      & ! area of leaf sunlit
               la_shade,    & ! area of leaf shaded
               lambda,      & !
               shade_psil,  & ! psil value for shaded leaf
               sun_psil,    & ! psil value for sunlit leaf
               tot_est_evap   ! total estimated evaporation

    if ( time%step .eq. 1 ) then
      gppt = 0.  ;  respt = 0.  ;  transt = 0.
    endif

    lai  = lafrac * stock_foliage / LMA   ! determine LAI distribution from Cf
    totn = avN * sum(lai)      ! total N
    Nla  = nfrac * totn        ! N per layer

    call boundary()

    call solar( time )

    call soil_processes( time , tot_est_evap )

    psis = weighted_SWP        ! weighted soil water potential

    do i = 1 , grid%canopy

       psil       = LWPstore(i) ! load layer's LWP from previous timestep
       co2amb     = coa         ! load current CO2 concentration
       psil_store = psil        ! store LWP 

       ! assign some local vars for convenience..
       frac_sun   = leaf_sun(i)
       la_sun     = frac_sun * lai(i)
       frac_shade = ( 1. - leaf_sun(i) )
       la_shade   = frac_shade * lai(i)  
    
       !-- first do sunlit leaf area --
       la = la_sun
       ! if there is active leaf area..
       if ( ( la .gt. 0. ) .and. ( tot_est_evap .gt. 0. ) )  &
            call leaf_processes( i , frac_sun , la , sunlit=.true. )

       sun_psil = psil   ! record LWP in sunlit leaves

       !-- now do shaded leaf area --
       psil = psil_store    ! reset LWP for shade leaf calculation
       la   = la_shade
       if ( ( la .gt. 0. ) .and. ( tot_est_evap .gt. 0. ) )  &
            call leaf_processes( i , frac_shade , la , sunlit=.false. )

       shade_psil  = psil  ! record LWP in shaded leaves

       psil        = frac_sun * sun_psil + frac_shade * shade_psil  ! calculate final LWP
       LWPstore(i) = psil                                           ! update store of LWP for next timestep

    enddo

    ! Calculate the total evapotranspiration during this time period..
    lambda  = 1000. * ( 2501.0 - 2.364 * temp )                            ! (J kg-1) latent heat of vapourisation
    totet   = 0.001 * time%seconds_per_step * transt(time%step) / lambda   ! (m3 m-2 t-1, m t-1)
    totevap = totevap + 1000. * totet                                      ! (mm)

    GPP = gppt(time%step) * time%seconds_per_step * 1e-6 * 12.  ! (g C) assimilated per time step

    call alloc_carbon_to_pools

  end subroutine timestep_calcs
  !
  !----------------------------------------------------------------------
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !----------------------------------------------------------------------
  !
  !
  !----------------------------------------------------------------------
  !
  subroutine boundary()

    ! determine boundary layer conductances at each canopy layer !

    use gv_clim,               only: atmos_press, gbw, temp_bot, temp_top, wind_spd
    use gv_scale_declarations, only: grid
    use gv_veg,                only: dimen, layer_height, tower_height

    implicit none

    ! local variables..
    integer :: i
    real    :: alpha, Dwv, mult, roughl, store, tempx, thick, u(grid%canopy), xx

    ! initialise..
    alpha  = 4.
    roughl = 0.075 * tower_height  ! tower height is canopy_height+4
    gbw    = 0.
    xx     = 0.
    store  = 1.0
 
    ! calculate..
    do i = 1 , grid%canopy
      xx    = xx + 1.
      tempx = temp_top - 0.1111 * ( xx - 1 ) * ( temp_top - temp_bot )
      Dwv   = .0000242 * ( ( ( tempx + 273.15 ) / 293.15 )**1.75 ) * 101300. / atmos_press    !Jones p 51
      ! windspeed in layer - no decay in windspeed in this open stand
      !             u(i,t)=windsp(t)
      mult   = exp( alpha * ( layer_height( i ) / tower_height - 1. ) )
      u( i ) = wind_spd * mult
      thick  = 0.004 * ( dimen / u( i ) )**0.5    ! boundary layer thickness
      ! conductance to water vapour (m s-1 - not mol m-2 s-1 as in Amthor p.333)
      ! i.e. remove P/RT
      gbw( i ) = Dwv / thick
    enddo

  end subroutine boundary
  !
  !----------------------------------------------------------------------
  !
  subroutine leaf_processes( i , leaf_fraction , leaf_area , sunlit )

    use gv_clim,               only: gbw, par_ratio, temp_bot, temp_top, wdbot, wdtop
    use gv_meteo,              only: gbb, la, nit, par, rad, temp, wdef
    use gv_scale_declarations, only: grid, time
    use gv_veg,                only: gppt, lai, Nla, respt, transt
    use leaf,                  only: assimilate
    use light,                 only: checkpar, longem, nirabs, parabs_shade, parabs_sun

    implicit none

    ! arguments..
    integer,intent(in) :: i  ! index
    real,intent(in)    :: leaf_area,     & !
                          leaf_fraction    ! fraction of leaf exposed (to sun/shade, depending on call)
    logical,intent(in) :: sunlit

    ! local variables..
    real    :: agr,   & !
               gsm,   & !
               modet, & !
               res      !

    nit = leaf_fraction * Nla(i)

    !  leaf area has all beam rad plus frac of diffuse rad..
    if ( sunlit ) then
      par = ( parabs_sun(i) + parabs_shade(i) * leaf_fraction ) / leaf_area
    else
      par = parabs_shade(i) * leaf_fraction / leaf_area
    endif

    ! net radiation = shortwave + longwave radiation
    rad = 0.001 * ( nirabs(i) * leaf_fraction / leaf_area     &
                     + longem(i) * leaf_fraction / leaf_area  &
                       + par / par_ratio )

    temp = temp_top - 0.1111 * ( i - 1 ) * ( temp_top - temp_bot )

    wdef = wdtop - 0.1111 * ( i - 1 ) * ( wdtop - wdbot )

    gbb  = gbw(i)                    ! boundary layer

    call set_leaf( i , leaf_fraction )    ! see below

    call assimilate( i , time , modet , agr , res , gsm )    ! see leaf.f90

    gppt(time%step)   = gppt(time%step)   + agr     ! umol CO2 assimilation m-2 ground area s-1

    respt(time%step)  = respt(time%step)  + res     ! umol CO2 respiration m-2 ground area s-1

    transt(time%step) = transt(time%step) + modet   ! evapotranspiration in W m-2 ground area

    checkpar(i) = checkpar(i) + par * la

  end subroutine leaf_processes
  !
  !----------------------------------------------------------------------
  !
  subroutine set_leaf( clayer , frac )

    ! sets Farquhar parameters and hydraulic parameters for each leaf layer

    use gv_metab, only: ht, layer_capac, rn, rplant, rsoil, vcm, vjm
    use gv_meteo, only: la, nit
    use gv_veg,   only: canopy_soil_resistance, capac, conductivity, gplant, kappac, kappaj, layer_height

    implicit none

    ! arguments..
    integer, intent(in) :: clayer
    real,    intent(in) :: frac

    ! local variables..
    real, parameter :: propresp = 1.0 ! proportional respiration constant based on N content
                                      !  and assumed temperature base of 10 oC.

    vcm = kappac * nit / la
    ! metabolic rates are umol/m2/s - so we assume here that nit is N/m2 (it's actually N/clayer)
    vjm = kappaj * nit / la
    rn = 0.105 * propresp  ! respiration constant umol CO2/g N at 10 deg C
    rn = rn / la           ! convert to resp per m2

    ht = layer_height( clayer )

    ! plant hydraulic resistances are determined by the amount of leaf area in the sunlit or shaded fraction
    if ( conductivity .eq. 1 ) then
       rplant = ht / ( gplant * la )  ! MPa s mmol-1 (per layer)
    else
       rplant = 1. / ( gplant * la )  ! conductance is constant with height
    endif
    layer_capac = capac * la

    rsoil = canopy_soil_resistance( clayer )  ! soil+root resistance of each canopy layer

    ! now correct rsoil according to the fraction of total leaf area in this run
    rsoil = rsoil / frac

  end subroutine set_leaf
  !
  !----------------------------------------------------------------------
  !
end module canopy
!
!------------------------------------------------------------------------
!
