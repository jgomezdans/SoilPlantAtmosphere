! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module soil_functions

  !! SPA_DALEC version, incorporates calculation of field capacity   !!
  !! Crank Nicholson model to solve temperature flux in soil profile !!
  !! new latent energy flux model (based on SWP of surface layer)    !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: soil_processes

  ! physical constants used in module..
  ! (if want to duplicate them elsewhere then perhaps declare higher up)
  real,parameter :: albedo    = 0.2,      & ! proportion of SW reflected by surface (%)
                    cpair     = 1012.,    & ! J kg-1 K-1
                    emiss     = 0.96,     & ! emissivity
                    grav      = 9.8067,   & ! acceleration due to gravity, m s-1
                    vonkarman = 0.41e0,   & ! von Karman's constant
                    Vw        = 18.05e-6    ! partial molal volume of water, m3 mol-1 at 20C

  save

contains
  !
  !-------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_processes( time , totestevap )

    ! > subroutine summary? < !

    use gv_clim,                only: atmos_press, ppt, temp_bot, temp_top, wdbot, wetev, wind_spd
    use gv_hourscale,           only: dayevap, freeze, evap_store, hour, hourppt, hourpress, hourrad, hourrnet, hourtemp, hourtime,&
                                      hourts, hourvpd, hourwind, overflow, Qc, Qe, Qh, Qn, surface_watermm, totet, underflow
    use gv_irradiance_sunshade, only: soilnet
    use gv_scale_declarations,  only: grid, time_holder
    use gv_snow_info,           only: snow_watermm
    use gv_soil_structure,      only: conduc, pptgain, soil_temp, soil_temp_nplus1, thickness, waterfrac, watergain, watericemm, &
                                      waterloss
    use gv_veg,                 only: ess, lai, modrnet, soiletmm
    use math_tools,             only: zbrent
    use soil_air,               only: slayer, soil_conductivity, soil_porosity, soil_resistance, soil_water_potential, &
                                      water_uptake_layer
    use spa_io,                 only: handle_output

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    real, intent(out)   :: totestevap

    ! local variables..
    integer :: i , p
    real    :: checkdiff, checkerror, cw, delw, fluxsum, lambda, pw, t1, &
               t2, top4, ts, waterchange, xacc

    !--- main calculations begin below ---

    pw         = sum(watericemm)
    xacc       = 0.0001
    Qh         = 0.
    Qe         = 0.
    Qn         = 0.
    Qc         = 0.
    hour       = time%step
    i          = time%step
    hourtemp   = temp_bot + freeze      !convert from oC to K
    hourvpd    = 0.1 * wdbot * hourtemp / 217.
    hourwind   = wind_spd
    hourrad    = soilnet
    hourpress  = atmos_press
    hourtime   = time%daytime
    hourppt    = ppt
    hourrnet   = modrnet
    surface_watermm = 0.      ! assume that water infiltrates within a timestep
    lambda     = 1000. * ( 2501.0 - 2.364 * ( hourtemp - 273.15 ) )  ! latent heat of vapourisation, J kg-1
    waterloss  = 0
    watergain  = 0
    pptgain    = 0

    ! Re-calculate soil water..
    do slayer = 1 , grid%core   ! loop over layers..
      conduc(slayer) = soil_conductivity( waterfrac(slayer) )
    enddo

    call soil_porosity
    call soil_water_potential
    call soil_resistance

    if ( ( ppt .gt. 0. ) .and. ( temp_top .lt. 0. ) ) then    ! add snowpack
       snow_watermm = snow_watermm + ppt
       ppt   = 0.
    end if
    if ( (snow_watermm .gt. 0. ) .and. ( temp_top .gt. 1 ) ) then    ! melt snow
       snow_watermm    = snow_watermm * 0.9    ! decay rate
       surface_watermm = snow_watermm * 0.1    ! how much water reaches the soil surface from melting snow?
    end if
    if ( snow_watermm .lt. 0.01 ) snow_watermm = 0.        ! remove tiny snow

    if ( temp_top .gt. 0 ) then    ! it is not snowing    
       call canopy_balance        ! interception of rain and canopy evap
       wetev(i) = evap_store     ! mm t-1
    endif

    t1 = hourtemp - 50.
    t2 = hourtemp + 50.   ! set temp bounds based on air temp
    TS = zbrent( 'soil_processes:energy' , energy , t1 , t2 , xacc )    ! calculate surface temp based on energy balance

    if ( sum( lai ) .gt. 0. ) then  ! leaves on trees so transpiration is on
      call water_uptake_layer( totestevap )
    endif

    hourts = ts
    call water_fluxes( time )   ! gains and losses
    call water_thermal          ! apply thermal corrections and implement water fluxes

    soil_temp_nplus1       = 0
    soil_temp(1)           = ts
    soil_temp_nplus1(1)    = ts
    soil_temp_nplus1(grid%core) = soil_temp(grid%core)

    call crank_nicholson        ! soil temperature profile
    call thaw_depth(ts)         ! thaw depths and icefractions

    ess(i)        = dayevap          ! (Wm-2)
    lambda        = 1000. * ( 2501.0 - 2.364 * ( ts - 273.15 ) )      ! latent heat of vapourisation (J kg-1)
    soiletmm( i ) = ess( i ) / lambda * time%seconds_per_step  ! convert to mm t-1


    top4 = 0. !water stored in top 4 layers
    do p = 1 , 4
      top4 = top4 + waterfrac( p ) * thickness( p )
    end do

    waterchange = 1e3 * ( sum(pptgain) + sum(watergain) - watergain(grid%core) ) &
                 - 1e3 * ( sum(waterloss) - waterloss(grid%core) )
    fluxsum    = soiletmm(i) + 1e3 * ( overflow + totet + underflow ) - surface_watermm
    checkerror = waterchange + fluxsum
    cw         = sum( watericemm )
    delw       = pw - cw
    checkdiff  = delw - fluxsum 

    ! write output..
    call handle_output( 5 , output_data=(/ ts, waterchange, &
                           fluxsum, checkerror, delw, checkdiff /) )

    if ( time%step .eq. time%steps_per_day )  &
        call handle_output( 6 , output_data=(/ ts, waterchange, &
                             fluxsum, checkerror, delw, checkdiff /) )

  end subroutine soil_processes
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine canopy_balance

    ! > subroutine summary? < !

    use gv_hourscale,      only: canopy_store, surface_watermm
    use gv_soil_structure, only: max_storage
    use math_tools,        only: dxsav, kmax, ode_int

    implicit none

    ! local variables..
    integer,parameter :: nvar = 2
    integer           :: nbad, nok
    real              :: eps, h1, hmin, x1, x2, ystart(nvar)

    eps  = 1.0e-4
    h1   = .001
    hmin = 0.0
    kmax = 100
    x1   = 1.
    x2   = 2.
    dxsav = ( x2 - x1 ) / 20.0
    ! initial conditions
    if ( canopy_store .lt. 1e-5 * max_storage ) canopy_store = 0.   ! empty store if it's tiny
    ystart( 1 ) = canopy_store  
    ystart( 2 ) = surface_watermm

    call ode_int( 'canopy_balance:canopy_water_store' , ystart , nvar , x1 , &
                   x2 , eps , h1 , hmin , nok , nbad , canopy_water_store )

    canopy_store    = ystart( 1 )
    surface_watermm = ystart( 2 )

  end subroutine canopy_balance
  !
  !----------------------------------------------------------------------
  !
  subroutine canopy_water_store( time_dummy, y , dydt )

    ! determines canopy water storage and evaporation, !
    ! and water reaching soil surface.                 !

    use gv_hourscale,          only: evap_store, hourppt
    use gv_scale_declarations, only: nmax, time
    use gv_soil_structure,     only: max_storage, through_fall
    use log_tools

    implicit none

    ! arguments..
    real,intent(in)  :: y(nmax)
    real,intent(in)  :: time_dummy ! dummy argument, provided for ode_int
    real,intent(out) :: dydt(nmax)

    ! local variables..
    real :: a, add_ground, add_store, b, drain_store, potential_evap, ratio

    add_store  = ( 1. - through_fall ) * hourppt ! rate of input of water to canopy storage per min
    add_ground = through_fall * hourppt          ! rate of input of water to ground
    ! function gives potential evaporation rate..
    potential_evap = wetted_surface_evap()
    ! rate of evaporation from storage NB store cannot exceed max storage..
    ratio      = min( 1. , y( 1 ) / max_storage )
    evap_store = potential_evap * ratio
    b = 3.7
    a = log( 0.002 ) - b * max_storage
    if ( y( 1 ) .gt. max_storage ) then
      ! rate of drainage from store (mm per timestep (in units of minutes))
      drain_store = exp( a + b * y(1) ) * ( time%seconds_per_step / 60. )
    else
      drain_store = 0.
    endif

    if ( ( y(2) .gt. 0. ) .and. ( y(2) .lt. 7e-37 ) ) then
      write(message,*) 'problem in Canopy_water_store: y(2) (',y(2),') is too small.'
      call write_log( message , msg_warning , __FILE__ , __LINE__ )
    endif

    dydt(1) = add_store - drain_store - evap_store  ! change in canopy storage
    dydt(2) = add_ground + drain_store              ! addition to soilwater
    if ( dydt(1) .gt. 100. ) then
      write(message,*) 'Problem in Canopy_water_store: rate of change of y (',dydt(1), &
                       ') is too large to be physically plausible.'
      call write_log( message  , msg_warning , __FILE__ , __LINE__ )
    endif

  end subroutine canopy_water_store
  !
  !----------------------------------------------------------------------
  !
  subroutine crank_nicholson()
 
    ! finite difference PDE solver for soil temperature profile !
 
    use gv_scale_declarations, only: grid, time
    use gv_soil_structure,     only: soil_temp, soil_temp_nplus1, thermal, thickness

    implicit none

    ! local variables..
    integer :: i
    real    :: beta, D, error, max_error, old_value, tdiffuse

    ! --calculations begin below--
    MAX_ERROR = 0.0000005
    beta      = 0.5

    ! Continue looping until error is sufficiently reduced
    ! (exit condition at end of loop)
    do
      error = 0.  ! reset error

      ! (Loop over all x-dimension nodes, except first and last)
      do i = 2 , grid%core - 1 

        thermal = thermal_conductivity( i )

        ! Walbroeck thermal conductivity, W m-1 K-1 is converted to J m-1 K-1 timestep-1...
        tdiffuse = time%seconds_per_step * thermal / heatcap( i )

        D = tdiffuse / ( thickness( i ) * thickness(i) )

        old_value = soil_temp_nplus1( i )          ! Store value of previous iteration

        ! Calculate the temperature at the new time step using an implicit method..
        soil_temp_nplus1( i ) = ( D / ( 1 + 2 * BETA * D ) ) &
                             * ( BETA * ( soil_temp_nplus1( i + 1 ) + soil_temp_nplus1( i - 1 ) ) &
                             +  ( 1 - BETA ) * ( soil_temp( i + 1 ) - 2 * soil_temp( i ) + soil_temp( i - 1 ) ) ) &
                             +  soil_temp( i ) / ( 1 + 2 * BETA * D )

        error = error + abs( old_value - soil_temp_nplus1( i ) )  ! Calculate the error

      enddo

      ! exit if total error sufficiently small..
      if ( error .le. max_error ) exit

    enddo

    ! Set the values at time n equal to the values at time n+1 for the next time step..

    ! (Loop over all x-dimension nodes, except first and last)
    do i = 2 , grid%core - 1  
      soil_temp( i ) = soil_temp_nplus1( i )
    enddo

  end subroutine crank_nicholson
  !
  !----------------------------------------------------------------------
  !
  real function energy( ts )

    ! Determines the surface energy balance. !
    ! (see Hinzmann et al. 1998)             !

    use gv_hourscale,          only: freeze, hourrad, hourtemp, Qc, Qe, Qh, Qn
    use gv_scale_declarations, only: boltz
    use gv_soil_structure,     only: soil_temp, thickness, waterfrac

    implicit none

    ! arguments..
    real,intent(in) :: ts

    ! local variables..
    real ::  downwelling_rad, gah, rho, upwelling_rad

    ! Sensible heat flux..
    rho = 353.0 / hourtemp                ! density of air kg m-3 (t-dependent)
    call exchange_coefficient( gah )      ! conductance to heat, m s-1

    Qh  = cpair * rho * gah * ( hourtemp - ts )

    ! Latent energy flux..
    ! (no evaporation if surface is frozen or dried-out)
    if ( ( soil_temp(1) .le. freeze ) .or. ( waterfrac(1) .le. 0. ) ) then
      Qe = 0.
    else
      Qe = qe_flux( ts )
    endif

    ! Net radiation (emitted LW varies with surface temp)..
    ! The brackets stop the enormous ts value from swamping
    ! the much smaller emissivity and boltzman values.
    upwelling_rad = ( 1d0 * emiss * boltz ) * ts**4
    downwelling_rad =  hourrad
    Qn = downwelling_rad - upwelling_rad

    ! (Hillel) Thermal conductivity of top layer..
    Qc = -thermal_conductivity(1) * ( ts - soil_temp(2) ) / ( 0.5 * thickness(1) )

    ! Energy balance..
    energy = Qh + Qe + Qn + Qc

  end function energy
  !
  !----------------------------------------------------------------------
  !
  subroutine exchange_coefficient( exco_out )

    ! heat or vapour exchange coefficient/boundary layer !
    ! cond, m s-1.   See Mat William's ref 1028          !

    use gv_hourscale, only: hourwind
    use gv_veg,       only: canopy_height, tower_height

    implicit none

    ! arguments..
    real, intent(out) :: exco_out

    ! local variables..
    real :: answer, log_frac, numerator

    ! Boundary layer conductance at ground level.
    ! Substitute lower altitude than canopy_height
    ! Heat exchange coefficient from Hinzmann
    ! 0.13 * canopy_height gives roughness length

    numerator = 1d0 * hourwind * ( 1d0 * vonkarman )**2
    log_frac  = log( tower_height / ( 0.13d0 * canopy_height ) )
    answer    = ( 1d0 * numerator ) / ( 1d0 * log_frac )**2
    exco_out  = 1d0 * answer

  end subroutine exchange_coefficient
  !
  !----------------------------------------------------------------------
  !
  real function heatcap(i)

    ! Walbroeck: impacts of freeze/thaw on energy fluxes !
    ! heat capacity, J m-3 K-1                           !

    use gv_soil_structure, only: iceprop, mineralfrac, organicfrac, soil_temp, thickness, waterfrac

    implicit none

    ! arguments..
    integer,intent(in) :: i

    ! local variables..
    real    :: delt, lhf, lw, volhc

    lhf  = 334000.     ! latent heat of fusion of ice, J kg-1
    delt = 1.0         ! temperature range over which freezing occurs

    volhc = 2e6 * mineralfrac(i)                               &
           + 2.5e6 * organicfrac(i)                            &
            + 4.2e6 * ( waterfrac(i) * ( 1. - iceprop(i) ) )   &
             + 1.9e6 * waterfrac(i) * iceprop(i)
              ! (J m-3 K-1), Hillel 1980 p. 294

    ! soil_temp is in kelvin, so convert to celcius and test if it lies in freezing range (-1<x<0)..
    if ( ( (soil_temp(i)-273.15) .le. 0. ) .and. ( (soil_temp(i)-273.15) .gt. -delt ) ) then
       lw      = 1000. * waterfrac(i) * ( thickness(i) / 0.1 )       ! liquid water content (kg m-3 soil)
       heatcap = volhc + lhf * lw / delt
    else
       heatcap = volhc
    endif

  end function heatcap
  !
  !----------------------------------------------------------------------
  !
  subroutine infiltrate

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    use gv_hourscale,          only: runoff, surface_watermm, overflow
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: porosity, pptgain, thickness, waterfrac, watergain, waterloss

    implicit none

    integer :: i
    real    :: add, wdiff

    add = surface_watermm * 0.001
    pptgain = 0.
    do i = 1 , grid%soil
       wdiff=max(0.,(porosity(i)-waterfrac(i))*thickness(i)-watergain(i)+waterloss(i))
       if(add.gt.wdiff)then
          pptgain(i)=wdiff
          add=add-wdiff
       else
          pptgain(i)=add
          add=0.
       endif
       if(add.le.0.)exit
    enddo
    if(add.gt.0.)then
       overflow=add
    else
       overflow=0.
    endif
    runoff=runoff+overflow

  end subroutine infiltrate
  !
  !----------------------------------------------------------------------
  !
  real function qe_flux( ts )

    ! latent energy loss from soil surface !

    use gv_hourscale,      only: freeze, gaw, gws, hourpress, hourtemp, hourvpd
    use gv_meteo,          only: Rcon
    use gv_soil_structure, only: drythick, porosity, SWP

    implicit none

    ! arguments..
    real, intent(in) :: ts

    ! local variables..
    real :: diff, ea, esat, esurf, lambda, por, rho, tort

    tort = 2.5              ! tortuosity
    por  = porosity(1)      ! porosity of surface layer

    lambda = 1000. * ( 2501.0 - 2.364 * ( ts - 273.15 ) )  ! latent heat of vapourisation (J kg-1)
    call exchange_coefficient( gaw )                       ! determine boundary layer conductance (m s-1)
    rho  = 353.0 / hourtemp                                ! density of air (kg m-3) (t-dependent)

    diff = 24.2e-6 * ( ts / 293.2 )**1.75                  ! (m2 s-1) diffusion coefficient for water

    ! saturation vapour pressure of air (kPa)...
    esat = 0.1 * exp( 1.80956664 + ( 17.2693882 * hourtemp - 4717.306081 ) / ( hourtemp - 35.86 ) )
    ea   = esat - hourvpd         ! vapour pressure of air
    ! saturation vapour pressure (kPA) at surface - assume saturation...
    esat  = 0.1 * exp( 1.80956664 + ( 17.2693882 * ts - 4717.306081 ) / ( ts - 35.86 ) )
    ! vapour pressure in soil airspace (kPa), dependent on soil water potential - Jones p.110. Vw=partial molal volume of water...
    esurf = esat * exp( 1e6 * SWP(1) * Vw / ( Rcon * Ts ) )
    ! soil conductance to water vapour diffusion (m s-1)...
    gws = por * diff / ( tort * drythick )

    qe_flux = lambda * rho * 0.622 / ( 0.001 * hourpress ) * ( ea - esurf ) / ( 1. / gaw + 1. / gws )

    if ( ts .lt. freeze ) qe_flux = 0.  ! no evaporation if surface is frozen

  end function qe_flux
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_balance( soil_layer )

    ! integrator for soil gravitational drainage !

    use gv_soil_structure, only: field_capacity, iceprop, porosity, thickness, waterfrac, watergain, waterloss
    use log_tools
    use math_tools,        only: dxsav, kmax, ode_int
    use soil_air,          only: drainlayer, liquid, slayer, soilpor, unsat

    implicit none

    ! arguments..
    integer,intent(in) :: soil_layer

    ! local variables..
    integer,parameter :: nvar = 1
    integer           :: nbad, nok
    real              :: change, eps, h1, hmin, newwf, x1, x2, ystart(nvar)

    ! --calculations begin below--
    eps        = 1.0e-4
    h1         = .001
    hmin       = 0.0
    kmax       = 100
    x1         = 1.
    x2         = 2.
    dxsav      = ( x2 - x1 ) / 20.0
    soilpor    = porosity( soil_layer )
    liquid     = waterfrac( soil_layer ) * ( 1. - iceprop( soil_layer ) )     ! liquid fraction
    drainlayer = field_capacity( soil_layer )
    ! unsaturated volume of layer below (m3 m-2)..
    unsat      = max( 0. , ( porosity( soil_layer+1 ) - waterfrac( soil_layer+1 ) ) &
                          * thickness( soil_layer+1 ) / thickness( soil_layer )     )
    slayer     = soil_layer

    ! initial conditions
    if ( ( liquid .gt. 0. ) .and. ( waterfrac( soil_layer ) .gt. drainlayer ) ) then
      ! there is liquid water..
      ystart(1) = waterfrac( soil_layer )    ! total layer   
      call ode_int( 'soil_balance:soil_water_store' , ystart , nvar , x1 , &
                     x2 , eps , h1 , hmin , nok , nbad , soil_water_store )
      newwf  = ystart(1)
      ! convert from waterfraction to absolute amount..
      change = ( waterfrac( soil_layer ) - newwf ) * thickness( soil_layer )
      watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
      waterloss(    soil_layer  ) = waterloss(   soil_layer   ) + change
    endif

    if ( waterloss( soil_layer ) .lt. 0. ) then
       call write_log( 'waterloss probem in soil_balance' , msg_error , __FILE__ , __LINE__ )
    endif

  end subroutine soil_balance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_store( time_dummy , y , dydt )

    ! determines gravitational water drainage !

    use gv_scale_declarations, only: nmax, time
    use soil_air,              only: drainlayer, liquid, soil_conductivity, unsat

    implicit none

    ! arguments..
    real,intent(in)  :: y(nmax)
    real,intent(in)  :: time_dummy ! dummy argument, provided for ode_int
    real,intent(out) :: dydt(nmax)

    ! local variables..
    real    :: drainage

    drainage = soil_conductivity( y(1) ) * time%seconds_per_step

    if ( y(1) .le. drainlayer ) then  ! gravitational drainage above field_capacity 
      drainage = 0.
    endif
    if ( drainage .gt. liquid ) then  ! ice does not drain
       drainage = liquid
    endif
    if ( drainage .gt. unsat ) then   ! layer below cannot accept more water than unsat
       drainage = unsat
    endif

    dydt(1) = -drainage               ! waterloss from this layer

  end subroutine soil_water_store
  !
  !----------------------------------------------------------------------
  !
  subroutine thaw_depth(ts)

    ! determines layer ice fraction, and depth of thaw !

    use gv_hourscale,          only: freeze
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: iceprop, soil_temp, thickness

    implicit none

    ! arguments..
    real,intent(in) :: ts

    ! local variables..
    integer :: k, liquid(0:grid%soil), numthaw
    real    :: depthtotop(grid%soil), middepth(grid%soil), &
               root, split, thaw(grid%soil), topt, topthick

    thaw       = -9999.
    numthaw    = 1          ! There may be two or more thaw points, so thaw is an array
    depthtotop = 0.         ! Records depth to start of each layer
    iceprop    = 0.
    if ( ts .gt. freeze ) then
      liquid(0) = 1
    else
      liquid(0) = 0
    endif
    do k = 1 , grid%soil   ! Check if layer centres are frozen
      if ( soil_temp(k) .gt. freeze ) then
        liquid(k) = 1
      else
        liquid(k) = 0
      endif
    enddo
    do k = 1 , grid%soil    ! Locate thaw depth
      if ( k .eq. 1 ) then  ! For layer 1, upper boundary is surface temperature
        topthick      = 0.
        middepth(k)   = 0.5 * thickness(k)
        depthtotop(k) = 0.
      else
        topthick        = thickness( k-1 )
        depthtotop( k ) = sum( depthtotop) + thickness( k-1 ) ! Increment
        middepth( k )   = depthtotop( k ) + 0.5 * thickness( k )
      endif
      if ( liquid( k-1 ) .ne. liquid( k ) ) then
        ! Change of state - locate thaw
        if ( k .eq. 1 ) then  !for layer 1, upper boundary is surface temperature
          topt = ts
        else
          topt = soil_temp( k-1 )
        endif
        ! Location of thaw relative to layer midpoint..
        root = ( freeze - soil_temp(k) ) * ( 0.5 * thickness(k) + 0.5 * topthick ) &
                     / ( topt -soil_temp(k) )
        thaw( numthaw ) = middepth( k ) - root  ! Determine thaw depth

        ! Now ice fractions
        if ( thaw( numthaw ) .gt. depthtotop( k ) ) then 
          ! Thaw is in layer k
          ! fraction of top half of layer thawed/unthawed..
          split = ( thaw( numthaw ) - depthtotop( k ) ) / ( 0.5 * thickness( k ) )
          if ( soil_temp(k) .lt. freeze ) then
            iceprop( k ) = iceprop( k ) + 0.5 * ( 1. - split )
          else
            iceprop( k ) = iceprop( k ) + 0.5 * split
            if ( k .gt. 1 ) iceprop( k-1 ) = iceprop( k-1 ) + 0.5  ! bottom half of k-1 is frozen
          endif
        else
          ! Thaw is in layer k-1
          ! fraction of bottom half of layer-1 thawed/unthawed..
          split = ( depthtotop( k ) - thaw( numthaw ) ) / ( 0.5 * thickness( k - 1 ) )
          if ( soil_temp( k-1 ) .lt. freeze ) then
            iceprop( k-1 ) = iceprop( k-1 ) + 0.5 * ( 1. - split )
          else
            iceprop( k-1 ) = iceprop( k-1 ) + 0.5 * split
            iceprop( k   ) = iceprop(  k  ) + 0.5  ! top half of layer k is frozen
          endif
        endif
        numthaw = numthaw + 1       ! next thaw has separate storage location
      else
        ! No change of state
        if ( liquid( k-1 ) + liquid( k ) .eq. 2. ) then
          ! Both water..
          iceprop( k ) = iceprop( k )
          if ( k .gt. 1 ) iceprop( k-1 ) = iceprop( k-1 )
        else
          ! Both ice..
          iceprop( k ) = iceprop( k ) + 0.5
          if ( k .gt. 1 ) iceprop( k-1 ) = iceprop( k-1 ) + 0.5
        endif
      endif
    enddo

  end subroutine thaw_depth
  !
  !----------------------------------------------------------------------
  !
  real function thermal_conductivity( i )

    ! thermal conductivity W m-1 K-1  !
    ! Hillel p.295                    !

    use gv_soil_structure, only: iceprop

    implicit none

    ! arguments..
    integer,intent(in) :: i

    if ( i .lt. 4. ) then
      thermal_conductivity = 0.34 + ( 0.58 - 0.34 ) * iceprop(i)
    else if ( ( i .ge. 4. ) .and. ( i .lt. 10. ) ) then
      thermal_conductivity = 1.5  + ( 2.05 -  1.5 ) * iceprop(i)
    else
      thermal_conductivity = 1.69 + ( 3.63 - 1.69 ) * iceprop(i)
    endif

  end function thermal_conductivity
  !
  !----------------------------------------------------------------------
  !
  subroutine water_fluxes( time )

    ! waterfrac is m3 m-3, soilwp is MPa !

    use gv_hourscale,          only: freeze, hourts, Qe, surface_watermm, totet, unintercepted
    use gv_scale_declarations, only: grid, time_holder
    use gv_soil_structure,     only: abovebelow, drythick, fraction_uptake, rooted_layers, snow, thickness, watergain, waterloss
    use log_tools
    use soil_air,              only: soil_resistance, soil_water_potential

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    integer            :: i, rr
    real               :: lambda

    snow    = 0
    lambda  = 1000. * ( 2501.0 - 2.364 * ( hourts - freeze ) )  ! latent heat of vapourisation (J kg-1)
    call wetting_layers( time )
    ! From which layer is evap water withdrawn?..
    if ( drythick .lt. thickness( 1 ) ) then
       rr = 1     ! The dry zone does not extend beneath the top layer
    else
       rr = 2     ! The dry zone does extend beneath the top layer
    endif
    if ( Qe .lt. 0. ) then
      ! Evaporation..
      waterloss( rr ) = waterloss( rr ) + 0.001 * ( -Qe / lambda ) * time%seconds_per_step  ! (t m-2 t-1, m t-1)
    else
      ! Dew formation..
      watergain( 1 )  = watergain( 1 )  + 0.001 * (  Qe / lambda ) * time%seconds_per_step  ! (t m-2 t-1, m t-1)
    endif

    totet = max( 0. , totet )
    do i = 1 , rooted_layers            !water loss from each layer 
      waterloss( i ) = waterloss( i ) + totet * fraction_uptake( i ) * abovebelow
    enddo

    do i = 1 , grid%soil
      call soil_balance( i )
      if ( waterloss( i ) .lt. 0. ) then
         write(message,*)'trouble with water loss = ',waterloss(i),' in layer ',i
         call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
      endif
    enddo
    call infiltrate

    ! Find SWP & soil resistance without updating waterfrac yet (do that in waterthermal)
    call soil_water_potential
    call soil_resistance

    if ( time%step .eq. 1 ) then
      unintercepted = surface_watermm
    else
      unintercepted = unintercepted + surface_watermm  ! determine total water reaching ground for day summation
    endif

  end subroutine water_fluxes
  !
  !----------------------------------------------------------------------
  !
  subroutine water_thermal

    ! redistribute heat according to !
    ! water movement in the soil     !

    use gv_hourscale,          only: dayevap, discharge, hourtemp, hourts, Qe, underflow
    use gv_scale_declarations, only: grid, time
    use gv_soil_structure,     only: iceprop, mineralfrac, organicfrac, pptgain, soil_temp, thickness, &
                                     waterfrac, watergain, watericemm, waterloss
    use log_tools

    implicit none

    ! local variables..
    integer :: i
    real    :: evap, heat, heatgain, heatloss, icecontent, lambda, newheat, volhc, watercontent

    do i = 1 , grid%soil    
      volhc = 2e6 * mineralfrac(i) + 2.5e6 * organicfrac(i) + &
               4.2e6 * ( waterfrac(i) * ( 1. - iceprop(i) ) ) + &
                1.9e6 * waterfrac( i ) * iceprop( i )    ! volumetric heat capacity of soil layer (J m-3 K-1)
      heat           = soil_temp(i)*volhc*thickness(i)        ! heat content of the layer (J m-2)
      ! water loss in m * heat capacity of water * temp => J m-2
      heatloss       = waterloss( i ) * 4.2e6 * soil_temp( i )
      ! water gain in m * heat capacity of water * temp => J m-2
      heatgain       = watergain( i ) * 4.2e6 * soil_temp( i ) + pptgain( i ) * 4.2e6 * hourtemp
      watercontent   = ( waterfrac(i) * ( 1. - iceprop(i) ) ) * thickness(i)    ! Liquid t or m of water
      icecontent     = ( waterfrac(i) *     iceprop(i)      ) * thickness(i)    ! Solid t or m of water
      ! now shift water content..
      watercontent   = max( 0. , watercontent + watergain(i) + pptgain( i ) - waterloss( i ) ) 
      waterfrac( i ) = ( watercontent + icecontent ) / thickness( i )
      if ( ( watercontent + icecontent ) .eq. 0 ) then
        iceprop(i) = 0.
        call write_log("water_thermal: iceprop is zero!" , msg_warning , __FILE__ , __LINE__ )
      else
        iceprop(i) = icecontent / ( watercontent + icecontent )
      endif
      if ( heatgain + heatloss .ne. 0. )then        ! net redistribution of heat
        newheat     = heat - heatloss + heatgain
        volhc       = 2e6 * mineralfrac(i) + 2.5e6 * organicfrac(i)   &
                     + 4.2e6 * ( waterfrac(i) * ( 1. - iceprop(i) ) ) &
                      + 1.9e6 * waterfrac(i) * iceprop(i)    ! volumetric heat capacity of soil layer, J m-3 K-1
        soil_temp(i) = newheat / ( volhc * thickness(i) )
      endif
      watericemm(i) = 1e3 * waterfrac(i) * thickness(i)    !mm of water in layer
    enddo
    lambda = 1000. * ( 2501.0 - 2.364 * ( hourts - 273.15 ) )    ! latent heat of vapourisation, J kg-1
    evap   = 0.001 * ( -Qe / lambda) * time%seconds_per_step    !m t-1

    dayevap   = -Qe    ! w m-2
    underflow = watergain( grid%core )
    discharge = discharge + watergain( grid%core ) * 1e3

  end subroutine water_thermal
  !
  !--------------------------------------------------------------------
  !
  real function wetted_surface_evap()

    ! evaporation from wetted surfaces (mm t-1) !

    use gv_hourscale,          only: freeze, hourrnet, hourtemp, hourwind, hourvpd
    use gv_scale_declarations, only: time
    use gv_veg,                only: canopy_height
    use log_tools

    implicit none

    ! local variables..
    real    :: d, htemp, lambda, psych, ra, rho, s, slope, z, z0

    ! boundary layer constants..
    d  = 0.75 * canopy_height
    z0 = 0.1 * canopy_height
    z  = canopy_height + 2.
    htemp  = hourtemp - freeze ! convert Kelvin to Celcius

    rho    = 353.0 / hourtemp                          ! density of air (g m-3) (t-dependent)
    ra     = ( 1. / vonkarman**2 * hourwind ) * ( log( z - d ) / z0 )**2 ! aerodynamic resistance
    s      = 6.1078 * 17.269 * 237.3 * exp( 17.269 * htemp / ( 237.3 + htemp ) )
    slope  = 100. * ( s / ( 237.3 + htemp )**2 )       ! slope of saturation vapour pressure curve (t-dependent)
    lambda = 1000. * ( 2501.0 - 2.364 * htemp )        ! latent heat of vapourisation
    psych  = 100. * ( 0.646 * exp( 0.00097 * htemp ) ) ! psychrometer constant
    if ( slope * hourrnet .gt. rho * cpair * hourvpd / ra ) then
      ! Penman-Monteith equation (kg m-2 s-1)
      wetted_surface_evap = ( slope * hourrnet + rho * cpair * hourvpd / ra ) / ( lambda * ( slope + psych ) )
    else
      ! Actually dewfall occurs
      wetted_surface_evap = 0.
    endif

    wetted_surface_evap = wetted_surface_evap * time%seconds_per_step ! convert from kg m-2 s-1 to mm t-1

    if ( wetted_surface_evap .lt. -1000 ) then
      call write_log( 'Problem in wetted_surface_evap' , msg_warning , __FILE__ , __LINE__ )
    endif

  end function wetted_surface_evap
  !
  !----------------------------------------------------------------------
  !
  subroutine wetting_layers( time )

    ! surface wetting and drying determines !
    ! thickness of dry layer and thus Qe    !

    use gv_hourscale,          only: freeze, hourtemp, Qe, surface_watermm
    use gv_scale_declarations, only: time_holder
    use gv_soil_structure,     only: drythick, porosity, snow, thickness, wettingbot, wettingtop
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    integer :: ar1(1), ar2, ar1b
    real    :: airspace, diff, dmin, lambda, netc

    dmin = 0.001
    ! latent heat of vapourisation (J kg-1)..
    lambda  = 1000. * ( 2501.0 - 2.364 * ( hourtemp - freeze ) ) 

    airspace = porosity( 1 )
    ! Soil LE should be withdrawn from the wetting layer with the smallest depth..
    ar1  = minloc( wettingbot , MASK = wettingbot .gt. 0. )
    ar1b = sum( ar1 )   ! convert ar1 to scalar, and make sure always greater than or equal to 1
    if ( ar1b .eq. 0. ) then
      ! There is no water left in soils!!
      call write_log( "There is no water left in any soil layer!" , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Calulate the net change in wetting in the top zone..
    netc = ( 0.001 * Qe / lambda * time%seconds_per_step ) / airspace + &
              ( surface_watermm * 0.001 + snow ) / airspace   !m
    if ( netc .gt. 0. ) then      ! Wetting
      ! resaturate the layer if top is dry and recharge is greater than drythick
      if ( ( netc .gt. wettingtop( ar1b ) ) .and. ( wettingtop( ar1b ) .gt. 0. ) ) then
        diff = netc - wettingtop( ar1b )  ! extra water to deepen wetting layer
        wettingtop( ar1b ) = 0.
        if ( ar1b .gt. 1 ) then
          ! Not in primary layer (primary layer can't extend deeper)..
          wettingbot( ar1b ) = wettingbot( ar1b ) + diff
        endif
        drythick = dmin
      else
        if ( wettingtop( ar1b ) .eq. 0. ) then
          ! surface is already wet, so extend depth of this wet zone..
          if ( ar1b .gt. 1 ) then
            ! not in primary layer (primary layer can't extend deeper)..
            wettingbot( ar1b ) = wettingbot( ar1b ) + netc          
            if ( wettingbot(ar1b) .ge. wettingtop(ar1b-1) ) then
              ! Layers are conterminous..
              wettingtop( ar1b - 1 ) = wettingtop(ar1b)
              wettingtop(   ar1b   ) = 0.     ! remove layer
              wettingbot(   ar1b   ) = 0.
            endif
          endif
        else
          ! Create a new wetting zone..
          wettingtop( ar1b + 1 ) = 0.
          wettingbot( ar1b + 1 ) = netc
        endif
        drythick = dmin
      endif
    else    ! drying
      wettingtop( ar1b ) = wettingtop( ar1b ) - netc         ! Drying increases the wettingtop depth
      if ( wettingtop( ar1b ) .ge. wettingbot( ar1b ) ) then ! Wetting layer is dried out.
        diff = wettingtop( ar1b ) - wettingbot( ar1b )       ! How much more drying is there?                 
        wettingtop( ar1b ) = 0.
        wettingbot( ar1b ) = 0.
        ar2 = ar1b - 1
        if ( ar2 .gt. 0 ) then  ! Move to deeper wetting layer
          wettingtop( ar2 ) = wettingtop( ar2 ) + diff    !dry out deeper layer
          drythick = max( dmin , wettingtop( ar2 ) )
        else    !no deeper layer
          drythick = thickness( 1 )   !layer 1 is dry
        endif
      else
        drythick = max( dmin , wettingtop( ar1b ) )
      endif
    endif
    if ( drythick .eq. 0. ) then
      call write_log( 'Problem in drythick' , msg_warning , __FILE__ , __LINE__ )
    endif

  end subroutine wetting_layers
  !
  !----------------------------------------------------------------------
  !
end module soil_functions
!
!------------------------------------------------------------------------
!
