! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module allocate_carbon

  !! Carbon water model for data assimilation !!
  !! part of the DALEC model.                 !!
  !!                                          !!
  !! 'N' - N-limitation added to v3           !!
  !! Cfmax - max foliar carbon (LAI)          !!
  !! temperature controlled phenology         !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: alloc_carbon_to_pools, roots

  ! variables local to this module..

  real :: turnover_foliage_switch,  & ! Switch to indicate turnover of..foliage pool is on/off
          turnover_labile_switch      ! ..labile pool is on/off

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine alloc_carbon_to_pools

    ! Predictor adjusts calculations depending !
    ! on whether evergreen or deciduous run.   !

    use gv_clim,               only: avtemp, daytempsum, gdd, max_fol, mint, temp_top
    use gv_scale_declarations, only: time, user_opts
    use gv_soil_structure,     only: resp_rate_temp_coeff, root_biomass 
    use gv_veg
    use spa_io,                only: handle_output

    implicit none

    ! local variables..
    real    :: alloc_to_resp_auto,          & ! amount of carbon to allocate to autotrophic respiration pool
               npp,                         & ! net primary production
               resp_cost_labile_to_foliage, & ! respiratory cost of moving carbon..from labile to foliage pools
               resp_cost_foliage_to_labile, & ! ..from foliage to labile pools
               resp_rate,                   & ! rate of respiration at given temperature
               ts_length                      ! length of timestep in hours

    ! length of time step in hours..
    ts_length = time%seconds_per_step / 3600.

    ! Carbon fluxes..
    resp_rate = 0.5 * exp( resp_rate_temp_coeff * temp_top )

    if ( user_opts%veg_is_deciduous ) then

       if ( time%day .le. 100 ) then
          gdd     = 0.    ! growing degree day summation starts after day 100
          max_fol = 1.    ! switch to determine if max foliage has been attained
       endif
       if ( time%step .eq. 1 ) then   ! reset at start of each day
          daytempsum = 0.
          mint       = 100
       endif
       ! time switches
       mint = min( temp_top , mint )
       daytempsum = daytempsum + temp_top
       if ( time%step .eq. time%steps_per_day ) then  ! determine average temp
          avtemp = daytempsum / real(time%steps_per_day)
          gdd    = gdd + avtemp          ! growing degree day heat sum 
          call handle_output( 2 , time ) ! write this to output
       endif
       if ( gdd .lt. GDD_thresh ) then ! Winter
          turnover_foliage_switch = 1.
          turnover_labile_switch  = 0.
       else                            ! (not winter..)
          if ( max_fol .eq. 1 ) then   ! Spring
             turnover_labile_switch  = 1.
             turnover_foliage_switch = 0.
          else                         ! Summer
             turnover_labile_switch  = 0.             
             turnover_foliage_switch = 0.
          endif
       endif
       if ( (stock_foliage .ge. max_stock_foliage ) .or. ( time%day .ge. 200 ) ) then
          max_fol = 0.
          turnover_labile_switch  = 0.
       endif
       if ( ( time%day .ge. 200 ) .and. ( mint .lt. min_Temp_thresh ) ) then ! drop leaves - N hemisphere
          turnover_foliage_switch = 1.
       endif

       ! calculate losses due to respiratory cost of moving carbon..
       ! ..from the labile pool to the vegetation pool
       resp_cost_labile_to_foliage = turnover_rate_labile  * stock_labile  * resp_cost_labile_trans * turnover_labile_switch * &
                                         resp_rate * ts_length
       ! ..from foliage to the labile pool during leaf death
       resp_cost_foliage_to_labile = turnover_rate_foliage * stock_foliage * resp_cost_labile_trans * turnover_foliage_switch * &
                                         resp_rate * ts_length * ( 1. - frac_leafLoss_litter )

    endif ! ends decid/evergreen if block.

    !-- Now update all the Carbon pools --

    ! autotrophic respiration..
    resp_auto = stock_resp_auto * turnover_rate_resp_auto * ts_length

    ! calculate how much to allocate to the autotrophic-respiration pool..
    if (user_opts%veg_is_deciduous) then
       alloc_to_resp_auto = GPP * frac_GPP_resp_auto + resp_cost_labile_to_foliage + resp_cost_foliage_to_labile
    else
       alloc_to_resp_auto = GPP * frac_GPP_resp_auto
    endif

    npp = ( 1 - frac_GPP_resp_auto ) * GPP        ! npp remaining

    ! Allocate to foliage and calculate how much NPP remains after leaf growth..
    if ( user_opts%veg_is_deciduous ) then
       alloc_to_foliage = min( max_stock_foliage - stock_foliage , npp * frac_alloc_foliage ) * &
                            turnover_labile_switch + alloc_from_labile
       npp              = npp - min( max_stock_foliage - stock_foliage , npp * frac_alloc_foliage ) * turnover_labile_switch  
    else
       alloc_to_foliage = npp * frac_alloc_foliage
       npp              = npp - alloc_to_foliage
    endif

    ! allocate to roots/plant-structure..
    alloc_to_roots      = npp * frac_alloc_roots
    alloc_to_stem  = npp * ( 1 - frac_alloc_roots )

    ! calculate losses due to litter..
    if (user_opts%veg_is_deciduous) then
       litterfall_foliage = turnover_rate_foliage * stock_foliage * ts_length * frac_leafLoss_litter * turnover_foliage_switch
    else
       litterfall_foliage = turnover_rate_foliage * stock_foliage * ts_length
    endif
    litterfall_roots     = turnover_rate_roots     * stock_roots     * ts_length
    litterfall_stem = turnover_rate_stem * stock_stem * ts_length

    ! heterotrophic respiration..
    resp_h_litter        = mineralisation_rate_litter        * stock_litter        * ts_length * resp_rate
    resp_h_soilOrgMatter = mineralisation_rate_soilOrgMatter * stock_soilOrgMatter * ts_length * resp_rate

    ! decomposition..
    decomposition = decomposition_rate * stock_litter * ts_length * resp_rate

    ! allocation to/from labile..
    if (user_opts%veg_is_deciduous) then

       alloc_to_labile = ( 1. - frac_leafLoss_litter ) * turnover_rate_foliage * stock_foliage * &
                           ( 1. - resp_cost_labile_trans ) * turnover_foliage_switch * resp_rate * ts_length

       alloc_from_labile = turnover_rate_labile * stock_labile * (1. - resp_cost_labile_trans) * &
                                turnover_labile_switch * resp_rate * ts_length
    else
      alloc_to_labile   = 0.
      alloc_from_labile = 0.
    endif

    ! update the pools..
    if ( user_opts%veg_is_deciduous ) then
       ! note that alloc-from-labile is not included as its value (from the previous timestep) is already within alloc-to-foliage
       stock_foliage = stock_foliage + alloc_to_foliage - litterfall_foliage - alloc_to_labile - resp_cost_foliage_to_labile
       stock_labile  = stock_labile  + alloc_to_labile  - alloc_from_labile  - resp_cost_labile_to_foliage
    else
       stock_foliage = stock_foliage + alloc_to_foliage - litterfall_foliage
       stock_labile  = 0.
    endif
    stock_litter        = stock_litter        + litterfall_foliage + litterfall_roots     - resp_h_litter - decomposition
    stock_resp_auto     = stock_resp_auto     + alloc_to_resp_auto - resp_auto
    stock_roots         = stock_roots         + alloc_to_roots     - litterfall_roots
    stock_soilOrgMatter = stock_soilOrgMatter + decomposition      - resp_h_soilOrgMatter + litterfall_stem
    stock_stem          = stock_stem          + alloc_to_stem      - litterfall_stem

  end subroutine alloc_carbon_to_pools
  !
  !----------------------------------------------------------------------
  !
  subroutine roots

    ! determines root distribution given root_biomass !

    use gv_scale_declarations, only: pi
    use gv_soil_structure,     only: layer_depth, max_depth, root_biomass, rooted_layers, root_length, &
                                     root_mass, root_radius, root_reach, root_k, surf_biomass, thickness
    use gv_veg,                only: stock_roots
    use math_tools,            only: zbrent

    implicit none

    ! local variables..
    integer        :: i
    real           :: cumdepth, curr, depth, mult, preb, prev, root_depth, root_cross_sec_area, slpa, x1, x2, xx(1)
    real,parameter :: xacc = 0.0001, & !
              root_density = 0.5e6     ! root density (g biomass m-3 root). Set as a constant value

    depth = 0.
    preb  = 0.

    ! convert from gC to g biomass
    root_biomass = stock_roots * 2.

    ! always provide a minimum root biomass
    root_biomass = max( 5. , root_biomass )

    root_cross_sec_area = pi * root_radius**2    ! root X-sectional area (m2)

    root_depth = max_depth * root_biomass / ( root_k + root_biomass )    ! rmass=fine root C

    i = size( shape( layer_depth ) )  
    xx = minloc( layer_depth , MASK=layer_depth .gt. root_depth )
    rooted_layers = int( xx(1) )
    root_reach = layer_depth( rooted_layers )    ! to what soil layer do the roots pentrate?

    ! ensure 50% of root mass is in top 25% of rooted layers 
    ! see pipochronosequence.xls for derivation of this relationship
    mult = min( 1 / thickness(1) , max( 2.0 , 11. * exp( -0.006 * root_biomass ) ) )

    ! assume surface root density is related to total root biomass by a scalar dependent on root biomass
    surf_biomass = root_biomass * mult

    if ( rooted_layers .gt. 1 ) then
      x1 = 0.1
      x2 = 10.
      ! determine slope of root distribution given rooting depth 
      ! and ratio of root mass to surface root density
      slpa = zbrent( 'roots:root_dist' , root_dist , x1 , x2 , xacc )
      prev = 1. / slpa
      cumdepth = 0.
      do i = 1 , rooted_layers
        cumdepth       = cumdepth + thickness(i)
        curr           = 1 / slpa * exp( -slpa * cumdepth )
        ! root_mass is g biomass, i.e. twice the C content
        root_mass(i)   = ( prev - curr ) * surf_biomass
        root_length(i) = root_mass( i ) / ( root_density * root_cross_sec_area )    ! (m m-3 soil)
        prev           = curr
      enddo
    else
       root_mass( 1 ) = root_biomass
    endif

  end subroutine roots
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  real function root_dist( xin )    ! parameters

    use gv_soil_structure

    implicit none

    ! arguments..
    real, intent(in) :: xin       ! slope parameter

    ! local variables..
    real    :: one, two

    one = ( 1. - exp( -xin * rooted_layers * thickness(1) ) ) / xin
    two = root_biomass / surf_biomass

    ! see pipo chronosequence.xls for this relationship
    root_dist = ( 1. - exp( -xin * root_reach ) ) / xin - root_biomass / surf_biomass

  end function root_dist
  !
  !----------------------------------------------------------------------
  !
end module allocate_carbon
!
!------------------------------------------------------------------------
!
