! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module light

  !!> description needed <!!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: solar
  ! Variables..
  public :: checkpar, leaf_sun, longem, nirabs, parabs_shade, parabs_sun


  real :: fdiff,    & ! ?
          skyabsp,  & ! ?
          soilabsn, & ! ?
          soilabsp, & ! ?
          soilabsl    ! ?

  real,dimension(:),allocatable :: &
         checkpar,     & ! check PAR radiation balances
         leaf_sun,     & ! fraction of leaf in a given canopy layer exposed to sun
         longem,       & ! long wave emission from each canopy layer
         nirabs,       & ! NIR abosrbed by each canopy layer
         parabs_shade, & ! PAR absorbed by each canopy layer in sun exposed leaves
         parabs_sun      ! PAR absorbed by each canopy layer in shade covered leaves

  real,parameter :: nirrefl   = 0.43,  & ! ?
                    nirtrans  = 0.26,  & ! ?
                    parrefl   = 0.11,  & ! Baldocchi
                    partrans  = 0.16,  & ! ?
                    soilpar   = 0.033, & ! ?
                    soilnir   = 0.023, & ! ?
                    spherical = 2.       ! spherical leaf angle distribution has a value 2.

  save

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine solar( time )

    use gv_clim,                only: par_ratio, par_top, sw_rad, temp_bot, temp_top 
    use gv_irradiance_sunshade, only: check, soilnet
    use gv_scale_declarations,  only: boltz, grid, pi, time_holder
    use gv_veg,                 only: lai, lat, modrnet  
    use spa_io,                 only: handle_output

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    integer :: counter, i
    real    :: ab, beam, clump, dayd, dec, diff, em, estfdiff, ff, kbm, lnet, long,       &
               nirbeam, nirdiff, period, PIb, PId, radnet, rtime, skyabsl, skyabsn, soem, &
               soilt, suml, sumn, sump, sun, SWR, totlong, totnir, totpar, xhour

    real,dimension(grid%canopy) :: abslong, abspar_shade, abspar_sun, absnir, sunfrac, temp, totparabs

    ! include an extra value each side bracketing grid%canopy..
    real,dimension(0:grid%canopy+1) :: downlong, downnir, downpar, uplong, upnir, uppar


    ! Check stuff is allocated before we get started..
    if (.not. allocated(    checkpar) )  allocate(    checkpar(grid%canopy) )
    if (.not. allocated(    leaf_sun) )  allocate(    leaf_sun(grid%canopy) )
    if (.not. allocated(      longem) )  allocate(      longem(grid%canopy) )
    if (.not. allocated(      nirabs) )  allocate(      nirabs(grid%canopy) )
    if (.not. allocated(parabs_shade) )  allocate(parabs_shade(grid%canopy) )
    if (.not. allocated(  parabs_sun) )  allocate(  parabs_sun(grid%canopy) )

    ! calculations..
    period    = real(time%steps_per_day)
    clump     = 1.    ! clumping factor - =1 means unclumped, random distribution
    ff        = 1.0   ! correction factor for stems etc on attenuation, currently unused and untested, set at 1.0.
    dayd      = real( time%day )
    totparabs = 0.
    rtime     = real( time%step )

    ! solar characteristics..
    dec = -23.4 * cos( ( 360. * ( dayd + 10. ) / 365. ) * pi / 180. ) * pi / 180.
    ! diffuse fraction of radiation (Amthor, 1994)

    ! Estimate diffuse fraction of radiation from measured vs maximum radiation..
    estfdiff = frac_diffuse_rad( dec , lat , par_top , rtime )
    ! day is divided into 48 1/2hour-long periods
    fdiff = estfdiff

    ! reset diffuse arrays
    downnir  = 0.  ;  upnir    = 0.  ;  downpar  = 0.
    uppar    = 0.  ;  downlong = 0.  ;  uplong   = 0.
    absnir = 0. ; abspar_sun = 0. ; abspar_shade = 0. ; abslong = 0.
    do i = 1 , grid%canopy
      temp( i ) = temp_top - 0.1111 * (i-1) * ( temp_top - temp_bot )
    enddo

    ! detemine light extinction coefficent from sun angle..
    xhour = pi / 180. * 15. * ( rtime - 0.5 * period ) * 24. / period

    ! sin(beta) - solar geometry affects attenuation
    sun = cos(lat) * cos(xhour) * cos(dec) + sin(lat) * sin(dec)

    ! determine extinction coefficient for direct beam radiation..
    ! (but only for light levels above a threshold)
    if ( sun .ge. 0.06 ) then
      kbm = 1. / ( spherical * sun )
    else
      ! At low sun angles extinction coefficient is zero..
      ! ( Only occurs near sunrise and sunset. Without this )
      ! ( correction we get very unrealistic estimates.     )
      sun = 0.
      kbm = 0.0
    endif

    if ( sun .gt. 0. ) then
      ! If sun has risen then partition par between beam and diffuse
      beam = ( 1. - fdiff ) * par_top    ! light attentuation(PPFD) initialise layer above canopy
      diff =    fdiff       * par_top
    else 
      beam = 0.
      diff = par_top
    endif

    ! ??
    SWR     = sw_rad
    PIb     = 0.48 * ( 1. - fdiff ) * SWR
    PId     = 0.5 * SWR - PIb
    nirbeam = ( 1. - fdiff ) * SWR - PIb
    nirdiff = fdiff * SWR - PId

    ! determine longwave incident from sky (Jones p.27) kW m-2 long=.350..
    long = 0.001 * boltz * ( temp_top + 273.15 - 20. )**4

    totnir  = nirbeam + nirdiff
    totpar  = beam + diff
    totlong = long

    ! reset (module) values to zero
    soilabsn = 0.  ;  soilabsp = 0.  ;  skyabsp = 0. ;  skyabsn = 0.
    soilabsl = 0.  ;  skyabsl  = 0.  ;  counter = 0  ;  em = 0.  ;  ab = 0.

    ! start the attenuation routines..
    do while ( counter .lt. 3 )

      ! multiple passes through the canopy - 3 is enough to account for every photon in general
      if ( totpar .gt. 1. ) then    ! if the sun is up, then call the routines....

        ! firsly calculate PAR attenuation - only do sunlit and shaded calculations for PAR..
        call attenuate_PAR( ff , lai , beam , diff , parrefl , partrans , soilpar , &
                            kbm , clump , abspar_sun , abspar_shade , uppar ,       &
                            downpar , soilabsp , skyabsp , sunfrac , sump           )

       ! set leaf_sun after only the first pass, when beam radiation is incident
        if ( counter .eq. 0. ) then
          do i = 1 , grid%canopy
            leaf_sun( i ) = sunfrac( i )
          enddo
        endif

        ! next do Near Infra Red..
        call attenuate_NIR( ff , lai , nirbeam , nirdiff , nirrefl , nirtrans , soilnir, kbm , &
                            clump , absnir , upnir , downnir , soilabsn , skyabsn , sumn       )

      else
        abspar_sun = 0.  ;  abspar_shade = 0.
        absnir     = 0.  ;  sunfrac      = 0.
      endif

      soilt = temp_top    ! use air temp from current timestep

      ! finally calculate the longwave..
      call longwave( ff , lai , long , suml , abslong , uplong , downlong , soilabsl , &
                     skyabsl , counter , temp , soilt , em , soem , ab , clump         )

      ! reset everything..
      beam = 0.  ;  diff   = 0.  ;  nirbeam = 0.  ;  nirdiff = 0.
      long = 0.  ;  radnet = 0.  ;  lnet    = 0.

      ! increment the counter..
      counter = counter + 1 

    enddo

    ! calculate output for use in other routines
    do i = 1 , grid%canopy
        parabs_sun(i) =   abspar_sun(grid%canopy+1-i)    ! beam par absorbed in each layer
      parabs_shade(i) = abspar_shade(grid%canopy+1-i)    ! diffuse par absorbed in each layer

      ! total shortwave radiation for heating = par + nir , per m2 ground area per layer
      nirabs(i)    =  absnir(grid%canopy+1-i)
      longem(i)    = abslong(grid%canopy+1-i) * 1000.
      lnet         = lnet   + longem(i)
      radnet       = radnet + nirabs(i) 
      totparabs(i) = abspar_sun(grid%canopy+1-i) + abspar_shade(grid%canopy+1-i)
    enddo

    !  Radiation absorbed by soil (W.m-2); no long emission, isothermal or otherwise
    soilnet = soilabsn + 1000. * soilabsl + soilabsp / par_ratio

    check = sum(abspar_sun) + sum(abspar_shade) + soilabsp + skyabsp
    if ( par_top .gt. 0 ) then
      call handle_output( 4 , output_data =  (/ soilnet, soilabsn, soilabsl,  &
                        soilabsp, par_top, fdiff, check, skyabsp, abspar_sun, &
                        parabs_sun, abspar_shade, parabs_shade, leaf_sun /)   )
    endif
    modRnet  = SWR - skyabsn - skyabsp / par_ratio - 1e3 * ( skyabsl - totlong )
    checkpar = 0.

  end subroutine solar
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine attenuate_NIR( ff , lai , bm , df , refl , trans , reflsoil , kbm , clump , &
                            absorb , updf , downdf , soilrad , skyrad , sum_rad )

    use gv_scale_declarations, only: grid

    implicit none

    ! arguments..
    real,intent(in)    :: bm, clump, df, ff, kbm, lai(grid%canopy), refl, reflsoil, trans
    real,intent(inout) :: absorb(grid%canopy), downdf(0:grid%canopy+1), skyrad, soilrad, updf(0:grid%canopy+1)
    real,intent(out)   :: sum_rad

    ! local variables..
    integer :: i
    real :: decay, intbm, intdf, kdf, leafrad
    real    :: beam(0:grid%canopy+1)


    ! calculations..
      beam(grid%canopy+1) = bm
    downdf(grid%canopy+1) = df

    ! diffuse radiation approximated by beta=30 degrees
    kdf = 0.5

    do i = grid%canopy , 1 , -1

       ! determine unintercepted radiation
       beam(i) = beam(i+1) * exp( -kbm * clump * lai(grid%canopy+1-i) )

       ! and thus intercepted radiation
       intbm = beam(i+1) - beam(i)

       ! now for diffuse
       decay = exp( -kdf * clump * lai(grid%canopy+1-i) )

       downdf(i) = downdf(i) + downdf(i+1) * decay

       intdf = downdf(i+1) * ( 1. - decay )

       ! correct for transmittance (trans)
       absorb(i) = absorb(i) + intbm * ( 1. - trans - refl ) / ff

       ! correct for reflectance (refl)
       absorb(i) = absorb(i) + intdf * ( 1. - trans - refl ) / ff

       ! add transmitted beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + trans * ( intbm + intdf )

       ! add reflected beam & diffuse to upward diffuse
       updf(i) = updf(i) + refl * ( intbm + intdf )

    enddo

    ! reflectance by soil of radiation that penetrates all vegetation
    updf(0) = reflsoil * ( beam(1) + downdf(1) )
    soilrad = soilrad + ( 1. - reflsoil ) * ( beam(1) + downdf(1) )

    ! reset the downwelling diffuse radiation..
    downdf = 0.

    do i = 1 , grid%canopy

       ! now return upwards through the canopy, dealing with refelcted radiation
       ! unintercepted
       decay = exp( -kdf * clump * lai(grid%canopy+1-i) )

       updf(i) = updf(i) + updf(i-1) * decay

       ! intercepted
       intdf = updf(i-1) * ( 1. - decay )

       ! absorbed determined from transmittance/reflectance
       absorb(i) = absorb(i) + intdf * ( 1. - trans - refl ) / ff

       ! add reflected beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + refl * intdf

       ! add transmitted beam & diffuse to upward diffuse
       updf(i)   = updf(i) + trans * intdf
       updf(i-1) = 0.

    enddo

    skyrad = skyrad + updf(grid%canopy)

    updf(grid%canopy) = 0.

    leafrad = sum( absorb(1:grid%canopy) )

    sum_rad = soilrad + skyrad + leafrad

  end subroutine attenuate_NIR
  !
  !----------------------------------------------------------------------
  !
  subroutine attenuate_PAR( ff , lai , bm , df , refl , trans , reflsoil ,  &
                            kbm , clump , absorb_sun, absorb_shade , updf , &
                            downdf , soilrad , skyrad , sunfrac , sum_rad   )

    use gv_scale_declarations, only: grid

    implicit none

    ! arguments..
    real,intent(in)    :: clump, bm, df, ff, kbm, lai(grid%canopy), refl, reflsoil, trans
    real,intent(inout) :: absorb_sun(grid%canopy), absorb_shade(grid%canopy), downdf(0:grid%canopy+1), &
                          skyrad, soilrad, updf(0:grid%canopy+1)
    real,intent(out)   :: sum_rad, sunfrac(grid%canopy)

    ! local variables..
    integer :: i
    real    :: beam(0:grid%canopy+1), cumlai, decay, intbm, intdf, kdf, leafrad, suncum, sunla(grid%canopy), sunprev


    ! calculations..
    beam(grid%canopy+1) = bm
    downdf(grid%canopy+1) = df

    ! initial values..
    cumlai = 0.  ;  sunprev = 0.

    ! For ease, we approximate diffuse radiation by beta=30 degrees,
    ! which gives an extinction coefficient of..
    kdf = 0.5

    ! Calculate how the radiation is diffused as it travels down
    ! through the canopy layers..
    do i = grid%canopy , 1 , -1

      ! If there is sunlight and some extinction in this layer then..
      if ( ( kbm .gt. 0. ) .and. ( bm .gt. 0. ) ) then

        ! first calculate cumulative leaf area..
        cumlai = cumlai + lai( grid%canopy+1 - i )

        ! then total sunlit for cumlai. ( kbm = 1/(2sinBeta) )
        suncum = ( 1 - EXP( -kbm * cumlai ) ) / kbm

        ! sunlit area in this layer is cumulative sunlit minus previous
        sunla(grid%canopy+1-i) = suncum - sunprev

        sunprev = suncum     ! save previous

        ! determine sunlight fraction
        if ( sunla( grid%canopy+1 - i ) .gt. 0. ) then
          sunfrac( grid%canopy+1 - i ) = sunla( grid%canopy+1 - i ) / lai( grid%canopy+1 - i )
        else
          sunfrac( grid%canopy+1 - i ) = 0.
        endif

        beam( i ) = bm    ! Beam radiation on sunlit leaves is constant

        ! Intercepted radiation dI = Io.k.Ls
        ! (dI=intercepted rad, Io=downwelling beam, Ls =sunlit leaf area, k=extinction)
        intbm = bm * kbm * sunla( grid%canopy+1 - i )

      else

        ! no sunlight..
        intbm      = 0. ! intercepted radiation
        beam(i)    = 0. !
        sunfrac(i) = 0. ! fraction of sunlit leaf

      endif

      ! now for diffuse
      decay     = exp( -kdf * clump * lai(grid%canopy+1-i) )            !attenuation factor
      downdf(i) = downdf(i) + downdf(i+1) * decay        !diffuse radiation that passes through layer
      intdf     = downdf(i+1) * ( 1. - decay )                !intercepted diffuse radiation in each layer

      ! Absorption of direct beam (correct interception for transmittance and reflectance)..
      absorb_sun(i)   = absorb_sun(i)   + intbm * ( 1. - trans - refl ) / ff

      ! Absorption of diffuse..
      absorb_shade(i) = absorb_shade(i) + intdf * ( 1. - trans - refl ) / ff

      ! add transmitted beam & diffuse to downward diffuse
      downdf( i ) = downdf( i ) + trans * ( intbm + intdf )

      ! add reflected beam & diffuse to upward diffuse
      updf( i )   = updf( i ) + refl * ( intbm + intdf )

    enddo

    if ( bm .gt. 0. ) then
      ! Direct beam radiation reaching soil surface is directly calculated by this eqn...
      beam(1) = ( 1. / kbm - suncum ) * bm * kbm
    endif

    ! reflectance by soil of radiation that penetrates all vegetation
    updf(0) = reflsoil * ( beam(1) + downdf(1) )
    soilrad = soilrad + ( 1. - reflsoil ) * ( beam(1) + downdf(1) )

    ! reset downwelling diffuse radiation..
    downdf = 0.

    do i = 1 , grid%canopy

      ! now return upwards through the canopy, dealing with reflected radiation
      ! unintercepted
      decay = exp( -kdf * clump * lai(grid%canopy+1-i) )
      updf(i) = updf(i) + updf(i-1) * decay

      ! intercepted
      intdf = updf( i - 1 ) * ( 1. - decay )

      ! absorbed determined from transmittance/reflectance
      absorb_shade(i) = absorb_shade(i) + intdf * ( 1. - trans - refl ) / ff

      ! add reflected beam & diffuse to downward diffuse
      downdf(i) = downdf(i) + refl * intdf

      ! add transmitted beam & diffuse to upward diffuse
      updf(i)   = updf(i) + trans * intdf
      updf(i-1) = 0.

    enddo

    skyrad  = skyrad + updf(grid%canopy)
    updf(grid%canopy) = 0.
    leafrad = 0.
    leafrad = sum( absorb_sun(1:grid%canopy) ) + sum( absorb_shade(1:grid%canopy) )
    sum_rad = soilrad + skyrad + leafrad

  end subroutine attenuate_PAR
  !
  !----------------------------------------------------------------------
  !
  real function frac_diffuse_rad( dec , lat , parsteps , rtime )
    
    ! Determines the ratio of actual-to-potential radiation (ksko)   !
    ! for the day, from the total potential daily PAR, and then uses !
    ! a relationship from Erbs et al (1982) to estimate fraction of  !
    ! incoming radiation that is diffuse.                            !

    use gv_scale_declarations, only: pi, time

    implicit none

    ! arguments..
    real,intent(in) :: dec, lat, parsteps, rtime

    ! local variables..
    real :: ff, hourangle, ko, ks, ksko, So

    ! calculations..
    So        = 1360.                ! solar constant (Wm-2)
    ks        = parsteps / 2.3       ! (MJ m-2 d-1)
    hourangle = pi / 180. * 15. * ( rtime - 0.5 * real(time%steps_per_day) ) * 24. / real(time%steps_per_day)

    ! extraterrestrial radiation (MJ m-2)
    ko        = So * ( sin( dec ) * sin( lat ) + cos( lat ) * cos( dec ) * cos( hourangle ) )
    ksko      = ks / ko  ! hourly ksko

    ! Diffuse fraction
    if (ksko .lt. 0.) then
      ff = 1.0
    else if ( ksko .lt. 0.22 ) then
      ff = 1.0 - 0.09 * ksko
    else if ( ksko .lt. 0.8 ) then
      ff = 0.9511                  &
          - 0.1604 * ksko          &
           + 4.388 * ksko**2       &
            - 16.638 * ksko**3     &
             + 12.336 * ksko**4
    else 
      ff = 0.165
    endif
    frac_diffuse_rad = ff 

  end function frac_diffuse_rad
  !
  !----------------------------------------------------------------------
  !
  subroutine longwave( ff, lai, df, sum_rad, absorb, updf, downdf, soilrad, skyrad, count, & 
                       Ta, Ts, totemit, soilem, totab, clump)

    use gv_scale_declarations, only: boltz, grid

    implicit none

    ! arguments..
    integer,intent(in) :: count
    real,intent(in)    :: clump, df, ff, lai(grid%canopy), Ta(grid%canopy), Ts
    real,intent(inout) :: absorb(grid%canopy), downdf(0:grid%canopy+1), skyrad, soilrad, totemit, updf(0:grid%canopy+1)
    real,intent(out)   :: soilem, sum_rad, totab

    ! local variables..
    integer :: i
    real    :: decay, emiss, eup, edown, intdf, kdf, leafrad


    emiss = 0.96
    downdf(grid%canopy+1) = df

    ! diffuse radiation approximated by beta=30 degrees
    kdf = 0.5
    do i = grid%canopy , 1 , -1

       decay = exp( -kdf * clump * lai(grid%canopy+1-i) )

       if ( count .eq. 0 ) then

         ! longwave radiative heat loss from top side of leaf (KW m-2)
         ! The '*(1-decay)' corrects emissions into the 1-D, vertical
         eup = 0.001 * emiss * boltz * ( Ta(grid%canopy+1-i) + 273.15 )**4 * ( 1. - decay )

         totemit = totemit + 2. * eup * lai(grid%canopy+1-i)

       else

         eup = 0.

       endif

       ! and bottom side..
       edown     = eup
       downdf(i) = downdf(i) + downdf(i+1) * decay + edown
       intdf     = downdf(i+1) * ( 1. - decay )

       ! correct for transmittance (trans) and reflectance (refl) & PAI/LAI ratio
       absorb(i) = absorb(i) + intdf * emiss / ff - eup - edown

       ! add transmitted diffuse to downward diffuse
       downdf(i) = downdf(i) + 0.5 * ( 1. - emiss ) * intdf

       ! add reflected diffuse to upward diffuse
       updf(i) = updf(i) + 0.5 * ( 1. - emiss ) * intdf + eup
       totab = totab + intdf * emiss
    enddo

    ! reflectance by soil of radiation that penetrates all vegetation
    soilem = emiss * boltz * ( Ts + 273.15 )**4 * 0.001
    if ( count .eq. 0 ) then
       updf(0) = 0.04 * downdf(1) + soilem
    else
       updf(0) = 0.04 * downdf(1)
    endif
    soilrad = soilrad + 0.96 * downdf(1)

    ! reset downwelling diffuse radiation..
    downdf = 0.

    do i = 1 , grid%canopy

       ! now return upwards through the canopy, dealing with refelcted radiation
       ! unintercepted
       decay   = exp( -kdf * clump * lai(grid%canopy+1-i) )

       updf(i) = updf(i) + updf(i-1) * decay

       ! intercepted
       intdf = updf(i-1) * ( 1. - decay )

       ! absorbed determined from transmittance/reflectance
       absorb(i) = absorb(i) + intdf * emiss / ff

       ! add reflected beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + 0.02 * intdf

       ! add transmitted beam & diffuse to upward diffuse
       updf(i)   = updf(i) + 0.02 * intdf
       updf(i-1) = 0.
       totab     = totab + intdf * 0.96
    enddo

    skyrad = skyrad + updf(grid%canopy)
    updf(grid%canopy) = 0.
    leafrad = sum( absorb(1:grid%canopy) )
    sum_rad = soilrad + skyrad + leafrad

  end subroutine longwave
  !
  !----------------------------------------------------------------------
  !
end module light
!
!------------------------------------------------------------------------
!
