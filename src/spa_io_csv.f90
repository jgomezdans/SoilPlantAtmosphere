! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_io_csv

  !! contains routines specific to reading/writing csv files. !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: close_output_csv, open_output_csv, read_met_csv, read_soil_csv, read_veg_csv, &
            write_assimilate_output_csv, write_output_csv, write_predictor_output_csv,    &
            write_solar_output_csv, write_soilday_output_csv


  ! to keep track of whether write-routine already called this timestep..
  real,allocatable :: last_call(:)

  ! Used in the write_csv_output s/r..
  real,dimension(:),allocatable :: lemod, mmmod, neemod, timeflux, wetle
  real :: sum_alloc_to_foliage, sum_alloc_to_labile, sum_alloc_to_roots, sum_alloc_to_stem,              &
          sum_alloc_from_labile, sum_litterfall_foliage, sum_litterfall_stem, sum_litterfall_roots,      &
          sum_stock_foliage, sum_stock_stem, sum_stock_roots, sum_stock_litter, sum_stock_soilOrgMatter, &
          sum_stock_labile, sum_stock_resp_auto, sum_resp_auto, sum_resp_h_litter,                       &
          sum_resp_h_soilOrgMatter, sum_decomposition, sum_GPP

  integer,parameter :: input_met_unit    = 20, &  ! Boundary condition
                       input_soil_unit   = 21, &  ! Initial
                       input_veg_unit    = 22, &  !    conditions
                       budget_unit       = 23, &  ! Various output..
                       daily_unit        = 24, &
                       energy_unit       = 25, &
                       fluxes_unit       = 26, &
                       gdd_unit          = 27, &
                       hourly_unit       = 28, &
                       ice_prop_unit     = 29, &
                       solar1_unit       = 30, &
                       solar2_unit       = 31, &
                       solar3_unit       = 32, &
                       solar4_unit       = 33, &
                       soil_status_unit  = 34, &
                       soil_temp_unit    = 35, &
                       soil_water_unit   = 36, &
                       stocks_unit       = 37, &
                       up_frac_unit      = 38, &
                       water_fluxes_unit = 39
            ! units 40-50 are reserved for 'layer_[i].csv' files
            ! units 51 & 52 are used in spa_restart.

  save

contains
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! public procedures first..
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine close_output_csv

     use gv_scale_declarations, only: grid
     use gv_veg,                only: lafrac

     implicit none

     ! local variables..
     integer :: i

     close( unit=daily_unit        , status='keep' )
     close( unit=gdd_unit          , status='keep' )
     close( unit=budget_unit       , status='keep' )
     close( unit=energy_unit       , status='keep' )
     close( unit=fluxes_unit       , status='keep' )
     close( unit=hourly_unit       , status='keep' )
     close( unit=ice_prop_unit     , status='keep' )
     do i = 1 , grid%canopy
       if (lafrac(i) .ne. 0.) close(unit=40+i,status='keep')
     enddo
     close( unit=solar1_unit       , status='keep' )
     close( unit=solar2_unit       , status='keep' )
     close( unit=solar3_unit       , status='keep' )
     close( unit=solar4_unit       , status='keep' )
     close( unit=soil_status_unit  , status='keep' )
     close( unit=soil_temp_unit    , status='keep' )
     close( unit=soil_water_unit   , status='keep' )
     close( unit=stocks_unit       , status='keep' )
     close( unit=water_fluxes_unit , status='keep' )

  end subroutine close_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine open_output_csv( outdir , decid_flag )

    ! This opens the standard output files !
    !  and writes a header to each.        !

    use gv_scale_declarations, only: fname_length, grid
    use gv_veg,                only: lafrac

    implicit none

    ! arguments..
    character(len=*),intent(in) :: outdir
    logical,intent(in)          :: decid_flag

    ! local variables..
    character(fname_length) :: filename
    integer                 :: i

    ! Open all output files (alphabetical order!)..
    call open_file( trim(outdir)//'budget.csv', budget_unit , header='delC,neesum' , daily=.true. )

    call open_file( trim(outdir)//'energy.csv', energy_unit , header = &
                     'Qh,Qe,Qn,Qc,airt,surfacet,soilt2,drythick' )

    call open_file( trim(outdir)//'daily.csv', daily_unit , header = &
         'daymm,modle,modess,modwet,SWP,rplant,CSRes,lsc,totalflux' , daily = .true. )

    if (decid_flag) then
      call open_file( trim(outdir)//'gdd.csv', gdd_unit , header='time,gdd,mint'  , daily = .true. )
    endif

    call open_file( trim(outdir)//'hourly.csv', hourly_unit , header = &
                      'gpp,ra,rh1,rh2,le,trans,ess,wetle,gpp,nee' )

    call open_file( trim(outdir)//'iceprop.csv', ice_prop_unit , header=&
            'i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15' )

    do i = 1 , grid%canopy
      if ( lafrac(i) .ne. 0. ) then 
        write(filename,"(A,i0.2,A)")"layer_",i,".csv"
        call open_file( trim(outdir)//trim(filename), 40+i , header = &
               'gs,agr,res,psil,ci,etr,tempdf,wdef,par,rad,cca,la,cistar,an' )
      endif
    enddo

    call open_file( trim(outdir)//'solar_part1.csv', solar1_unit , header=&
                   'soilnet,soilabsn,soilabls,soilabsp,par_top' )

    call open_file( trim(outdir)//'solar_part2.csv', solar2_unit , header=&
       'sum(abspar_sun),(1.-fdiff)*par_top,(parabs_sun(i),i=1,10)' )

    call open_file( trim(outdir)//'solar_part3.csv', solar3_unit , header=&
     'sum(abspar_shade),fdiff*par_top,(parabs_shade(i),i=1,10),skyabsp' )

    call open_file( trim(outdir)//'solar_part4.csv', solar4_unit , header=&
                           'fdiff,(leaf_sun(i),i=1,10),check' )

    call open_file( trim(outdir)//'soil_status.csv', soil_status_unit , header=&
                   'of,totet,uf,sin,wtr,wch,flux,chk,w+1,w-1,delw,'    &
                 //'chkdff,ppt+,w+,w-,w+c,w-c,snow weight,snow height' )

    call open_file( trim(outdir)//'soil_temp.csv', soil_temp_unit , header=&
       'soil-temp(layer 1),soil-temp(layer 2),soil-temp(layer 3),'     &
     //'soil-temp(layer 4),soil-temp(layer 5),soil-temp(layer 6),'     &
     //'soil-temp(layer 7),soil-temp(layer 8),soil-temp(layer 9),'     &
     //'soil-temp(layer 10),soil-temp(layer 11),soil-temp(layer 12),'  &
     //'soil-temp(layer 13),soil-temp(layer 14),soil-temp(layer 15)'   )

    call open_file( trim(outdir)//'soilwater.csv',  soil_water_unit , header=&
     'w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w_swp' )

    call open_file( trim(outdir)//'upfrac.csv',     up_frac_unit , header=&
     'up1,up2,up3,up4,up5,up6,up7,up8,up9,up10,up11,up12,up13,up14,up15' )

    call open_file( trim(outdir)//'waterfluxes.csv', water_fluxes_unit , &
                    header='modess,runoff,ppt,trans,delw,disc,cwtr,check,diff,modwet,unint,canst,snow_watermm' , &
                    daily=.true. )

    ! Dependent upon plant-func-type..
    if (decid_flag) then

      call open_file( trim(outdir)//'fluxes.csv', fluxes_unit , header=&
        'Ra,Af,Aw,Ar,Lf,Lw,Lr,Rh1,Rh2,D,G,Atolab,Afromlab,neesum,sumresp,daymm' , daily=.true. )

      call open_file( trim(outdir)//'stocks.csv', stocks_unit , header=&
        'stock_foliage,stock_stem,stock_roots,stock_litter,stock_soilOrgMatter,stock_labile,stock_resp_auto' , daily=.true. )

    else

      call open_file( trim(outdir)//'fluxes.csv', fluxes_unit , header=&
        'Ra,Af,Aw,Ar,Lf,Lw,Lr,Rh1,Rh2,D,G,neesum,sumresp,daymm' , daily=.true. )

      call open_file( trim(outdir)//'stocks.csv', stocks_unit , header=&
        'stock_foliage,stock_stem,stock_roots,stock_litter,stock_soilOrgMatter,stock_resp_auto' , daily=.true. )

    endif

  end subroutine open_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_met_csv( filename , met )

    ! open the met driver (csv) file, and load the data from it !

    use gv_scale_declarations, only: met_drivers, time, user_opts
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    type(met_drivers)           :: met

    ! local variables..
    character(len=200) :: line
    integer            :: i, num_lines, status
    real               :: dummy

    ! open file
    call open_file( filename , input_met_unit , readonly=.true. )

    ! get header out of the way..
    read(unit=input_met_unit,fmt=*) line

    ! count the number of remaining lines in the file..
    status = 0    ;     num_lines = 0
    do
      read(input_met_unit,fmt=*,iostat=status) line
      if ( status .ne. 0. ) exit
      num_lines = num_lines + 1
    enddo

    ! make the variables the correct size..
    allocate( met%ambient_co2( num_lines ) , &
                      met%par( num_lines ) , &
                   met%precip( num_lines ) , &
                   met%sw_rad( num_lines ) , &
             met%sfc_pressure( num_lines ) , &
                     met%temp( num_lines ) , &
                      met%vpd( num_lines ) , &
                 met%wind_spd( num_lines )   )

    ! Go back to start of file..
    rewind( input_met_unit )

    ! get header out of the way again..
    read(unit=input_met_unit,fmt=*) line

    ! and now load all data..
    do i = 1 , num_lines
      read(unit=input_met_unit,fmt=*)  dummy, met%temp(i),      &
                          met%ambient_co2(i), met%wind_spd(i),  &
                          met%sw_rad(i),      met%vpd(i),       &
                          met%par(i),         met%precip(i)
    enddo

    ! Check sanity of temperature..
    if ( ( maxval(met%temp) .gt. 150. ) .or. ( minval(met%temp) .lt. -100 ) ) then
      write(message,*) "Temperature appears to be in units other than Celcius or Kelvin."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Is precip in volume-per-timestep (mm), or rate (mm/sec)?
    if ( user_opts%precip_is_rate ) then
      call write_log("Converting precipitation from rate to volume by "&
                   //"multiplying by nos of seconds per timestep.")
      met%precip = met%precip * time%seconds_per_step
    endif

    ! Surface pressure not provided..
    met%sfc_pressure = met%const_sfc_pressure

    ! In case user doesn't want to use co2 in the file..
    if ( .not. user_opts%use_co2_from_met_file ) met%ambient_co2 = met%const_ambient_co2

    ! In case user doesn't want to use the par in the file..
    if ( .not. user_opts%use_par_from_met_file )  met%par = 2.3 * met%sw_rad 

    ! Adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2 ) met%wind_spd = 0.2

  end subroutine read_met_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_soil_csv( filename , decid_flag )

    ! open and read the contents of the soil.csv file. !
 
    use gv_hydrol,            only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations,only: grid
    use gv_snow_info,         only: snowheight, snowweight
    use gv_soil_structure,    only: iceprop, layer_depth, mineralfrac, organicfrac, resp_rate_temp_coeff, &
                                    root_radius, thickness, waterfrac, soil_temp

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    logical,intent(in)          :: decid_flag

    ! local variables..
    character(len=200) :: header
    integer            :: i

    ! Open the soil parameters file..
    call open_file( filename , input_soil_unit , readonly=.true. )

    ! read soil properties..
    read(unit=input_soil_unit,fmt=*)header
    read(unit=input_soil_unit,fmt=*)header,(thickness(i),  i=1,grid%core)
    read(unit=input_soil_unit,fmt=*)header,(layer_depth(i),i=1,grid%core)
    read(unit=input_soil_unit,fmt=*)header,(organicfrac(i),i=1,grid%core)
    read(unit=input_soil_unit,fmt=*)header,(mineralfrac(i),i=1,grid%core)
    read(unit=input_soil_unit,fmt=*)header,(waterfrac(i),  i=1,grid%core)
    read(unit=input_soil_unit,fmt=*)header,(soil_temp(i),  i=1,grid%core)
    read(unit=input_soil_unit,fmt=*)header,(iceprop(i),    i=1,grid%core)
    if (decid_flag) then
      ! convert to array if sand content varies in all layers
      read(unit=input_soil_unit,fmt=*)header,soil_frac_sand(1)
      ! set parameter constant for all soil layers
      soil_frac_sand=soil_frac_sand(1)
      ! convert to array if clay content varies in all layers
      read(unit=input_soil_unit,fmt=*)header,soil_frac_clay(1)
      ! set parameter constant for all soil layers
      soil_frac_clay=soil_frac_clay(1)
    else
      read(unit=input_soil_unit,fmt=*)header,(soil_frac_sand(i),i=1,grid%core)
      read(unit=input_soil_unit,fmt=*)header,(soil_frac_clay(i),i=1,grid%core)
    endif
    read(unit=input_soil_unit,fmt=*)header,root_radius
    read(unit=input_soil_unit,fmt=*)header,snowweight
    read(unit=input_soil_unit,fmt=*)header,snowheight
    read(unit=input_soil_unit,fmt=*)header,resp_rate_temp_coeff

  end subroutine read_soil_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_veg_csv(filename,decid_flag)

    ! open and read the contents of the veg.csv file. !

    use gv_scale_declarations, only: grid, pi
    use gv_soil_structure,     only: max_depth,max_storage,root_k,rootresist,through_fall
    use gv_veg
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    logical,intent(in)          :: decid_flag

    ! local variables..
    character(len=200) :: header, dummy
    integer            :: i

    ! Open the plant physiological parameters
    call open_file( filename , input_veg_unit , readonly=.true. )

    !read vegetation parameters
    read(unit=input_veg_unit,fmt=*)header,(lafrac(i),i=1,grid%canopy)       ! LA fraction
    read(unit=input_veg_unit,fmt=*)header,(nfrac(i),i=1,grid%canopy)        ! nitrogen fraction in each layer
    read(unit=input_veg_unit,fmt=*)header,(layer_height(i),i=1,grid%canopy) ! layer heights
    read(unit=input_veg_unit,fmt=*)header,avN                         ! average foliar N, gN m-2 leaf area
    read(unit=input_veg_unit,fmt=*)header,gplant                      ! conductivity
    read(unit=input_veg_unit,fmt=*)header,minlwp                      ! critical LWP
    read(unit=input_veg_unit,fmt=*)header,iota                        ! stomatal efficiency
    read(unit=input_veg_unit,fmt=*)header,capac                       ! leaf capacitance    
    read(unit=input_veg_unit,fmt=*)header,lat                         ! latitude
    lat = lat * pi / 180.                                             !  -> convert from degrees to radians
    read(unit=input_veg_unit,fmt=*)header,dummy                       !
    read(unit=input_veg_unit,fmt=*)header,dimen                       ! leaf size
    read(unit=input_veg_unit,fmt=*)header,rootresist                  ! root resistivity
    read(unit=input_veg_unit,fmt=*)header,tower_height                ! height of measurement tower
    read(unit=input_veg_unit,fmt=*)header,conductivity                ! Does conductance vary with stem length? 0=NO, 1=YES
    read(unit=input_veg_unit,fmt=*)header,kappac                      ! ??
    read(unit=input_veg_unit,fmt=*)header,kappaj                      ! ??
    read(unit=input_veg_unit,fmt=*)header,LMA                         ! leaf mass per area
    read(unit=input_veg_unit,fmt=*)header,max_depth                   ! max rooting depth
    read(unit=input_veg_unit,fmt=*)header,root_k                      ! root mass for reaching 50% max depth
    read(unit=input_veg_unit,fmt=*)header,decomposition_rate          ! decomposition rate
    read(unit=input_veg_unit,fmt=*)header,frac_GPP_resp_auto          ! fraction of gpp respired
    read(unit=input_veg_unit,fmt=*)header,frac_alloc_foliage          ! npp allocated to foliage
    read(unit=input_veg_unit,fmt=*)header,frac_alloc_roots            ! remaining npp allocated to fine roots
    read(unit=input_veg_unit,fmt=*)header,turnover_rate_foliage       ! turnover rate of foliage
    read(unit=input_veg_unit,fmt=*)header,turnover_rate_stem          ! turnover rate of stem pool
    read(unit=input_veg_unit,fmt=*)header,turnover_rate_roots         ! turnover rate of fine roots
    read(unit=input_veg_unit,fmt=*)header,mineralisation_rate_litter        ! mineralisation rate of litter
    read(unit=input_veg_unit,fmt=*)header,mineralisation_rate_soilOrgMatter ! mineralisation rate of SOM/CWD
    if (decid_flag) then
      read(unit=input_veg_unit,fmt=*)header,GDD_thresh                ! GDD threshhold
      read(unit=input_veg_unit,fmt=*)header,min_Temp_thresh           ! minimum temperature threshold
      read(unit=input_veg_unit,fmt=*)header,frac_leafLoss_litter      ! fraction of leaf loss to litter
      read(unit=input_veg_unit,fmt=*)header,turnover_rate_labile      ! turnover rate of labile pool
      read(unit=input_veg_unit,fmt=*)header,resp_cost_labile_trans    ! respiratory cost of labile transfers
      read(unit=input_veg_unit,fmt=*)header,max_stock_foliage         ! maximum carbon allowed to be stored in foliage pool
    endif
    read(unit=input_veg_unit,fmt=*)header,turnover_rate_resp_auto     ! turnover rate of autotrophic respiration pool
    read(unit=input_veg_unit,fmt=*)header,stock_foliage               ! carbon stock of..foliage pool
    read(unit=input_veg_unit,fmt=*)header,stock_stem                  ! ..stem pool
    read(unit=input_veg_unit,fmt=*)header,stock_roots                 ! ..roots pool
    read(unit=input_veg_unit,fmt=*)header,stock_litter                ! ..litter pool
    read(unit=input_veg_unit,fmt=*)header,stock_soilOrgMatter         ! ..soil organic matter pool
    if (decid_flag) then
      read(unit=input_veg_unit,fmt=*)header,stock_labile              ! ..labile pool
    endif
    read(unit=input_veg_unit,fmt=*)header,stock_resp_auto             ! ..autotrophic respiration pool
    read(unit=input_veg_unit,fmt=*)header,through_fall
    read(unit=input_veg_unit,fmt=*)header,max_storage

  end subroutine read_veg_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_assimilate_output_csv( time , output )

    ! moved writes from assimilate (leaf.f90) to here !

    use gv_clim,               only: coa
    use gv_metab,              only: an, ci
    use gv_meteo,              only: gbb, la, par, psil, rad, temp, wdef
    use gv_scale_declarations, only: grid, time_holder
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    real,             intent(in) :: output(:)

    ! local variables..
    integer                 :: output_layer

    ! If first time through then we need to allocate..
    if ( .not. allocated(last_call) ) then
      allocate(last_call(grid%canopy))
      last_call = 0.
    endif

    output_layer = int( output(1) )        ! convert back from real to integer
    output_layer = max( output_layer , 1 ) ! make sure index is at least 1

    ! Check we have not been called already this time step!
    if ( last_call(output_layer) .lt. time%daytime ) then

      call write_to_file( 40+output_layer , (/ output(2:4), psil, ci, output(5), &
                                               output(6)-temp, wdef, par, rad,   &
                                               coa, la, an, gbb*output(7)       /) )
 
      last_call(output_layer) = time%daytime

    else

      write(message,*)'write_assimilate_output_csv has already been called this timestep: ignoring additional call'
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )

    endif

  end subroutine write_assimilate_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_output_csv( time , decid_flag , prevwater )

    ! write to output files. HARD-CODED. !

    use gv_clim,                only: daypar, dayppt, temp_bot, temp_top, wetev
    use gv_hourscale,           only: canopy_store, discharge, runoff, unintercepted
    use gv_irradiance_sunshade, only: check
    use gv_scale_declarations,  only: grid, met, time_holder
    use gv_snow_info,           only: snow_watermm
    use gv_soil_structure,      only: watericemm, weighted_swp
    use gv_veg

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    logical,intent(in)           :: decid_flag
    real,intent(inout)           :: prevwater

    ! local variables..
    integer :: i
    real    :: check2, currentC, currentwater, daymm, dayrad, daytrans, &
               delC, delwater, lambda, lsc, modess, modess2, modess3 ,  &
               modle, modwet, modwetle, neesum, rplant, sumresp, totalflux
    real,dimension(time%steps_per_day) :: posess

    ! if first time called..
    if ( .not. allocated(    lemod ) ) allocate(    lemod(time%steps_per_day) )
    if ( .not. allocated(    mmmod ) ) allocate(    mmmod(time%steps_per_day) )
    if ( .not. allocated(   neemod ) ) allocate(   neemod(time%steps_per_day) )
    if ( .not. allocated( timeflux ) ) allocate( timeflux(time%steps_per_day) )
    if ( .not. allocated(    wetle ) ) allocate(    wetle(time%steps_per_day) )

    ! calculations...
    if ( time%step .eq. 1 ) then
      ! reset sums
      sum_alloc_to_foliage = 0. ; sum_alloc_to_labile = 0. ; sum_alloc_from_labile  = 0.
      sum_alloc_to_roots   = 0. ; sum_alloc_to_stem = 0.
      sum_decomposition = 0.
      sum_litterfall_foliage = 0. ; sum_litterfall_roots = 0. ; sum_litterfall_stem = 0.
      sum_resp_auto = 0.       ; sum_resp_h_litter = 0. ; sum_resp_h_soilOrgMatter = 0.
      sum_stock_foliage   = 0. ; sum_stock_labile = 0.  ; sum_stock_litter = 0.
      sum_stock_resp_auto = 0. ; sum_stock_roots  = 0.  ; sum_stock_soilOrgMatter = 0.
      sum_stock_stem = 0.
      timeflux = 0.
    endif

    lambda              = 1000. * ( 2501.0 - 2.364 * temp_top )                  ! latent heat of vapourisation (J kg-1)
    wetle( time%step )  = lambda * wetev( time%step ) / time%seconds_per_step    ! wet leaves evap convert (mm t-1 to Wm-2)
    ! estimate NEE (umol m-2 s-1)...
    neemod( time%step ) = ( resp_auto + resp_h_litter + resp_h_soilOrgMatter - GPP ) / time%seconds_per_step * 1e6 / 12.
    ! sum modelled canopy & modelled soil LE, and wet canopy evap..
    lemod( time%step )  = transt( time%step ) + ess( time%step ) + wetle( time%step )
    ! lemod x length of timestep..
    mmmod( time%step )  = lemod( time%step) * ( time%seconds_per_step / lambda ) 

    ! sum flux from each tree layer..
    do i = 1 , grid%canopy
      timeflux( time%step ) = timeflux( time%step ) + flux( time%step , i )   ! (mmol m-2 GA s-1)
    enddo

    call write_to_file( hourly_unit , (/ GPP, resp_auto, resp_h_litter, resp_h_soilOrgMatter, lemod(time%step), transt(time%step), &
                        ess(time%step), wetle(time%step), gppt(time%step), neemod(time%step) /) )

    sum_alloc_to_foliage     = sum_alloc_to_foliage     + alloc_to_foliage
    sum_alloc_to_stem        = sum_alloc_to_stem        + alloc_to_stem
    sum_alloc_to_roots       = sum_alloc_to_roots       + alloc_to_roots
    sum_alloc_to_labile      = sum_alloc_to_labile      + alloc_to_labile
    sum_alloc_from_labile    = sum_alloc_from_labile    + alloc_from_labile
    sum_decomposition        = sum_decomposition        + decomposition
    sum_GPP                  = sum_GPP                  + GPP
    sum_litterfall_foliage   = sum_litterfall_foliage   + litterfall_foliage
    sum_litterfall_stem      = sum_litterfall_stem      + litterfall_stem
    sum_litterfall_roots     = sum_litterfall_roots     + litterfall_roots
    sum_resp_auto            = sum_resp_auto            + resp_auto
    sum_resp_h_litter        = sum_resp_h_litter        + resp_h_litter
    sum_resp_h_soilOrgMatter = sum_resp_h_soilOrgMatter + resp_h_soilOrgMatter
    sum_stock_foliage        = sum_stock_foliage        + stock_foliage
    sum_stock_labile         = sum_stock_labile         + stock_labile
    sum_stock_litter         = sum_stock_litter         + stock_litter
    sum_stock_roots          = sum_stock_roots          + stock_roots
    sum_stock_soilOrgMatter  = sum_stock_soilOrgMatter  + stock_soilOrgMatter
    sum_stock_stem           = sum_stock_stem           + stock_stem
    sum_stock_resp_auto      = sum_stock_resp_auto      + stock_resp_auto

    if ( time%step .eq. time%steps_per_day ) then
      modle     = sum(lemod) * time%seconds_per_step * 1e-6      ! modelled ecosystem water loss to atmosphere (MJ m-2 d-1)
      daytrans  = sum(transt) * time%seconds_per_step * 1e-6     ! modelled canopy transpiration (MJ m-2 d-1)
      modess    = sum(soiletmm)                                  ! modelled soil evaporation (mm d-1)
      modess2   = sum(ess) * time%seconds_per_step * 1e-6        ! modelled soil evaporation (MJ m-2 d-1)
      posess    = max( 0. , ess )                                ! remove dew component
      modess3   = sum(posess) * time%seconds_per_step * 1e-6     ! modelled soil evaporation (MJ m-2 d-1)
      dayrad    = daypar * 1e-6 / 2.3                            ! measured irradiance (MJ m-2 d-1)
      modwet    = sum(wetev)                                     ! modelled evap from wetted canopy (mm d-1)
      modwetle  = sum(wetle) * time%seconds_per_step * 1e-6      ! modelled evap from wet canopy (MJ m-2 d-1)
      daymm     = sum(mmmod)                                     ! total ET (mm d-1)
      totalflux = sum(timeflux) * time%seconds_per_step * 0.001 * 18 * 0.001  ! total flux  (mm d-1 or kg m-2 ga d-1)
      sumresp   = resp_auto + resp_h_litter + resp_h_soilOrgMatter
      neesum    = sumresp - GPP

      ! determine leaf specific conductance for 2nd canopy layer
      if ( lai(2) .gt. 0 ) then
        if ( conductivity .eq. 1 ) then
          rplant = layer_height(2) / ( gplant * lai(2) ) ! MPa s mmol-1 layer
        else
          rplant = 1. / ( gplant * lai(2) )         ! Conductance is constant with height
        endif
        lsc = ( 1. / ( rplant + canopy_soil_resistance(2) ) ) / lai(2)
      else
        rplant = -999.
        lsc    = -999.
      endif

      call write_to_file( daily_unit,                                            &
                          (/ daymm, modle, modess, modwet, weighted_SWP, rplant, &
                             canopy_soil_resistance(2), lsc, totalflux /) ,      &
                          daily = .true.                                          )
      if (decid_flag) then

        ! deciduous stocks..
        call write_to_file( stocks_unit,                                               &
                            (/ stock_foliage, stock_stem, stock_roots, stock_litter,   &
                               stock_soilOrgMatter, stock_labile, stock_resp_auto /) , &
                            daily = .true.                                              )

        ! fluxes, daily sums..
        call write_to_file( fluxes_unit ,                                                                       &
                            (/ sum_resp_auto, sum_alloc_to_foliage, sum_alloc_to_stem, sum_alloc_to_roots,      &
                               sum_litterfall_foliage, sum_litterfall_stem, sum_litterfall_roots,               &
                               sum_resp_h_litter, sum_resp_h_soilOrgMatter, sum_decomposition, sum_GPP,         &
                               sum_alloc_to_labile, sum_alloc_from_labile, neesum, sumresp, daymm          /) , &
                            daily = .true.                                                                       )
        

        currentC     = stock_foliage + stock_stem + stock_roots + stock_litter + stock_soilOrgMatter + &
                            stock_labile + stock_resp_auto

      else

        ! evergreen stocks..
        call write_to_file( stocks_unit ,                                                 &
                            (/ stock_foliage, stock_stem, stock_roots, stock_litter,      &
                               stock_soilOrgMatter, stock_resp_auto                  /) , &
                             daily = .true.                                                )

        ! fluxes, daily sums..
        call write_to_file( fluxes_unit ,                                                                    &
                             (/ sum_resp_auto, sum_alloc_to_foliage, sum_alloc_to_stem, sum_alloc_to_roots,  &
                                sum_litterfall_foliage, sum_litterfall_stem, sum_litterfall_roots,           &
                                sum_resp_h_litter, sum_resp_h_soilOrgMatter, sum_decomposition, sum_GPP,     &
                                neesum, sumresp, daymm /) ,                                                  &
                            daily = .true.                                                                    )

        currentC     = stock_foliage + stock_stem + stock_roots + stock_litter + stock_soilOrgMatter + stock_resp_auto

      endif

      currentwater = sum( watericemm )
      delwater     = currentwater - prevwater    ! change in soil water storage (mm)
      check2       = -1. * ( modess + runoff * 1e3 + totevap + discharge - unintercepted ) ! total water fluxes
      !check and delwater should be the same magnitude

      call write_to_file( water_fluxes_unit ,                                       &
                          (/ modess, runoff*1e3, dayppt, totevap, delwater,         &
                             discharge, currentwater, check2, check-delwater ,      &
                             modwet, unintercepted, canopy_store, snow_watermm /) , &
                          daily = .true.                                             )

      delC      = currentC - prevC
      prevwater = currentwater
      prevC     = currentC

      call write_to_file( budget_unit, (/delC, neesum/) , daily = .true. )

    endif

  end subroutine write_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_predictor_output_csv( output )

    ! moved writes from predictor (allocate_carbon.f90) to here !

    implicit none

    ! arguments..
    real,intent(in) :: output(:)

    call write_to_file( gdd_unit , output , daily = .true. )

  end subroutine write_predictor_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_solar_output_csv( output )

    ! moved writes from solar (light.f90) to here !

    use gv_clim,               only: par_top
    use gv_scale_declarations, only: grid

    implicit none

    ! arguments..
    real, intent(in) :: output(:)

    ! local variables..
    real,dimension(grid%canopy) :: abspar_shade_out, abspar_sun_out, &
                                   leaf_sun_out, parabs_shade_out, parabs_sun_out

    ! output =  (/ soilnet, soilabsn, soilabsl, soilabsp, par_top, fdiff, check, skyabsp,  &
    !                      abspar_sun(grid%canopy),  &
    !                      parabs_sun(grid%canopy),  &
    !                    abspar_shade(grid%canopy),  &
    !                    parabs_shade(grid%canopy),  &
    !                        leaf_sun(grid%canopy)   /)

    abspar_sun_out   = output( 9               : 9+  grid%canopy-1 )
    parabs_sun_out   = output( 9+  grid%canopy : 9+2*grid%canopy-1 )
    abspar_shade_out = output( 9+2*grid%canopy : 9+3*grid%canopy-1 )
    parabs_shade_out = output( 9+3*grid%canopy : 9+4*grid%canopy-1 )
    leaf_sun_out     = output( 9+4*grid%canopy : 9+5*grid%canopy-1 )

    call write_to_file( solar1_unit , output(1:5) )

    call write_to_file( solar2_unit , (/ sum(abspar_sun_out), (1.-output(6))*output(5), parabs_sun_out /) )

    call write_to_file( solar3_unit , (/ sum(abspar_shade_out), output(6)*output(5) , parabs_shade_out , output(8) /) )

    call write_to_file( solar4_unit , (/ output(6), leaf_sun_out, output(7) /) )

  end subroutine write_solar_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_soilday_output_csv( flag , output )

    ! moved writes from soilday (soil_functions.f90) to here !

    use gv_clim,               only: temp_bot, temp_top
    use gv_hourscale,          only: freeze, hourtemp, overflow, Qc, Qe, Qh, Qn, surface_watermm, &
                                     totet, underflow
    use gv_scale_declarations, only: grid
    use gv_snow_info,          only: snowheight, snowweight
    use gv_soil_structure,     only: drythick, fraction_uptake, iceprop, pptgain, soil_temp, waterfrac, &
                                     watergain, watericemm, waterloss, weighted_SWP
    use log_tools

    implicit none

    ! arguments..
    integer,intent(in) :: flag
    real,   intent(in) :: output(:)

    select case (flag)
    case (1)

      call write_to_file( energy_unit , (/ Qh, Qe, Qn, Qc, hourtemp-freeze, &
                                          output(1)-freeze, soil_temp(2)-freeze, drythick*1e3/) )

      call write_to_file( ice_prop_unit , iceprop(1:15) )

      call write_to_file( soil_status_unit, (/ 1e3*overflow, 1e3*totet, 1e3*underflow,              &
                                            surface_watermm, sum(watericemm), output(2:4),         &
                                            1e3*watergain(1), 1e3*waterloss(1), output(5:6),       &
                                            sum(pptgain), sum(watergain), sum(waterloss) ,         &
                                            watergain(grid%core), waterloss(grid%core), snowweight, snowheight /) )

      call write_to_file( soil_temp_unit , (/soil_temp(1:15)-273.15, temp_bot/) )

      call write_to_file( soil_water_unit , (/ waterfrac(1:15), weighted_swp /) )

    case (2)

      call write_to_file( up_frac_unit , fraction_uptake(1:15) )

    case default

      write(message,*)"flag supplied to handle_output:",flag," was not recognised"
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )

    end select

  end subroutine write_soilday_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine open_file( filename , unitnos , readonly , header , daily )
 
    ! This s/r is a wrapper to open(). It performs !
    ! sanity checks when opening files.  Control   !
    ! of write-permission is available through the !
    ! optional readonly argument, and if header is !
    ! provided it will be written to the file.     !

    use log_tools

    implicit none

    ! arguments..
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: unitnos
    logical,optional, intent(in) :: readonly, daily
    character(len=*),optional, &
                      intent(in) :: header

    ! local variables..
    integer            :: ios
    logical            :: print_day_msg, read_mode
    character(len=6)   :: header_fmt
    character(len=20)  :: file_status, write_status
    character(len=200) :: full_header

    ios = 0

    ! determine whether reading or writing..
    ! (assume that it is okay to write)
    read_mode    = .false.
    file_status  = 'replace'
    write_status = 'readwrite'
    if ( present( readonly ) ) then
      if ( readonly ) then
        read_mode    = .true.
        file_status  = 'old'
        write_status = 'read'
      endif
    endif

    ! determine if output data will be daily-avg'd..
    ! (assume that msg is not wanted)
    print_day_msg = .false.
    if ( present( daily ) ) print_day_msg = daily

    ! open file..
    open( unit=unitnos , file=trim(filename) , iostat=ios ,    &
          status=trim(file_status) , action=trim(write_status) )

    ! check opened file okay..
    if ( ios .ne. 0 ) then
      write(message,*)"Problem opening file : ",trim(filename)
      call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! write header if provided..
    if ( present( header ) ) then
      if ( read_mode ) then
        write(message,*)"You cannot write a header to a file that "&
                      //"you have opened in read-only mode!"
        call write_log( message , msg_warning , __FILE__ , __LINE__ )
      else
        if ( print_day_msg ) then
          write(unitnos,*)'!NB. This file contains data that is averaged or summed over daily periods.!'
        endif
        full_header='Time (days),'//trim(header)
        write(header_fmt,'(A2,I3.3,A1)')'(A',len(trim(full_header)),')'
        write(unitnos,header_fmt)trim(full_header)
      endif
    endif

  end subroutine open_file
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_to_file( unit_nos , out_data , daily )

    use gv_scale_declarations, only: time
    use log_tools

    implicit none

    ! arguments..
    integer,intent(in)          :: unit_nos
    real,   intent(in)          :: out_data(:)
    logical,optional,intent(in) :: daily

    ! local variables..
    logical           :: daily_data, file_open
    integer           :: n
    character(len=7)  :: write_permit
    character(len=20) :: write_fmt

    if (present(daily)) then
      daily_data = daily
    else 
      daily_data = .false.
    endif

    ! check file is open and can be written to..
    inquire( unit=unit_nos, opened=file_open, write=write_permit )

    if ( .not. file_open ) then

      ! Unit-nos is not for a file that is open!!
      write(message,*)"Unit_nos:",unit_nos," is not connected to an open file!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    else if ( trim(write_permit) .ne. 'YES' ) then

      ! Unit-nos is for a file in read-only setting!!
      write(message,*)"Unit number:",unit_nos," is not open in write-mode!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    else

      ! Okay to do the write..
      n = size( out_data )
      if ( daily_data ) then
        write(write_fmt,'(a,i0,a)') "(i0.3," , n ,'(",",f0.7))'
        write( unit=unit_nos , fmt=trim(write_fmt) )  time%day , out_data
      else
        write(write_fmt,'(a,i0,a)') "(f0.2," , n ,'(",",f0.7))'
        write( unit=unit_nos , fmt=trim(write_fmt) )  time%daytime , out_data
      endif
    endif

  end subroutine write_to_file
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
end module spa_io_csv
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
