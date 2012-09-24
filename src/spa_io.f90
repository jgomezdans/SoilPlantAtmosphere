! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_io

  !! coordinates calling of all io routines, including    !!
  !! reading user-config and opening/reading input files. !!

  use gv_scale_declarations, only: fname_length
#ifdef USE_NETCDF
  use spa_io_netcdf,         only: nc_met_file
#endif

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: handle_output, start_spa, update_met_drivers

#ifdef USE_NETCDF
  ! Not sure this is the right place to put this..
  type(nc_met_file),pointer :: met_nc
#endif

  ! variables.. 
  character(fname_length) :: spa_config_filename = "SPA.config" ! default value
  integer                 :: met_slice = 0,    &
                             prev_step_of_day

  save
  
contains
  !
  !-------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine handle_output( flag , time , iwater , output_data )

    ! wrapper to the various possible writes. !

    use gv_clim,               only: gdd, mint
    use gv_scale_declarations, only: time_holder, user_opts
    use log_tools
    use spa_io_csv
! #ifdef USE_NETCDF
!     use spa_io_netcdf
! #endif

    implicit none

    integer,intent(in)      :: flag
    ! flag takes the following values...
    ! 0   = open output files
    ! 1   = perform a standard write
    ! 2-7 = specific writes for subroutines within SPA
    ! 8 = spare
    ! 9   = close output files
    type(time_holder),optional,intent(in) :: time
    real,optional,intent(inout)           :: iwater
    real,dimension(:),optional            :: output_data

    if ( ( flag .ge. 3 ) .and. ( flag .le. 6) &
         .and. ( .not. present(output_data) ) ) then
      write(message,*)"When handle_output is called with flag=",flag,&
                         " you MUST supply output_data."
      call write_log( trim(message), msg_fatal, __FILE__ , __LINE__ )
    endif

    select case (flag)
    case (0) !open output files and write headers..
      if ( user_opts%std_csv_output ) &
            call open_output_csv( user_opts%output_directory, user_opts%veg_is_deciduous )
    case (1)
      if ( user_opts%std_csv_output ) then
        if ( present( iwater ) ) then
          call write_output_csv( time, user_opts%veg_is_deciduous , iwater )
       else
          write(message,*)"When handle_output is called with flag=1 you MUST supply iwater."
          call write_log( trim(message), msg_fatal, __FILE__ , __LINE__ )
        endif
      endif
    case (2)
      call write_predictor_output_csv( (/gdd, mint/) )
    case (3)
      call write_assimilate_output_csv( time, output_data )
    case (4)
      call write_solar_output_csv( output_data )
    case (5)
      call write_soilday_output_csv( 1 , output_data )
    case (6)
      call write_soilday_output_csv( 2 , output_data )
    case (9)
      if (user_opts%std_csv_output)   call close_output_csv
    case default
      write(message,*)"flag supplied to handle_output:",flag," was not recognised"
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end select

  end subroutine handle_output
  !
  !----------------------------------------------------------------------
  !
  subroutine start_spa( initialwater )

    ! Read the user config, then call the relevant routines !
    ! to open files and do any initialisation required.     !

    use config_tools,          only: ConfigSection
    use gv_scale_declarations, only: met, user_opts
    use log_tools
    use spa_cmd_line,          only: parse_spa_cmd_line
    use spa_config,            only: read_user_config, update_parameters_from_user_config
    use spa_initialise
    use spa_io_csv
#ifdef USE_NETCDF
    use spa_io_netcdf
    use spa_restart
#endif

    implicit none

    ! arguments.. 
    real,intent(out) :: initialwater

    ! local variables..
    type(ConfigSection), pointer :: section_names

    ! Find out the config filename..
    call parse_spa_cmd_line( spa_config_filename )

    ! Read user's config file..
#ifdef USE_NETCDF
    allocate( met_nc, met_nc%header )
    call read_user_config( spa_config_filename , user_opts , section_names , met_nc )
#else
    call read_user_config( spa_config_filename , user_opts , section_names )
#endif

    ! Open the meteorology file..(but don't load any data from it yet)
    if ( user_opts%met_file_is_nc ) then
#ifdef USE_NETCDF
       ! load basic data associated with file...
       call write_log( "opening the NC meteorology input file" , msg_info , __FILE__ , __LINE__ )
       call open_met_nc( met_nc%header , met_nc%time , met_nc%lat , met_nc%lon )

       ! ! Check latitude and longitude desired by user are inside bounds of this file..
       ! call write_log( "checking it contains the desired lat-lon" , msg_info , __FILE__ , __LINE__ )
       ! call check_LatLon( met_nc%lat%values, met_nc%lon%values, met_nc%grid )

       ! Load all the meteorological data into the pointer..
       call load_met_nc_data( met_nc )
       call write_log_div
#else
       write(message,*) "SPA was NOT compiled with the Netcdf library, so cannot open the NC " &
                      //"meteorology input file.  See the Makefile for the switch to change this." 
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
#endif
    else

       ! Load all the meteorological data into the pointer..
       call write_log( "opening the CSV meteorology input file" , msg_info , __FILE__ , __LINE__ )
       call read_met_csv( user_opts%met_filename , met )
       call write_log_div

    endif

    if ( .not. user_opts%loop_met_driver ) then
      ! check the met-driver file has enough timeslices to perform the run..
      call write_log( "checking met file contains enough time-slices for run length" , msg_info , __FILE__ , __LINE__ )
      call check_enough_timeslices( met )
    endif

#ifdef USE_NETCDF
    ! Setup lists of SPA variables/dimensions
    ! (useful in various input/output files)..
    call setup_spa_lists
#endif

    ! How are we starting model?..
    if ( user_opts%load_restart ) then

#ifdef USE_NETCDF
      ! Restarting from a previous dump..
      call write_log( "Starting SPA from restart file." , msg_info , __FILE__ , __LINE__ )
      call restart_model( user_opts%restart_filename )
#else
       write(message,*) "SPA was NOT compiled with the Netcdf library, so cannot create or read" &
                      //" from restart files.  See the Makefile for the switch to change this." 
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
#endif

    else

      ! Standard initialisation..

      ! Read the soils file..    
      call write_log( "opening the soils file" , msg_info , __FILE__ , __LINE__ )
      call read_soil_csv( user_opts%soils_filename, user_opts%veg_is_deciduous )

      ! Read the veg file..
      call write_log( "opening the veg file" , msg_info , __FILE__ , __LINE__ )
      call read_veg_csv( user_opts%veg_filename , user_opts%veg_is_deciduous )

      ! Initialise the soils..
      call write_log( "Initialising the soils" , msg_info , __FILE__ , __LINE__ )
      call initialise_soils( initialwater )

      ! Initialise the veg.. (On day 1 set leaf WP in each layer)
      call write_log( "Initialising the veg" , msg_info , __FILE__ , __LINE__ )
      call initialise_veg

    endif

    ! Make sure parameters read in from the user-config overwrite any in the files/default values..
    call write_log( "Reading user-defined params from user-config" , msg_info , __FILE__ , __LINE__ )
    call update_parameters_from_user_config( section_names )

    ! Open output files..
    call write_log( "Opening output files" , msg_info , __FILE__ , __LINE__ )
    call handle_output( 0 )

  end subroutine start_spa
  !
  !----------------------------------------------------------------------
  !
  subroutine update_met_drivers( time )

    ! wrapper to call appropriate nc/csv routine !
    ! for reading next slice of meteorology data !

    use gv_clim,               only: atmos_press, coa, daypar, dayppt, par_top, ppt, sw_rad, temp_bot, &
                                     temp_top, vpd_bot, vpd_top, wdbot, wdtop, wetev, wind_spd
    use gv_hourscale,          only: discharge, runoff
    use gv_scale_declarations, only: met, time_holder, user_opts
    use gv_veg,                only: flux, gppt, respt, totass, totevap, totres, transt
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time 

    ! check we haven't been called already this step...
    if ( time%step .ne. prev_step_of_day ) then

      ! iterate forward..
      met_slice = met_slice + 1

      if ( user_opts%loop_met_driver ) then
        ! if exceeding length of data then reset to start of data..
        if ( met_slice .gt. size(met%temp,1) ) met_slice = 1
      endif

      ! On first hour of day, set counters to zero..
      if ( time%step .eq. 1 ) then
        daypar  = 0. ; dayppt  = 0. ; discharge = 0. ; flux   = 0.
        gppt    = 0. ; respt   = 0. ; runoff    = 0. ; totass = 0.
        totevap = 0. ; totres  = 0. ; transt    = 0. ; wetev  = 0. 
      endif

      ! load this time-step's driver into each var..
      atmos_press = met%sfc_pressure( met_slice )
      coa         =  met%ambient_co2( met_slice )
      par_top     =          met%par( met_slice )
      ppt         =       met%precip( met_slice )
      sw_rad      =       met%sw_rad( met_slice )
      temp_bot    =         met%temp( met_slice )
      temp_top    = temp_bot
      vpd_bot     =          met%vpd( met_slice )
      vpd_top     = vpd_bot
      wind_spd    =     met%wind_spd( met_slice )

      ! Calculate the absolute water deficit (g m-3)
      wdtop = vpd_top * 217. / ( 0.1 * ( temp_top + 273.4 ) )
      wdbot = vpd_bot * 217. / ( 0.1 * ( temp_bot + 273.4 ) )

      !--- The following are only used for output files.. ---

      ! keep a (daily) running total of the precip..
      dayppt = dayppt + ppt

      ! sum the daily energy (PAR multiplied by nos of seconds per timestep)..
      daypar = daypar + par_top * time%seconds_per_step

      prev_step_of_day = time%step

      call write_log( "Loaded met data for step" , msg_info , __FILE__ , __LINE__ )

    endif

  end subroutine update_met_drivers
  !
  !----------------------------------------------------------------------
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !----------------------------------------------------------------------
  !
  subroutine check_enough_timeslices( met )

    ! Check that there are enough timeslices in the driver  !
    ! file to run the model for the number of days the user !
    ! wants to run for.                                     !

    use gv_scale_declarations, only: met_drivers, time
    use log_tools

    implicit none

    ! arguments..
    type(met_drivers),intent(in) :: met

    ! local variables..
    integer :: nos_of_timeslices

    ! nos of available timeslices = length of time dimension.
    ! nos of required timeslices  = user_run_length in days * nos of steps per day (we read at every step)

    nos_of_timeslices = size( met%temp )

    if ( nos_of_timeslices .lt. time%nos_of_years * time%days_in_year * time%steps_per_day ) then
      write(message,*)"There are not enough timeslices in the input file to run the model "//&
                      "for the number of days declared in the config file (default is 365)" 
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

  end subroutine check_enough_timeslices
  !
  !----------------------------------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine check_LatLon( file_lats , file_lons , file_grid )

    ! Check the lat/lon defined in the user-config falls within !
    ! that of the input file. If they are, retrieve the indices !
    ! of the points that encompass the user's desired location. !

    use gv_scale_declarations, only: grid
    use log_tools
    use netcdf_tools
    use spa_io_netcdf,         only: gridcalc

    implicit none

    ! arguments..
    real,dimension(:),intent(in) :: file_lats, file_lons
    type(gridcalc),pointer       :: file_grid

    ! local variables..

    logical :: status 
    integer :: tmp(1), grid_index
    real    :: spacing, value

    status = is_value_within_dim( file_lats, grid%latitude )
    if ( .not. status ) then
      write(message,*) "Requested latitude lies outside of bounds of input file latitudes!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    status = is_value_within_dim( file_lons, grid%longitude )
    if ( .not. status ) then
      write(message,*)"Requested longitude lies outside of bounds of input file longitudes!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! If we are still going, the required lat/lon must lie within the
    ! bounds of the dimensions, so find where..

    ! Check if the user has picked a point on the longitude grid...
    ! Find the point closest..
    tmp = minloc( file_lons - grid%longitude ) 
    grid_index = tmp(1)
    value = file_lons( grid_index )
    spacing = file_lons(2) - file_lons(1) ! assume constant grid spacing
    if ( value .eq. grid%longitude ) then
       ! set this as our index..
       file_grid%indices(1) = grid_index
       file_grid%no_bilinear = .true.
    else
       ! if not, find the nearest point just less than it...
       tmp = minloc( abs( file_lons - grid%longitude + 0.5*spacing ) )
       file_grid%indices(1) = tmp(1)
    endif

    ! Check if the user has picked a point on the latitude grid...
    tmp = minloc( file_lats - grid%latitude )
    grid_index = tmp(1)
    value = file_lats( grid_index )
    spacing = file_lats(2) - file_lats(1) ! assume constant grid spacing
    if ( ( value - grid%latitude ) .eq. grid%latitude ) then
       file_grid%indices(2) = grid_index
    else
       ! we only avoid bi-linear if both are true...
       file_grid%no_bilinear = .false.
       tmp = minloc( abs( file_lats - grid%latitude + 0.5*spacing ) )
       file_grid%indices(2) = tmp(1)
    endif

  end subroutine check_LatLon
  !
  !----------------------------------------------------------------------
  !
  logical function is_value_within_dim( dim , point )

    ! check whether a point lies within the bounds of a dimension !

    ! arguments..
    real,intent(in) :: dim(:), point

    ! check if value lies within it...
    if ( ( minval(dim) .lt. point ) .and. ( maxval(dim) .gt. point ) ) then
       is_value_within_dim = .true.
    else 
       is_value_within_dim = .false.
    endif

  end function is_value_within_dim
  !
  !----------------------------------------------------------------------
  !
  subroutine setup_spa_lists()

    ! Sets up a common structure which can be used by !
    ! all of SPA's input/output routines, containing  !
    ! information (and pointers to) SPA variables.    !

    use gv_clim,               only: daytempsum, gdd, max_fol
    use gv_hourscale,          only: canopy_store, hourts
    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid, time
    use gv_snow_info,          only: snowheight, snowweight
    use gv_soil_structure       ! loads!!
    use gv_veg                  ! loads!!
    use linked_lists,       only: append, new_item, spa_dims, spa_vars, spaitem, print_item
    use log_tools
    use netcdf_tools

    implicit none

    ! Build a list of the dimensions in spa..
    if ( .not. associated( spa_dims ) ) allocate( spa_dims )
    call append( spa_dims , new_item( 'grid%core'          , grid%core          ) )
    call append( spa_dims , new_item( 'nos_canopy_layers'  , grid%canopy        ) )
    call append( spa_dims , new_item( 'nos_soil_layers'    , grid%soil          ) )
    call append( spa_dims , new_item( 'nos_wetting_laters' , grid%wetting       ) )
    call append( spa_dims , new_item( 'steps_per_day'      , time%steps_per_day ) )
    call append( spa_dims , new_item( 'time (years)'       , time%year          ) )

    ! Build a list of the variables in spa..
    if ( .not. associated( spa_vars ) ) allocate( spa_vars )  
    call append( spa_vars , new_item( 'alloc_from_labile'                 , alloc_from_labile                 ) )
    call append( spa_vars , new_item( 'alloc_to_foliage'                  , alloc_to_foliage                  ) )
    call append( spa_vars , new_item( 'alloc_to_labile'                   , alloc_to_labile                   ) )
    call append( spa_vars , new_item( 'alloc_to_roots'                    , alloc_to_roots                    ) )
    call append( spa_vars , new_item( 'alloc_to_stem'                , alloc_to_stem                ) )
    call append( spa_vars , new_item( 'avN'                               , avn                               ) )
    call append( spa_vars , new_item( 'canopy_height'                     , canopy_height                     ) )
    call append( spa_vars , new_item( 'canopy_store'                      , canopy_store                      ) )
    call append( spa_vars , new_item( 'capac'                             , capac                             ) )
    call append( spa_vars , new_item( 'conductivity'                      , conductivity                      ) )
    call append( spa_vars , new_item( 'daytempsum'                        , daytempsum                        ) )
    call append( spa_vars , new_item( 'decomposition'                     , decomposition                     ) )
    call append( spa_vars , new_item( 'decomposition_rate'                , decomposition_rate                ) )
    call append( spa_vars , new_item( 'dimen'                             , dimen                             ) )
    call append( spa_vars , new_item( 'drythick'                          , drythick                          ) )
    call append( spa_vars , new_item( 'frac_alloc_foliage'                , frac_alloc_foliage                ) )
    call append( spa_vars , new_item( 'frac_alloc_roots'                  , frac_alloc_roots                  ) )
    call append( spa_vars , new_item( 'frac_GPP_resp_auto'                , frac_GPP_resp_auto                ) )
    call append( spa_vars , new_item( 'frac_leafLoss_litter'              , frac_leafLoss_litter              ) )
    call append( spa_vars , new_item( 'GDD'                               , GDD                               ) )
    call append( spa_vars , new_item( 'GDD_thresh'                        , GDD_thresh                        ) )
    call append( spa_vars , new_item( 'GPP'                               , GPP                               ) )
    call append( spa_vars , new_item( 'gplant'                            , gplant                            ) )
    call append( spa_vars , new_item( 'hourts'                            , hourts                            ) )
    call append( spa_vars , new_item( 'iota'                              , iota                              ) )
    call append( spa_vars , new_item( 'kappaC'                            , kappaC                            ) )
    call append( spa_vars , new_item( 'kappaJ'                            , kappaJ                            ) )
    call append( spa_vars , new_item( 'lat'                               , lat                               ) )
    call append( spa_vars , new_item( 'litterfall_foliage'                , litterfall_foliage                ) )
    call append( spa_vars , new_item( 'litterfall_roots'                  , litterfall_roots                  ) )
    call append( spa_vars , new_item( 'litterfall_stem'              , litterfall_stem              ) )
    call append( spa_vars , new_item( 'LMA'                               , LMA                               ) )
    call append( spa_vars , new_item( 'max_depth'                         , max_depth                         ) )
    call append( spa_vars , new_item( 'max_fol'                           , max_fol                           ) )
    call append( spa_vars , new_item( 'max_stock_foliage'                 , max_stock_foliage                 ) )
    call append( spa_vars , new_item( 'max_storage'                       , max_storage                       ) )
    call append( spa_vars , new_item( 'minlwp'                            , minlwp                            ) )
    call append( spa_vars , new_item( 'min_Temp_thresh'                   , min_Temp_thresh                   ) )
    call append( spa_vars , new_item( 'resp_cost_labile_trans'            , resp_cost_labile_trans            ) )
    call append( spa_vars , new_item( 'resp_rate_temp_coeff'              , resp_rate_temp_coeff              ) )
    call append( spa_vars , new_item( 'resp_auto'                         , resp_auto                         ) )
    call append( spa_vars , new_item( 'resp_h_litter'                     , resp_h_litter                     ) )
    call append( spa_vars , new_item( 'resp_h_soilOrgMatter'              , resp_h_soilOrgMatter              ) )
    call append( spa_vars , new_item( 'rooted_layers'                     , rooted_layers                     ) )
    call append( spa_vars , new_item( 'root_k'                            , root_k                            ) )
    call append( spa_vars , new_item( 'root_radius'                       , root_radius                       ) )
    call append( spa_vars , new_item( 'rootresist'                        , rootresist                        ) )
    call append( spa_vars , new_item( 'snowheight'                        , snowheight                        ) )
    call append( spa_vars , new_item( 'snowweight'                        , snowweight                        ) )
    call append( spa_vars , new_item( 'stock_foliage'                     , stock_foliage                     ) )
    call append( spa_vars , new_item( 'stock_labile'                      , stock_labile                      ) )
    call append( spa_vars , new_item( 'stock_litter'                      , stock_litter                      ) )
    call append( spa_vars , new_item( 'stock_soilOrgMatter'               , stock_soilOrgMatter               ) )
    call append( spa_vars , new_item( 'stock_resp_auto'                   , stock_resp_auto                   ) )
    call append( spa_vars , new_item( 'stock_roots'                       , stock_roots                       ) )
    call append( spa_vars , new_item( 'stock_stem'                   , stock_stem                   ) )
    call append( spa_vars , new_item( 'through_fall'                      , through_fall                      ) )
    call append( spa_vars , new_item( 'tower_height'                      , tower_height                      ) )
    call append( spa_vars , new_item( 'turnover_rate_foliage'             , turnover_rate_foliage             ) )
    call append( spa_vars , new_item( 'turnover_rate_labile'              , turnover_rate_labile              ) )
    call append( spa_vars , new_item( 'mineralisation_rate_litter'        , mineralisation_rate_litter        ) )
    call append( spa_vars , new_item( 'turnover_rate_resp_auto'           , turnover_rate_resp_auto           ) )
    call append( spa_vars , new_item( 'turnover_rate_roots'               , turnover_rate_roots               ) )
    call append( spa_vars , new_item( 'mineralisation_rate_soilOrgMatter' , mineralisation_rate_soilOrgMatter ) )
    call append( spa_vars , new_item( 'turnover_rate_stem'           , turnover_rate_stem           ) )

    ! vector valued items..
    call append( spa_vars , new_item( 'conduc'         , conduc         , 'grid%core'    ) )
    call append( spa_vars , new_item( 'iceprop'        , iceprop        , 'grid%core'    ) )
    call append( spa_vars , new_item( 'lafrac'         , lafrac         , 'grid%canopy'  ) )
    call append( spa_vars , new_item( 'lai'            , lai            , 'grid%canopy'  ) )
    call append( spa_vars , new_item( 'layer_depth'    , layer_depth    , 'grid%core'    ) )
    call append( spa_vars , new_item( 'layer_height'   , layer_height   , 'grid%canopy'  ) )
    call append( spa_vars , new_item( 'LWPstore'       , LWPstore       , 'grid%canopy'  ) )
    call append( spa_vars , new_item( 'mineralfrac'    , mineralfrac    , 'grid%core'    ) )
    call append( spa_vars , new_item( 'nfrac'          , nfrac          , 'grid%canopy'  ) )
    call append( spa_vars , new_item( 'organicfrac'    , organicfrac    , 'grid%core'    ) )
    call append( spa_vars , new_item( 'root_length'    , root_length    , 'grid%core'    ) )
    call append( spa_vars , new_item( 'root_mass'      , root_mass      , 'grid%core'    ) )
    call append( spa_vars , new_item( 'thickness'      , thickness      , 'grid%core'    ) )
    call append( spa_vars , new_item( 'soil_frac_clay' , soil_frac_clay , 'grid%core'    ) )
    call append( spa_vars , new_item( 'soil_frac_sand' , soil_frac_sand , 'grid%core'    ) )
    call append( spa_vars , new_item( 'soil_temp'      , soil_temp      , 'grid%core'    ) )
    call append( spa_vars , new_item( 'waterfrac'      , waterfrac      , 'grid%core'    ) )
    call append( spa_vars , new_item( 'watericemm'     , watericemm     , 'grid%soil'    ) )
    call append( spa_vars , new_item( 'wettingbot'     , wettingbot     , 'grid%wetting' ) )
    call append( spa_vars , new_item( 'wettingtop'     , wettingtop     , 'grid%wetting' ) )
      
  end subroutine setup_spa_lists
#endif
  !
  !----------------------------------------------------------------------
  !
end module spa_io
!
!------------------------------------------------------------------------
!
