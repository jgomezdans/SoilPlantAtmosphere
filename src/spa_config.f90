! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_config

  !! This module reads the user-provided  !!
  !! runtime-configuration file for SPA.  !!

  use gv_scale_declarations, only: fname_length
#ifdef USE_NETCDF
  use spa_io_netcdf ! for the nc_met/nc_soils/nc_veg type definitions
#endif

  implicit none


  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: read_user_config, update_parameters_from_user_config

  interface check_allocated
    module procedure check_allocated_1d, check_allocated_2d
  end interface check_allocated

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine read_user_config( config_filename , config_options , section_names , met_nc )
#else
  subroutine read_user_config( config_filename , config_options , section_names )
#endif

    ! reads the user config, looking for the filenames !
    ! output dir and SPA configuration flags.          !

    use config_tools
    use gv_scale_declarations, only: fname_length, grid, met, time, user_config_holder
    use log_tools
#ifdef USE_NETCDF
    use netcdf_tools
#endif

    ! arguments..
    character(fname_length),intent(in) :: config_filename !(input)  where to look for the config file
    type(user_config_holder)           :: config_options  !(output) holds the user's config choices.
    type(ConfigSection), pointer       :: section_names   !(output) structure holding sections of the configuration file   
#ifdef USE_NETCDF
    type(nc_met_file),pointer          :: met_nc          !(output)
#endif

    ! local variables..
    type(ConfigSection), pointer :: section
    integer                      :: ios
    character(len=100)           :: outdirtestfile

    ! Get the section names..
    call ConfigRead( config_filename, section_names )

    ! print the user's config to the logfile
    call write_log( "The contents of the configuration file are.." , msg_info , __FILE__ , __LINE__ )
    call PrintConfig( section_names, get_logunit() )
    call write_log_div

    ! Get the General options...
    call GetSection(section_names,section,'Options')
    if ( associated(section) ) then
      call GetValue(section,'veg_is_deciduous'  ,     config_options%veg_is_deciduous )
      if (config_options%veg_is_deciduous) then
        call write_log( "Deciduous simulation." , msg_info , __FILE__ , __LINE__ )
      else
        call write_log( "Evergreen simulation." , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'start_from_restart',     config_options%load_restart    )
      call GetValue(section,'layers_in_canopy',       grid%canopy         )
      ! sanity check..
      if ( grid%canopy .lt. 1 ) then
        call write_log( "Must be at least one canopy layer!" , msg_fatal , __FILE__ , __LINE__ )
      elseif ( grid%canopy .gt. 20 ) then
        call write_log( "Use of more than 20 canopy layers is not recommended!" , msg_warning , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'layers_in_soil',         grid%soil           )
      ! sanity check..
      if ( grid%soil .lt. 3 ) then
        call write_log( "Must be at least 3 soil layers!" , msg_fatal , __FILE__ , __LINE__ )
      elseif ( grid%soil .gt. 40 ) then
        call write_log( "Use of more than 40 soil layers is not recommended!" , msg_warning , __FILE__ , __LINE__ )
      endif
      grid%core = grid%soil + 1
      call GetValue(section,'latitude_deg_north',     grid%latitude                  )
      call GetValue(section,'longitude_deg_east',     grid%longitude                 )
      call GetValue(section,'number_of_years',        time%nos_of_years              )
      call GetValue(section,'days_per_year',          time%days_in_year              )
      call GetValue(section,'steps_per_day',          time%steps_per_day             )
      ! sanity check..
      if ( time%steps_per_day .lt. 24 ) then
        call write_log( "Use of less than 24 steps per day is not recommended (resolving the dirunal cycle is important!)" , &
                         msg_warning , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'loop_met_driver',       config_options%loop_met_driver           )
      call GetValue(section,'log_messages_to_screen:', config_options%print_msg_at_severity )
      call set_msg_level( config_options%print_msg_at_severity )
    else
      write(message,*)"Could not find [Options] section in the config file.  Cannot continue without it!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Get the initialisation input files...
    call GetSection(section_names,section,'Initialisation')
    if ( associated(section) ) then
      if ( config_options%load_restart ) then
        ! Restart model..
        call GetValue(section,'restart',config_options%restart_filename )
        if (config_options%restart_filename=='') then
          write(message,*)'Restart filename must be supplied in the user-config file.'
          call write_log( message , msg_fatal, __FILE__ , __LINE__ )
        endif
      else
        ! Standard initialisation..
        ! soils
        call GetValue(section,'soils',config_options%soils_filename)
        if (config_options%soils_filename=='') then
          write(message,*)'Soils filename must be supplied in the user-config file.'
          call write_log( message , msg_fatal, __FILE__ , __LINE__ )
        endif
        ! vegetation
        call GetValue(section,'vegetation',config_options%veg_filename)
        if (config_options%veg_filename=='') then
          write(message,*)'Vegetation filename must be supplied in the user-config file.'
          call write_log( message , msg_fatal , __FILE__ , __LINE__ )
        endif
      endif
    else
      write(message,*)"Could not find [Initialisation] section in the config file."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Meteorology driver
    call GetSection( section_names , section , 'Meteorology driver' )
    if ( associated(section) ) then
      call GetValue( section , 'filename' , config_options%met_filename )
      if (config_options%met_filename=='') then
        write(message,*) "Meteorology filename must be supplied in the user-config file."
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      call GetValue( section , 'file_is_netcdf' , config_options%met_file_is_nc )
      call GetValue( section , 'file_has_co2' ,  config_options%use_co2_from_met_file  )
      if ( .not. config_options%use_co2_from_met_file ) then 
        call write_log( "Using a constant CO2 concentration." , msg_info , __FILE__ , __LINE__ )
        call GetValue( section , 'constant_CO2_concentration' , met%const_ambient_co2 )
      endif
      call GetValue( section , 'file_has_par' ,  config_options%use_par_from_met_file  )
      if ( .not. config_options%use_par_from_met_file ) then
        call write_log( "PAR will be calculated from SW-radiation." , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue( section , 'file_has_sfc_pressure' , config_options%met_file_has_sfc_press )
#ifdef USE_NETCDF
      if ( config_options%met_file_is_nc ) then
        allocate( met_nc%header )
        met_nc%header%name = config_options%met_filename
        call get_met_nc_varnames_from_config( section , met_nc )
      end if
#endif
      call GetValue( section , 'precip_is_rate' , config_options%precip_is_rate )
    else
      write(message,*)"Could not find [Meteorology driver] section in the config file."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Get the general output detail..
    call GetSection(section_names,section,'Output')
    if (associated(section)) then
      call GetValue( section , 'directory'               , config_options%output_directory        )
      ! check the directory exists by trying to create a file in it..
      outdirtestfile = trim(config_options%output_directory)//'spa_test_dir_exists'
      open( iostat=ios , unit=999 , file=outdirtestfile , status='new' )
      close( 999 , status='delete' )
      if ( ios .ne. 0 ) then
        write(message,*)"Cannot find output directory: ",trim(config_options%output_directory)
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      call GetValue( section , 'standard_csv_output'     , config_options%std_csv_output          )
      call GetValue( section , 'write_restart_file'      , config_options%make_restart_file       )
      call GetValue( section , 'restart_write_frequency' , config_options%restart_write_frequency )
    end if

    call allocate_arrays

  end subroutine read_user_config
  !
  !----------------------------------------------------------------------------
  !
  subroutine update_parameters_from_user_config( section_names )

    ! make sure any parameters provided by the user overwrite !
    ! the values taken from the input files.                  !

    use config_tools
    use gv_scale_declarations, only: met, user_opts
    use gv_snow_info,          only: snowheight, snowweight
    use gv_soil_structure,     only: root_radius

    implicit none

    ! arguments..
    type(ConfigSection), pointer :: section_names ! holds structure of the configuration file   

    ! local variables..
    type(ConfigSection), pointer :: section

    ! For each section, we look for possible parameters.
    ! GetValue only changes the third argument IF it
    ! finds something.

    call GetSection( section_names , section , 'Met Parameters' )
    if (associated(section)) then
    endif

    call GetSection(section_names,section,'Soil Parameters')
    if (associated(section)) then
       call GetValue(section,'snowweight',snowweight)
       call GetValue(section,'snowheight',snowheight)
    endif

    call GetSection(section_names,section,'Veg Parameters')
    if (associated(section)) then
       call GetValue(section,'root_radius',root_radius)
    endif

  end subroutine update_parameters_from_user_config
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine allocate_arrays( )

    ! The size of some arrays depends on grid-size, which !
    ! can be set at runtime. These arrays must then be    !
    ! allocated at the start of a run, rather than being  !
    ! hard-coded. This routine allocates them.            !

    use gv_clim,               only: gbw, wetev
    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid, time
    use gv_soil_structure
    use gv_veg

    implicit none

    ! clim..
    call check_allocated( gbw   , grid%canopy )
    call check_allocated( wetev , time%steps_per_day )

    ! hydrol..
    call check_allocated( soil_frac_clay , grid%core )
    call check_allocated( soil_frac_sand , grid%core )

    ! soil..
    if ( .not. allocated(cond1) ) then
      ! (it's fair to assume that if one hasn't been, none of them have)..
      allocate( cond1(grid%core), cond2(grid%core), cond3(grid%core), potA(grid%core), potB(grid%core) )
    endif
    if ( .not. allocated(           conduc ) )  allocate(           conduc( grid%core ) )
    if ( .not. allocated(   field_capacity ) )  allocate(   field_capacity( grid%core ) )
    if ( .not. allocated(  fraction_uptake ) )  allocate(  fraction_uptake( grid%core ) )
    if ( .not. allocated(          iceprop ) )  allocate(          iceprop( grid%core ) )
    if ( .not. allocated(      layer_depth ) )  allocate(      layer_depth( grid%core ) )
    if ( .not. allocated(      mineralfrac ) )  allocate(      mineralfrac( grid%core ) )
    if ( .not. allocated(      organicfrac ) )  allocate(      organicfrac( grid%core ) )
    if ( .not. allocated(         porosity ) )  allocate(         porosity( grid%core ) )
    if ( .not. allocated(          pptgain ) )  allocate(          pptgain( grid%core ) )
    if ( .not. allocated(      root_length ) )  allocate(      root_length( grid%core ) )
    if ( .not. allocated(        root_mass ) )  allocate(        root_mass( grid%core ) )
    if ( .not. allocated(        soil_temp ) )  allocate(        soil_temp( grid%core ) )
    if ( .not. allocated( soil_temp_nplus1 ) )  allocate( soil_temp_nplus1( grid%core ) )
    if ( .not. allocated(            soilR ) )  allocate(            soilR( grid%core ) )
    if ( .not. allocated(           soilR1 ) )  allocate(           soilR1( grid%core ) )
    if ( .not. allocated(           soilR2 ) )  allocate(           soilR2( grid%core ) )
    if ( .not. allocated(              SWP ) )  allocate(              SWP( grid%core ) )
    if ( .not. allocated(        thickness ) )  allocate(        thickness( grid%core ) )
    if ( .not. allocated(        waterfrac ) )  allocate(        waterfrac( grid%core ) )
    if ( .not. allocated(        watergain ) )  allocate(        watergain( grid%core ) )
    if ( .not. allocated(       watericemm ) )  allocate(       watericemm( grid%soil ) )
    if ( .not. allocated(        waterloss ) )  allocate(        waterloss( grid%core ) )
    if ( .not. allocated( wettingbot ) )  then
      allocate( wettingbot( grid%wetting ) ) 
      allocate( wettingtop( grid%wetting ) )
    endif

    ! veg..
    if ( .not. allocated( canopy_soil_resistance ) )  allocate( canopy_soil_resistance( grid%canopy )    )
    if ( .not. allocated(                    ess ) )  allocate( ess( time%steps_per_day )                )
    if ( .not. allocated(                   gppt ) )  allocate( gppt( time%steps_per_day )      )
    if ( .not. allocated(                 lafrac ) )  allocate( lafrac( grid%canopy )           )
    if ( .not. allocated(                    lai ) )  allocate( lai( grid%canopy )              )
    if ( .not. allocated(           layer_height ) )  allocate( layer_height( grid%core )       )
    if ( .not. allocated(               LWPstore ) )  allocate( LWPstore( grid%canopy )         )
    if ( .not. allocated(                  nfrac ) )  allocate( nfrac( grid%canopy )            )
    if ( .not. allocated(                    nla ) )  allocate( nla( grid%canopy )              )
    if ( .not. allocated(                  respt ) )  allocate( respt( time%steps_per_day )     )
    if ( .not. allocated(               soiletmm ) )  allocate( soiletmm( time%steps_per_day )  )
    if ( .not. allocated(                 transt ) )  allocate( transt( time%steps_per_day )    )

    call check_allocated( flux , varsize=(/time%steps_per_day,grid%canopy/) )

  end subroutine allocate_arrays
  !
  !----------------------------------------------------------------------
  !
  subroutine check_allocated_1d( var, varsize )

    use log_tools

    implicit none

    ! arguments..
    real,allocatable,intent(inout) :: var(:)
    integer,intent(in) :: varsize

    if ( .not. allocated( var ) ) then
      ! allocate it..
      allocate( var( varsize ) )
    else
      ! check it is actually the desired size..
       if ( size(var) .ne. varsize ) then
         write(message,*)"Variable not allocated to same size as stated varsize!"
         call write_log( message , msg_error , __FILE__ , __LINE__ )
       endif
     endif

  end subroutine
  !
  !----------------------------------------------------------------------
  !
  subroutine check_allocated_2d( var, varsize )

    use log_tools

    implicit none

    ! arguments..
    real,allocatable,intent(inout)  :: var(:,:)
    integer,dimension(:),intent(in) :: varsize

    ! local variables..
    integer :: varshape(2)

    if ( .not. allocated( var ) ) then
      ! allocate it..
      allocate( var( varsize(1),varsize(2) ) )
    else
      ! check it is actually the desired size..
       varshape = shape(var)
       if ( ( varshape(1) .ne. varsize(1) ) .or. ( varshape(2) .ne. varsize(2) ) ) then
         write(message,*)"Variable not allocated to same size as stated varsize!"
         call write_log( message , msg_error , __FILE__ , __LINE__ )
       endif
     endif

  end subroutine check_allocated_2d
  !
  !----------------------------------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine get_met_nc_varnames_from_config( section , met_nc )

    ! If the input met file is netcdf we need to check the !
    !  [Meteorology driver] section for details of what    !
    !  the variables are called in the met input file.     !

    use config_tools
    use gv_scale_declarations, only: met, user_opts
    use log_tools
    use spa_io_netcdf,         only: nc_met_file

    implicit none

    ! arguments..
    type(ConfigSection),pointer :: section
    type(nc_met_file),pointer   :: met_nc

    ! define default names for dimensions, and then go see if the user
    ! has given different ones...
    ! (I'm sure there should be a better way of doing this...)
    allocate( met_nc%time   ) ; met_nc%time%name = 'time'
    allocate( met_nc%lat    ) ; met_nc%lat%name  = 'lat'
    allocate( met_nc%lon    ) ; met_nc%lon%name  = 'lon'
    call GetValue( section, 'time',      met_nc%time%name   )
    call GetValue( section, 'latitude',  met_nc%lat%name    )
    call GetValue( section, 'longitude', met_nc%lon%name    )

    ! define default names for variables.. nd then go see if the user
    ! has given different ones...
    ! (I'm sure there should be a better way of doing this...)
    allocate( met_nc%coa    ) ; met_nc%coa%name    = 'co2'
    allocate( met_nc%par    ) ; met_nc%par%name    = 'par'
    allocate( met_nc%ppt    ) ; met_nc%ppt%name    = 'ppt'
    allocate( met_nc%sat    ) ; met_nc%sat%name    = 'temp'
    allocate( met_nc%sfc_p  ) ; met_nc%sfc_p%name  = 'pressure'
    allocate( met_nc%swrad  ) ; met_nc%swrad%name  = 'sw'
    allocate( met_nc%vpd    ) ; met_nc%vpd%name    = 'vpd'
    allocate( met_nc%windsp ) ; met_nc%windsp%name = 'windsp'
    call GetValue( section, 'air_temperature',         met_nc%sat%name      )
    call GetValue( section, 'temp_in_kelvin',          met_nc%sat_in_kelvin )
    if ( user_opts%use_co2_from_met_file ) then
      call GetValue( section, 'carbon_dioxide',        met_nc%coa%name      )
    else
      met_nc%coa%name = ''
    endif
    if ( user_opts%use_par_from_met_file ) then
      call GetValue( section, 'photo_active_rad',       met_nc%par%name     )
    else
      met_nc%par%name = ''
    endif
    call GetValue( section, 'precipitation',           met_nc%ppt%name      )
    if ( user_opts%met_file_has_sfc_press ) then
      call GetValue( section, 'surface_pressure',      met_nc%sfc_p%name    )
    else
      met_nc%sfc_p%name = ''
      call GetValue(section,'constant_sfc_pressure',   met%const_sfc_pressure )
    endif
    call GetValue( section, 'sw_radiation',            met_nc%swrad%name    )
    call GetValue( section, 'vapour_pressure_deficit', met_nc%vpd%name      )
    call GetValue( section, 'wind_speed',              met_nc%windsp%name   )

    call write_log("The variable names SPA is expecting to find in met nc file are:", &
                     msg_info, __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%time%name),' ',trim(met_nc%lat%name),' ',trim(met_nc%lon%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%coa%name),' ',trim(met_nc%par%name),' ',  &
                     trim(met_nc%ppt%name),' ',trim(met_nc%sat%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%sfc_p%name),' ',trim(met_nc%swrad%name),' ', &
                     trim(met_nc%vpd%name),' ',trim(met_nc%windsp%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )

  end subroutine get_met_nc_varnames_from_config
#endif
  !
  !----------------------------------------------------------------------
  !
end module spa_config
!
!----------------------------------------------------------------------
!
