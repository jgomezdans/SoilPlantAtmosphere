! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module spa_io_netcdf

  !! This module declares derived types describing SPA-specific  !!
  !! input/output files, and procedures for reading those files. !!

  use gv_scale_declarations, only: fname_length
  use netcdf_tools,          only: nc_header, nc_dimension, nc_variable

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: load_met_nc_data, open_met_nc
  ! Variables..
  public :: gridcalc, nc_met_file


  ! Required for finding particular spot on grid that we will use..
  type gridcalc
    integer,dimension(2) :: indices  = (1,1)   ! bottom-left corner of box of grid-pts
    ! that encompass user-desired location.
    logical           :: no_bilinear = .false. ! bilinear flag (assume false=>do calcs)
  end type gridcalc

  ! Input Meteorolgy file structure..
  type nc_met_file
    type(nc_header),pointer    :: header => null() ! information about the file, such as name & dims
    type(nc_dimension),pointer :: time   => null() ! time
    type(nc_dimension),pointer :: lat    => null() ! latitude
    type(nc_dimension),pointer :: lon    => null() ! longitude
    type(nc_variable),pointer  :: sat    => null() ! surface air temperature
    logical                    :: sat_in_kelvin = .False. ! assume SAT in Celcius
    type(nc_variable),pointer  :: coa    => null() ! carbon dioxide atmospheric concentration
    type(nc_variable),pointer  :: par    => null() ! photosyntheticaly active radiation
    type(nc_variable),pointer  :: ppt    => null() ! precipitation
    type(nc_variable),pointer  :: swrad  => null() ! shortwave radiation
    type(nc_variable),pointer  :: sfc_p  => null() ! surface atmospheric pressure
    type(nc_variable),pointer  :: vpd    => null() ! vapour pressure deficit
    type(nc_variable),pointer  :: windsp => null() ! wind speed
    type(gridcalc),pointer     :: grid   => null() ! point on grid that we will use
  end type nc_met_file

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine load_met_nc_data( ncfile )

    ! Read the met input file and load data to pointers !
    ! Unlike the csv file, which is read line-by-line as!
    ! it is needed, the netcdf file contents are read   !
    ! all in one go.                                    !

    use gv_scale_declarations, only: met, time, user_opts
    use log_tools
    use netcdf_tools,          only: get_nc_var

    implicit none

    ! arguments..
    type(nc_met_file)  :: ncfile

    ! local variables..
    integer :: time_length

    time_length = size(ncfile%time%values)

    ! Carbon dioxide..
    allocate( met%ambient_co2(time_length) )
    if ( user_opts%use_co2_from_met_file ) then   
      call write_log("Try to load Ambient CO2 concentration from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%coa )
      met%ambient_co2 = reshape( ncfile%coa%data%grid_3d , (/ time_length /) )
      call write_log("..successful.")
    else
      ! use a default value..
      met%ambient_co2 = met%const_ambient_co2
      call write_log("Ambient CO2 concentration set to default value (330ppm)")
    endif

    ! Precipitation
    allocate( met%precip(time_length) )
    call write_log("Try to load Precipitation from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%ppt )
    call write_log ("..successful.")
    met%precip = reshape( ncfile%ppt%data%grid_3d , (/ time_length /) )
    if ( user_opts%precip_is_rate ) then
      call write_log("Converting precipitation from rate to volume by "&
                   //"multiplying by nos of seconds per timestep.")
      met%precip = met%precip * time%seconds_per_step
    endif

    ! Short wave radiation
    allocate( met%sw_rad(time_length) )
    call write_log("Try to load SW radiation from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%swrad )
    call write_log("..successful.")
    met%sw_rad = reshape( ncfile%swrad%data%grid_3d  , (/ time_length /) )

    ! Surface air temperature
    allocate( met%temp(time_length) )
    call write_log("Try to load Surface air temperature from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%sat )
    ! Adjust to Kelvin, if needed, and also check that units seem sane..
    if ( ncfile%sat_in_kelvin ) then
      if ( minval( ncfile%sat%data%grid_3d ) .lt. 100. ) then
        write(message,*)"Units of temperature appear to be wrong. You "//&
             "specified Kelvin, but minimum is below 100K; is this correct??"
        call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
      endif
      ncfile%sat%data%grid_3d = ncfile%sat%data%grid_3d - 273.15 ! convert T to Celcius
    else
      if ( maxval( ncfile%sat%data%grid_3d ) .gt. 150. ) then
        write(message,*)"Units of temperature appear to be wrong. You "//&
             "specified Celcius, but maximum is above 150C; is this correct??"
        call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
      endif
    endif
    call write_log("..successful.")
    met%temp = reshape( ncfile%sat%data%grid_3d , (/ time_length /) )

    ! Surface pressure
    allocate( met%sfc_pressure(time_length) )
    if ( user_opts%met_file_has_sfc_press ) then
      call write_log("Try to load Surface pressure from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%sfc_p )
      call write_log("..successful.")
      met%sfc_pressure = reshape( ncfile%sfc_p%data%grid_3d , (/ time_length /) )
    else
      ! use a default value..
      met%sfc_pressure = met%const_sfc_pressure
      call write_log("Surface pressure set to default value (990mb)")
    endif

    ! Vapour pressure deficit
    allocate( met%vpd(time_length) )
    call write_log("Try to load vapour pressure deficit from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%vpd )
    call write_log("..successful.")
    met%vpd = reshape( ncfile%vpd%data%grid_3d , (/ time_length /) )

    ! Wind speed
    allocate( met%wind_spd(time_length) )
    call write_log("Try to load wind speed from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%windsp )
    call write_log("..successful.")
    met%wind_spd = reshape( ncfile%windsp%data%grid_3d , (/ time_length /) )
    ! adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2 )  met%wind_spd = 0.2

    ! Photosynthetically active radiation
    allocate( met%par(time_length) )
    if ( user_opts%use_par_from_met_file ) then
      call write_log("Try to load photosynthetically active radiation from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%par )
      call write_log("..successful.")
      met%par = reshape( ncfile%par%data%grid_3d , (/ time_length /) )
    else
      ! PAR isn't available in most climate models, so
      ! instead we use a loose approximation...
      met%par = 2.3 * met%sw_rad
    endif

  end subroutine load_met_nc_data
  !
  !----------------------------------------------------------------------
  !
  subroutine open_met_nc( header , time , lat , lon )

    ! Open a netcdf input file, and fill out the header with !
    ! information on handles to the file, and the dimensions !
    ! of time, latitude and longitude.                       !

    use netcdf_tools, only: get_dim_info, nc_header, nc_dimension, open_nc_file

    implicit none

    ! arguments..
    type(nc_header),pointer    :: header
    type(nc_dimension),pointer :: time, lat, lon

    ! open netCDF file..
    call open_nc_file( header )

    ! populate each dimension with their dim & var handle
    !  ids, their length and their actual values..
    call get_dim_info( header%id , time%name , time%dim_id , &
                                      time%var_id , time%values )
    call get_dim_info( header%id, lat%name , lat%dim_id , &
                                      lat%var_id , lat%values )
    call get_dim_info( header%id, lon%name , lon%dim_id , &
                                      lon%var_id , lon%values  )

  end subroutine open_met_nc
  !
  !----------------------------------------------------------------------
  !
end module spa_io_netcdf
!
!------------------------------------------------------------------------
!
