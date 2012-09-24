! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


program main_spa

  !! Soil Plant Atmosphere model coupled to DALEC !!

  use allocate_carbon,       only: roots
  use canopy,                only: timestep_calcs
  use gv_scale_declarations, only: time
  use log_tools
  use spa_io,                only: handle_output, start_spa, update_met_drivers

  implicit none

  integer            :: yr, day, step
  real               :: iwater

  ! start logging..
  call open_log( unit=100 , fname="spa.log" )

  ! read user config, open files & initialise (spa_io.f90)
  call start_spa( iwater )

  ! for each year..
  do yr = 1 , time%nos_of_years

     call increment_time( 'year' , time )

     ! for each day...
     do day = 1 , time%days_in_year

        call increment_time( 'day' , time)

        ! Update root structure..
        call write_log( 'Updating root structure.' , msg_info , __FILE__ , __LINE__ )
        call roots

        ! for each sub-daily slice...
        do step = 1 , time%steps_per_day

           call increment_time( 'step' , time )

           ! get met driver for step (spa_io.f90)
           call write_log( 'Get next chunk of met data' , msg_info , __FILE__ , __LINE__ )
           call update_met_drivers( time )

           ! transfer carbon and water (canopy.f90)
           call write_log('Entering the timestep calculations' , msg_info , __FILE__ , __LINE__ )
           call timestep_calcs( time )

           ! write output if needed (spa_io.f90)
           call write_log('Dealing with any output (if needed)' , msg_info , __FILE__ , __LINE__ )
           call handle_output( 1 , time , iwater=iwater )

           call write_log_div

        enddo

        call write_log('              ' , msg_info , __FILE__ , __LINE__ )

        call write_log_div

     enddo

     ! End of year, sync the log..
     call sync_log

#ifdef USE_NETCDF
     ! check if it is time to write a restart file...
     call check_restart
#endif

  enddo

#ifdef USE_NETCDF
  ! Check that, if restarts were wanted, we wrote at least one!..
  call check_restart( end_of_run = .True. )
#endif

  ! finish logging..
  call close_log( 0 )

contains
  !
  !-------------------------------------------
  !
  subroutine increment_time( type , time )

    ! Increment parts of the time-holder variable. !
  
    use gv_scale_declarations, only: time_holder
    use log_tools,             only: write_log

    implicit none

    ! arguments..
    character(*),intent(in)         :: type
    type(time_holder),intent(inout) :: time

    select case (type)
    case ('year')
      ! next year, reset day & step..
      time%year = time%year + 1
      time%day  = 0
      time%step = 0
      write(message,*)'year is: ',time%year
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('day')
      ! next day, reset step...
      time%day  = time%day + 1
      time%step = 0
      write(message,*)'day is: ',time%day
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('step')
      ! next step.
      ! Also update the count of total-steps, and the daytime..
      time%step = time%step + 1
      write(message,*)'step-of-day is: ',time%step
      call write_log( trim(message), msg_info  , __FILE__ , __LINE__ )
      time%steps_count = time%steps_count + 1
      write(message,*)'number of steps so far is: ',time%steps_count
      call write_log( trim(message), msg_info , __FILE__ , __LINE__ )

      ! Re-calculate the current time in units of days..
      ! Consider step 1 to be at 00:00 on the day, and the
      ! last step of the day to be at or before 23.59...
      time%daytime = ( time%year -1 ) * time%days_in_year    &
                      + time%day                                  &
                       + ( real(time%step) - 1 ) / real(time%steps_per_day)
      write(message,*)'daytime is: ',time%daytime
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case default
      write(message,*)'Do not recognise type of update: ',type
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    endselect

  end subroutine increment_time
  !
  !-------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine check_restart( end_of_run )

    ! Check if it is time to write a restart file !

    use gv_scale_declarations, only: time, user_opts
    use log_tools
    use spa_restart

    implicit none

    ! arguments..
    logical, intent(in), optional :: end_of_run

    ! local variables..
    logical :: final_call

    ! If user didn't want dumps, then exit the s/r now..
    if ( .not. user_opts%make_restart_file ) return

    if (present(end_of_run)) then
      final_call = end_of_run
    else
      final_call = .false.
    endif

    if ( final_call ) then
      ! check that we wrote at least one restart dump.
      if ( prev_restart_write .eq. 0. ) then
        ! Not written yet. Do write now..
        call write_restart_file( time , user_opts%veg_is_deciduous )
        ! And tell user..
        write(message,*) "The run completed before any restart files" &
              //" were written. A restart file has been written now."
        call write_log( message , msg_warning , __FILE__ , __LINE__ )
      end if      
    else
      ! Check if it is time to write a restart file yet..
      if ( time%year .eq. prev_restart_write + user_opts%restart_write_frequency ) then
         call write_restart_file( time , user_opts%veg_is_deciduous )
         prev_restart_write = time%year
      end if
    end if

  end subroutine check_restart
#endif
  !
  !-------------------------------------------
  !
end program main_spa
