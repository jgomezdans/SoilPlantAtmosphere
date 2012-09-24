! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_cmd_line

  !! Routines to handle command line input. !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: parse_spa_cmd_line

  integer :: cmd_line_len = 5000

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine parse_spa_cmd_line( config_filename )

    ! Looks through the command line, and does appropriate   !
    ! things with arguments etc. Returns name of config-file !

    use gv_scale_declarations, only: fname_length
    use log_tools

    implicit none

    ! arguments..
    character(len=fname_length),intent(out) :: config_filename 

    ! local variables..
    integer                     :: i, nos_of_args
    character(len=fname_length) :: argument, this_program

    ! Get the program name..
    call get_command_argument( 0 , this_program )

    ! Get the number of arguments..
    nos_of_args = command_argument_count()

    ! if user hasn't supplied any arguments at all, then ask for name..
    if ( nos_of_args .gt. 0 ) then

       ! loop over command line arguments
       i = 0
       do
          i = i + 1
          if ( i .gt. nos_of_args ) exit

          call get_command_argument( i , argument )

          ! check if it is an option
          if ( argument(1:1) .eq. '-' ) then
             select case ( trim( argument ) )
             case ('-h','--help') ! print help and exit..

                call print_cmd_line_help( this_program )

             case ('-v','--version') ! print help and exit..

                call print_spa_version( this_program )

             case default ! print help and exit..

                write(*,*) 'Unknown option ',trim( argument )
                call print_cmd_line_help( this_program )

             end select
          else
             ! it's not an option - assume it's our config file...
             call get_command_argument( 1 , config_filename )
             call write_log( "Will use runtime parameters from : "//trim(config_filename) )
          end if
       end do

    else ! user didn't supply any arguments, so...

       write(*,*)
       write(*,*) 'Enter name of configuration file to be read (type "exit" to quit):'
       read(*,'(a)') config_filename
       if ( trim(config_filename) .eq. "exit" ) then
         write(*,*) 'Exiting now.'
         call write_log( "User requested run stop." , msg_info , __FILE__ , __LINE__ )
         call close_log( 0 )
         stop
       else
         write(*,*) 'Thank you, continuing.'
       endif

    end if

  end subroutine parse_spa_cmd_line
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine print_cmd_line_help( this_program )

    ! Print a help statement !

    implicit none

    character(len=*),intent(in) :: this_program

    write(*,*) 'Usage: ',trim(this_program),' [options] configname'
    write(*,*) 'where [options] are'
    write(*,*) '  -h|--help       print this message'
    write(*,*) '  -v|--version    print the version number'

    stop "Stopped."

  end subroutine print_cmd_Line_help
  !
  !----------------------------------------------------------------------
  !
  subroutine print_spa_version( this_program )

    ! Print a help statement !

    use gv_scale_declarations, only: model_sci_version, model_bug_version 

    implicit none

    ! argument..
    character(len=*),intent(in) :: this_program

    ! local variable..
    character(len=6) :: version

    write(version,'(I1,".",I1,".",I0.1)') floor(model_sci_version/10.),&
                           modulo(model_sci_version,10),model_bug_version
    write(*,*) "This is version ",trim(version)," of ",trim(this_program),"."

    stop "Stopped."

  end subroutine print_spa_version
  !
  !----------------------------------------------------------------------
  !
end module spa_cmd_line
!
!------------------------------------------------------------------------
!
