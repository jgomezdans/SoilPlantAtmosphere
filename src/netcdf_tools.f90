! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module netcdf_tools

  !! This module declares generic derived types and procedures    !!
  !! that allow us to more easily manipulate the data associated  !!
  !! with netcdf files.                                           !!
  !!                                                              !!
  !! The derived types allow us to declare a general 'type', such !!
  !! as nc_dimension, with a bunch of properties, such as 'name', !!
  !! and then build more complex things, such as the 'nc_met_file'!!
  !! which has a nc_header (containing several dimensions) and a  !!
  !! series of 'nc_variables..                                    !!

  use gv_scale_declarations, only: fname_length

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: change_mode, close_nc_output, create_coord_dim, &
            create_dim, create_var, get_dim_info,           &
            get_global_att, get_nc_var, open_nc_file,       &
            put_global_att, put_var_nc
  ! Variables..
  public :: nc_dimension, nc_header, nc_variable


  ! A generic dimension...
  type nc_dimension
    character(fname_length)   :: name     ! name of dimension
    integer                   :: dim_id = -999    ! dim id
    integer                   :: var_id = -999    ! var id
    real,dimension(:),pointer :: values => null() ! values of dimension
  end type nc_dimension

  ! Handles the various possibilities of a variable's value.
  type nc_data
    real                          :: scalar  = 0.
    real,dimension(:),pointer     :: vector  => null()
    real,dimension(:,:),pointer   :: grid_2d => null()
    real,dimension(:,:,:),pointer :: grid_3d => null()
  end type nc_data

  ! A generic variable
  ! Only one of the data types will be actually be used, controlled
  ! by the ndims.
  type nc_variable
    character(fname_length) :: name = ''       ! name of variable
    integer                 :: id   = -999     ! var id
    real                    :: mdi  = -1e32    ! missing data indicator
    integer                 :: ndims = -1.     ! nos of dims
    type(nc_data),pointer   :: data => null()  ! handles various sizes of data
  end type nc_variable

  ! A typical file header (comprising file info and dimensions)..
  type nc_header
     character(fname_length)    :: name        =  ''     ! file name
     integer                    :: id          =  -999   ! file handle
     logical                    :: define_mode = .False. ! indicates read/write (false) or define (true)
     logical                    :: exists      = .False. ! did the file exist before opening?
  end type nc_header

  ! --- interfaces to handle over-loading of s/r names ---
  interface get_dim_info
    module procedure get_dim_info, get_dim_info2
  end interface

  interface get_global_att
    module procedure get_int_global_att, get_real_global_att
  end interface

  interface get_nc_var
    module procedure get_nc_var, get_nc_var2
  end interface

  interface put_global_att
    module procedure put_int_global_att, put_real_global_att
  end interface

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine change_mode( nc_handle , define )

    ! Take a given ncdf handle and puts the file into define !
    ! (True) or write (False) mode.  As part of this, it     !
    ! first check what mode the file is currently in.        !

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    type(nc_header) :: nc_handle
    logical         :: define    ! True => want file in define mode
    ! False => want file in read/write mode

    ! local variables..
    integer :: status

    ! Only do something if the file's current mode
    !  differs from what is wanted...
    if ( nc_handle%define_mode .neqv. define ) then

       if ( define ) then ! assume file is in read/write mode..

          status = nf90_redef( nc_handle%id )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          nc_handle%define_mode = .True.
          call write_log( "Changed file to define mode" , msg_info , __FILE__ , __LINE__ )

       else ! assume file is in define mode..

          status = nf90_enddef( nc_handle%id )
          if ( status .ne. 0) then
            call handle_nc_error( __FILE__ , status , __LINE__ )
          endif
          nc_handle%define_mode = .False.
          call write_log( "Changed file to read/write mode" , msg_info , __FILE__ , __LINE__ )

       endif

    endif

  end subroutine change_mode
  !
  !----------------------------------------------------------------------
  !
  subroutine close_nc_output( header )

    ! close output file !

    use netcdf,   only: nf90_close
    use log_tools

    implicit none

    ! arguments..
    type(nc_header), pointer :: header

    ! local variables..
    integer :: status

    call write_log_div
    write(message,*) 'Closing output file '//trim( header%name )
    call write_log( message , msg_info )

    ! close the file..
    status = nf90_close( header%id )
    call handle_nc_error( __FILE__ , status , __LINE__ )
    
  end subroutine close_nc_output
  !
  !-------------------------------------------------------------------
  !
  subroutine create_coord_dim( file_id , dim_name , dim_length , dim_id , nc_type , var_id )

    ! calls appropriate other routines to make !
    ! a coordinate dimension.                  !

    implicit none

    ! arguments..
    integer,         intent(in)  :: dim_length, file_id, nc_type
    character(len=*),intent(in)  :: dim_name
    integer,         intent(out) :: dim_id, var_id

    call create_dim( file_id , dim_name, dim_length , dim_id )
    call create_var( file_id , dim_name , (/ dim_id /)  , nc_type , var_id )

  end subroutine create_coord_dim
  !
  !----------------------------------------------------------------------
  !
  subroutine create_dim( file_id , dim_name , dim_length , dim_id )

    ! Creates a new dimension in a file. !

    use linked_lists
    use netcdf

    ! arguments..
    integer,         intent(in)  :: file_id, dim_length
    character(len=*),intent(in)  :: dim_name
    integer,         intent(out) :: dim_id

    ! local variables..
    integer :: status

    status = nf90_def_dim( file_id , dim_name , dim_length , dim_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine create_dim
  !
  !----------------------------------------------------------------------
  !
  subroutine create_var( file_id , var_name , dim_ids , nc_type, var_id , from_file_id, from_var_id )

    ! Create the structure for a variable in a netcdf file. !
    ! If given another file's id and variable id it will    !
    !  copy the "units" and "long-name" attributes of that  !
    !  variable to this one.                                !
    ! nc_type is either 1 (integer) or 2 (real).            !

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    integer,intent(in)          :: file_id, dim_ids(:) , nc_type
    character(len=*),intent(in) :: var_name
    integer,intent(in),optional :: from_file_id, from_var_id
    integer,intent(out)         :: var_id

    ! local variables..
    integer :: status

    ! create variable..
    select case (nc_type)
    case (1) ! integer
      status = nf90_def_var( file_id , var_name , nf90_int , dim_ids , var_id )
      call handle_nc_error( __FILE__ , status , __LINE__ )
    case (2) ! float
      status = nf90_def_var( file_id , var_name , nf90_float , dim_ids , var_id )
      call handle_nc_error( __FILE__ , status , __LINE__ )
    case default
      write(message,*) "This type of variable is not recognised, only integer" &
                    // " or float are recognised for writing to netcdf files."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end select

    if ( present(from_file_id) .and. present(from_var_id) ) then
       ! copy units attribute..
       status = nf90_copy_att( from_file_id, from_var_id, 'units', file_id, var_id)
       call handle_nc_error( __FILE__ , status , __LINE__ )

       ! copy long-name attribute..
       status = nf90_copy_att( from_file_id, from_var_id, 'long_name', file_id, var_id)
       call handle_nc_error( __FILE__ , status , __LINE__ )
    endif

  end subroutine create_var
  !
  !-------------------------------------------------------------------
  !
  subroutine get_dim_info( file_id , name, dim_id, var_id, data )

    ! Returns the dim-id,var-id and size of !
    !  a given dimension in a given file.   !

    use netcdf

    ! arguments..
    integer,         intent(in)  :: file_id
    character(len=*),intent(in)  :: name
    integer,         intent(out) :: dim_id, var_id
    real,            pointer     :: data(:)

    ! local variables..
    integer                 :: dim_length, status
    real,allocatable,target :: dummy(:)

    ! get the dimension id..
    status = nf90_inq_dimid( file_id , trim(name) , dim_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! ..the variable id..
    status = nf90_inq_varid( file_id , trim(name) , var_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! the dimensions size..
    status = nf90_inquire_dimension( file_id, dim_id, len=dim_length )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! and finally the dimension data..
    allocate( dummy(dim_length) )
    status = nf90_get_var( file_id , var_id , dummy )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    data => dummy

  end subroutine get_dim_info
  !
  !-------------------------------------------------------------------
  !
  subroutine get_dim_info2( file_id , dim )

    ! Returns the dim-id,var-id and size of !
    !  a given dimension in a given file.   !

    use linked_lists, only: spaitem
    use log_tools
    use netcdf

    ! arguments..
    integer,   intent(in) :: file_id
    type(spaitem),pointer :: dim

    ! local variables..
    integer :: dim_length, nc_type, status
    integer(kind=1),allocatable :: short(:)
    real(kind=8),   allocatable :: double(:)

    ! get the dimension id..
    status = nf90_inq_dimid( file_id , dim%name , dim%dim_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! the dimension's size (in case we need it)..
    status = nf90_inquire_dimension( file_id, dim%dim_id, len=dim_length )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! ..the variable id (if it has one)..
    status = nf90_inq_varid( file_id , dim%name , dim%var_id )
    if ( status .eq. nf90_noerr ) then
      ! coordinate dimension, continue...

      ! what data type is the dimension data (int or real)...
      status = nf90_inquire_variable( file_id , dim%var_id , xtype=nc_type )
      call handle_nc_error( __FILE__ , status , __LINE__ )

      ! Load the data..
      select case (nc_type)
      case (nf90_short)
        ! short-form (16-bit) integer, convert to normal..
        allocate( short(dim_length) )
        status = nf90_get_var( file_id , dim%var_id , short )
        call handle_nc_error( __FILE__ , status , __LINE__ )
        dim%int_vector = int( short , kind=4 )
      case (nf90_int)
        status = nf90_get_var( file_id , dim%var_id , dim%int_vector )
      case (nf90_float)
        status = nf90_get_var( file_id , dim%var_id , dim%real_vector )
        call handle_nc_error( __FILE__ , status , __LINE__ )
      case (nf90_double)
        allocate( double(dim_length) )
        status = nf90_get_var( file_id , dim%var_id , double )
        call handle_nc_error( __FILE__ , status , __LINE__ )
        dim%real_vector = real( double , kind=4 )
      case default ! i.e. types nf90_byte or nf90_char
        call write_log( "dim should not be byte or char!" , msg_fatal , &
                                                  __FILE__ , __LINE__ )
      end select
      call handle_nc_error( __FILE__ , status , __LINE__ )

    else if ( status .eq. -49 ) then
      ! -49 is the errorcode for variable-not-found, which implies
      ! that this dimension is _not_ a coordinate dimension.
      ! Thus it must be an integer scalar..
      dim%data_type   = 1
      allocate(dim%int_scalar)
      dim%int_scalar  = dim_length
    else
      ! some other problem occured..
      call handle_nc_error( __FILE__ , status , __LINE__ )
    endif

  end subroutine get_dim_info2
  !
  !-------------------------------------------------------------------
  !
  subroutine get_nc_var( file_id, var_handle )

    ! for the given var (assumes %name is already put in)
    ! it (1) checks to see if name is in file,
    ! (2) gets dim-ids for var from file
    ! (3) allocates var%data to be correct size..

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    integer :: file_id
    type(nc_variable) :: var_handle

    ! local variables..
    integer                      :: i,n_dims,status,datatype
    integer,allocatable          :: dim_id_list(:), dim_length(:)
    character(len=nf90_max_name),allocatable :: dim_name_list(:)
    real,allocatable             :: data_4d(:,:,:,:)
    real(kind(0.d0)),allocatable :: dble_data_3d(:,:,:), dble_data_4d(:,:,:,:)

    ! get the variable's id..
    status = nf90_inq_varid( file_id , var_handle%name , var_handle%id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! get the number of dim's that variable has..
    status = nf90_inquire_variable( file_id , var_handle%id , xtype=datatype , ndims=n_dims )
    call handle_nc_error(  __FILE__ , status , __LINE__ )

    ! get the id of the dimensions associated with a particular variable..
    allocate( dim_id_list( n_dims ) )
    status = nf90_inquire_variable( file_id , var_handle%id , dimids=dim_id_list )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! and the sizes of those dims..
    allocate( dim_length( n_dims ) , dim_name_list( n_dims ) )
    do i = 1 , n_dims
      status = nf90_inquire_dimension( file_id , dim_id_list(i) , len=dim_length(i) , name=dim_name_list(i) )
      call handle_nc_error( __FILE__ , status , __LINE__ )
    enddo

    ! now get the data, and trim a dimension if it happens to be a 4-d object..
    allocate(var_handle%data)
    select case (n_dims)
    case (3)
       call write_log("Data variable is 3-d" , msg_info , __FILE__ , __LINE__ )
       if ( associated(var_handle%data%grid_3d ) ) deallocate(var_handle%data%grid_3d)
       allocate( var_handle%data%grid_3d( dim_length(1) , dim_length(2) , dim_length(3) ) )
       ! check in case input-data is a double (we will Not pass along!)..
       if ( datatype .eq. nf90_double ) then
          allocate( dble_data_3d( dim_length(1) , dim_length(2) , dim_length(3) ) )
          status = nf90_get_var( file_id , var_handle%id , dble_data_3d )
          call handle_nc_error(  __FILE__ , status , __LINE__ )
          var_handle%data%grid_3d = real( dble_data_3d )
       else
          status = nf90_get_var( file_id , var_handle%id , var_handle%data%grid_3d )
          call handle_nc_error(  __FILE__ , status , __LINE__ )
       endif
    case (4)
       write(message,*)"Data variable is 4-d - need to remove a dimension!"
       call write_log( message , msg_info , __FILE__ , __LINE__ )
       allocate( data_4d( dim_length(1) , dim_length(2) , dim_length(3) , dim_length(4) ) )
       ! check in case input-data is a double (we will Not pass along!)..
       if ( datatype .eq. nf90_double ) then
          allocate( dble_data_4d( dim_length(1) , dim_length(2) , dim_length(3) , dim_length(4) ) )
          status = nf90_get_var( file_id , var_handle%id , dble_data_4d )
          call handle_nc_error( __FILE__ ,  status , __LINE__ )
          data_4d = real( dble_data_4d , kind(data_4d) )
       else
          status = nf90_get_var( file_id , var_handle%id , data_4d )
          call handle_nc_error( __FILE__ ,  status , __LINE__ )
       endif
       ! now check which, if any, of the dimensions can be concatenated (ie are of size 1)
       if ( associated(var_handle%data%grid_3d) ) deallocate(var_handle%data%grid_3d)
       if ( dim_length(1) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(1))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
          allocate( var_handle%data%grid_3d( dim_length(2), dim_length(3), dim_length(4) ) )
          var_handle%data%grid_3d = data_4d( 1 , 1:dim_length(2) , 1:dim_length(3) , 1:dim_length(4) )
       else if ( dim_length(2) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(2))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
          allocate( var_handle%data%grid_3d( dim_length(1), dim_length(3), dim_length(4) ) )
          var_handle%data%grid_3d = data_4d( 1:dim_length(1) , 1 , 1:dim_length(3) , 1:dim_length(4) )
       else if ( dim_length(3) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(3))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
          allocate( var_handle%data%grid_3d( dim_length(1), dim_length(2), dim_length(4) ) )
          var_handle%data%grid_3d = data_4d( 1:dim_length(1) , 1:dim_length(2) , 1 , 1:dim_length(4) )
       else if ( dim_length(4) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(4))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
          allocate( var_handle%data%grid_3d( dim_length(1), dim_length(2), dim_length(3) ) )
          var_handle%data%grid_3d = data_4d( 1:dim_length(1) , 1:dim_length(2) , 1:dim_length(3) , 1 )
       endif
    case default ! n-dims case
       write(message,*) "I don't know what to do with a variable with this many dimensions:",n_dims
       call write_log( trim(message) , msg_fatal, __FILE__ , __LINE__ )
    end select

    if ( allocated(  dim_id_list  ) )  deallocate(  dim_id_list  )
    if ( allocated(   dim_length  ) )  deallocate(   dim_length  )
    if ( allocated( dim_name_list ) )  deallocate( dim_name_list )
    if ( allocated(    data_4d    ) )  deallocate(    data_4d    )

  end subroutine get_nc_var
  !
  !-------------------------------------------------------------------
  !
  subroutine get_nc_var2( file_id , var )

    !>   <!

    use linked_lists, only: spaitem
    use log_tools
    use netcdf

    implicit none

    ! arguments..
    integer,intent(in)    :: file_id
    type(spaitem),pointer :: var

    ! local variables..
    integer                     :: i, n_dims, nc_type, status
    integer,allocatable         :: dim_id_list(:), dim_length(:)
    character(len=nf90_max_name),&
                    allocatable :: dim_name_list(:)
    integer(kind=1),allocatable :: short(:)
    real(kind=8),   allocatable :: double(:)

    ! get the variable's id..
    status = nf90_inq_varid( file_id , var%name , var%var_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! check how many dim's the variable in the file has..
    status = nf90_inquire_variable( file_id , var%var_id , ndims=n_dims )
    call handle_nc_error(  __FILE__ , status , __LINE__ )

    ! is this the same as the spa-variable expects?
    if ( n_dims .ne. var%ndims ) then
       write(message,*) "The number of dimensions in the file does not match" &
            //" those expected for the variable we've been asked to load."
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    if ( n_dims .gt. 0 ) then
       ! get the id of the dimensions associated with a particular variable..
       allocate( dim_id_list( n_dims ) )
       status = nf90_inquire_variable( file_id , var%var_id , dimids=dim_id_list )
       call handle_nc_error( __FILE__ , status , __LINE__ )

       ! and the sizes of those dims..
       allocate( dim_length( n_dims ) , dim_name_list( n_dims ) )
       do i = 1 , n_dims
          status = nf90_inquire_dimension( file_id , dim_id_list(i) , &
               len=dim_length(i) , name=dim_name_list(i) )
          call handle_nc_error( __FILE__ , status , __LINE__ )
       enddo

!!! Check that the variable's number of elements (ie dim-lengths) 
!!!  in the restart-file match the number expected.

    endif


    ! We should set the datatype based on what's in the file.
    ! ie read the data, (by inquiring into the file) and then 
    ! declare the appropriate datatype and appropriate pointer..

    ! what data type is the dimension data (int or real)...
    status = nf90_inquire_variable( file_id , var%var_id , xtype=nc_type )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! Knowing the data-type and the number of
    !  dimensions we can now load the data..
    select case (n_dims)
    case (0)        
       ! n_dims == 0 => scalar
       select case (nc_type)
       case (nf90_short)
          ! short-form (16-bit) integer, convert to normal..
          allocate( short(1) )
          status = nf90_get_var( file_id , var%var_id , short )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          var%data_type = 1
          allocate(var%int_scalar)
          var%int_scalar = int( short(1) , kind=4 )

       case (nf90_int)
          var%data_type = 1
          allocate(var%int_scalar)
          status = nf90_get_var( file_id , var%var_id , var%int_scalar )

       case (nf90_float)
          var%data_type = 3
          allocate(var%real_scalar)
          status = nf90_get_var( file_id , var%var_id , var%real_scalar )
          call handle_nc_error( __FILE__ , status , __LINE__ )

       case (nf90_double)
          allocate( double(1) )
          status = nf90_get_var( file_id , var%var_id , double )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          var%data_type = 3
          allocate(var%real_scalar)
          var%real_scalar = real( double(1) , kind=4 )

       case default ! i.e. types nf90_byte or nf90_char
          call write_log( "dim should not be byte or char!" , msg_fatal , &
               __FILE__ , __LINE__ )
       end select
       call handle_nc_error( __FILE__ , status , __LINE__ )
    case (1)
       ! n_dims == 1 => vector
       select case (nc_type)
       case (nf90_short)
          ! short-form (16-bit) integer, convert to normal..
          allocate( short( dim_length(1) ) )
          status = nf90_get_var( file_id , var%var_id , short )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          var%data_type = 2
          allocate( var%int_vector( dim_length(1) ) )
          var%int_vector = int( short , kind=4 )

       case (nf90_int)
          var%data_type = 2
          allocate( var%int_vector( dim_length(1) ) )
          status = nf90_get_var( file_id , var%var_id , var%int_vector )

       case (nf90_float)
          var%data_type = 4
          allocate( var%real_vector( dim_length(1) ) )
          status = nf90_get_var( file_id , var%var_id , var%real_vector )
          call handle_nc_error( __FILE__ , status , __LINE__ )

       case (nf90_double)
          allocate( double( dim_length(1) ) )
          status = nf90_get_var( file_id , var%var_id , double )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          var%data_type = 4
          allocate( var%real_vector( dim_length(1) ) )
          var%real_vector = real( double , kind=4 )

       case default ! i.e. types nf90_byte or nf90_char
          call write_log( "dim should not be byte or char!" , msg_fatal , &
               __FILE__ , __LINE__ )
       end select
       call handle_nc_error( __FILE__ , status , __LINE__ )
    case default
       ! n_dims >  1 => (for now) ignore
       write(message,*)"cannot yet handle more than 1 dimension!" 
       call write_log( message, msg_fatal , __FILE__ , __LINE__ )
    end select

    call handle_nc_error(  __FILE__ , status , __LINE__ )

    if ( allocated(  dim_id_list  ) )  deallocate(  dim_id_list  )
    if ( allocated(   dim_length  ) )  deallocate(   dim_length  )
    if ( allocated( dim_name_list ) )  deallocate( dim_name_list )

  end subroutine get_nc_var2
  !
  !-------------------------------------------------------------------
  !
  subroutine handle_nc_error( file , status , line )

    ! handle netCDF error !

    use netcdf
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: file    ! name of f90 file error occured in
    integer,         intent(in) :: status  ! netCDF return value
    integer,optional,intent(in) :: line    ! line number error occured at

    ! local variables..
    integer :: linenos

    if ( .not. present(line) ) then
       linenos = -1
    else
       linenos = line
    endif

    if ( status .ne. nf90_noerr ) then
       call write_log( nf90_strerror(status) , type=msg_fatal , file=file , line=linenos )
    end if

  end subroutine handle_nc_error
  !
  !-------------------------------------------------------------------
  !
  subroutine open_nc_file( file_info , for_output , append )

    ! opens a netcdf file for reading/writing. !

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    type(nc_header),pointer     :: file_info     ! Structure; contains filename, returned with extra info
    logical,intent(in),optional :: for_output, & ! is file to be opened in write mode?
                                   append        ! if file exists should we overwrite, or append to?

    ! local variables..
    integer :: status
    logical :: write_mode, append_mode

    ! Open in write mode? (if not given, assume read-only)..
    if (.not. present( for_output ) ) then
       write_mode = .False.
    else
       write_mode = for_output
       ! Overwrite file if it already exists? (if not given, assume appending)..
       if ( .not. present(append) ) then
         append_mode = .True.
       else
         append_mode = append
       endif
    endif

    if ( write_mode ) then
      ! We want to write to this file..

      ! First off, check if file already exists...
      inquire( file=trim(file_info%name) , exist=file_info%exists ,iostat=status )
      if ( status .ne. 0 ) call write_log( "Could not get information about nc file" , &
                                              msg_fatal , __FILE__ , __LINE__ )

      if ( file_info%exists ) then
        ! It already exists..
        if ( append_mode ) then
          ! ..and we want to append to it...
          status = nf90_open( trim(file_info%name) , nf90_write , file_info%id )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          file_info%define_mode = .False.
          write(message,*)"Opened existing file "//trim(file_info%name)//" for appending to."
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
        else
          ! ..and we want to overwrite it..
          status = nf90_create( trim(file_info%name) , nf90_clobber , file_info%id )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          file_info%define_mode = .True.
          write(message,*)"Overwriting existing file "//trim(file_info%name)//"."
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
        endif
      else
         ! It doesn't exist yet..
         status = nf90_create( trim(file_info%name) , nf90_noclobber , file_info%id )
         call handle_nc_error( __FILE__ , status , __LINE__ )
         file_info%define_mode = .True.
         write(message,*)"Created new file "//trim(file_info%name)//" for writing to."
         call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
      endif
    else
      ! We only want to read from it..
      status = nf90_open( trim(file_info%name) , nf90_nowrite , file_info%id )
      call handle_nc_error( __FILE__ , status , __LINE__ )
      file_info%define_mode = .False.
      write(message,*)"Opened existing file "//trim(file_info%name)//" for reading from."
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
      file_info%exists = .True.
    endif

  end subroutine open_nc_file
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module 
  !
  !----------------------------------------------------------------------
  !
  subroutine copy_global_atts( from_nc , to_nc )

    ! Copy the global attributes from one file to another. !

    use netcdf

    implicit none

    ! arguments..
    type(nc_header),intent(in) :: from_nc, to_nc

    ! local variables..
    integer :: i, nos_gl_atts, status
    character(len=nf90_max_name) :: name

    ! get number of global attributes..
    status = nf90_inquire( from_nc%id, nAttributes=nos_gl_atts )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    do i = 1 , nos_gl_atts

       ! get the name associated with each attribute..
       status = nf90_inq_attname( from_nc%id, nf90_global, i, name)
       call handle_nc_error( __FILE__ , status , __LINE__ )

       ! specify a copy of this attribute by name..
       status = nf90_copy_att( from_nc%id, nf90_global, trim(name), to_nc%id, nf90_global )
       call handle_nc_error( __FILE__ , status , __LINE__ )

    enddo

  end subroutine copy_global_atts
  !
  !----------------------------------------------------------------------
  !
  subroutine get_int_global_att( file_id , att_name , data )

    use netcdf

    implicit none

    ! arguments..
    integer,            intent( in) :: file_id
    character(len=*),   intent( in) :: att_name
    integer,allocatable,intent(out) :: data(:)

    ! local variables..
    integer :: attsize, status

    status = nf90_inquire_attribute( file_id , NF90_GLOBAL , att_name , len=attsize )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    allocate( data(attsize) )
    status = nf90_get_att( file_id , NF90_GLOBAL , att_name , data )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine get_int_global_att
  !
  !----------------------------------------------------------------------
  !
  subroutine get_real_global_att( file_id , att_name , data )

    use netcdf

    implicit none

    ! arguments..
    integer,         intent( in) :: file_id
    character(len=*),intent( in) :: att_name
    real,allocatable,intent(out) :: data(:)

    ! local variables..
    integer :: attsize, status

    status = nf90_inquire_attribute( file_id , NF90_GLOBAL , att_name , len=attsize )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    allocate( data(attsize) )
    status = nf90_get_att( file_id , NF90_GLOBAL , att_name , data )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine get_real_global_att
  !
  !----------------------------------------------------------------------
  !
  subroutine put_int_global_att( file_id , att_name , data )

    use netcdf, only: nf90_global, nf90_put_att

    implicit none

    ! arguments..
    integer,         intent(in) :: file_id, data(:)
    character(len=*),intent(in) :: att_name

    ! local variables..
    integer :: status

    status = nf90_put_att( file_id , NF90_GLOBAL , att_name , data )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine put_int_global_att
  !
  !----------------------------------------------------------------------
  !
  subroutine put_real_global_att( file_id , att_name , data )

    use netcdf, only: nf90_global, nf90_put_att

    implicit none

    ! arguments..
    integer,         intent(in) :: file_id
    character(len=*),intent(in) :: att_name
    real,            intent(in) :: data(:)

    ! local variables..
    integer :: status

    status = nf90_put_att( file_id , NF90_GLOBAL , att_name , data )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine put_real_global_att
  !
  !----------------------------------------------------------------------
  !
  subroutine put_var_nc( file_id , var )

    ! Write a variable/coordinate-dimension's !
    ! data to a netcdf file.                  !

    use linked_lists, only: spaitem
    use log_tools
    use netcdf,       only: nf90_put_var

    implicit none

    ! arguments..
    integer,   intent(in) :: file_id   ! handle to nc-file to write to
    type(spaitem),pointer :: var

    ! local variables..
    integer :: status

    select case ( var%data_type )
    case (1)
      status = nf90_put_var( file_id , var%var_id , var%int_scalar )
    case (2)
      status = nf90_put_var( file_id , var%var_id , var%int_vector )
    case (3)
      status = nf90_put_var( file_id , var%var_id , var%real_scalar )
    case (4)
      status = nf90_put_var( file_id , var%var_id , var%real_vector )
    case (5)
      status = nf90_put_var( file_id , var%var_id , var%real_grid_2d )
    case (6)
      status = nf90_put_var( file_id , var%var_id , var%real_grid_3d )
    case default
      write(message,*) "Data-type: ",var%data_type," not recognised."
      call write_log( message, msg_fatal , __FILE__ , __LINE__ )
    end select

    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine put_var_nc
  !
  !----------------------------------------------------------------------
  !
end module netcdf_tools
!
!------------------------------------------------------------------------
!
