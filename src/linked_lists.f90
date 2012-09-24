! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module linked_lists

  !! This module declares the basic forms of  !!
  !! linked-lists that are used to keep track !!
  !! of variabes, which amongst other things  !!
  !! makes it easier to process SPA's IO.     !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: spalist, spaitem, spa_vars, spa_dims
  ! Variables..
  public :: append, dims_from_var_list, get_item, new_item, print_item

  ! An extension for holding SPA variables, with space for netcdf-bits
  ! as well.  Note we use the same construct to build SPA dimensions,
  ! but they can only have datatypes 1,2 or 4 (ie be an integer value,
  ! such as the nos-of-soil-layers, or be a vector, such as time (and
  ! time%year => int_vec, whilst time%daytime => real_vec).
  integer,parameter :: str_length = 35
  type spaitem
    ! basic bits..
    character(len=str_length) :: name =''
    character(len=str_length) :: dim_name = ''
    ! data info..
    integer         :: data_type = 0 ! 1=>int_scalar, 2=> int_vector, etc.
    integer,pointer :: int_scalar     => null()
    integer,pointer :: int_vector(:)   => null()
    real,   pointer :: real_scalar      => null()
    real,   pointer :: real_vector(:)    => null()
    real,   pointer :: real_grid_2d(:,:)  => null()
    real,   pointer :: real_grid_3d(:,:,:) => null()
    ! netcdf bits..
    integer         :: dim_id = -999   ! dim id
    integer         :: var_id = -999   ! var id
    real            :: mdi  = -1e32    ! missing data indicator
    integer         :: ndims = -1      ! nos of dims
    ! navigation bits..
    type(spaitem),pointer :: next_ptr => null()
    type(spaitem),pointer :: prev_ptr => null()
    type(spalist),pointer :: up_ptr   => null()
  end type spaitem
  type spalist
    type(spaitem),pointer :: first_ptr => null()
  end type spalist

  ! an interface to handle all the different possible type/shapes
  !  of variable..
  ! (only the interface is public, the individual routines are private)
  interface new_item
    module procedure new_item_int, new_item_int_vec, new_item_real, &
                   new_item_real_vec, new_item_real_2d, new_item_real_3d
  end interface

  ! The base lists that SPA's variables and dimensions are defined in
  ! (The lists are filled in spa_io, but are declared here in)
  ! (order to make them accessible to most other SPA modules.)
  type(spalist),pointer,save :: spa_vars , spa_dims

contains
  !
  !----------------------------------------------------------------------
  !
  ! PUBLIC PROCEDURES FIRST..
  !
  !----------------------------------------------------------------------
  !
  subroutine append( list , item )

    ! Append an item to a list. !

    implicit none

    ! arguments..
    type(spalist),intent(inout), target :: list
    type(spaitem),               target :: item

    ! local variables..
    type(spaitem), pointer :: last

    if (associated(item%up_ptr)) call remove(item)
    item%up_ptr => list

    if (associated(list%first_ptr)) then
      last => list%first_ptr%prev_ptr
      last%next_ptr => item
      item%prev_ptr => last
      list%first_ptr%prev_ptr => item
    else
      list%first_ptr => item
      item%prev_ptr => item
    end if

  end subroutine append
  !
  !----------------------------------------------------------------------
  !
  subroutine dims_from_var_list( dim_list , var_list , vars_dim_list )

    ! Given a list of available dimensions, and a list of   !
    ! desired variables, this s/r will construct a new list !
    ! of the dimensions required for those variables.       !
    ! (as long the requested dimension exists within the    !
    ! available dimensions list).                           !

    use log_tools

    implicit none

    ! argument..
    type(spalist),pointer :: dim_list, var_list ! input
    type(spalist),pointer :: vars_dim_list      ! output

    ! local variables..
    character(len=500 )   :: dims_so_far = repeat('  ',250 )
    type(spaitem),pointer :: p
    integer               :: status

    p => var_list%first_ptr
    do
      ! are we at the end of the list...
      if ( .not. associated(p) ) exit
      ! if not, check whether item has a dimension..
      if ( p%dim_name .ne. '' ) then
        ! if yes, then check if we already
        !  have this dim in our list..
        status = index( trim(dims_so_far) , trim(p%dim_name) )
        if (status .eq. 0 ) then

          ! dimension not in list, add it to list..
          call append( vars_dim_list , get_item( dim_list , p%dim_name ) )
          ! and keep track of which dims we've got..
          if ( len(trim(dims_so_far)) .eq. 0 ) then
            dims_so_far = trim(adjustl(p%dim_name))
          else
            write(dims_so_far,*) trim(dims_so_far)//' '&
                                     //trim(adjustl(p%dim_name))
          endif
        endif
      endif
      p => p%next_ptr
    end do

    ! Sanity check..
    if ( len(trim(dims_so_far)) .eq. 0 ) then
      write(message,*)"Warning! The provided var-list required " &
                    //"no dimensions - is that expected?"
      call write_log( message , msg_warning , __FILE__ , __LINE__ )
    endif

  end subroutine dims_from_var_list
  !
  !----------------------------------------------------------------------
  !
  function get_item( list , item_name )

    ! returns a copy of an item in a list but without the !
    ! list references (ie. up/next/prev ptrs are reset)   !

    implicit none

    ! arguments..
    type(spalist),intent(in),target :: list
    character(*), intent(in)     :: item_name

    ! function result..
    type(spaitem),pointer :: get_item

    ! local variables...
    type(spaitem), pointer :: p

    p => find_item( list , trim(item_name) )

    ! if find_item failed, then exit early..
    if ( .not. associated( p ) ) return

    ! copy the values from one item to another
    allocate(get_item)
    get_item => copy( p )

    ! last bit of initialisatin is to set it to loop on itself..
    get_item%prev_ptr => get_item
    get_item%up_ptr   => null()

  end function get_item
  !
  !----------------------------------------------------------------------
  !
  subroutine print_item( item )

    ! s/r to print out info about an item (for debugging) !

    implicit none

    ! argument..
    type(spaitem),pointer :: item

    if ( .not. associated(item) ) then
       print*,"   This item is null!"
       return
    endif

    print*,"   Item's name       :   ",item%name
    if ( item%dim_name .ne. '' ) &
       print*,"   Item's dim name   :   ",item%dim_name
    print '(A,I1)',"    Item's data_type  :   ",item%data_type
    select case( item%data_type)
    case (1)
       if (associated(item%int_scalar)) &
            print*,"   Item's int_scalar :",item%int_scalar
    case (2)
       if (associated(item%int_vector)) &
            print*,"   Item's int_vector :",item%int_vector
    case (3)
       if (associated(item%real_scalar)) &
            print*,"   Item's real_scalar:",item%real_scalar
    case (4)
       if (associated(item%real_vector)) &
            print*,"   Item's real_vector:",item%real_vector
    case default
       return
    end select
    print*,' '

  end subroutine print_item
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module.
  !
  !----------------------------------------------------------------------
  !
  subroutine check_name_length( var_name )

    ! make sure that the name we want to write will actually fit !

    use log_tools

    implicit none

    ! arguments..
    character(*) :: var_name

    ! check the name we want to write is going to fit!..
    if ( len(var_name) .gt. str_length ) then
      write(message,*) "Variable name: ",new_line('x'),&
                       "   ",var_name,new_line('x'),&
                       "   has more than",str_length," characters!"
      call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    endif

  end subroutine check_name_length
  !
  !----------------------------------------------------------------------
  !
  function copy( item )

    ! Returns a copy of just the  !
    ! data-contents of a variable !

    implicit none

    ! arguments..
    type(spaitem),target :: item

    ! function result..
    type(spaitem), pointer :: copy

    allocate(copy)

    ! copy over the name,dims and flag..
    copy%name       = item%name
    if ( item%dim_name .ne. '' ) &
      copy%dim_name = item%dim_name
    copy%data_type  = item%data_type
    copy%dim_id     = item%dim_id
    copy%var_id     = item%var_id
    copy%mdi        = item%mdi
    copy%ndims      = item%ndims

    ! copy over the appropriate part of the data..
    select case ( copy%data_type )
    case (1)
      copy%int_scalar  => item%int_scalar
    case (2)
      copy%int_vector  => item%int_vector
    case (3)
      copy%real_scalar => item%real_scalar
    case (4)
      copy%real_vector => item%real_vector
    case (5)
      copy%real_grid_2d => item%real_grid_2d
    case (6)
      copy%real_grid_3d => item%real_grid_3d
    case default
      return
    end select

  end function copy
  !
  !----------------------------------------------------------------------
  !
  function find_item( list , name )

    ! finds where in a list is the item !
    ! whose name matches that given.    !

    use log_tools

    implicit none

    ! arguments..
    type(spalist), target      :: list
    character(*),intent(in) :: name

    ! function result..
    type(spaitem), pointer :: find_item

    ! check we have a list to start with..
    if ( .not. associated(list%first_ptr) ) return

    ! start at the beginning of the list..
    find_item => list%first_ptr
    do
      ! check if this item has the right name..
      if ( find_item%name .eq. name ) then
        ! ..if yes, then we're done..
        return
      else
        ! ..if no, we have to check if there is a
        !   next item in the list to move onto..
        if ( .not. associated( find_item%next_ptr ) ) then
          ! ..if there isn't, we clear find_item, print an
          !   an error message and exit.
          nullify(find_item)
          write(message,*) "Could not find ",name," in list."
          call write_log( message , msg_warning , __FILE__ , __LINE__ )
          return
        endif
      endif
      find_item => find_item%next_ptr
    enddo

  end function find_item
  !
  !----------------------------------------------------------------------
  !
  ! The next four functions will return a new item, with the
  !  correct data_type number set and data-type associated. 
  ! We use overloading of the name 'new_item' in order to avoid
  !  having to remember which of the four funcs we need
  ! (Which routine to use is based upon the supplied arguments)
  !
  function new_item_int( var_name , var )

    ! see above !

    implicit none

    ! arguments..
    character(*),  target,intent(in) :: var_name
    integer,       target,intent(in) :: var

    ! function result..
    type(spaitem),pointer :: new_item_int

    call check_name_length( var_name )

    allocate(new_item_int)
    new_item_int%name        =  var_name
    new_item_int%data_type   =  1
    new_item_int%ndims       =  0
    new_item_int%int_scalar  => var
    new_item_int%prev_ptr    => new_item_int

  end function new_item_int
  !
  !----------------------------------------------------------------------
  !
  function new_item_int_vec( var_name , var , dim_name )

    ! see above !

    implicit none

    ! arguments..
    character(*),target,intent(in) :: var_name
    integer,     target,intent(in) :: var(:)
    character(*),       intent(in) :: dim_name

    ! function result..
    type(spaitem),pointer :: new_item_int_vec

    call check_name_length( var_name )

    allocate(new_item_int_vec)
    new_item_int_vec%name       = var_name
    new_item_int_vec%data_type  = 2
    new_item_int_vec%int_vector => var
    new_item_int_vec%dim_name   = dim_name
    new_item_int_vec%ndims      = 1
    new_item_int_vec%prev_ptr   => new_item_int_vec

  end function new_item_int_vec
  !
  !----------------------------------------------------------------------
  !
  function new_item_real( var_name , var )

    ! see above !

    implicit none

    ! arguments..
    character(*),         target,intent(in) :: var_name
    real,                 target,intent(in) :: var

    ! function result..
    type(spaitem),pointer :: new_item_real

    call check_name_length( var_name )

    allocate(new_item_real)
    new_item_real%name        =  var_name
    new_item_real%data_type   =  3
    new_item_real%real_scalar => var
    new_item_real%ndims       =  0
    new_item_real%prev_ptr    => new_item_real

  end function new_item_real
  !
  !----------------------------------------------------------------------
  !
  function new_item_real_vec( var_name , var , dim_name )

    ! see above !

    implicit none

    ! arguments..
    character(*),target,intent(in) :: var_name
    real,        target,intent(in) :: var(:)
    character(*),target,intent(in) :: dim_name

    ! function result..
    type(spaitem),pointer :: new_item_real_vec

    call check_name_length( var_name )

    allocate(new_item_real_vec)
    new_item_real_vec%name        =  var_name
    new_item_real_vec%data_type   =  4
    new_item_real_vec%real_vector => var
    new_item_real_vec%dim_name    =  dim_name
    new_item_real_vec%ndims       =  1
    new_item_real_vec%prev_ptr    => new_item_real_vec

  end function new_item_real_vec
  !
  !----------------------------------------------------------------------
  !
  function new_item_real_2d( var_name , var , dim_name )

    ! see above !

    implicit none

    ! arguments..
    character(*),target,intent(in) :: var_name
    real,        target,intent(in) :: var(:,:)
    character(*),target,intent(in) :: dim_name

    ! function result..
    type(spaitem),pointer :: new_item_real_2d

    call check_name_length( var_name )

    allocate(new_item_real_2d)
    new_item_real_2d%name         =  var_name
    new_item_real_2d%data_type    =  5
    new_item_real_2d%real_grid_2d => var
    new_item_real_2d%dim_name     =  dim_name
    new_item_real_2d%ndims        =  2
    new_item_real_2d%prev_ptr     => new_item_real_2d

  end function new_item_real_2d
  !
  !----------------------------------------------------------------------
  !
  function new_item_real_3d( var_name , var , dim_name )

    ! see above !

    implicit none

    ! arguments..
    character(*),target,intent(in) :: var_name
    real,        target,intent(in) :: var(:,:,:)
    character(*),target,intent(in) :: dim_name

    ! function result..
    type(spaitem),pointer :: new_item_real_3d

    call check_name_length( var_name )

    allocate(new_item_real_3d)
    new_item_real_3d%name         =  var_name
    new_item_real_3d%data_type    =  6
    new_item_real_3d%real_grid_3d => var
    new_item_real_3d%dim_name     =  dim_name
    new_item_real_3d%ndims        =  3
    new_item_real_3d%prev_ptr     => new_item_real_3d

  end function new_item_real_3d
  !
  !----------------------------------------------------------------------
  !
  subroutine remove( item )

    ! Remove an item from a list, but keep it and its value. !

    implicit none

    ! arguments..
    type(spaitem), intent(inout), target :: item

    ! local variables..
    type(spalist), pointer :: list

    list => item%up_ptr
    if (associated(list)) then
      if (associated(item%prev_ptr, item)) then
        ! Single item in list.
        nullify(list%first_ptr)
      else if (.not.associated(item%next_ptr)) then
        ! Last item in list.
        list%first_ptr%prev_ptr => item%prev_ptr
        nullify(item%prev_ptr%next_ptr)
        item%prev_ptr => item
      else if (associated(list%first_ptr, item)) then
        ! First item in list.
        list%first_ptr => item%next_ptr          ! first = next.
        item%prev_ptr%prev_ptr => item%next_ptr  ! last%prev = item%next.
        item%next_ptr%prev_ptr => item%prev_ptr  ! next%prev = last.
      end if
    end if
    nullify(item%up_ptr)

  end subroutine remove
  !
  !----------------------------------------------------------------------
  !
end module linked_lists
!
!------------------------------------------------------------------------
!
