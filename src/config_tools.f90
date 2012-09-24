! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module config_tools

  !! This module provides the tools required to be able to read a !!
  !!  configuration file.                                         !!
  !! The reading of a given file will depend on the model, e.g.   !!
  !!  SPA requires a meteorology driver file, and needs to be     !!
  !!  told certain things about that file.  This is specific to   !!
  !!  SPA and so is controlled elsewhere (within spa_config.f90)  !!
  !!                                                              !!
  !! The code for these tools is taken from the Glimmer project,  !!
  !!  which itself took them from the Generic Mapping Tools       !!
  !!  project and converted them to Fortran 90.                   !!

  use gv_scale_declarations, only: fname_length

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: ConfigRead, GetSection, GetValue, PrintConfig
  ! Variables..
  public :: ConfigSection


  integer, parameter ::       dp = kind(1.0d0), & !
                        namelen  = 50,          & !
                        valuelen = 200,         & !
                        linelen  = 250            ! 

  type ConfigValue
     character(len=namelen)     :: name = ''
     character(len=valuelen)    :: value
     type(ConfigValue), pointer :: next=>null()
  end type ConfigValue

  type ConfigSection
     character(len=namelen) :: name = ''
     logical :: used = .false.
     type(ConfigValue), pointer :: values=>null()
     type(ConfigSection), pointer :: next=>null()
  end type ConfigSection

  type ConfigData
     ! This type exists so that we can have
     ! arrays of config data, since f90 doesn't
     ! allow arrays of pointers
     type(ConfigSection), pointer :: config=>null()
  end type ConfigData

  interface GetValue
     module procedure GetValueDouble, GetValueReal, GetValueInt, GetValueChar, GetValueLogical, &
          GetValueDoubleArray, GetValueRealArray, GetValueIntArray, GetValueCharArray
  end interface

  interface ConfigSetValue
     module procedure ConfigSetValueData, ConfigSetValueSec
  end interface

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine ConfigRead(fname,config)
    ! read configuration file
    use log_tools
    implicit none
    character(len=*), intent(in) :: fname
    ! name of configuration file
    type(ConfigSection), pointer :: config

    ! local variables..
    type(ConfigSection), pointer :: this_section
    type(ConfigValue), pointer ::  this_value
    logical there
    integer unit,ios,linenr
    character(len=linelen) :: line

    inquire (exist=there,file=fname)
    if (.not.there) then
       write(message,*) 'ConfigRead: Cannot open configuration file '//trim(fname)
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    unit=99
    open(unit,file=trim(fname),status='old')
    ios=0
    linenr=0
    config=>null()
    this_section=>null()
    do while(ios.eq.0)
       read(unit,fmt='(a250)',iostat=ios) line
       line = adjustl(line)
       if (ios.ne.0) then
          exit
       end if
       if (.not.(line(1:1).eq.'!' .or. line(1:1).eq.'#' .or. line(1:1).eq.';' .or. line(1:1).eq.' ')) then
          ! handle comments
          if (line(1:1).eq.'[') then
             ! new section
             call handle_section(linenr,line,this_section)
             this_value=>null()
             if (.not.associated(config)) then
                ! this is the first section in config file
                config=>this_section
             end if
          else
             ! handle value
             if (.not.associated(this_section)) then
                write(message,*) 'ConfigRead: No section defined yet'
                call write_log( message , msg_error , __FILE__ , __LINE__ )
                write(message,*) trim(adjustl(fname)), linenr
                call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
             end if
             call handle_value(linenr,line,this_value)
             if (.not.associated(this_section%values)) then
                this_section%values => this_value
             end if
          end if
       end if
       linenr = linenr + 1
    end do
    close(unit)
  end subroutine ConfigRead
  !
  !----------------------------------------------------------------------
  !
  subroutine GetSection(config,found,name)
    ! Find and return section with name
    implicit none
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: found
    character(len=*),intent(in) :: name

    found=>config
    do while(associated(found))
       if (name.eq.trim(found%name)) then
          found%used = .true.
          return
       end if
       found=>found%next
    end do
  end subroutine GetSection
  !
  !----------------------------------------------------------------------
  !
  subroutine PrintConfig(config,unit)
    implicit none
    type(ConfigSection), pointer :: config
    integer,optional :: unit

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val

    sec=>config
    do while( associated(sec) )
       if ( present(unit) ) then
          write(unit,*) sec%name
       else
          write(*,*) sec%name
       endif
       val=>sec%values
       do while( associated(val) )
          if ( present(unit) ) then
             write(unit,*) '  ',trim(val%name),' == ', trim(val%value)
          else
             write(*,*) '  ',trim(val%name),' == ', trim(val%value)
          endif
          val=>val%next
       end do
       if ( present(unit) ) then
          write(unit,*)
       else
          write(*,*)
       endif
       sec=>sec%next
    end do
  end subroutine PrintConfig
  !
  !----------------------------------------------------------------------
  !
  ! Private procedures below, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine CheckSections(config)
    ! traverse linked list and check that all sections have been used
    use log_tools
    implicit none
    type(ConfigSection), pointer :: config

    ! local variables..
    type(ConfigSection), pointer :: cf

    cf=>config
    do 
      if (associated(cf)) exit 
      if (.not.cf%used) then
        write(message,*)'CheckSections: Unused section: '//trim(cf%name)
        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
      end if
      cf=>cf%next
    end do
  end subroutine CheckSections
  !
  !----------------------------------------------------------------------
  !
  subroutine ConfigAsString(config,string)

    implicit none

    type(ConfigSection), pointer :: config
    character(*),intent(out) :: string

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val
    character, parameter :: endline = achar(10)

    string=''

    sec=>config
    do while(associated(sec))
       string=trim(string)//'['//trim(sec%name)//']'//trim(endline)
       val=>sec%values
       do while(associated(val))
          string=trim(string)//trim(val%name)//': '//trim(val%value)//trim(endline)
          val=>val%next
       end do
       sec=>sec%next
    end do
  end subroutine ConfigAsString
  !
  !----------------------------------------------------------------------
  !
  logical function ConfigSectionHasTag(section,tag)

    type(ConfigSection), pointer :: section
    character(*) :: tag
    character(200) :: testtag

    ConfigSectionHasTag=.false.
    if (ConfigSectionHasValue(section,'tag',testtag)) then
       if (trim(tag)==trim(testtag)) then
          ConfigSectionHasTag=.true.
       end if
    end if

  end function ConfigSectionhasTag
  !
  !----------------------------------------------------------------------
  !
  logical function ConfigSectionHasValue(section,valname,val)

    type(ConfigSection), pointer :: section
    type(ConfigValue), pointer :: thisval
    character(*) :: valname,val

    ConfigSectionHasValue=.false.
    val=''

    if (.not.associated(section)) return

    thisval=>section%values
    do
       if (.not.associated(thisval)) exit
       if (trim(valname)==trim(thisval%name)) then
          val=trim(thisval%value)
          ConfigSectionHasValue=.true.
          exit
       else
          thisval=>thisval%next
       end if
    end do

  end function ConfigSectionHasValue
  !
  !----------------------------------------------------------------------
  !
  subroutine ConfigSetValueData(config,secname,valname,value,tag)
    ! Either overwrite a given key-value pair,
    ! or create a new one

    type(ConfigData) :: config
    character(*) :: secname,valname,value
    character(*),optional :: tag

    call ConfigSetValueSec(config%config,secname,valname,value,tag)

  end subroutine ConfigSetValueData
  !
  !----------------------------------------------------------------------
  !
  subroutine ConfigSetValueSec(config,secname,valname,value,tag)
    ! Either overwrite a given key-value pair,
    ! or create a new one
    ! tag is a label that you can give to a particular config section
    ! allowing the identification of the right section when
    ! several with the same name are present (e.g. [CF output])

    type(ConfigSection), pointer :: config
    character(*) :: secname,valname,value
    character(*),optional :: tag
    type(ConfigSection), pointer :: found
    type(ConfigSection), pointer :: newsec
    type(ConfigValue), pointer :: val
    type(ConfigValue), pointer :: newval
    type(ConfigValue), pointer :: newtag
    logical :: tagflag

    ! Find or create correct section

    if (.not.associated(config)) allocate(config)

    found=>config
    do
       if (associated(found)) then
          if (present(tag)) then
             tagflag=ConfigSectionHasTag(found,tag)
          else
             tagflag=.true.
          end if
          if ((trim(secname)==trim(found%name)).and.tagflag) then
             exit
          else
             if (associated(found%next)) then
                found=>found%next
             else
                allocate(newsec)
                found%next=>newsec
                found=>found%next
                found%name=trim(secname)
                if (present(tag)) then
                   allocate(newtag)
                   newtag%name='tag'
                   newtag%value=trim(tag)
                   found%values=>newtag
                end if
                exit
             end if
          end if
       else
          exit
       end if
    end do

    ! Add or create key-value pair

    if (.not.associated(found%values)) then
       allocate(newval)
       found%values=>newval
       found%values%name=valname
       found%values%value=value
    else
       val=>found%values
       do
          if (trim(valname)==trim(val%name)) then
             val%value=value
             exit
          else
             if (associated(val%next)) then
                val=>val%next
             else
                allocate(newval)
                val%next=>newval
                val%next%name=valname
                val%next%value=value
                exit
             end if
          end if
       end do
    end if

  end subroutine ConfigSetValueSec
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueDoubleArray(section,name,val,numval)
    ! get real array value
    use log_tools
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(kind=dp), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables..
    character(len=valuelen) :: value,tmp
    real(kind=dp), dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do
    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return

10  call write_log('GetValueDoubleArray: Array error in config file - check syntax' , &
                      msg_fatal , __FILE__ , __LINE__ )

  end subroutine GetValueDoubleArray
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueRealArray(section,name,val,numval)
    ! get real array value
    use log_tools
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real, pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables..
    character(len=valuelen) :: value,tmp
    real, dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return

10  call write_log('GetValueRealArray: Array error in config file - check syntax' , &
                     msg_fatal , __FILE__ , __LINE__ )

  end subroutine GetValueRealArray
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueIntArray(section,name,val,numval)
    ! get integer array value
    use log_tools
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    integer, pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables..
    character(len=valuelen) :: value,tmp
    integer, dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return

10  call write_log('GetValueIntArray: Array error in config file - check syntax' , &
                    msg_fatal , __FILE__ , __LINE__ )

  end subroutine GetValueIntArray
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueCharArray(section,name,val,numval)
    ! get character array value
    use log_tools
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    character(len=80), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables..
    character(len=valuelen) :: value
    character(80), dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tempval(i)=value(1:ind-1)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return

  end subroutine GetValueCharArray
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueReal(section,name,val)
    ! get real value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real :: val

    ! local variables..
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueReal
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueDouble(section,name,val)
    ! get double value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(kind=dp) :: val

    ! local variables..
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueDouble
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueInt(section,name,val)
    ! get integer value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    integer :: val

    ! local variables..
    character(len=valuelen) :: value
    integer temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueInt
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueChar(section,name,val)
    ! get character value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    character(len=*) :: val

    type(ConfigValue), pointer :: value

    value=>section%values
    do while(associated(value))
       if (name.eq.trim(value%name)) then
          val = value%value
          return
       end if
       value=>value%next
    end do
  end subroutine GetValueChar
  !
  !----------------------------------------------------------------------
  !
  subroutine GetValueLogical(section,name,val)
    ! get logical value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    logical :: val

    ! local variables..
    character(len=valuelen) :: value
    integer itemp
    logical ltemp
    integer ios

    value=''
    call GetValueChar(section,name,value)
    read(value,*,iostat=ios) itemp
    if (ios==0) then
       val = ( itemp .eq. 1 )
    end if

    read(value,*,iostat=ios) ltemp
    if (ios==0) then
       val = ltemp
    end if

  end subroutine GetValueLogical
  !
  !----------------------------------------------------------------------
  !
  subroutine handle_section(linenr,line,section)
    use log_tools
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigSection), pointer :: section

    ! local variables..
    integer i

    do i=1,linelen
       if (line(i:i).eq.']') then
          exit
       end if
    end do
    if (line(i:i).ne.']') then
       write(message,*) 'handle_section: Cannot find end of section ',linenr
       call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    end if

    call InsertSection(trim(adjustl(line(2:i-1))),section)
  end subroutine handle_section
  !
  !----------------------------------------------------------------------
  !  
  subroutine handle_value(linenr,line,value)
    use log_tools
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigValue), pointer :: value

    ! local variables..
    integer i

    do i=1,linelen
       if (line(i:i).eq.'=' .or. line(i:i).eq.':') then
          exit
       end if
    end do
    if (.not.(line(i:i).eq.'=' .or. line(i:i).eq.':')) then
       write(message,*) 'handle_value: Cannot find = or : ',linenr
       call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    end if

    call InsertValue(trim(adjustl(line(:i-1))), trim(adjustl(line(i+1:))),value)
  end subroutine handle_value
  !
  !----------------------------------------------------------------------
  !
  subroutine InsertSection(name,section)
    ! add a new section
    implicit none
    character(len=*), intent(in) :: name
    type(ConfigSection), pointer :: section
    type(ConfigSection), pointer :: new_sec

    allocate(new_sec)
    new_sec%name = name

    if (associated(section)) then
       if (associated(section%next)) then
          new_sec%next => section%next
       end if
       section%next=>new_sec
    end if
    section=>new_sec
  end subroutine InsertSection
  !
  !----------------------------------------------------------------------
  !
  subroutine InsertValue(name,val,value)
    ! add a new value
    implicit none
    character(len=*), intent(in) :: name, val
    type(ConfigValue), pointer :: value
    type(ConfigValue), pointer :: new_value

    allocate(new_value)

    new_value%name = name
    new_value%value = val

    if(associated(value)) then
       if (associated(value%next)) then
          new_value%next => value%next
       end if
       value%next => new_value
    end if
    value=>new_value
  end subroutine InsertValue
  !
  !----------------------------------------------------------------------
  !
end module config_tools
!
!------------------------------------------------------------------------
!
