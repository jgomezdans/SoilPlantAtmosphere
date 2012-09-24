! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_restart

  !! > Module description < !!

  use linked_lists
  use netcdf_tools

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: restart_model, write_restart_file
  ! Variables..
  public :: prev_restart_write, year_restart

  integer,parameter   :: restart_in_unit    = 50 ! unit nos for input restart file
  integer,parameter   :: restart_out_unit   = 51 ! unit nos for output restart file
  integer             :: prev_restart_write = 0. ! Year of last write to restart file
  integer,allocatable :: year_restart(:)         ! Year restart file applies to

  ! linked-list of variables (and required dimensions) in restart file..
  type(spalist),pointer,save :: restart_var_list , restart_dim_list

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine restart_model( infile )

    ! Read the restart file and load the data !
    ! into the relevant variables.            !

    use gv_scale_declarations, only: grid, time
    use gv_veg
    use linked_lists,          only: spaitem
    use log_tools
    use netcdf_tools
    use soil_air

    implicit none

    ! arguments..
    character(len=*),intent(in) :: infile

    ! local variables..
    type(nc_header),pointer :: nc_in_header
    type(spaitem),pointer   :: p
    real                    :: dummy

    ! If restart-list doesn't exist yet, create it..
    if ( .not. associated(restart_var_list) ) call setup_restart_lists()

    ! open the (NetCDF-format) restart file..
    allocate( nc_in_header )
    nc_in_header%name = trim(infile)
    call open_nc_file( nc_in_header )

    ! check that the dimensions declared in the file
    ! match those currently in use in SPA..
    call check_dims_match( nc_in_header%id )

    p => restart_var_list%first_ptr
    do
      ! check if we have reached the end of the list yet..
      if ( .not. associated(p) ) exit
      ! ..if not, then get that var from the file..
      call get_restart_var( nc_in_header%id , p )
      p => p%next_ptr
    enddo

    ! We do need to check the parameters loaded from the config file
    ! match those in the restart file.
    ! So, how will we write them to the restart file?

    ! Allocate variables needed later..
    if ( .not. allocated(      ess ) ) allocate(      ess( time%steps_per_day )                     )
    if ( .not. allocated(     flux ) ) allocate(     flux( time%steps_per_day , grid%canopy ) )
    if ( .not. allocated(     gppt ) ) allocate(     gppt( time%steps_per_day )                     )
    if ( .not. allocated(    respt ) ) allocate(    respt( time%steps_per_day )                     )
    if ( .not. allocated( soiletmm ) ) allocate( soiletmm( time%steps_per_day )                     )
    if ( .not. allocated(   transt ) ) allocate(   transt( time%steps_per_day )                     )

    ! initialise some soil-related bits..
    call saxton_parameters()
    call water_retention()
    call soil_resistance()
    call soil_water_potential()
    call water_uptake_layer( dummy )

    ! Set the starting year to be that of the restart file..
    call get_global_att( nc_in_header%id , 'SPA year completed' , year_restart )
    time%year          = year_restart(1)
    prev_restart_write = time%year

  end subroutine restart_model
  !
  !----------------------------------------------------------------------
  !
  subroutine write_restart_file( time, plant_func_type )

    use gv_scale_declarations, only: fname_length,time_holder,user_opts
    use log_tools
    use netcdf_tools,          only: change_mode,close_nc_output,nc_header

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    logical, intent(in)          :: plant_func_type

    ! local variables..
    character(fname_length) :: restart_name
    type(nc_header),pointer :: nc_out => null()

    allocate(nc_out)

    ! If restart-list doesn't exist yet, create it..
    if ( .not. associated(restart_var_list) ) call setup_restart_lists()

    ! Create name of output file..
    if ( plant_func_type ) then
      ! deciduous..
      write(restart_name,'(a18,i0.2,a3)') "deciduous.restart.",time%year,".nc"
    else
      ! evergreen..
      write(restart_name,'(a18,i0.2,a3)') "evergreen.restart.",time%year,".nc"
    endif
    restart_name = trim(adjustl(user_opts%output_directory)) &
                     //trim(adjustl(restart_name))
    write(message,*)"Creating restart file: "//trim(restart_name)
    call write_log( message , msg_info , __FILE__ , __LINE__ )

    ! Open netcdf file for output..
    nc_out%name = trim(restart_name)
    call open_nc_file( nc_out , for_output=.True. , append=.False. )

    ! Set SPA-time as a global attribute..
    call put_global_att( nc_out%id ,  'SPA year completed' , (/ time%year /) )

    ! Create new dimensions..
    call create_nc_dims_from_list( nc_out%id , restart_dim_list )

    ! Create variables..
    call create_nc_vars_from_list( nc_out%id , restart_var_list , restart_dim_list )

    ! Switch netcdf file into write mode..
    call change_mode( nc_out , .False. )

    ! Write the data to the files..
    call write_restart_nc( nc_out%id , restart_dim_list , restart_var_list )

    ! Close the restart file..
    call close_nc_output( nc_out )

    ! Finished with the pointer..
    deallocate( nc_out )

    ! Restart the model from this dump..
    !  ( This means that whatever error was made in writing )
    !  ( the dump file will be incorporated into this run.  )
    call restart_model( restart_name )

  end subroutine write_restart_file
  !
  !----------------------------------------------------------------------
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !----------------------------------------------------------------------
  !
  subroutine check_dims_match( file_id )

    ! Check that the dimensions in a given netcdf match !
    ! those within SPA.  At the same time, populate the !
    ! the netcdf info (dim_id and var_id) part of the   ! 
    ! restart_dims_list structure.                      !

    use linked_lists
    use log_tools
    use netcdf_tools

    implicit none

    ! arguments..
    integer,intent(in) :: file_id

    ! local variables..
    type(spaitem),pointer    :: p, tmp

    ! Work through the restart dims list..
    p => restart_dim_list%first_ptr
    do
      ! finish when none left..
      if (.not. associated( p ) ) exit
      ! Work based on data type.  note vector => coord dim (unusual in SPA)
      allocate(tmp)
      tmp%name   = p%name
      call get_dim_info( file_id , tmp )
      select case ( p%data_type )
      case (1)
        if ( p%int_scalar .ne. tmp%int_scalar ) then
          write(message,*) "Dimension ",trim(p%name)," in input file does " &
                         //"not have same size as current SPA setup"
          call write_log( message , msg_fatal , __FILE__ , __LINE__ )
        end if
      case (2)
        ! check that the length of the dimension in the file
        !  matches the size we have defined it to be in the model.
        if ( size(p%int_vector) .ne. size(tmp%int_vector) ) then
          write(message,*) "Dimension ",trim(p%name)," in input file does " &
                         //"not have same size as current SPA setup"
          call write_log( message , msg_fatal , __FILE__ , __LINE__ )
        end if
      case (4)
        if ( size(p%real_vector) .ne. size(tmp%real_vector) ) then
          write(message,*) "Dimension ",trim(p%name)," in input file does " &
                         //"not have same size as current SPA setup"
          call write_log( message , msg_fatal , __FILE__ , __LINE__ )
        end if
      end select
      p => p%next_ptr
      deallocate(tmp)
    end do

  end subroutine check_dims_match
  !
  !----------------------------------------------------------------------
  !
  subroutine create_nc_dims_from_list( file_id , dimlist )

    ! Define dimensions in the output file (ie give !
    !  them properties, but no data).               !
    ! Dimensions can be scalar (ie simply a length) !
    !  or vector (a coordinate dimension).  If      !
    !  scalar they must be integer, but if vector   !
    !  their contents can be integer or real.       !

    use netcdf,       only: nf90_unlimited
    use netcdf_tools, only: create_dim, create_coord_dim

    implicit none

    ! arguments..
    integer,intent(in)    :: file_id
    type(spalist),pointer :: dimlist

    ! local variables..
    type(spaitem),pointer :: p

    p => dimlist%first_ptr
    do
      ! finish when none left..
      if (.not. associated( p ) ) exit
      select case (p%data_type)
      case (1)
        call create_dim( file_id , p%name, p%int_scalar , p%dim_id )
      case (2)
        if ( index( p%name , 'time' ) .ne. 0 ) then
          ! this dimension is related to 'time' in some way
          ! and so is (probably) going to be unlimited..
          call create_coord_dim( file_id , p%name, nf90_unlimited , &
                                           p%dim_id , 1 , p%var_id )
        else
          ! a 'normal' dimension..
          call create_coord_dim( file_id , p%name, size(p%int_vector) , &
                                               p%dim_id , 1 , p%var_id )
        endif
      case (4)
        if ( index( p%name , 'time' ) .ne. 0 ) then
          ! this dimension is related to 'time' in some way
          ! and so is (probably) going to be unlimited..
          call create_coord_dim( file_id , p%name, nf90_unlimited , &
                                           p%dim_id , 2 , p%var_id )
        else
          ! a 'normal' dimension..
          call create_coord_dim( file_id , p%name, size(p%real_vector) , &
                                                p%dim_id , 2 , p%var_id )
        endif
      end select
      p => p%next_ptr
    enddo

  end subroutine create_nc_dims_from_list
  !
  !----------------------------------------------------------------------
  !
  subroutine create_nc_vars_from_list( file_id , varlist , dimlist )

    ! Define variables in the output file (ie give !
    !  them properties, but no data).              !
    ! Where variables have a dimension (ie are not !
    !  scalar) then supply their dim-id.           !

    use linked_lists
    use log_tools
    use netcdf_tools

    implicit none

    ! arguments..
    integer,intent(in)    :: file_id
    type(spalist),pointer :: dimlist, varlist

    ! local variables..
    integer               :: nc_type
    type(spaitem),pointer :: d, p
    integer,allocatable   :: dimnos(:)

    p => varlist%first_ptr
    do
      ! finish when none left..
      if (.not. associated( p ) ) exit
      write(message,*) "Creating: ",p%name
      call write_log( message, msg_info )

      select case (p%data_type)
      case (1,2) ! integer
        nc_type = 1
      case (3,4,5,6) ! real
        nc_type = 2
      end select

      ! For each var-pointer we know if it has a dimension 
      !  or not because the dim_name will be non-blank..
      if ( p%dim_name .eq. '' ) then
        allocate(dimnos(0))
        call create_var( file_id , p%name , dimnos , nc_type , p%var_id )
        deallocate(dimnos)
      else
        ! find out what dim-id is associated with that dimension name..
        ! we need to be able to get an item from the list..
        d => get_item( dimlist , p%dim_name )
        call create_var( file_id , p%name , (/ d%dim_id /) , nc_type , p%var_id )
      endif
      p => p%next_ptr
    enddo

  end subroutine create_nc_vars_from_list
  !
  !----------------------------------------------------------------------
  !
  subroutine get_restart_var( fileid , restart_var )

    !>  <!

    use log_tools

    implicit none

    ! arguments..
    integer, intent(in)   :: fileid
    type(spaitem),pointer :: restart_var

    ! local variables
    type(spaitem),pointer :: tmp_var => null()

    allocate(tmp_var)
    call write_log( "Loading : "//trim(restart_var%name) )
    tmp_var%name  = restart_var%name
    tmp_var%ndims = restart_var%ndims
    call get_nc_var( fileid , tmp_var )

    ! Check that the variable's datatypes match those expected..
    if ( tmp_var%data_type .ne. restart_var%data_type ) then
      write(message,*) "The datatype of the variable in the file" &
                     //  " do not match those expected."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Load 
    select case( restart_var%data_type )
      case (1)
        restart_var%int_scalar  = tmp_var%int_scalar
      case (2)
        restart_var%int_vector  = tmp_var%int_vector
      case (3)
        restart_var%real_scalar = tmp_var%real_scalar
      case (4)
        restart_var%real_vector = tmp_var%real_vector
    end select

  end subroutine get_restart_var
  !
  !----------------------------------------------------------------------
  !
  subroutine setup_restart_lists()

    ! This subroutine creates a list of the variables that !
    ! are needed for a restart file from the big list of   !
    ! SPA variables.                                       !

    use linked_lists, only: append, dims_from_var_list, get_item, &
                            print_item, spa_dims, spa_vars, spaitem
    use log_tools

    implicit none

    ! Check that the main SPA lists (from which we subset) have
    !  already been setup..
    if ( ( .not. associated(spa_vars) ) .or. & 
                                  ( .not. associated(spa_dims ) ) ) then
      write(message,*) "This routine requires that setup_spa_lists "&
                     //"be run at some point prior to the call."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Build the list of restart variables..
    if ( .not. associated(restart_var_list)) allocate( restart_var_list )
    call append( restart_var_list , get_item( spa_vars , 'avN'                               ) )
    call append( restart_var_list , get_item( spa_vars , 'alloc_to_foliage'                  ) )
    call append( restart_var_list , get_item( spa_vars , 'alloc_to_labile'                   ) )
    call append( restart_var_list , get_item( spa_vars , 'alloc_from_labile'                 ) )
    call append( restart_var_list , get_item( spa_vars , 'alloc_to_roots'                    ) )
    call append( restart_var_list , get_item( spa_vars , 'alloc_to_stem'                ) )
    call append( restart_var_list , get_item( spa_vars , 'canopy_height'                     ) )
    call append( restart_var_list , get_item( spa_vars , 'canopy_store'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'capac'                             ) )
    call append( restart_var_list , get_item( spa_vars , 'conduc'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'conductivity'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'daytempsum'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'decomposition'                     ) )
    call append( restart_var_list , get_item( spa_vars , 'decomposition_rate'                ) )
    call append( restart_var_list , get_item( spa_vars , 'dimen'                             ) )
    call append( restart_var_list , get_item( spa_vars , 'drythick'                          ) )
    call append( restart_var_list , get_item( spa_vars , 'frac_alloc_foliage'                ) )
    call append( restart_var_list , get_item( spa_vars , 'frac_alloc_roots'                  ) )
    call append( restart_var_list , get_item( spa_vars , 'frac_GPP_resp_auto'                ) )
    call append( restart_var_list , get_item( spa_vars , 'frac_leafLoss_litter'              ) )
    call append( restart_var_list , get_item( spa_vars , 'GDD'                               ) )
    call append( restart_var_list , get_item( spa_vars , 'GDD_thresh'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'gplant'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'GPP'                               ) )
    call append( restart_var_list , get_item( spa_vars , 'hourts'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'iceprop'                           ) )
    call append( restart_var_list , get_item( spa_vars , 'iota'                              ) )
    call append( restart_var_list , get_item( spa_vars , 'kappaC'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'kappaJ'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'lafrac'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'lai'                               ) )
    call append( restart_var_list , get_item( spa_vars , 'lat'                               ) )
    call append( restart_var_list , get_item( spa_vars , 'layer_depth'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'layer_height'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'litterfall_foliage'                ) )
    call append( restart_var_list , get_item( spa_vars , 'litterfall_stem'              ) )
    call append( restart_var_list , get_item( spa_vars , 'litterfall_roots'                  ) )
    call append( restart_var_list , get_item( spa_vars , 'LMA'                               ) )
    call append( restart_var_list , get_item( spa_vars , 'LWPstore'                          ) )
    call append( restart_var_list , get_item( spa_vars , 'max_depth'                         ) )
    call append( restart_var_list , get_item( spa_vars , 'max_fol'                           ) )
    call append( restart_var_list , get_item( spa_vars , 'max_stock_foliage'                 ) )
    call append( restart_var_list , get_item( spa_vars , 'max_storage'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'mineralfrac'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'mineralisation_rate_soilOrgMatter' ) )
    call append( restart_var_list , get_item( spa_vars , 'minlwp'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'min_Temp_thresh'                   ) )
    call append( restart_var_list , get_item( spa_vars , 'nfrac'                             ) )
    call append( restart_var_list , get_item( spa_vars , 'organicfrac'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'resp_auto'                         ) )
    call append( restart_var_list , get_item( spa_vars , 'resp_cost_labile_trans'            ) )
    call append( restart_var_list , get_item( spa_vars , 'resp_h_litter'                     ) )
    call append( restart_var_list , get_item( spa_vars , 'resp_h_soilOrgMatter'              ) )
    call append( restart_var_list , get_item( spa_vars , 'resp_rate_temp_coeff'              ) )
    call append( restart_var_list , get_item( spa_vars , 'rooted_layers'                     ) )
    call append( restart_var_list , get_item( spa_vars , 'root_length'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'root_k'                            ) )
    call append( restart_var_list , get_item( spa_vars , 'root_mass'                         ) )
    call append( restart_var_list , get_item( spa_vars , 'root_radius'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'rootresist'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'snowheight'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'snowweight'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'soil_frac_clay'                    ) )
    call append( restart_var_list , get_item( spa_vars , 'soil_frac_sand'                    ) )
    call append( restart_var_list , get_item( spa_vars , 'soil_temp'                         ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_foliage'                     ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_labile'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_litter'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_resp_auto'                   ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_roots'                       ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_soilOrgMatter'               ) )
    call append( restart_var_list , get_item( spa_vars , 'stock_stem'                   ) )
    call append( restart_var_list , get_item( spa_vars , 'thickness'                         ) )
    call append( restart_var_list , get_item( spa_vars , 'through_fall'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'tower_height'                      ) )
    call append( restart_var_list , get_item( spa_vars , 'turnover_rate_foliage'             ) )
    call append( restart_var_list , get_item( spa_vars , 'turnover_rate_labile'              ) )
    call append( restart_var_list , get_item( spa_vars , 'mineralisation_rate_litter'        ) )
    call append( restart_var_list , get_item( spa_vars , 'turnover_rate_resp_auto'           ) )
    call append( restart_var_list , get_item( spa_vars , 'turnover_rate_roots'               ) )
    call append( restart_var_list , get_item( spa_vars , 'turnover_rate_stem'           ) )
    call append( restart_var_list , get_item( spa_vars , 'waterfrac'                         ) )
    call append( restart_var_list , get_item( spa_vars , 'watericemm'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'wettingbot'                        ) )
    call append( restart_var_list , get_item( spa_vars , 'wettingtop'                        ) )

    ! For a netcdf file we also need to have a list of the
    !  dimensions that are associated with these variables..
    if ( .not. associated( restart_dim_list ) ) allocate( restart_dim_list )
    call dims_from_var_list( spa_dims , restart_var_list , restart_dim_list )

  end subroutine setup_restart_lists
  !
  !----------------------------------------------------------------------
  !
  subroutine write_restart_nc( fileid , dimlist , varlist )

    ! This writes the data to the coordinate-dimensions ! 
    !  and variables.                                   !

    use linked_lists
    use log_tools
    use netcdf_tools

    implicit none

    ! arguments..
    integer,intent(in)    :: fileid
    type(spalist),pointer :: dimlist, varlist

    ! local variables..
    type(spaitem),pointer :: p

    p => dimlist%first_ptr
    do
      if (.not. associated(p)) exit
      ! only put data in for datatypes of coordinate dims..
      select case(p%data_type)
      case (2,4) ! integer/real vector
        call write_log( "Writing  : "//trim(p%name ) )
        call put_var_nc( fileid , p )
      end select
      p => p%next_ptr
    enddo

    p => varlist%first_ptr
    do
      if (.not. associated(p)) exit
      call write_log( "Writing  : "//trim(p%name ) )
      call put_var_nc( fileid , p )
      p => p%next_ptr
    enddo

  end subroutine write_restart_nc
  !
  !----------------------------------------------------------------------
  !
end module spa_restart
