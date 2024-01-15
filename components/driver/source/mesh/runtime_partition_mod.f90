!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Functions/Subroutines specific to code path that partitions
!>        global_mesh_type objects at runtime.(Best endeavours)
module runtime_partition_mod

  use constants_mod,           only: i_def, r_def, str_def, l_def, &
                                     str_max_filename
  use global_mesh_mod,         only: global_mesh_type
  use log_mod,                 only: log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR,   &
                                     LOG_LEVEL_INFO
  use local_mesh_mod,          only: local_mesh_type
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type
  use ncdf_quad_mod,           only: ncdf_quad_type

  use partition_mod, only: partition_type,                 &
                           partitioner_interface,          &
                           partitioner_cubedsphere_serial, &
                           partitioner_cubedsphere,        &
                           partitioner_planar

  use local_mesh_collection_mod,  only: local_mesh_collection
  use global_mesh_collection_mod, only: global_mesh_collection

  !------------------------------
  ! Configuration modules
  !------------------------------
  use base_mesh_config_mod, only: GEOMETRY_SPHERICAL,      &
                                  GEOMETRY_PLANAR,         &
                                  TOPOLOGY_FULLY_PERIODIC, &
                                  TOPOLOGY_NON_PERIODIC

  use partitioning_config_mod, only: PANEL_DECOMPOSITION_AUTO,   &
                                     PANEL_DECOMPOSITION_ROW,    &
                                     PANEL_DECOMPOSITION_COLUMN, &
                                     PANEL_DECOMPOSITION_CUSTOM

  implicit none

  private
  public :: get_partition_parameters
  public :: create_local_mesh
  public :: create_local_mesh_maps

contains

!> @brief Sets common partition parameters to be applied to global meshes.
!>
!> @param[in]   total_ranks      Total number of MPI ranks in this job
!> @param[out]  xproc            Number of ranks in mesh panel x-direction
!> @param[out]  yproc            Number of ranks in mesh panel y-direction
!> @param[out]  partitioner_ptr  Mesh partitioning strategy
subroutine get_partition_parameters( configuration, total_ranks, &
                                     xproc, yproc, partitioner_ptr )

  implicit none

  type(namelist_collection_type), intent(in) :: configuration

  integer(i_def), intent(in)  :: total_ranks
  integer(i_def), intent(out) :: xproc
  integer(i_def), intent(out) :: yproc

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  ! Locals
  integer(i_def) :: ranks_per_panel
  integer(i_def) :: start_factor
  integer(i_def) :: end_factor
  integer(i_def) :: fact_count
  logical(l_def) :: found_factors

  character(len=str_def) :: domain_desc

  integer(i_def), parameter :: max_factor_iters = 10000

  type(namelist_type), pointer :: base_mesh_nml    => null()
  type(namelist_type), pointer :: partitioning_nml => null()

  integer(i_def) :: panel_decomposition
  integer(i_def) :: panel_xproc
  integer(i_def) :: panel_yproc
  integer(i_def) :: geometry
  integer(i_def) :: topology


  !============================================================================
  ! 0.0 Extract configuration variables
  !============================================================================
  base_mesh_nml    => configuration%get_namelist('base_mesh')
  partitioning_nml => configuration%get_namelist('partitioning')

  call base_mesh_nml%get_value( 'geometry', geometry )
  call base_mesh_nml%get_value( 'topology', topology )
  call partitioning_nml%get_value( 'panel_xproc', panel_xproc )
  call partitioning_nml%get_value( 'panel_yproc', panel_yproc )
  call partitioning_nml%get_value( 'panel_decomposition', panel_decomposition )

  partitioning_nml => null()
  base_mesh_nml    => null()

  !============================================================================
  ! 0.1 Initalise variables
  !============================================================================
  partitioner_ptr => null()


  ! 1.0 Setup the partitioning strategy
  !===================================================================
  if ( geometry == geometry_spherical .and. &
       topology == topology_fully_periodic ) then

    ! Assume that we have a cubed sphere (and not a global lon-lat mesh)
    if (total_ranks == 1 .or. mod(total_ranks,6) == 0) then

      ranks_per_panel = total_ranks/6
      domain_desc = "6x"

      if (total_ranks == 1) then
        ! Serial run job
        ranks_per_panel = 1
        partitioner_ptr => partitioner_cubedsphere_serial
        call log_event( "Using serial cubed sphere partitioner", &
                        LOG_LEVEL_INFO )

      else
        ! Paralled run job
        partitioner_ptr => partitioner_cubedsphere
        call log_event( "Using parallel cubed sphere partitioner", &
                        LOG_LEVEL_INFO )
      end if

    else
      call log_event( "Total number of processors must be 1 (serial) "// &
                      "or a multiple of 6 for a cubed-sphere domain.",   &
                      LOG_LEVEL_ERROR )
    end if

  else ! Planar/LAM mesh

    ranks_per_panel = total_ranks
    domain_desc = ""

    partitioner_ptr => partitioner_planar
    call log_event( "Using planar mesh partitioner ", &
                    LOG_LEVEL_INFO )
  end if

  ! 2.0 Setup Panel decomposition
  !===================================================================
  select case ( panel_decomposition )

  case( PANEL_DECOMPOSITION_AUTO )

    ! For automatic partitioning, attempt to partition into the squarest
    ! possible partitions by finding the two factors of ranks_per_panel
    ! that are closest to sqrt(ranks_per_panel).
    start_factor  = nint(sqrt(real(ranks_per_panel, kind=r_def)), kind=i_def)
    end_factor    = max(1,(start_factor-max_factor_iters))
    found_factors = .false.

    do fact_count=start_factor, end_factor, -1
      if (mod(ranks_per_panel,fact_count) == 0) then
        found_factors = .true.
        exit
      end if
    end do

    if (found_factors) then
      xproc = fact_count
      yproc = ranks_per_panel/fact_count
    else
      call log_event( "Could not automatically partition domain.", &
                      LOG_LEVEL_ERROR )
    end if

  case( PANEL_DECOMPOSITION_ROW )
    xproc = ranks_per_panel
    yproc = 1

  case( PANEL_DECOMPOSITION_COLUMN )
    xproc = 1
    yproc = ranks_per_panel

  case( PANEL_DECOMPOSITION_CUSTOM )
    ! Use the values provided from the partitioning namelist
    xproc = panel_xproc
    yproc = panel_yproc

    if (xproc*yproc /= ranks_per_panel) then
      call log_event( "The values of panel_xproc and panel_yproc "// &
                      "are inconsistent with the total number of "// &
                      "processors available.", LOG_LEVEL_ERROR )
    end if

  case default

    call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

  end select

  if (total_ranks > 1) then
    write(log_scratch_space, '(2(A,I0))' ) &
        'Panel decomposition: ', xproc, 'x', yproc
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end if

end subroutine get_partition_parameters


!> @brief  Loads the given list of global meshes names, partitions them
!!         and creates local meshes from them.
!>
!> @param[in]  input_mesh_file    Input file to load meshes from.
!> @param[in]  mesh_names[:]      Array of requested mesh names to load
!!                                from the mesh input file
!> @param[in]  local_rank         Number of the local MPI rank
!> @param[in]  total_ranks        Total number of MPI ranks in this job
!> @param[in]  xproc              Number of ranks in mesh panel x-direction
!> @param[in]  yproc              Number of ranks in mesh panel y-direction
!> @param[in]  stencil_depth      Depth of cells outside the base cell
!!                                of stencil.
!> @param[in]  partitioner_ptr    Mesh partitioning strategy
subroutine create_local_mesh( mesh_names,              &
                              local_rank, total_ranks, &
                              xproc, yproc,            &
                              stencil_depth,           &
                              partitioner_ptr )

  implicit none

  character(len=str_def), intent(in) :: mesh_names(:)

  integer(i_def), intent(in) :: local_rank
  integer(i_def), intent(in) :: total_ranks
  integer(i_def), intent(in) :: xproc
  integer(i_def), intent(in) :: yproc
  integer(i_def), intent(in) :: stencil_depth

  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr

  type(global_mesh_type), pointer :: global_mesh_ptr
  type(partition_type)            :: partition
  type(local_mesh_type)           :: local_mesh

  integer(i_def) :: local_mesh_id, i

  do i=1, size(mesh_names)

    global_mesh_ptr => global_mesh_collection%get_global_mesh( mesh_names(i) )

    ! Create partition
    partition = partition_type( global_mesh_ptr, &
                                partitioner_ptr, &
                                xproc, yproc,    &
                                stencil_depth,   &
                                local_rank, total_ranks )
    ! Create local_mesh
    call local_mesh%initialise( global_mesh_ptr, partition )

    ! Make sure the local_mesh cell owner lookup is correct
    ! (Can only be done when the code is running on its full set of MPI tasks)
    call local_mesh%init_cell_owner()
    local_mesh_id = local_mesh_collection%add_new_local_mesh( local_mesh )

  end do

end subroutine create_local_mesh


!> @brief    Creates the local mesh intergrid maps from available global
!!           intergrid maps from file.
!> @details  Global meshes which have been read into the model's global mesh
!!           collection will have a list of target mesh names. These target mesh
!!           names (if any) indicate the valid intergrid maps available in the
!!           mesh file. This routine will read in the appropriate intergrid
!!           maps convert the LiD-LiD map them to the appropriate local
!!           mesh object.
!!
!> @param[in]  input_mesh_file  Input file to load mesh maps from.
subroutine create_local_mesh_maps( input_mesh_file )

  implicit none

  character(len=str_max_filename) :: input_mesh_file

  type(ncdf_quad_type) :: file_handler

  character(str_def), allocatable :: source_mesh_names(:)
  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: gid_mesh_map(:,:,:)
  integer(i_def),     allocatable :: lid_mesh_map(:,:,:)

  integer(i_def) :: i, j, n, x, y
  integer(i_def) :: n_meshes

  type(global_mesh_type), pointer :: source_global_mesh => null()

  type(local_mesh_type), pointer :: source_local_mesh => null()
  type(local_mesh_type), pointer :: target_local_mesh => null()

  integer(i_def) :: ntarget_per_source_cell_x, ntarget_per_source_cell_y
  integer(i_def) :: ncells
  integer(i_def) :: target_local_mesh_id

  ! Read in the maps for each global mesh
  !=================================================================
  call file_handler%file_open(trim(input_mesh_file))

  allocate( source_mesh_names, &
            source=global_mesh_collection%get_mesh_names() )
  n_meshes = global_mesh_collection%n_meshes()

  ! Loop over every source mesh
  do i=1, n_meshes

    ! Get the global and local source mesh
    source_global_mesh => &
        global_mesh_collection%get_global_mesh( source_mesh_names(i) )
    source_local_mesh => &
        local_mesh_collection%get_local_mesh( source_mesh_names(i) )
    call source_global_mesh%get_target_mesh_names( target_mesh_names )

    if (allocated(target_mesh_names)) then

      ! Loop over each target mesh
      do j=1, size(target_mesh_names)
        target_local_mesh => &
           local_mesh_collection%get_local_mesh( target_mesh_names(j) )

        if ( associated(target_local_mesh) ) then

          ! Read in the global mesh map
          call file_handler%read_map( source_mesh_names(i), &
                                      target_mesh_names(j), &
                                      gid_mesh_map )

          ! Create the local mesh map
          ntarget_per_source_cell_x = size(gid_mesh_map, 1)
          ntarget_per_source_cell_y = size(gid_mesh_map, 2)
          ncells = source_local_mesh%get_num_cells_in_layer()
          allocate( lid_mesh_map( ntarget_per_source_cell_x, &
                                  ntarget_per_source_cell_y, &
                                  ncells ) )

          ! Convert global cell IDs in the global mesh map
          ! into local cell IDs in a local mesh map
          do x=1, ntarget_per_source_cell_x
            do y=1, ntarget_per_source_cell_y
              do n=1, ncells
                lid_mesh_map(x, y, n) = target_local_mesh%get_lid_from_gid( &
                    gid_mesh_map(x, y, source_local_mesh%get_gid_from_lid(n)) )
              end do
            end do
          end do

          ! Put the local mesh map in the local mesh
          target_local_mesh_id = target_local_mesh%get_id()
          call source_local_mesh%add_local_mesh_map( target_local_mesh_id, &
                                                     lid_mesh_map )

          if ( allocated(gid_mesh_map) ) deallocate( gid_mesh_map )
          if ( allocated(lid_mesh_map) ) deallocate( lid_mesh_map )

        end if

      end do

      if ( allocated( target_mesh_names ) ) then
        deallocate( target_mesh_names )
      end if
    end if
  end do

  if ( allocated( source_mesh_names ) ) then
    deallocate( source_mesh_names)
  end if

  call file_handler%file_close()

  return
end subroutine create_local_mesh_maps

end module runtime_partition_mod
