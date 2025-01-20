!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Functions/Subroutines specific to code path that partitions
!>        global_mesh_type objects at runtime.
module runtime_partition_mod

  use constants_mod,           only: i_def, r_def, str_def, l_def, &
                                     str_max_filename
  use global_mesh_mod,         only: global_mesh_type
  use log_mod,                 only: log_event,         &
                                     log_scratch_space, &
                                     log_level_error,   &
                                     log_level_debug
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

  implicit none

  private
  public :: get_partition_parameters
  public :: create_local_mesh
  public :: create_local_mesh_maps

  integer, public, parameter :: mesh_cubedsphere = 34
  integer, public, parameter :: mesh_planar      = 28

  integer, public, parameter :: decomp_auto   = 30
  integer, public, parameter :: decomp_row    = 68
  integer, public, parameter :: decomp_column = 12

  interface get_partition_parameters
    procedure get_partition_parameters_auto
    procedure get_partition_parameters_custom
  end interface get_partition_parameters

contains

!> @brief Determines the common parameters for partitioning global mesh domains.
!>
!> This routine supports only 'Planar' or 'Cubed-Sphere` mesh types specified
!> using this module's enumeration parameters, [mesh_cubedsphere|mesh_planar].
!>
!> Individual partitions are restricted to a "retangular" mesh panel.
!> Decomposition of partitions across a given panel are specified using this
!> module's enumeration parameters, [decomp_auto|decomp_row|decomp_column].
!>
!> @param[in]   mesh_selection      Choice of supported mesh
!> @param[in]   panel_decomposition Choice of panel decomposition
!> @param[in]   total_ranks         Total number of MPI ranks
!> @param[out]  panel_xproc         Ranks in mesh panel local x-direction
!> @param[out]  panel_yproc         Ranks in mesh panel local y-direction
!> @param[out]  partitioner_ptr     Mesh partitioning strategy
!=============================================================================
subroutine get_partition_parameters_auto( mesh_selection,      &
                                          panel_decomposition, &
                                          total_ranks,         &
                                          panel_xproc,         &
                                          panel_yproc,         &
                                          partitioner_ptr )

  implicit none

  integer, intent(in) :: mesh_selection
  integer, intent(in) :: panel_decomposition

  integer(i_def), intent(in)  :: total_ranks

  integer(i_def), intent(out) :: panel_xproc
  integer(i_def), intent(out) :: panel_yproc

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  ! Locals
  integer(i_def) :: ranks_per_panel
  integer(i_def) :: start_factor
  integer(i_def) :: end_factor
  integer(i_def) :: fact_count
  logical(l_def) :: found_factors

  integer(i_def), parameter :: max_factor_iters = 10000

  !===========================================================================
  ! Partition strategy
  !===========================================================================
  call get_partition_strategy( mesh_selection, total_ranks, &
                               ranks_per_panel, partitioner_ptr )


  !===========================================================================
  ! Panel decomposition
  !===========================================================================
  select case (panel_decomposition)

  case (decomp_auto)

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
      panel_xproc = fact_count
      panel_yproc = ranks_per_panel/fact_count
    else
      call log_event( "Could not automatically partition domain.", &
                      log_level_error )
    end if

  case (decomp_row)
    panel_xproc = ranks_per_panel
    panel_yproc = 1

  case (decomp_column)
    panel_xproc = 1
    panel_yproc = ranks_per_panel

  case default
    call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", log_level_error )

  end select

  if (total_ranks > 1) then
    write(log_scratch_space, '(2(A,I0))' ) &
        'Panel decomposition: ', panel_xproc, 'x', panel_yproc
    call log_event( log_scratch_space, log_level_debug )
  end if

end subroutine get_partition_parameters_auto


!> @brief Determines the common parameters for partitioning global mesh domains.
!>
!> This routine supports only 'Planar' or 'Cubed-Sphere` mesh types specified
!> using this module's enumeration parameters, [mesh_cubedsphere|mesh_planar].
!>
!> @param[in]   mesh_selection     Choice of supported mesh
!> @param[in]   total_ranks        Total number of MPI ranks
!> @param[in]   panel_xproc        Ranks mesh panel local x-direction
!> @param[in]   panel_yproc        Ranks mesh panel local x-direction
!> @param[out]  partitioner_ptr    Mesh partitioning strategy
!=============================================================================
subroutine get_partition_parameters_custom( mesh_selection, total_ranks, &
                                            panel_xproc, panel_yproc,    &
                                            partitioner_ptr )

  implicit none

  integer, intent(in) :: mesh_selection

  integer(i_def), intent(in) :: total_ranks
  integer(i_def), intent(in) :: panel_xproc
  integer(i_def), intent(in) :: panel_yproc

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  ! Locals
  integer(i_def) :: ranks_per_panel

  !===========================================================================
  ! Partition strategy
  !===========================================================================
  call get_partition_strategy( mesh_selection, total_ranks, &
                               ranks_per_panel, partitioner_ptr )

  !===========================================================================
  ! Panel decomposition
  !===========================================================================
  if (panel_xproc*panel_yproc /= ranks_per_panel) then
     call log_event( "The values of panel_xproc and panel_yproc "// &
                     "are inconsistent with the total number of "// &
                     "processors available.", log_level_error )
  end if

  if (total_ranks > 1) then
    write(log_scratch_space, '(2(A,I0))' ) &
        'Panel decomposition: ', panel_xproc, 'x', panel_yproc
    call log_event( log_scratch_space, log_level_debug )
  end if

end subroutine get_partition_parameters_custom


!> @brief Determines the partition stratedy to used for supported meshes
!>
!> This routine supports only 'Planar' or 'Cubed-Sphere` mesh types specified
!> using this module's enumeration parameters, [mesh_cubedsphere|mesh_planar].
!>
!> @param[in]   mesh_selection     Choice of supported mesh
!> @param[in]   total_ranks        Total number of MPI ranks per mesh
!> @param[out]  ranks_per_panel    Ranks per mesh panel
!> @param[out]  partitioner_ptr    Mesh partitioning strategy
!=============================================================================
subroutine get_partition_strategy( mesh_selection, total_ranks, &
                                   ranks_per_panel, partitioner_ptr )

  implicit none

  integer, intent(in) :: mesh_selection

  integer(i_def), intent(in)  :: total_ranks
  integer(i_def), intent(out) :: ranks_per_panel

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  partitioner_ptr => null()

  ! Determine the partitioning strategy
  !===================================================================
  select case (mesh_selection)

  case (mesh_cubedsphere)

    if (total_ranks == 1 .or. mod(total_ranks,6) == 0) then

      ranks_per_panel = total_ranks/6

      if (total_ranks == 1) then
        ! Serial run job
        ranks_per_panel = 1
        partitioner_ptr => partitioner_cubedsphere_serial
        call log_event( "Using serial cubed sphere partitioner", &
                        log_level_debug )

      else
        ! Paralled run job
        partitioner_ptr => partitioner_cubedsphere
        call log_event( "Using parallel cubed sphere partitioner", &
                        log_level_debug )
      end if

    else
      call log_event( "Total number of processors must be 1 (serial) "// &
                      "or a multiple of 6 for a cubed-sphere domain.",   &
                      log_level_error )
    end if

  case (mesh_planar)

    ranks_per_panel = total_ranks
    partitioner_ptr => partitioner_planar
    call log_event( "Using planar mesh partitioner ", &
                    log_level_debug )

  end select

end subroutine get_partition_strategy


!> @brief  Loads the given list of global meshes names, partitions them
!!         and creates local meshes from them.
!>
!> @param[in]  input_mesh_file        Input file to load meshes from.
!> @param[in]  mesh_names[:]          Array of requested mesh names to load
!!                                    from the mesh input file
!> @param[in]  local_rank             Number of the local MPI rank
!> @param[in]  total_ranks            Total number of MPI ranks in this job
!> @param[in]  xproc                  Number of ranks in mesh panel x-direction
!> @param[in]  yproc                  Number of ranks in mesh panel y-direction
!> @param[in]  generate_inner_haloes  Generate inner halo regions
!!                                    to overlap comms & compute
!> @param[in]  stencil_depth          Depth of cells outside the base cell
!!                                    of stencil.
!> @param[in]  partitioner_ptr        Mesh partitioning strategy
subroutine create_local_mesh( mesh_names,              &
                              local_rank, total_ranks, &
                              xproc, yproc,            &
                              stencil_depth,           &
                              generate_inner_haloes,   &
                              partitioner_ptr )

  implicit none

  character(len=str_def), intent(in) :: mesh_names(:)

  integer(i_def), intent(in) :: local_rank
  integer(i_def), intent(in) :: total_ranks
  integer(i_def), intent(in) :: xproc
  integer(i_def), intent(in) :: yproc
  integer(i_def), intent(in) :: stencil_depth

  logical(l_def), intent(in) :: generate_inner_haloes

  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr

  type(global_mesh_type), pointer :: global_mesh_ptr
  type(partition_type)            :: partition
  type(local_mesh_type)           :: local_mesh

  integer(i_def) :: local_mesh_id, i

  do i=1, size(mesh_names)

    global_mesh_ptr => global_mesh_collection%get_global_mesh( mesh_names(i) )

    ! Create partition
    partition = partition_type( global_mesh_ptr,       &
                                partitioner_ptr,       &
                                xproc, yproc,          &
                                stencil_depth,         &
                                generate_inner_haloes, &
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
