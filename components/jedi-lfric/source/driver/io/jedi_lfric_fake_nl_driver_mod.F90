!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the fake nl DA model
!>
module jedi_lfric_fake_nl_driver_mod

  use add_mesh_map_mod,         only: assign_mesh_maps
  use checksum_alg_mod,         only: checksum_alg
  use constants_mod,            only: str_def, i_def, l_def, r_def
  use create_mesh_mod,          only: create_mesh
  use driver_model_data_mod,    only: model_data_type
  use driver_time_mod,          only: init_time, get_calendar
  use driver_mesh_mod,          only: init_mesh
  use driver_fem_mod,           only: init_fem, final_fem
  use inventory_by_mesh_mod,    only: inventory_by_mesh_type
  use field_collection_mod,     only: field_collection_type
  use field_mod,                only: field_type
  use jedi_lfric_fake_nl_init_mod, &
                                only: create_da_model_data, &
                                      initialise_da_model_data
  use log_mod,                  only: log_event,          &
                                      log_scratch_space,  &
                                      LOG_LEVEL_ERROR
  use mesh_mod,                 only: mesh_type
  use mesh_collection_mod,      only: mesh_collection
  use extrusion_mod,            only: extrusion_type,         &
                                      uniform_extrusion_type, &
                                      TWOD
  use model_clock_mod,          only: model_clock_type
  use mpi_mod,                  only: mpi_type
  use namelist_collection_mod,  only: namelist_collection_type
  use namelist_mod,             only: namelist_type

  use jedi_lfric_increment_alg_mod, &
                                only: jedi_lfric_increment_alg
  !> @todo: Test code should not appear in the component
  !> @{
  use jedi_lfric_tests_config_mod, &
                                only: test_field
  !> @}
  use jedi_lfric_fake_nl_extrusion_mod, &
                                only: create_extrusion
  use jedi_lfric_fake_nl_init_files_mod, &
                                only: init_jedi_lfric_files

#ifdef USE_XIOS
  use driver_io_mod,            only: init_io, final_io, filelist_populator, &
                                      get_io_context
  use io_config_mod,            only: use_xios_io
  use io_context_mod,           only: io_context_type
  use lfric_xios_context_mod,   only: lfric_xios_context_type, advance
  use lfric_xios_write_mod,     only: write_state
#endif

  !------------------------------------
  ! Configuration modules
  !------------------------------------
  use base_mesh_config_mod, only: GEOMETRY_SPHERICAL, &
                                  GEOMETRY_PLANAR

  implicit none

  private
  ! To run local LFRic mini-app
  public initialise, step, finalise

  ! To be moved at a later date
  type( model_clock_type ), public, allocatable :: model_clock

  ! Coordinate field
  type(field_type), pointer, public :: chi(:)    => null()
  type(field_type), pointer, public :: panel_id  => null()
  type(mesh_type),  pointer, public :: mesh      => null()
  type(mesh_type),  pointer, public :: twod_mesh => null()
  type(inventory_by_mesh_type)      :: chi_inventory
  type(inventory_by_mesh_type)      :: panel_id_inventory

contains


  !> @brief Initialise the model mesh, fem, clock, and IO
  !>
  !> @param [in]    configuration Application configuration object
  !> @param [inout] program_name  The program name
  !> @param [inout] mpi           The mpi communicator
  subroutine initialise( configuration, program_name, mpi )

    implicit none

    type(namelist_collection_type) :: configuration

    character(len=*), intent(inout) :: program_name
    class(mpi_type),  intent(inout) :: mpi

    class(extrusion_type),        allocatable :: extrusion
    type(uniform_extrusion_type), allocatable :: extrusion_2d

    character(str_def)              :: base_mesh_names(1)
    character(str_def), allocatable :: twod_names(:)

    character(str_def) :: prime_mesh_name

    logical(l_def) :: apply_partition_check

    integer(i_def) :: geometry
    integer(i_def) :: stencil_depth
    real(r_def)    :: domain_bottom
    real(r_def)    :: domain_top
    real(r_def)    :: scaled_radius

    type(namelist_type), pointer :: base_mesh_nml => null()
    type(namelist_type), pointer :: extrusion_nml => null()
    type(namelist_type), pointer :: planet_nml    => null()

    integer(i_def) :: i
    integer(i_def), parameter :: one_layer = 1_i_def

    !--------------------------------------
    ! 0.0 Extract namelist variables
    !--------------------------------------
    base_mesh_nml => configuration%get_namelist('base_mesh')
    planet_nml    => configuration%get_namelist('planet')
    extrusion_nml => configuration%get_namelist('extrusion')

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call extrusion_nml%get_value( 'domain_top', domain_top )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    base_mesh_nml => null()
    planet_nml    => null()
    extrusion_nml => null()

    !--------------------------------------
    ! 1.0 Create the meshes
    !--------------------------------------
    base_mesh_names(1) = prime_mesh_name

    !--------------------------------------
    ! 1.1 Create the required extrusions
    !--------------------------------------
    select case (geometry)
    case (geometry_planar)
      domain_bottom = 0.0_r_def
    case (geometry_spherical)
      domain_bottom = scaled_radius
    case default
      call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
    end select

    allocate( extrusion, source=create_extrusion() )
    extrusion_2d = uniform_extrusion_type( domain_top,    &
                                           domain_bottom, &
                                           one_layer, TWOD )


    !-------------------------------------------------------------------------
    ! 1.2 Create the required meshes
    !-------------------------------------------------------------------------
    stencil_depth = 1
    apply_partition_check = .false.
    call init_mesh( configuration,       &
                    mpi%get_comm_rank(), &
                    mpi%get_comm_size(), &
                    base_mesh_names,     &
                    extrusion,           &
                    stencil_depth,       &
                    apply_partition_check )

    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps(twod_names)


    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)

    call init_time( model_clock )

#ifdef USE_XIOS
    if ( use_xios_io ) then
      call initialise_io( program_name, mpi, model_clock )
    end if
#endif

  end subroutine initialise


#ifdef USE_XIOS
  !> @brief Initialise the model IO
  !>
  !> @param [inout] program_name The program name
  !> @param [inout] mpi          The mpi communicator
  !> @param [inout] model_clock  The model clock
  subroutine initialise_io( program_name, mpi, model_clock )

    implicit none

    character(len=*),        intent(inout) :: program_name
    class(mpi_type),         intent(inout) :: mpi
    type(model_clock_type),  intent(inout) :: model_clock

    procedure(filelist_populator), pointer :: fl_populator => null()
    class(io_context_type),        pointer :: model_io_context => null()

    ! Initialise I/O context
    fl_populator => init_jedi_lfric_files
    call init_io( program_name, mpi%get_comm(), chi_inventory, panel_id_inventory, &
    model_clock, get_calendar(), populate_filelist=fl_populator )

    ! Do initial step
    model_io_context => get_io_context()
    if (model_clock%is_initialisation()) then
      select type (model_io_context)
      type is (lfric_xios_context_type)
        call advance(model_io_context, model_clock)
      end select
    end if

  end subroutine initialise_io
#endif


  !> @brief Performs a single timestep of the fake nl model
  !>
  !> @param [inout] model_data  The model data instance
  subroutine step( model_data )

    implicit none

    type(model_data_type), intent(inout) :: model_data

    type(field_collection_type), pointer :: depository => null()
    type(field_type),            pointer :: working_field => null()
#ifdef USE_XIOS
    class(io_context_type),      pointer :: model_io_context => null()
#endif

    depository => model_data%get_field_collection("depository")
    call depository%get_field( test_field, working_field )

#ifdef USE_XIOS
    if ( use_xios_io ) then
      ! Switch to the main model I/O context
      model_io_context => get_io_context()
      call model_io_context%set_current()
    end if
#endif

    call jedi_lfric_increment_alg( working_field )

  end subroutine step


  !> @brief Initialise the model IO
  !>
  !> @param [in]    program_name     The program name
  subroutine finalise( program_name )

    implicit none

    character(len=*),            intent(in)    :: program_name

#ifdef USE_XIOS
    if ( use_xios_io ) then
      call final_io()
    end if
#endif

    call final_fem()

  end subroutine finalise


end module jedi_lfric_fake_nl_driver_mod
