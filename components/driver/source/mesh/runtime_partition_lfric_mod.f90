!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Subroutines to for interfacing runtime partitioning using LFRic
!>        infrastructure objects.

module runtime_partition_lfric_mod

  use constants_mod,           only: i_def, l_def
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type
  use partition_mod,           only: partitioner_interface
  use runtime_partition_mod,   only: get_partition_parameters_args_only &
                                         => get_partition_parameters

  use partitioning_config_mod, only: panel_decomposition_auto,   &
                                     panel_decomposition_row,    &
                                     panel_decomposition_column, &
                                     panel_decomposition_custom

  implicit none

  private
  public :: get_partition_parameters

  interface get_partition_parameters
    procedure get_partition_parameters_cfg
    procedure get_partition_parameters_nml
  end interface get_partition_parameters

contains

subroutine get_partition_parameters_cfg( configuration,  &
                                         mesh_selection, &
                                         total_ranks,    &
                                         panel_xproc,    &
                                         panel_yproc,    &
                                         partitioner_ptr )

  implicit none

  type(namelist_collection_type), intent(in) :: configuration

  integer,        intent(in) :: mesh_selection
  integer(i_def), intent(in) :: total_ranks

  integer(i_def), intent(inout) :: panel_xproc
  integer(i_def), intent(inout) :: panel_yproc

  type(namelist_type), pointer :: partitioning

  procedure(partitioner_interface), intent(out), pointer :: partitioner_ptr

  partitioning => configuration%get_namelist('partitioning')

  call get_partition_parameters_nml( partitioning,     &
                                     mesh_selection,   &
                                     total_ranks,      &
                                     panel_xproc,      &
                                     panel_yproc,      &
                                     partitioner_ptr )


end subroutine get_partition_parameters_cfg

subroutine get_partition_parameters_nml( partitioning,   &
                                         mesh_selection, &
                                         total_ranks,    &
                                         panel_xproc,    &
                                         panel_yproc,    &
                                         partitioner_ptr )

  use runtime_partition_mod, only: decomp_auto, &
                                   decomp_row,  &
                                   decomp_column
  implicit none


  type(namelist_type), intent(in), pointer :: partitioning

  integer,        intent(in) :: mesh_selection
  integer(i_def), intent(in) :: total_ranks

  integer(i_def), intent(inout) :: panel_xproc
  integer(i_def), intent(inout) :: panel_yproc

  procedure(partitioner_interface), intent(out), pointer :: partitioner_ptr

  logical :: custom
  integer :: decomp
  integer :: panel_decomposition

  call partitioning%get_value( 'panel_decomposition', panel_decomposition )

  custom = .false.
  decomp = decomp_auto

  select case (panel_decomposition)

  case ( panel_decomposition_auto )
    decomp = decomp_auto

  case ( panel_decomposition_row )
    decomp = decomp_row

  case ( panel_decomposition_column )
    decomp = decomp_column

  case ( panel_decomposition_custom )
    custom = .true.
    call partitioning%get_value( 'panel_xproc', panel_xproc )
    call partitioning%get_value( 'panel_yproc', panel_yproc )

  case default

  end select


  if (custom) then

    call get_partition_parameters_args_only( mesh_selection, &
                                             total_ranks,    &
                                             panel_xproc,    &
                                             panel_yproc,    &
                                             partitioner_ptr )
  else
    call get_partition_parameters_args_only( mesh_selection, &
                                             decomp,         &
                                             total_ranks,    &
                                             panel_xproc,    &
                                             panel_yproc,    &
                                             partitioner_ptr )
  end if

end subroutine get_partition_parameters_nml


end module runtime_partition_lfric_mod
