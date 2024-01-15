!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------

!> @brief Module container for query functions related to science.
module query_mod

  use constants_mod,   only: i_def
  use global_mesh_mod, only: global_mesh_type
  use log_mod,         only: log_event, log_scratch_space, &
                             log_level_error

  implicit none

  private
  public  :: valid_for_global_model, &
             check_uniform_partitions


contains

!> @brief  Queries whether a global mesh is valid for a global model domain.
!>         A negative result does not imply that the mesh is valid for a
!>         regional model.
!> @param[in]  global_mesh
!> @return     answwer
function valid_for_global_model( global_mesh ) result ( answer )

   implicit none

   type(global_mesh_type), intent(in) :: global_mesh
   logical :: answer

   if ( global_mesh%is_topology_periodic() .and. &
        global_mesh%is_coord_sys_ll()      .and. &
        global_mesh%is_geometry_spherical() ) then

      ! Note: if these conditions where satisfied from the
      !       planar mesh generator then the model would be
      !       a torus, though as this is not supported, this
      !       is only .true. for a sphere.
      answer = .true.
  else
      answer = .false.
  end if

end function valid_for_global_model

!> @brief  Checks that a global mesh/partitioning strategy combination will
!>         produce equally sized partitions. This is an important check
!>         for related meshes which require the partition stratedy to
!>         to produce partitions over the same geographical location.
!>
!> @param[in]  global_mesh  The mesh to be checked
!> @param[in]  xproc        The number of partitions across a mesh panel
!>                          in the panel x-direction.
!> @param[in]  yproc        The number of partitions across a mesh panel
!>                          in the panel y-direction.

!==========================================================================
function check_uniform_partitions( global_mesh, xproc, yproc ) result (answer)

  use reference_element_mod, only : W, E

  implicit none

  type(global_mesh_type), intent(in), pointer :: global_mesh
  integer(i_def),         intent(in)          :: xproc, yproc

  integer(i_def) :: void_cell    ! Cell id that marks the cell as a cell outside of the partition.
  integer(i_def) :: w_cell       ! The id of a cell on the western edge of the domain
  integer(i_def) :: cell_next(4) ! The cells around the cell being queried
  integer(i_def) :: cell_next_e  ! The cell to the east of the cell being queried
  integer(i_def) :: panel_edge_ncells_x  ! number of cells across a panel of the input mesh in x-direction
  integer(i_def) :: panel_edge_ncells_y  ! number of cells across a panel of the input mesh in y-direction
  logical :: periodic_xy(2) ! Is mesh periodic in the x/y-axes

  logical :: answer

  if ( valid_for_global_model(global_mesh) ) then
    !================================================
    ! Treat as global model with panelled sphere mesh
    !================================================

    ! Generate the "C" number (as in Cnn) from global mesh data
    ! Assume the panels on the globe are square
    panel_edge_ncells_x = int(sqrt(real( global_mesh%get_ncells() / &
                                         global_mesh%get_npanels() )))
    panel_edge_ncells_y = panel_edge_ncells_x

  else
    !================================================
    ! Treat as regional model a single panel mesh.
    !================================================

    ! Planar or LAM mesh mesh might be non-square - so find the dimensions
    ! First find a cell on the west edge of the domain
    ! If periodic, cell id 1 can be used as mesh conectivity loops round
    w_cell = 1

    void_cell   = global_mesh%get_void_cell()
    periodic_xy = global_mesh%get_mesh_periodicity()
    if ( .not. periodic_xy(1) ) then
      ! If not periodic in E-W direction then walk West until you reach mesh
      ! edge defined by the void cell.
      call global_mesh%get_cell_next(w_cell,cell_next)
      do while (cell_next(W) /= void_cell)
        w_cell = cell_next(W)
        call global_mesh%get_cell_next(w_cell,cell_next)
      end do
    end if

    ! Work out number of cells in x and y directions
    panel_edge_ncells_x = 1
    ! Starting at the West edge of the mesh, walk East until you reach either
    ! the cell you started at (periodic) or a void cell (LAM)
    ! - this determines the number of cells in the x direction
    call global_mesh%get_cell_next(w_cell,cell_next)
    cell_next_e = cell_next(E)
    do while (cell_next_e /= w_cell .and. cell_next_e /= void_cell)
      panel_edge_ncells_x = panel_edge_ncells_x + 1
      call global_mesh%get_cell_next(cell_next_e, cell_next)
      cell_next_e = cell_next(E)
    end do
    ! Infer num_cells_y from the total domin size and num_cells_x
    panel_edge_ncells_y = global_mesh%get_ncells() / panel_edge_ncells_x

  end if

  if ( mod(panel_edge_ncells_x,xproc) /= 0 .or. &
       mod(panel_edge_ncells_y,yproc) /= 0 ) then
    answer = .false.
  else
    answer = .true.
  end if

end function check_uniform_partitions

end module query_mod
