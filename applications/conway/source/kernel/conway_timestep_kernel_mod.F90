!-----------------------------------------------------------------------------
! (C) Crown copyright 2015 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @todo Edit the kernel with tutorial tasks #3982


!> @brief  Create a next timestep state from the state at the current timestep
!>         following the rules of Conway's Game of Life
!>

module conway_timestep_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_INTEGER,      &
                                    GH_READ, GH_WRITE,         &
                                    CELL_COLUMN, STENCIL, REGION
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W3
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.

  type, public, extends(kernel_type) :: conway_timestep_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                 &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  W3, STENCIL(REGION)), &
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, W3)                   &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: conway_timestep_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: conway_timestep_code

contains

!> @brief Calculates horizontal Smagorinsky diffusion for a tracer variable
!! @param[in] nlayers Number of layers in the mesh
!! @param[in] current_state The state of the game at the current timestep
!! @param[in] map_w3_stencil_size Number of cells in the stencil at the base
!!                                of the column for Wtheta
!! @param[in] map_w3_stencil Array holding the dofmap for the stencil at the
!!                           base of the column for Wtheta
!! @param[inout] next_state The generated state of the game at the next timestep
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3 space
!! @param[in] undf_w3  Number of unique degrees of freedom for w3 space
!! @param[in] map_w3 Cell dofmap for w3 space
subroutine conway_timestep_code( nlayers,             &
                                 current_state,       &
                                 map_w3_stencil_size, &
                                 map_w3_stencil,      &
                                 next_state,          &
                                 ndf_w3,              &
                                 undf_w3,             &
                                 map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), intent(in) :: map_w3_stencil_size
  integer(kind=i_def), intent(in) :: map_w3_stencil(ndf_w3,map_w3_stencil_size)

  integer(kind=i_def), intent(inout) :: next_state(undf_w3)
  integer(kind=i_def), intent(in)    :: current_state(undf_w3)

  ! Internal variables
  integer(kind=i_def) :: i, k
  integer(kind=i_def) :: neighbours

  ! The layout of the cells in the stencil is:
  !
  !     ------------|
  !     |   |   |   |
  !     | 9 | 8 | 7 |
  !     |------------
  !     |   |   |   |
  !     | 2 | 1 | 6 |
  !     -------------
  !     |   |   |   |
  !     | 3 | 4 | 5 |
  !     -------------

  ! If the full stencil isn't available, we must be at the domain edge/corner.
  ! Simply set the new state to 0 for now, and exit the routine.
  if (map_w3_stencil_size < 9_i_def) then
    do k = 0, nlayers
      next_state(map_w3(1) + k) = 0_i_def
    end do
    return
  end if

  ! Calculate next cell state from current
  do k = 0, nlayers

    ! Calucate the number of "live" neighbours
    neighbours = 0_i_def
    do i = 2, 9
      neighbours = neighbours + current_state(map_w3_stencil(1,i) + k)
    end do

    ! Apply the rules of Conway's Game of Life
    select case(neighbours)
      case (0:1, 4:8)
        next_state(map_w3(1)+k) = 0_i_def
      case (2)
        ! no change
      case (3)
        next_state(map_w3(1)+k) = 1_i_def
    end select

  end do

end subroutine conway_timestep_code

end module conway_timestep_kernel_mod
