!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the index for switching rho_x/rho_y at cubed sphere edges
!> @details Calculates the start and end indices for when rho_x/rho_y should be
!!          switched at cubed sphere edges.
!!
!> @note This kernel only works with the lowest order elements.
!!

module cosmic_index_kernel_mod

  use argument_mod,       only : arg_type,                 &
                                 GH_FIELD, GH_REAL,        &
                                 CELL_COLUMN, GH_WRITE,    &
                                 GH_READ, GH_INTEGER,      &
                                 STENCIL, X1D, Y1D,        &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_tran, i_def, r_def
  use fs_continuity_mod,  only : W3, W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: cosmic_index_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                                      &
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),               & ! ix_start
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),               & ! ix_end
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),               & ! iy_start
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),               & ! iy_end
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(X1D)), & ! panel_id_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(Y1D))  & ! panel_id_y
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: cosmic_index_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: cosmic_index_code

contains

  !> @brief Compute the advective increment in x using PPM for the advective fluxes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] ix_start          The starting index for change in rho_x
  !> @param[in,out] ix_end            The end index for change in rho_x
  !> @param[in,out] iy_start          The starting index for change in rho_y
  !> @param[in,out] iy_end            The end index for change in rho_y
  !> @param[in]     panel_id_x        Panel ID of cell with stencil in x direction
  !> @param[in]     stencil_size_x    Local length of panel ID x stencil
  !> @param[in]     stencil_map_x     Dofmap for the panel ID x stencil
  !> @param[in]     panel_id_y        Panel ID of cell with stencil in y direction
  !> @param[in]     stencil_size_y    Local length of panel ID y stencil
  !> @param[in]     stencil_map_y     Dofmap for the panel ID y stencil
  !> @param[in]     ndf_wi            Number of degrees of freedom for index variables per cell
  !> @param[in]     undf_wi           Number of unique degrees of freedom for index variables
  !> @param[in]     map_wi            Map for index variables

  subroutine cosmic_index_code( nlayers,        &
                                ix_start,       &
                                ix_end,         &
                                iy_start,       &
                                iy_end,         &
                                panel_id_x,     &
                                stencil_size_x, &
                                stencil_map_x,  &
                                panel_id_y,     &
                                stencil_size_y, &
                                stencil_map_y,  &
                                ndf_wi,         &
                                undf_wi,        &
                                map_wi )

    use ffsl_cubed_sphere_edge_mod, only: get_index_rho_x, get_index_rho_y

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wi
    integer(kind=i_def), intent(in) :: ndf_wi
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_y

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_wi), intent(in) :: map_wi
    integer(kind=i_def), dimension(ndf_wi,stencil_size_x), intent(in) :: stencil_map_x
    integer(kind=i_def), dimension(ndf_wi,stencil_size_y), intent(in) :: stencil_map_y

    ! Arguments: Fields
    integer(kind=i_def), dimension(undf_wi), intent(inout) :: ix_start
    integer(kind=i_def), dimension(undf_wi), intent(inout) :: ix_end
    integer(kind=i_def), dimension(undf_wi), intent(inout) :: iy_start
    integer(kind=i_def), dimension(undf_wi), intent(inout) :: iy_end
    real(kind=r_def),    dimension(undf_wi), intent(in)    :: panel_id_x
    real(kind=r_def),    dimension(undf_wi), intent(in)    :: panel_id_y

    ! Local fields
    integer(kind=i_def) :: ipanel_x_local(1:stencil_size_x)
    integer(kind=i_def) :: ipanel_y_local(1:stencil_size_y)

    ! Stencils
    integer(kind=i_def) :: stencil_half_x
    integer(kind=i_def) :: stencil_half_y

    ! Indices
    integer(kind=i_def) :: jj

    ! Local stencils match those used in FFSL flux kernels
    ! e.g. | 1 | 2 | 3 | 4 | 5 | for extent 2

    ! Get the half point of the stencil
    stencil_half_x = (stencil_size_x + 1_i_def) / 2_i_def
    stencil_half_y = (stencil_size_y + 1_i_def) / 2_i_def

    ! Set up local panel ID
    do jj = 1, stencil_half_x
      ipanel_x_local(jj) = int(panel_id_x(stencil_map_x(1,stencil_half_x+1-jj)), i_def)
    end do
    do jj = stencil_half_x+1, stencil_size_x
      ipanel_x_local(jj) = int(panel_id_x(stencil_map_x(1,jj)), i_def)
    end do
    do jj = 1, stencil_half_y
      ipanel_y_local(jj) = int(panel_id_y(stencil_map_y(1,stencil_half_y+1-jj)), i_def)
    end do
    do jj = stencil_half_y+1, stencil_size_y
      ipanel_y_local(jj) = int(panel_id_y(stencil_map_y(1,jj)), i_def)
    end do

    ! Get the index for the direction change
    call get_index_rho_x(ix_start(map_wi(1)), ix_end(map_wi(1)), ipanel_y_local, stencil_size_y, stencil_half_y)
    call get_index_rho_y(iy_start(map_wi(1)), iy_end(map_wi(1)), ipanel_x_local, stencil_size_x, stencil_half_x)

  end subroutine cosmic_index_code

end module cosmic_index_kernel_mod
