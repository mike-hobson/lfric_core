!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Arbitary kernel for testing fields with a non-spatial dimension
!>
!> @details Sets the field equal to the position along the non-spatial dimension
!>

module non_spatial_kernel_mod

    use argument_mod,      only : arg_type, GH_FIELD, GH_REAL, GH_SCALAR, &
                                  GH_INTEGER, GH_READ, GH_WRITE, CELL_COLUMN, &
                                  ANY_DISCONTINUOUS_SPACE_1
    use constants_mod,     only : r_def, i_def
    use fs_continuity_mod, only : W3
    use kernel_mod,        only : kernel_type
    use log_mod,           only : log_scratch_space, log_event, LOG_LEVEL_INFO

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------

    type, public, extends(kernel_type) :: non_spatial_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                                 &
             arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
             arg_type(GH_SCALAR, GH_INTEGER, GH_READ) &
             /)
        integer :: operates_on = CELL_COLUMN
    contains
        procedure, nopass :: non_spatial_code
    end type

    !---------------------------------------------------------------------------
    ! Contained functions/subroutines
    !---------------------------------------------------------------------------
    public :: non_spatial_code

contains

    !> @brief Set dof equal to position along non-spatial dimension
    !> @details finds the maximum value in the neighbouring fields (stencil + fields above / below)
    !>          and based on that will increase the value of the current field by blend_percentage * difference
    !> @param[in] number_of_layers: the number of layers
    !> @param[in] field: target field to be worked on
    !> @param[in] length_nsd: The length of the first non-spatial dimension
    !> @param[in] degrees_of_freedom: degrees of freededom
    !> @param[in] unique_degrees_of_freedom: total number of degrees of freedom across the column
    !> @param[in] field_dof_map: dof map for bottom cell of field
    subroutine non_spatial_code(n_layers, &
            field, &
            length_nsd, &
            degrees_of_freedom, &
            unique_degrees_of_freedom, &
            field_dof_map)

        implicit none

        !> Arguments

        integer(kind = i_def), intent(in) :: n_layers
        integer(kind = i_def), intent(in) :: degrees_of_freedom
        integer(kind = i_def), intent(in) :: unique_degrees_of_freedom
        integer(kind = i_def), dimension(degrees_of_freedom), intent(in) :: field_dof_map

        integer(kind = i_def), intent(in) :: length_nsd
        real(kind = r_def), dimension(unique_degrees_of_freedom), intent(inout) :: field

        !> Internal Vars
        integer(kind = i_def) :: non_spatial_pos, layer, dof

        !> processing

        do non_spatial_pos = 0, length_nsd - 1
            do layer = 0, n_layers - 1
                do dof = 1, degrees_of_freedom
                    field(field_dof_map(dof) + (non_spatial_pos)*(n_layers) + layer) = non_spatial_pos + 1
                end do
            end do
        end do

    end subroutine non_spatial_code

end module non_spatial_kernel_mod
