!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Kernel which calculates the reciprocal of a field
!>        apart from when the field is zero in which case it is zero

module set_reciprocal_mod

use kernel_mod,    only : kernel_type
use argument_mod,  only : arg_type,            &
                          GH_REAL,GH_FIELD,    &
                          GH_WRITE, GH_READ,   &
                          ANY_DISCONTINUOUS_SPACE_1, &
                          CELL_COLUMN
use constants_mod, only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: set_reciprocal_type
  private
  type(arg_type) :: meta_args(2) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: set_reciprocal_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: set_reciprocal_code
contains

!> @brief Kernel which calculates the reciprocal of a field
!>        apart from when the field is zero in which case it is zero
!! @param[in]     nlayers          Number of layers
!! @param[in,out] field_reciprocal Field to write the reciprocal to
!! @param[in]     input_field      Input field to read from
!! @param[in]     ndf              Number of degrees of freedom per cell
!! @param[in]     undf             Total number of unique degrees of freedom
!! @param[in]     map              Dofmap for the cell at the base of the column

subroutine set_reciprocal_code(nlayers,                       &
                               field_reciprocal, input_field, &
                               ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def),                  intent(in)    :: nlayers
  integer(kind=i_def),                  intent(in)    :: ndf
  integer(kind=i_def),                  intent(in)    :: undf
  integer(kind=i_def), dimension(ndf),  intent(in)    :: map
  real(kind=r_def),    dimension(undf), intent(inout) :: field_reciprocal
  real(kind=r_def),    dimension(undf), intent(in)    :: input_field

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 0, nlayers - 1
    do df = 1, ndf
      if (abs(input_field(map(df) + k)) > tiny(1.0_r_def)) then
        field_reciprocal(map(df) + k) = 1.0_r_def/input_field(map(df) + k)
      else
        field_reciprocal(map(df) + k) = 0.0_r_def
      end if
    end do
  end do

end subroutine set_reciprocal_code

end module set_reciprocal_mod
