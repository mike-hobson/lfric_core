!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Creates wet-rho from dry-rho and the moisture fields
!>
module wetrho_kernel_mod

  use argument_mod,         only : arg_type,                         &
                                  GH_FIELD, GH_READ, GH_READWRITE,   &
                                  CELLS
  use fs_continuity_mod,    only : WTheta
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: wetrho_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/             &
        arg_type(GH_FIELD,   GH_READWRITE, WTheta), &
        arg_type(GH_FIELD*6, GH_READ,      WTheta)  &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::wetrho_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface wetrho_kernel_type
    module procedure wetrho_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public wetrho_code

contains

type(wetrho_kernel_type) function wetrho_kernel_constructor() result(self)
  implicit none
  return
end function wetrho_kernel_constructor

!> @brief Adds the effect of water vapour to the exner pressure
!! @param[in] nlayers  Number of layers
!! @param[in,out] rho  Density in wth space
!! @param[in] mr_v     Water vapour mixing ratio
!! @param[in] mr_cl    Liquid cloud mixing ratio
!! @param[in] mr_r     Rain mixing ratio
!! @param[in] mr_ci    Ice cloud mixing ratio
!! @param[in] mr_s     Snow mixing ratio
!! @param[in] mr_g     Graupel mixing ratio
!! @param[in] ndf_wth  Number of degrees of freedom per cell
!! @param[in] undf_wth Total number of degrees of freedom
!! @param[in] map_wth  Dofmap for the cell at the base of the column
subroutine wetrho_code(nlayers,                         &
                       rho,                             &
                       mr_v,mr_cl,mr_r,mr_ci,mr_s,mr_g, &
                       ndf_wth, undf_wth, map_wth)

  use constants_mod, only: r_def, i_def

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers 
  integer(kind=i_def),                     intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth

  real(kind=r_def), intent(inout), dimension(undf_wth) :: rho
  real(kind=r_def), intent(in),    dimension(undf_wth) :: mr_v, mr_cl, mr_r, &
                                                          mr_ci, mr_s, mr_g

  integer(kind=i_def) :: k

  do k = 0, nlayers

    rho(map_wth(1) + k) = rho(map_wth(1) + k) * ( 1.0_r_def +            &
                          mr_v(map_wth(1) + k) + mr_cl(map_wth(1) + k) + &
                          mr_r(map_wth(1) + k) + mr_ci(map_wth(1) + k) + &
                          mr_s(map_wth(1) + k) + mr_g(map_wth(1) + k) )

  end do

end subroutine wetrho_code

end module wetrho_kernel_mod
