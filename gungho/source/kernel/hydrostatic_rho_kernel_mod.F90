!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes LHS of Galerkin projection and solves equation in W3 space

module hydrostatic_rho_kernel_mod

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_READ, GH_WRITE,    &
                                       WTHETA, W3,                     &
                                       GH_BASIS, GH_DIFF_BASIS,        &
                                       CELLS
use constants_mod,              only : r_def
use planet_config_mod,          only : gravity, cp, rd, p_zero
use idealised_config_mod,       only : test
use kernel_mod,                 only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: hydrostatic_rho_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,  GH_WRITE,  W3),                             &
       arg_type(GH_FIELD,  GH_READ,   WTHETA),                         &
       arg_type(GH_FIELD,  GH_READ,   WTHETA),                         &
       arg_type(GH_FIELD,  GH_READ,   W3)                              &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::hydrostatic_rho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface hydrostatic_rho_kernel_type
   module procedure hydrostatic_rho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public hydrostatic_rho_code
contains

type(hydrostatic_rho_kernel_type) &
   function hydrostatic_rho_kernel_constructor() result(self)
  return
end function hydrostatic_rho_kernel_constructor

!> @brief Computes hydrostatic pressure
subroutine hydrostatic_rho_code(nlayers, rho, theta, height_wth, height_w3, &
                                ndf_w3, undf_w3, map_w3,                    &
                                ndf_wth, undf_wth, map_wth)

  !Arguments
  integer, intent(in) :: nlayers, ndf_w3, ndf_wth, undf_w3, undf_wth
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_wth), intent(in) :: map_wth
  real(kind=r_def), dimension(undf_w3), intent(inout) :: rho
  real(kind=r_def), dimension(undf_wth), intent(in) :: theta, height_wth
  real(kind=r_def), dimension(undf_w3), intent(in) :: height_w3

  !Internal variables
  integer                     :: k
  real(kind=r_def)            :: dz, theta_cell
  real(kind=r_def)            :: exner(nlayers)

  real(kind=r_def), parameter :: p_surf=100000.0

! exner at the first level
  dz = height_w3(map_w3(1))-height_wth(map_wth(1))
  exner(1) = (p_surf/p_zero)**(rd/cp) - gravity * dz /(cp*theta(map_wth(1)))

! exner on other levels
  do k = 1, nlayers-1
    dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
    exner(k+1) = exner(k) - gravity * dz/(cp*theta(map_wth(1)+k))
  end do

! density from eqn of state
  do k = 1, nlayers
    theta_cell = 0.5*(theta(map_wth(1)+k)+theta(map_wth(1)+k-1))
    rho(map_w3(1)+k-1) = (p_zero*exner(k)**(1.0/(rd/cp)-1.0))/(rd*theta_cell)
  end do

end subroutine hydrostatic_rho_code

end module hydrostatic_rho_kernel_mod
