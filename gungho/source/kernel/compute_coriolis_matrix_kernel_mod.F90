!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the Coriolis operator to apply the rotation vector Omega to
!>        the wind fields.
!> @details The form of Coriolis operator is:
!> \f[ <v,2\Omega \times v> \f] where v is the test/trial function for
!> the velocity space and Omega is the rotation vector of the domain.
!> The Coriolis terms will then be applied by multiplying this operator
!> by a wind field.
!>
!
module compute_coriolis_matrix_kernel_mod

use constants_mod,           only: i_def, r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,             &
                                   GH_OPERATOR, GH_FIELD,           &
                                   GH_READ, GH_WRITE,               &
                                   ANY_SPACE_9,                     &
                                   GH_BASIS, GH_DIFF_BASIS,         &
                                   CELLS, GH_QUADRATURE_XYoZ
use fs_continuity_mod,       only: W2

use coordinate_jacobian_mod, only: coordinate_jacobian
use base_mesh_config_mod,    only: geometry,           &
                                   geometry_spherical, &
                                   f_lat
use rotation_vector_mod,     only: rotation_vector_fplane,  &
                                   rotation_vector_sphere
use planet_config_mod,       only: scaled_omega
use cross_product_mod,       only: cross_product
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_coriolis_matrix_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W2, W2),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)                    &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_coriolis_matrix_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface compute_coriolis_matrix_kernel_type
   module procedure compute_coriolis_matrix_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_coriolis_matrix_code
contains

type(compute_coriolis_matrix_kernel_type) &
                        function compute_coriolis_matrix_constructor() result(self)
  implicit none
  return
end function compute_coriolis_matrix_constructor

!> @brief Compute the Coriolis operator to apply the rotation vector Omega to
!>        the wind fields.
!!
!! @param[in] cell     Identifying number of cell.
!! @param[in] nlayers  Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[inout] mm    Local stencil or Coriolis operator.
!! @param[in] chi1  Physical coordinate in the 1st dir.
!! @param[in] chi2  Physical coordinate in the 2nd dir.
!! @param[in] chi3  Physical coordinate in the 3rd dir.
!! @param[in] ndf   Degrees of freedom per cell.
!! @param[in] basis Vector basis functions evaluated at quadrature points.
!! @param[in] ndf_chi  Degrees of freedum per cell for chi field.
!! @param[in] undf_chi Unique degrees of freedum  for chi field.
!! @param[in] map_chi  Dofmap for the cell at the base of the column, for the
!!                     space on which the chi field lives
!! @param[in] basis_chi Basis functions evaluated at
!!                      quadrature points.
!! @param[in] diff_basis_chi Differential basis functions evaluated at
!!                           quadrature points.
!! @param[in] nqp_h    Number of horizontal quadrature points.
!! @param[in] nqp_v    Number of vertical quadrature points.
!! @param[in] wqp_h    Horizontal quadrature weights.
!! @param[in] wqp_v    Vertical quadrature weights.
subroutine compute_coriolis_matrix_code(cell, nlayers, ncell_3d,     &
                                        mm,                          &
                                        chi1, chi2, chi3,            &
                                        ndf, basis,                  &
                                        ndf_chi, undf_chi, map_chi,  &
                                        basis_chi, diff_basis_chi,   &
                                        nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: cell, nqp_h, nqp_v
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ncell_3d
  integer(kind=i_def), intent(in)    :: ndf

  real(kind=r_def), dimension(3,ndf,nqp_h,nqp_v), intent(in) :: basis

  integer(kind=i_def), intent(in)    :: ndf_chi
  integer(kind=i_def), intent(in)    :: undf_chi
  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  real(kind=r_def),    intent(inout) :: mm(ndf,ndf,ncell_3d)
  real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: chi1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi3(undf_chi)
  real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
  real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, ik
  integer(kind=i_def)                          :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(3,nqp_h,nqp_v)   :: rotation_vector
  real(kind=r_def), dimension(3)               :: omega_cross_u
  real(kind=r_def), dimension(3)               :: jac_u, jac_v


  ! Loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers

     ! Indirect the chi coord field here
     do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
     end do

    ! Calculate rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) and Jacobian
    if ( geometry == geometry_spherical ) then
      call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e,       &
                                  chi3_e, basis_chi, rotation_vector)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, scaled_omega, f_lat,           &
                                  rotation_vector)
    end if

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,    &
                             diff_basis_chi, jac, dj)

    ik = k + (cell-1)*nlayers
    mm(:,:,ik) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        do df2 = 1, ndf
          jac_u = matmul(jac(:,:,qp1,qp2),basis(:,df2,qp1,qp2))
          omega_cross_u = wqp_h(qp1)*wqp_v(qp2)                                &
                        *cross_product(rotation_vector(:,qp1,qp2),jac_u)       &
                        /dj(qp1,qp2)
          do df = 1, ndf
             jac_v = matmul(jac(:,:,qp1,qp2),basis(:,df,qp1,qp2))
             mm(df,df2,ik) = mm(df,df2,ik) - dot_product(jac_v,omega_cross_u)
          end do
        end do
      end do
    end do
  end do ! end of k loop

end subroutine compute_coriolis_matrix_code

end module compute_coriolis_matrix_kernel_mod
