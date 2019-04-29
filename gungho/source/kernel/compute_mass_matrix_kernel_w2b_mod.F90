!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Provides access to the members of the w2b_kernel class.
module compute_mass_matrix_kernel_w2b_mod

  use argument_mod,            only: arg_type, func_type,     &
                                     GH_OPERATOR, GH_FIELD,   &
                                     GH_READ, GH_WRITE,       &
                                     ANY_SPACE_9,             &
                                     GH_BASIS, GH_DIFF_BASIS, &
                                     CELLS, GH_QUADRATURE_XYoZ
  use constants_mod,           only: i_def, r_def
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use fs_continuity_mod,       only: W2broken
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w2b_type
    private
    type(arg_type) :: meta_args(2) = (/                       &
        arg_type(GH_OPERATOR, GH_WRITE, W2broken, W2broken),  &
        arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)          &
        /)
    type(func_type) :: meta_funcs(2) = (/         &
        func_type(ANY_SPACE_9, GH_DIFF_BASIS),    &
        func_type(W2broken, GH_BASIS)             &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_mass_matrix_w2b_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface compute_mass_matrix_kernel_w2b_type
    module procedure compute_mass_matrix_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public compute_mass_matrix_w2b_code

contains

  type(compute_mass_matrix_kernel_w2b_type) function compute_mass_matrix_constructor() &
     result(self)
    implicit none
    return
  end function compute_mass_matrix_constructor

  !> @brief Computes the mass matrix for the w2b space.
  !!
  !! @param[in] cell      Identifying number of cell.
  !! @param[in] nlayers   Number of layers.
  !! @param[in] ncell_3d  ncell*ndf
  !! @param[out] mm       Local stencil or mass matrix.
  !! @param[in] chi1      Physical coordinate in the 1st dir.
  !! @param[in] chi2      Physical coordinate in the 2nd dir.
  !! @param[in] chi3      Physical coordinate in the 3rd dir.
  !! @param[in] ndf_w2b   Degrees of freedom per cell.
  !! @param[in] basis_w2b Vector basis functions evaluated at quadrature points.
  !! @param[in] ndf_chi   Degrees of freedum per cell for chi field.
  !! @param[in] undf_chi  Unique degrees of freedum  for chi field.
  !! @param[in] map_chi   Dofmap for the cell at the base of the column, for the
  !!                      space on which the chi field lives.
  !! @param[in] diff_basis_chi Vector differential basis functions evaluated at
  !!                           quadrature points.
  !! @param[in] nqp_h     Number of horizontal quadrature points.
  !! @param[in] nqp_v     Number of vertical quadrature points.
  !! @param[in] wqp_h     Horizontal quadrature weights.
  !! @param[in] wqp_v     Vertical quadrature weights.
  subroutine compute_mass_matrix_w2b_code(cell, nlayers, ncell_3d,     &
                                          mm,                          &
                                          chi1, chi2, chi3,            &
                                          ndf_w2b, basis_w2b,          &
                                          ndf_chi, undf_chi, map_chi,  &
                                          diff_basis_chi,              &
                                          nqp_h, nqp_v, wqp_h, wqp_v)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: cell, nqp_h, nqp_v
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_3d
    integer(kind=i_def), intent(in)    :: ndf_w2b
    integer(kind=i_def), intent(in)    :: ndf_chi
    integer(kind=i_def), intent(in)    :: undf_chi
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)

    real(kind=r_def),    intent(out)   :: mm(ndf_w2b,ndf_w2b,ncell_3d)
    real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in)    :: chi1(undf_chi)
    real(kind=r_def),    intent(in)    :: chi2(undf_chi)
    real(kind=r_def),    intent(in)    :: chi3(undf_chi)
    real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
    real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)

    real(kind=r_def), dimension(3, ndf_w2b, nqp_h, nqp_v), intent(in) :: basis_w2b

    ! Internal variables
    integer(kind=i_def)                             :: df, df2, k, ik
    integer(kind=i_def)                             :: qp1, qp2

    real(kind=r_def), dimension(ndf_chi)            :: chi1_e, chi2_e, chi3_e
    real(kind=r_def)                                :: integrand
    real(kind=r_def), dimension(nqp_h, nqp_v)       :: dj
    real(kind=r_def), dimension(3, 3, nqp_h, nqp_v) :: jac

    ! Loop over layers: Start from 1 as in this loop k is not an offset
    do k = 1, nlayers
      ik = k + (cell - 1) * nlayers

      ! Indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
      end do

      ! Compute Jacobian
      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                               diff_basis_chi, jac, dj)

      ! Run over dof extent of W2Broken
      do df2 = 1, ndf_w2b
        do df = df2, ndf_w2b ! Exploiting symmetry of the operator
          ! Initialize
          mm(df, df2, ik) = 0.0_r_def
          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              ! Mass term: dot(w2b_basis_i, w2b_basis_j) * dx (mapped to reference space)
              integrand = wqp_h(qp1) * wqp_v(qp2)                                               &
                        * dot_product(matmul(jac(:, :, qp1, qp2), basis_w2b(:, df, qp1, qp2)),  &
                                      matmul(jac(:, :, qp1, qp2), basis_w2b(:, df2, qp1, qp2))) &
                        / dj(qp1, qp2)
              mm(df, df2, ik) = mm(df, df2, ik) + integrand
            end do
          end do
        end do
        do df = df2, 1, -1
          ! Mass matrix is symmetric
          mm(df, df2, ik) = mm(df2, df, ik)
        end do
      end do

    end do ! End of layer loop

  end subroutine compute_mass_matrix_w2b_code

end module compute_mass_matrix_kernel_w2b_mod
