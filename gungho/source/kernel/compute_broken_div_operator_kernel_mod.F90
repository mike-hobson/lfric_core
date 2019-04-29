!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the broken (cell-local) divergence operator 
module compute_broken_div_operator_kernel_mod

  use argument_mod,              only: arg_type, func_type,     &
                                       GH_OPERATOR, GH_FIELD,   &
                                       GH_READ, GH_WRITE,       &
                                       ANY_SPACE_1,             &
                                       GH_BASIS, GH_DIFF_BASIS, &
                                       CELLS, GH_QUADRATURE_XYoZ
  use constants_mod,             only: r_def, i_def
  use coordinate_jacobian_mod,   only: coordinate_jacobian
  use fs_continuity_mod,         only: W2broken, W3
  use finite_element_config_mod, only: rehabilitate
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_broken_div_operator_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                &
        arg_type(GH_OPERATOR, GH_WRITE, W3, W2broken), &
        arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_1)   &
        /)
    type(func_type) :: meta_funcs(3) = (/     &
        func_type(W3, GH_BASIS),              &
        func_type(W2broken, GH_DIFF_BASIS),   &
        func_type(ANY_SPACE_1, GH_DIFF_BASIS) &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_broken_div_operator_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface compute_broken_div_operator_kernel_type
    module procedure compute_broken_div_operator_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public compute_broken_div_operator_code

contains

  type(compute_broken_div_operator_kernel_type) function compute_broken_div_operator_constructor() &
     result(self)
    implicit none
    return
  end function compute_broken_div_operator_constructor

  !> @brief Computes the broken (cell-local) divergence operator 
  !! @param[in] cell     Cell number
  !! @param[in] nlayers  Number of layers.
  !! @param[in] ncell_3d ncell*ndf
  !! @param[out] broken_div Local stencil of the broken div operator.
  !! @param[in] chi1     Data array for chi in the 1st dir
  !! @param[in] chi2     Data array for chi in the 2nd dir
  !! @param[in] chi3     Data array for chi in the 3rd dir
  !! @param[in] ndf_w3   Number of degrees of freedom per cell.
  !! @param[in] basis_w3 Scalar basis functions
  !!                     evaluated at quadrature points.
  !! @param[in] ndf_w2b  Number of degrees of freedom per cell.
  !! @param[in] diff_basis_w2b Differential vector basis
  !!                     functions evaluated at quadrature points.
  !! @param[in] ndf_chi  Number of degrees of freedom per cell for chi
  !!                     field.
  !! @param[in] undf_chi Number of unique degrees of freedom  for chi
  !!                     field.
  !! @param[in] map_chi  Dofmap for the cell at the
  !!                     base of the column, for the space on which the chi field
  !!                     lives.
  !! @param[in] diff_basis_chi Vector differential
  !!                    basis functions evaluated at quadrature points.
  !! @param[in] nqp_h    Number of horizontal quadrature points
  !! @param[in] nqp_v    Number of vertical quadrature points
  !! @param[in] wqp_h    Horizontal quadrature weights
  !! @param[in] wqp_v    Vertical quadrature weights
  subroutine compute_broken_div_operator_code(cell, nlayers, ncell_3d,       &
                                              broken_div,                    &
                                              chi1, chi2, chi3,              &
                                              ndf_w3, basis_w3,              &
                                              ndf_w2b, diff_basis_w2b,       &
                                              ndf_chi, undf_chi,             &
                                              map_chi, diff_basis_chi,       &
                                              nqp_h, nqp_v, wqp_h, wqp_v)

    implicit none 

    ! Argument declarations
    integer(kind=i_def),                     intent(in) :: cell, nqp_h, nqp_v
    integer(kind=i_def),                     intent(in) :: nlayers
    integer(kind=i_def),                     intent(in) :: ncell_3d
    integer(kind=i_def),                     intent(in) :: ndf_w3, ndf_w2b
    integer(kind=i_def),                     intent(in) :: ndf_chi, undf_chi
    integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi

    real(kind=r_def), intent(in) :: diff_basis_chi(3, ndf_chi, nqp_h, nqp_v)
    real(kind=r_def), intent(in) :: basis_w3(1, ndf_w3, nqp_h, nqp_v)
    real(kind=r_def), intent(in) :: diff_basis_w2b(1, ndf_w2b, nqp_h, nqp_v)

    real(kind=r_def), dimension(ndf_w3, ndf_w2b, ncell_3d), intent(out)   :: broken_div
    real(kind=r_def), dimension(undf_chi),                  intent(in)    :: chi1
    real(kind=r_def), dimension(undf_chi),                  intent(in)    :: chi2
    real(kind=r_def), dimension(undf_chi),                  intent(in)    :: chi3
    real(kind=r_def), dimension(nqp_h),                     intent(in)    :: wqp_h
    real(kind=r_def), dimension(nqp_v),                     intent(in)    :: wqp_v

    ! Internal variables
    integer(kind=i_def)                             :: df, df2, df3, k, ik
    integer(kind=i_def)                             :: qp1, qp2
    real(kind=r_def), dimension(ndf_chi)            :: chi1_e, chi2_e, chi3_e
    real(kind=r_def)                                :: integrand
    real(kind=r_def), dimension(nqp_h, nqp_v)       :: dj
    real(kind=r_def), dimension(3, 3, nqp_h, nqp_v) :: jac

    ! Loop over layers
    do k = 0, nlayers - 1
      ik = k + 1 + (cell - 1) * nlayers

      ! Indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k)
        chi2_e(df) = chi2(map_chi(df) + k)
        chi3_e(df) = chi3(map_chi(df) + k)
      end do

      ! Compute Jacobian
      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                               diff_basis_chi, jac, dj)

      ! Run over dof extent of W2Broken
      do df2 = 1, ndf_w2b
        ! Run over dof extent of W3
        do df3 = 1, ndf_w3
          ! Initialize
          broken_div(df3, df2, ik) = 0.0_r_def

          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              ! Weak divergence operator: div(W2_basis) * W3_basis * dx
              if ( rehabilitate ) then
                ! With rehabilitation 
                !   divergence mapping is div(x) -> ! \hat{div}(\hat{x})
                integrand = wqp_h(qp1) * wqp_v(qp2)          &
                          * basis_w3(1, df3, qp1, qp2)       &
                          * diff_basis_w2b(1, df2, qp1, qp2) 
              else
                ! Without rehabilitation
                !   divergence mapping is div(x) -> ! \hat{div}(\hat{x})/det(J)
                integrand = wqp_h(qp1) * wqp_v(qp2)          &
                          * basis_w3(1, df3, qp1, qp2)       &
                          * diff_basis_w2b(1, df2, qp1, qp2) &
                          / dj(qp1, qp2) 
              end if
              broken_div(df3, df2, ik) = broken_div(df3, df2, ik) + integrand
            end do
          end do
        end do ! End of W3 dof loop
      end do ! End of W2Broken dof loop

    end do ! End of layer loop

  end subroutine compute_broken_div_operator_code

end module compute_broken_div_operator_kernel_mod
