!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes the initial theta field

!> @details The kernel computes initial theta perturbation field for theta in the space
!>          of horizontally discontinuous, vertically continuous polynomials

module initial_theta_kernel_mod

    use argument_mod,                  only: arg_type,  &
        GH_FIELD, GH_WRITE, GH_READ,                    &
        W0, ANY_SPACE_1, GH_BASIS,                      &
        GH_DIFF_BASIS,                                  &
        CELLS
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type
    use idealised_config_mod,          only: test

    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_theta_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                               &
            arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_1),                  &
            arg_type(GH_FIELD*3, GH_READ, W0)                             &
            /)
        integer :: iterates_over = CELLS

    contains
        procedure, nopass :: initial_theta_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! overload the default structure constructor for function space
    interface initial_theta_kernel_type
        module procedure initial_theta_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public initial_theta_code
contains

    type(initial_theta_kernel_type) function initial_theta_kernel_constructor() result(self)
        return
    end function initial_theta_kernel_constructor

    !> @brief Computes the initial theta field
    !! @param[in] nlayers Number of layers
    !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
    !! @param[in] undf_wtheta Number of total degrees of freedom for wtheta
    !! @param[in] map_wtheta Dofmap for the cell at the base of the column
    !! @param[inout] theta Potential temperature
    !! @param[in] ndf_w0 Number of degrees of freedom per cell
    !! @param[in] ndf_w0 Total number of degrees of freedom
    !! @param[in] undf_w0 Number of total degrees of freedom for w0
    !! @param[in] map_w0 Dofmap for the cell at the base of the column
    !! @param[in] w0_basis Basis functions evaluated at gaussian quadrature points
    !! @param[inout] chi_1 X component of the w0 coordinate field
    !! @param[inout] chi_2 Y component of the w0 coordinate field
    !! @param[inout] chi_3 Z component of the w0 coordinate field
    subroutine initial_theta_code(nlayers, ndf_wtheta, undf_wtheta, map_wtheta, theta, &
                                  ndf_w0, undf_w0, map_w0, w0_basis, chi_1, chi_2, chi_3)

        use analytic_temperature_profiles_mod, only : analytic_temperature

        implicit none

        !Arguments
        integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, ndf_w0, undf_wtheta, undf_w0
        integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
        integer(kind=i_def), dimension(ndf_w0), intent(in) :: map_w0
        real(kind=r_def), dimension(undf_wtheta),          intent(inout) :: theta
        real(kind=r_def), dimension(undf_w0),              intent(in)    :: chi_1, chi_2, chi_3
        real(kind=r_def), dimension(1,ndf_w0,ndf_wtheta),  intent(in)    :: w0_basis

        !Internal variables
        integer(kind=i_def)                 :: df, df0, k
        real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
        real(kind=r_def)                    :: x(3)

        ! compute the pointwise theta profile

        do k = 0, nlayers-1
          do df0 = 1, ndf_w0
            chi_1_e(df0) = chi_1( map_w0(df0) + k)
            chi_2_e(df0) = chi_2( map_w0(df0) + k)
            chi_3_e(df0) = chi_3( map_w0(df0) + k)
          end do

          do df = 1, ndf_wtheta
            x(:) = 0.0_r_def
            do df0 = 1, ndf_w0
              x(1) = x(1) + chi_1_e(df0)*w0_basis(1,df0,df)
              x(2) = x(2) + chi_2_e(df0)*w0_basis(1,df0,df)
              x(3) = x(3) + chi_3_e(df0)*w0_basis(1,df0,df)
            end do

            theta(map_wtheta(df) + k) = analytic_temperature(x, test)
          end do
        end do

    end subroutine initial_theta_code

end module initial_theta_kernel_mod
