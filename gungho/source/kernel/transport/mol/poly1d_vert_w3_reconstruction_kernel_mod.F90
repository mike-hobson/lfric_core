!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes vertical reconstructions through fitting a high order 1D
!!        upwind reconstruction.
!> @details Computes the reconstruction for a tracer field using a high order
!!          polynomial fit to the integrated tracer values. The stencil used
!!          for the polynomial is centred on the upwind cell.
!!          Near the boundaries the order of reconstruction may be reduced
!!          if there are not enough points to compute desired order.
!!          This method is only valid for the lowest order elements based on
!!          a reference cube.
module poly1d_vert_w3_reconstruction_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_INC, GH_READ, GH_BASIS,   &
                              CELL_COLUMN, GH_EVALUATOR,   &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              outward_normals_to_vertical_faces
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_w3_reconstruction_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,   W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
       /)
  type(func_type) :: meta_funcs(1) = (/                                 &
       func_type(W2, GH_BASIS)                                          &
       /)
  type(reference_element_data_type) :: meta_reference_element(1) = (/   &
       reference_element_data_type( outward_normals_to_vertical_faces ) &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: poly1d_vert_w3_reconstruction_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly1d_vert_w3_reconstruction_code

contains

!> @brief Computes the vertical reconstructions for a tracer.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] reconstruction Mass reconstruction field to compute
!> @param[in]     wind           Wind field
!> @param[in]     tracer         Tracer field
!> @param[in]     coeff          Array of polynomial coefficients for interpolation
!> @param[in]     ndata          Number of data points per dof location
!> @param[in]     global_order   Desired order of polynomial reconstruction
!> @param[in]     logspace       If true (=1), then perform interpolation in log space;
!!                               this should be a logical but this is not currently supported in
!!                               PSyclone, see Issue #1248
!!                               TODO #3010: the use of logical types is now supported
!!                               so this ticket should remove these integers
!> @param[in]     ndf_w2         Number of degrees of freedom per cell
!> @param[in]     undf_w2        Number of unique degrees of freedom for the reconstruction &
!!                               wind fields
!> @param[in]     map_w2         Dofmap for the cell at the base of the column
!> @param[in]     basis_w2       Basis function array evaluated at w2 nodes
!> @param[in]     ndf_w3         Number of degrees of freedom per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for the tracer field
!> @param[in]     map_w3         Cell dofmaps for the tracer space
!> @param[in]     ndf_c          Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c         Total number of degrees of freedom for the coeff space
!> @param[in]     map_c          Dofmap for the coeff space
!> @param[in]     nfaces_re_v    Number of vertical faces (used by PSyclone to size
!!                               coeff array)
!> @param[in]     outward_normals_to_vertical_faces Vector of normals to the
!!                                                  reference element vertical
!!                                                  "outward faces"
subroutine poly1d_vert_w3_reconstruction_code( nlayers,                           &
                                               reconstruction,                    &
                                               wind,                              &
                                               tracer,                            &
                                               coeff,                             &
                                               ndata,                             &
                                               global_order,                      &
                                               logspace,                          &
                                               ndf_w2,                            &
                                               undf_w2,                           &
                                               map_w2,                            &
                                               basis_w2,                          &
                                               ndf_w3,                            &
                                               undf_w3,                           &
                                               map_w3,                            &
                                               ndf_c,                             &
                                               undf_c,                            &
                                               map_c,                             &
                                               nfaces_re_v,                       &
                                               outward_normals_to_vertical_faces )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), intent(in)                    :: ndf_c
  integer(kind=i_def), intent(in)                    :: undf_c
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_re_v

  real(kind=r_def), dimension(undf_w2), intent(inout) :: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)    :: tracer
  ! ndata = (global_order+1)*nfaces_re_v
  ! (i, j, map(df) + k ) => i - 1 + (j-1)*nfaces_re_v + k*ndata + map_c(df)
  real(kind=r_def), dimension(undf_c),  intent(in)    :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  real(kind=r_def), intent(in)  :: outward_normals_to_vertical_faces(:,:)

  integer(kind=i_def), intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def)            :: k, df, ij, p, stencil, order, id, &
                                    m, ijkp, boundary_offset, vertical_order, &
                                    idx
  integer(kind=i_def)            :: vert_offset
  real(kind=r_def)               :: direction
  real(kind=r_def), dimension(2) :: v_dot_n
  real(kind=r_def)               :: polynomial_tracer

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  vert_offset = ndf_w2 - nfaces_re_v

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  do id = 1,nfaces_re_v
    df = id + vert_offset
    v_dot_n(id) =  dot_product(basis_w2(:,df,df), outward_normals_to_vertical_faces(:,id))
  end do

  ij = map_w3(1)

  ! Vertical reconstruction computation
  do k = 0, nlayers - 1
    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))

    boundary_offset = 0
    if ( order > 0 ) then
      if ( k == 0 )           boundary_offset =  1
      if ( k == nlayers - 1 ) boundary_offset = -1
    end if

    do id = 1,nfaces_re_v

      df = id + vert_offset
      ! Check if this is the upwind cell
      direction = wind(map_w2(df) + k )*v_dot_n(id)
      if ( direction >= 0.0_r_def ) then
        if (logspace == 1_i_def) then
          ! Interpolate log(tracer)
          ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
          !                 = tracer_1**c_1*tracer_2**c_2...
          ! Note that we further take the absolute value before raising to the
          ! fractional power. This code should only be used for a positive
          ! quantity, but adding in the abs ensures no errors are thrown
          ! if negative numbers are passed through in redundant calculations
          ! in the haloes
          polynomial_tracer = 1.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
            polynomial_tracer = polynomial_tracer &
                              * abs(tracer( ijkp ))**coeff( idx )
          end do
        else
          polynomial_tracer = 0.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
            polynomial_tracer = polynomial_tracer &
                              + tracer( ijkp )*coeff( idx )
          end do
        end if
        reconstruction(map_w2(df) + k ) = polynomial_tracer
      end if
    end do
  end do

  ! Boundary conditions
  ! When computing a reconstruction the wind = 0 so reconstruction is zero
  ! For averaging the wind = +/-1 so the resulting average is non-zero
  ! Make sure we always compute the boundary values
  order = min(vertical_order, 2)
  boundary_offset = 0

  ! The boundary values only have a cell on one side of them so the
  ! reconstruction above may not have performed. Therefore we ensure
  ! that the values at the boundaries are computed

  k = 0
  id = 1
  if ( order > 0 ) boundary_offset = 1
  df = id + vert_offset
  if (logspace == 1_i_def) then
    ! Interpolate log(tracer)
    ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
    !                 = tracer_1**c_1*tracer_2**c_2...
    polynomial_tracer = 1.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_tracer = polynomial_tracer &
                        * abs(tracer( ijkp ))**coeff( idx )
    end do
  else
    polynomial_tracer = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_tracer = polynomial_tracer &
                        + tracer( ijkp )*coeff( idx )
    end do
  end if
  reconstruction(map_w2(df) + k ) = polynomial_tracer

  k = nlayers-1
  id = 2
  if ( order > 0 ) boundary_offset =  -1
  df = id + vert_offset
  if (logspace == 1_i_def) then
    ! Interpolate log(tracer)
    ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
    !                 = tracer_1**c_1*tracer_2**c_2...
    polynomial_tracer = 1.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_tracer = polynomial_tracer &
                        * abs(tracer( ijkp ))**coeff( idx )
    end do
  else
    polynomial_tracer = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_tracer = polynomial_tracer &
                        + tracer( ijkp )*coeff( idx )
    end do
  end if
  reconstruction(map_w2(df) + k ) = polynomial_tracer

  deallocate( smap )

end subroutine poly1d_vert_w3_reconstruction_code

end module poly1d_vert_w3_reconstruction_kernel_mod
