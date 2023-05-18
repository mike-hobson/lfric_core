!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Receives a field from another coupled component
!> @details  Takes incoming field data from Oasis (sent from another coupled
!>           component), re-orders it and copies the data into an LFRic field
!>           object.
! Note that the word "component", in this module, refers to the model
! components that are being coupled together by Oasis.

module cpl_field_receive_mod

#ifdef MCT
  use mod_oasis,     only: oasis_get, oasis_recvd, oasis_recvout
#endif
  use field_mod,     only: field_type, field_proxy_type
  use constants_mod, only: i_def, r_def, l_def, imdi
  use log_mod,       only: log_event,       &
                           LOG_LEVEL_DEBUG, &
                           LOG_LEVEL_ERROR, &
                           log_scratch_space

  implicit none

  private

  public :: cpl_field_receive

  contains

  !> @brief Receives field from another coupled component
  !>
  !> @param [in]     rfield       Field to be received
  !> @param [in]     mtime        Current model time
  !> @param [out]    ldex         Field exchange flag
  !> @param [in,out] ldfail       Failure flag for receive operation
  !> @param [in]     icpl_size    Length of the coupling data being received
  !> @param [in]     slength      Max length of coupling field names
  !> @param [in]     slocal_index Index to sort data for sending
  !
  subroutine cpl_field_receive( rfield, mtime, ldex, ldfail, &
                                icpl_size, slength, slocal_index)
  implicit none

  type(field_type), intent(inout)              :: rfield
  integer(i_def),   intent(in)                 :: mtime
  logical(l_def),   intent(out)                :: ldex
  logical(l_def),   intent(inout)              :: ldfail
  integer(i_def),   intent(in)                 :: icpl_size
  integer(i_def),   intent(in)                 :: slength
  integer(i_def),   intent(in)                 :: slocal_index(:)

#ifdef MCT
  ! Number of data-levels
  integer(i_def)                               :: nlev
  ! Name of the verialbe to be sent
  character(len=slength)                       :: rname
  ! Proxy of the field
  type(field_proxy_type)                       :: rfield_proxy
  ! Oasis id for varialble or data level
  integer(i_def)                               :: rvar_id
  ! Error return by oasis_get
  integer(i_def)                               :: ierror
  ! Data received from OASIS
  real(r_def)                                  :: rdata(icpl_size)
  ! Data ordered according to LFRic convention
  real(r_def)                                  :: wdata(icpl_size)
  ! Index over coupling data length
  integer(i_def)                               :: i
  ! Index over data levels
  integer(i_def)                               :: k

  rname        = trim(adjustl(rfield%get_name()))
  rfield_proxy = rfield%get_proxy()
  nlev         = rfield_proxy%vspace%get_ndata()

  ldex = .false.

  do k = 1, nlev
    rvar_id = rfield%get_cpl_id(k)
    if (rvar_id /= imdi) then
      call oasis_get(rvar_id, mtime, rdata(:), ierror)
      if (ierror == oasis_recvd .or. ierror == oasis_recvout) then
        do i = 1, icpl_size
          wdata(slocal_index(i)) = rdata(i)
        enddo
        rfield_proxy%data(k:icpl_size*nlev:nlev) = wdata(:)
        ldex = .true.
        write(log_scratch_space, '(3A)' ) "cpl_field_receive: field ", &
                           trim(rname), " received"
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      else
        write(log_scratch_space, '(3A)' ) "cpl_field_receive: field ", &
                           trim(rname), " NOT exchanged on this timestep"
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      endif
    else
      ldfail = .true.
      write(log_scratch_space, '(3A)' ) "PROBLEM cpl_field_receive: field ", &
                                         trim(rname), " cpl_id NOT set"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
  enddo

  if (.not. ldfail) call rfield_proxy%set_dirty()
#else
  ldex = .false.
  write(log_scratch_space, '(A)' ) &
               "cpl_field_receive: to use OASIS cpp directive MCT must be set"
  call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif
  end subroutine cpl_field_receive

end module cpl_field_receive_mod
