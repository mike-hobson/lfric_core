!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Container for time axes.
!>
module gungho_time_axes_mod

  use linked_list_mod, only : linked_list_type

  implicit none

  private

  !> Collection of time axes.
  !>
  type, public :: gungho_time_axes_type

    private

    !> Time varying ancillaries time axis.
    type(linked_list_type), public :: ancil_times_list

    !> Time varying LBC time axis.
    type(linked_list_type), public :: lbc_times_list

    !> Time varying linearisation state time axis.
    !>
    !> @todo Is this part of the linear model?
    !>
    type(linked_list_type), public :: ls_times_list

  end type gungho_time_axes_type

end module gungho_time_axes_mod
