!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Conway knows what configuration it needs.
!>
module conway_mod

  implicit none

  private

  character(*), public, parameter ::                        &
      conway_required_namelists(5) =  [ 'base_mesh     ', &
                                          'extrusion     ', &
                                          'finite_element', &
                                          'partitioning  ', &
                                          'planet        ']

end module conway_mod
