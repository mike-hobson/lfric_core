!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page catalyst ParaView Catalyst visualisation mini app
!> Test program for demonstrating in-situ visualisation with ParaView
!> Catalyst.
!>
!> @brief Main program used to run an example for in-situ visualisation.

program catalyst_demo

  use cli_mod,             only: get_initial_filename
  use catalyst_demo_driver_mod, only: initialise, run, finalise

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

end program catalyst_demo
