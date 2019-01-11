!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Init functionality for gungho model
!> @details Creates and initialises prognostic fields, also runtime_constants
!>          specific to the model

module init_gungho_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type, &
                                             write_diag_interface, &
                                             checkpoint_interface, &
                                             restart_interface
  use finite_element_config_mod,      only : element_order, &
                                             vorticity_in_w1
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta
  use function_space_collection_mod , only : function_space_collection
  use init_prognostic_fields_alg_mod, only : init_prognostic_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO
  use mr_indices_mod,                 only : nummr
  use moist_dyn_mod,                  only : num_moist_factors
  use runtime_constants_mod,          only : create_runtime_constants
  use io_mod,                         only : xios_write_field_node, &
                                             xios_write_field_face, &
                                             checkpoint_xios,   &
                                             checkpoint_netcdf, &
                                             restart_netcdf,    &
                                             restart_xios
  use io_config_mod,                  only : use_xios_io,     &
                                             write_diag,      &
                                             checkpoint_read, &
                                             checkpoint_write

  implicit none


contains

  !>@brief Initialise the gungho model
  !> @param[in] mesh_id Identifier of the mesh
  !> @param[in,out] chi Spatial coordinates
  !> @param[in,out] u Wind field
  !> @param[in,out] rho Density field
  !> @param[in,out] theta Potential temperature field
  !> @param[in,out] exner Exner pressure field
  !> @param[in,out] mr Moisture mixing ratios
  !> @param[in,out] moist_dyn Auxilliary fields for moist dynamics
  !> @param[in,out] xi Vorticity
  subroutine init_gungho( mesh_id, chi, u, rho, theta, exner, mr, &
                          moist_dyn, xi )

    implicit none

    integer(i_def), intent(in)               :: mesh_id
    ! Prognostic fields
    type( field_type ), intent(inout)        :: u, rho, theta, exner
    type( field_type ), intent(inout)        :: mr(nummr)
    ! Diagnostic fields
    type( field_type ), intent(inout)        :: xi
    type( field_type ), intent(inout)        :: moist_dyn(num_moist_factors)
    ! Coordinate fields
    type( field_type ), intent(in)           :: chi(:)

    integer(i_def)                           :: imr

    procedure(write_diag_interface), pointer :: tmp_write_diag_ptr
    procedure(checkpoint_interface), pointer :: tmp_checkpoint_ptr
    procedure(restart_interface), pointer    :: tmp_restart_ptr

    call log_event( 'GungHo: initialisation...', LOG_LEVEL_INFO )

    ! Create prognostic fields
    theta = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, Wtheta) )
    if ( vorticity_in_w1 ) then
      xi    = field_type( vector_space = &
                          function_space_collection%get_fs(mesh_id, element_order, W1) )
    else
      xi    = field_type( vector_space = &
                          function_space_collection%get_fs(mesh_id, element_order, W2) )
    end if
    u     = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, W2) )
    rho   = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, W3) )
    exner = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, W3) )

    do imr = 1,nummr
      mr(imr) = field_type( vector_space = &
      function_space_collection%get_fs(mesh_id, element_order, theta%which_function_space()) )
    end do

    ! Auxilliary fields holding moisture-dependent factors for dynamics
    do imr = 1, num_moist_factors
      moist_dyn(imr) = field_type( vector_space = &
      function_space_collection%get_fs(mesh_id, element_order, theta%which_function_space()) )
    end do

    ! Set I/O behaviours for diagnostic output

    if (write_diag .and. use_xios_io) then

       ! Set diagnostic output handlers
 
       ! Face domain

       tmp_write_diag_ptr => xios_write_field_face

       ! Vector fields that are projected to scalar components
       call xi%set_write_diag_behaviour(tmp_write_diag_ptr)
       call u%set_write_diag_behaviour(tmp_write_diag_ptr)

       ! Scalar fields
       call rho%set_write_diag_behaviour(tmp_write_diag_ptr)
       call exner%set_write_diag_behaviour(tmp_write_diag_ptr)

       ! Theta is a special case as it can be on face (if function space is WTheta)
       ! or node (if function space is W0)
       if (theta%which_function_space() == Wtheta) then

         call theta%set_write_diag_behaviour(tmp_write_diag_ptr)

       else

        tmp_write_diag_ptr => xios_write_field_node

        call theta%set_write_diag_behaviour(tmp_write_diag_ptr)

       end if

       ! Moisture uses the same type of field write as Theta

       call theta%get_write_diag_behaviour(tmp_write_diag_ptr)

       do imr = 1,nummr
         call mr(imr)%set_write_diag_behaviour(tmp_write_diag_ptr)
       end do

    end if 

    if ( checkpoint_write .or. checkpoint_read) then

      if ( use_xios_io ) then

        ! Use XIOS for checkpoint / restart

        tmp_checkpoint_ptr => checkpoint_xios
        tmp_restart_ptr => restart_xios

        call log_event( 'GungHo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )

      else

        ! Use old checkpoint and restart methods

        tmp_checkpoint_ptr => checkpoint_netcdf
        tmp_restart_ptr => restart_netcdf

       call log_event( 'GungHo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )

      end if

      call xi%set_checkpoint_behaviour(tmp_checkpoint_ptr)
      call u%set_checkpoint_behaviour(tmp_checkpoint_ptr)
      call rho%set_checkpoint_behaviour(tmp_checkpoint_ptr)
      call theta%set_checkpoint_behaviour(tmp_checkpoint_ptr)
      call exner%set_checkpoint_behaviour(tmp_checkpoint_ptr)

      call xi%set_restart_behaviour(tmp_restart_ptr)
      call u%set_restart_behaviour(tmp_restart_ptr)
      call rho%set_restart_behaviour(tmp_restart_ptr)
      call theta%set_restart_behaviour(tmp_restart_ptr)
      call exner%set_restart_behaviour(tmp_restart_ptr)

      do imr = 1,nummr
        call mr(imr)%set_checkpoint_behaviour(tmp_checkpoint_ptr)
        call mr(imr)%set_restart_behaviour(tmp_restart_ptr)
      end do

    end if
    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, chi)

    ! Initialise prognostic fields
    call init_prognostic_fields_alg( u, rho, theta, exner, mr, moist_dyn, &
                                     xi)

    nullify( tmp_write_diag_ptr, tmp_checkpoint_ptr, tmp_restart_ptr )

    call log_event( 'Gungho initialised', LOG_LEVEL_INFO )

  end subroutine init_gungho
end module init_gungho_mod
