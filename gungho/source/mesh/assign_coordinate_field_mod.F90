!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module to assign the values of the coordinates of the mesh to a field
module assign_coordinate_field_mod

  use base_mesh_config_mod, only : geometry, &
                                   geometry_spherical
  use constants_mod,        only : r_def, i_def, i_native, l_def
  use log_mod,              only : log_event, LOG_LEVEL_ERROR
  use planet_config_mod,    only : scaled_radius
  use mesh_collection_mod,  only : mesh_collection
  use coord_transform_mod,  only : xyz2llr, identify_panel, &
                                   xyz2alphabetar, alphabetar2xyz
  use finite_element_config_mod, only: spherical_coord_system, &
                                       spherical_coord_system_xyz, &
                                       spherical_coord_system_abh, &
                                       spherical_coord_system_llh
  implicit none

  private

  public :: assign_coordinate_field
  ! Make this public only for unit-testing
  public :: assign_coordinate_xyz

contains

!> @brief Subroutine which assigns the values of the coordinates of the mesh
!! to a field
!> @details An array of size 3 for the type field is passed in to be populated.
!! The field proxy is used to break encapsulation and access the function space
!! and the data attributes of the field so that its values can be assigned.
!! calls two subroutines, get_cell_coords from the mesh generator and then
!! assign_coordinate on a column by column basis
!! @param[in,out] chi      Model coordinate array of size 3 (x,y,z) of fields
!! @param[in]     panel_id Field giving the ID of mesh panels
!! @param[in]     mesh_id  Id of mesh on which this field is attached
!! @param[in]     spherical_coords Logical labelling whether this is the
!!                                 spherical coordinate field
  subroutine assign_coordinate_field(chi, panel_id, mesh_id, spherical_coords)

    use field_mod,             only: field_type, field_proxy_type
    use reference_element_mod, only: reference_element_type
    use mesh_mod,              only: mesh_type
    use mesh_constructor_helper_functions_mod, &
                               only: domain_size_type
    implicit none

    type( field_type ),  intent( inout )        :: chi(3)
    type( field_type ),  intent( inout )        :: panel_id
    integer(kind=i_def), intent( in )           :: mesh_id
    logical(kind=l_def), intent( in ), optional :: spherical_coords

    integer(kind=i_def),                pointer :: map(:,:)          => null()
    integer(kind=i_def),                pointer :: map_pid(:,:)      => null()
    real(kind=r_def),                   pointer :: dof_coords(:,:)   => null()
    type(mesh_type),                    pointer :: mesh              => null()
    class(reference_element_type),      pointer :: reference_element => null()

    type(field_proxy_type) :: chi_proxy(3)
    type(field_proxy_type) :: panel_id_proxy
    type(domain_size_type) :: domain_size

    real(kind=r_def), allocatable :: column_coords(:,:,:)
    real(kind=r_def), allocatable :: dz(:)  ! dz(nlayers) array
    real(kind=r_def), allocatable :: vertex_coords(:,:)

    integer(kind=i_def) :: cell
    integer(kind=i_def) :: undf, ndf, nlayers
    integer(kind=i_def) :: undf_pid, ndf_pid, nlayers_pid
    integer(kind=i_def) :: nverts

    integer(kind=i_native) :: alloc_error
    integer(kind=i_def)    :: depth
    logical(kind=l_def)    :: spherical_coords_flag


    ! Break encapsulation and get the proxy.
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    undf    = chi_proxy(1)%vspace%get_undf()
    ndf     = chi_proxy(1)%vspace%get_ndf()
    nlayers = chi_proxy(1)%vspace%get_nlayers()

    panel_id_proxy = panel_id%get_proxy()
    undf_pid = panel_id_proxy%vspace%get_undf()
    ndf_pid  = panel_id_proxy%vspace%get_ndf()
    nlayers_pid = panel_id_proxy%vspace%get_nlayers()

    map => chi_proxy(1)%vspace%get_whole_dofmap()
    map_pid => panel_id_proxy%vspace%get_whole_dofmap()

    allocate ( dz(nlayers), STAT = alloc_error )
    if ( alloc_error /= 0 ) then
      call log_event( " assign_coordinate_field: Unable to allocate "// &
                      "local array dz(nlayers) ", LOG_LEVEL_ERROR )
    end if

    mesh => mesh_collection%get_mesh( mesh_id )
    call mesh%get_dz(dz)

    reference_element => mesh%get_reference_element()
    call reference_element%get_vertex_coordinates( vertex_coords )
    nverts = reference_element%get_number_vertices()

    allocate( column_coords(3,nverts,nlayers ) )
    dof_coords => chi_proxy(1)%vspace%get_nodes( )

    domain_size =  mesh%get_domain_size()

    ! Set spherical coordinates flag
    spherical_coords_flag = .false.
    if (present(spherical_coords)) spherical_coords_flag = spherical_coords

    if ( (.not. spherical_coords_flag) .or. &
          ( spherical_coord_system == spherical_coord_system_xyz ) ) then

      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call mesh%get_column_coords(cell,column_coords)

        call calc_panel_id( nlayers_pid, nverts,   &
                            ndf_pid, undf_pid,     &
                            map_pid(:,cell),       &
                            column_coords(:,:,1),  &
                            panel_id_proxy%data    )

        call assign_coordinate_xyz( nlayers,                 &
                                    ndf,                     &
                                    nverts,                  &
                                    undf,                    &
                                    map(:,cell),             &
                                    dz,                      &
                                    chi_proxy(1)%data,       &
                                    chi_proxy(2)%data,       &
                                    chi_proxy(3)%data,       &
                                    column_coords,           &
                                    dof_coords,              &
                                    vertex_coords,           &
                                    domain_size%maximum%x,   &
                                    domain_size%minimum%y,   &
                                    panel_id_proxy%data,     &
                                    ndf_pid,                 &
                                    undf_pid,                &
                                    map_pid(:,cell)          )
      end do

    else

      if ( spherical_coord_system /= spherical_coord_system_abh ) then
        call log_event('Your spherical coordinate system has not been implemented yet', LOG_LEVEL_ERROR)
      end if

      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call mesh%get_column_coords(cell,column_coords)

        call calc_panel_id( nlayers_pid, nverts,   &
                            ndf_pid, undf_pid,     &
                            map_pid(:,cell),       &
                            column_coords(:,:,1),  &
                            panel_id_proxy%data    )

        call assign_coordinate_sph( nlayers,                 &
                                    ndf,                     &
                                    nverts,                  &
                                    undf,                    &
                                    map(:,cell),             &
                                    chi_proxy(1)%data,       &
                                    chi_proxy(2)%data,       &
                                    chi_proxy(3)%data,       &
                                    column_coords,           &
                                    dof_coords,              &
                                    vertex_coords,           &
                                    panel_id_proxy%data,     &
                                    ndf_pid,                 &
                                    undf_pid,                &
                                    map_pid(:,cell)          )
      end do

    end if

    ! As we have correctly set the chi fields into their full halos,
    ! mark their halos as clean, out to the full halo depth
    ! This is necessary so that subsequent kernel calls don't try to
    ! halo_swap the Wchi field which is read-only
    depth = mesh%get_halo_depth()
    call chi_proxy(1)%set_clean(depth)
    call chi_proxy(2)%set_clean(depth)
    call chi_proxy(3)%set_clean(depth)

    deallocate ( dz, column_coords, vertex_coords )

  end subroutine assign_coordinate_field

!> @brief Assigns the cubed sphere panel ID values to the panel_id field
!> @details A scalar field is passed in and all values are assigned to
!! be the panel IDs which are calculated from the coordinates. For planar
!! geometry the ID is just 1 everywhere.
!! @param[in]     nlayers   Number of layers for the panel_id field
!! @param[in]     nverts    Number of reference element vertices
!! @param[in]     ndf_pid   Number of DoFs per cell for the panel_id field
!! @param[in]     undf_pid  Universal number of DoFs for the panel_id field
!! @param[in]     map_pid   DoF map for the panel_id field
!! @param[in]     column_base_coords Coordinates for vertices at base of column
!! @param[in,out] panel_id  Field (to be calculated) with the ID of cubed sphere panels
  subroutine calc_panel_id( nlayers,            &
                            nverts,             &
                            ndf_pid,            &
                            undf_pid,           &
                            map_pid,            &
                            column_base_coords, &
                            panel_id            )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers, nverts, ndf_pid, undf_pid
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)

    real(kind=r_def), intent(in)  :: column_base_coords(3,nverts,1)
    real(kind=r_def), intent(out) :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: vert, panel

    real(kind=r_def) :: interp_weight, x, y, z

    ! Set all panel_ids to 1 for non-spherical geometry
    if ( geometry /= geometry_spherical ) then
      panel_id(map_pid(1):map_pid(1)+nlayers-1) = 1.0_r_def
    else
      ! Find coordinates at centre of cell at bottom of column
      x = 0.0_r_def
      y = 0.0_r_def
      z = 0.0_r_def
      interp_weight = 1.0_r_def/nverts

      do vert = 1,nverts
        ! evaluate at cell centre
        x = x + interp_weight*column_base_coords(1,vert,1)
        y = y + interp_weight*column_base_coords(2,vert,1)
        z = z + interp_weight*column_base_coords(3,vert,1)
      end do

      panel = identify_panel(x, y, z)
      panel_id(map_pid(1):map_pid(1)+nlayers-1) = real(panel,r_def)
    end if

  end subroutine calc_panel_id

!> @brief Determines and assigns the (X,Y,Z) coordinates for a single column
!! @param[in]  nlayers       integer: loop bound
!! @param[in]  ndf           integer: array size and loop bound
!! @param[in]  nverts        integer: array size and loop bound
!! @param[in]  undf          integer: array size and loop bound
!! @param[in]  map           integer array: indirection map
!! @param[in]  dz            Mesh layer thickness
!! @param[out] chi_1         real array: size undf x coord
!! @param[out] chi_2         real array: size undf y coord
!! @param[out] chi_3         real array: size undf z coord
!! @param[in]  column_coords real array: (3,nverts,nlayers)
!! @param[in]  chi_hat_node  real array: (3,ndf)
!! @param[in]  chi_hat_vert  real array: (nverts,3)
!! @param[in]  domain_x      real: domain extent in x direction for planar mesh
!! @param[in]  domain_y      real: domain extent in y direction for planar mesh
!! @param[in]  panel_id      real: field giving ID of mesh panel
!! @param[in]  ndf_pid       integer: number of DoFs per cell for panel_id space
!! @param[in]  undf_pid      integer: number of universal DoFs for panel_id space
!! @param[in]  map_pid       integer: DoF map for panel_id space
  subroutine assign_coordinate_xyz( nlayers,       &
                                    ndf,           &
                                    nverts,        &
                                    undf,          &
                                    map,           &
                                    dz,            &
                                    chi_1,         &
                                    chi_2,         &
                                    chi_3,         &
                                    column_coords, &
                                    chi_hat_node,  &
                                    chi_hat_vert,  &
                                    domain_x,      &
                                    domain_y,      &
                                    panel_id,      &
                                    ndf_pid,       &
                                    undf_pid,      &
                                    map_pid        )

    use reference_element_mod, only: SWB, SEB, NEB, NWB, SWT, SET, NET, NWT

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nlayers, ndf, nverts, undf
    integer(kind=i_def), intent(in)  :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)  :: map(ndf), map_pid(ndf_pid)
    real(kind=r_def),    intent(in)  :: dz(nlayers)
    real(kind=r_def),    intent(out) :: chi_1(undf), chi_2(undf), chi_3(undf)
    real(kind=r_def),    intent(in)  :: column_coords(3,nverts,nlayers)
    real(kind=r_def),    intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)
    real(kind=r_def),    intent(in)  :: domain_x, domain_y
    real(kind=r_def),    intent(in)  :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: k, df, dfk, vert, panel

    real(kind=r_def)    :: interp_weight, x, y, z, radius_correction
    real(kind=r_def)    :: alpha, beta, radius
    real(kind=r_def)    :: v_x, v_y, v_z, v_a, v_b, v_r
    real(kind=r_def)    :: vertex_local_coords(3,nverts)

    radius_correction = 1.0_r_def
    panel = int(panel_id(map_pid(1)))

    ! Compute the representation of the coordinate field
    do k = 0, nlayers-1
      vertex_local_coords(:,:) = column_coords(:,:,k+1)
      if (  geometry /= geometry_spherical ) then
        ! Check if point cell is on right or bottom boundary,
        ! assumes a monotonic coordinate field
        if ( column_coords(1,SEB,k+1) < column_coords(1,SWB,k+1)) then
        ! On x boundary
          vertex_local_coords(1,SEB) = domain_x
          vertex_local_coords(1,NEB) = domain_x
          vertex_local_coords(1,SET) = domain_x
          vertex_local_coords(1,NET) = domain_x
        end if
        if ( column_coords(2,SWB,k+1) > column_coords(2,NWB,k+1)) then
        ! On y boundary
          vertex_local_coords(2,SWB) = domain_y
          vertex_local_coords(2,SEB) = domain_y
          vertex_local_coords(2,SWT) = domain_y
          vertex_local_coords(2,SET) = domain_y
        end if
      end if

      if ( spherical_coord_system == spherical_coord_system_abh ) then
        do df = 1, ndf
          dfk = map(df)+k
          ! Compute interpolation weights
          alpha = 0.0_r_def
          beta = 0.0_r_def
          radius = 0.0_r_def
          do vert = 1,nverts
             interp_weight = &
                   (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))
             v_x = vertex_local_coords(1,vert)
             v_y = vertex_local_coords(2,vert)
             v_z = vertex_local_coords(3,vert)

             call xyz2alphabetar(v_x,v_y,v_z,panel,v_a,v_b,v_r)
             alpha = alpha + interp_weight*v_a
             beta = beta + interp_weight*v_b
             radius = radius + interp_weight*v_r
          end do

          call alphabetar2xyz(alpha, beta, radius, panel, &
                              chi_1(dfk), chi_2(dfk), chi_3(dfk))

        end do
      else

        do df = 1, ndf
          ! Compute interpolation weights
          x = 0.0_r_def
          y = 0.0_r_def
          z = 0.0_r_def
          do vert = 1,nverts
             interp_weight = &
                   (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))

             x = x + interp_weight*vertex_local_coords(1,vert)
             y = y + interp_weight*vertex_local_coords(2,vert)
             z = z + interp_weight*vertex_local_coords(3,vert)
          end do
          ! For spherical domains we need to project x,y,z back onto
          ! spherical shells
          if ( geometry == geometry_spherical ) then
             radius_correction = scaled_radius + &
                                 sum(dz(1:k)) + chi_hat_node(3,df)*dz(k+1)
             radius_correction = radius_correction/sqrt(x*x + y*y + z*z)
          end if
          dfk = map(df)+k
          chi_1(dfk) = x*radius_correction
          chi_2(dfk) = y*radius_correction
          chi_3(dfk) = z*radius_correction

        end do
      end if ! Spherical coord system

    end do

  end subroutine assign_coordinate_xyz

!> @brief Determines and assigns the (X,Y,Z) coordinates for a single column
!! @param[in]  nlayers       integer: loop bound
!! @param[in]  ndf           integer: array size and loop bound
!! @param[in]  nverts        integer: array size and loop bound
!! @param[in]  undf          integer: array size and loop bound
!! @param[in]  map           integer array: indirection map
!! @param[out] chi_1         real array: size undf x coord
!! @param[out] chi_2         real array: size undf y coord
!! @param[out] chi_3         real array: size undf z coord
!! @param[in]  column_coords real array: (3,nverts,nlayers)
!! @param[in]  chi_hat_node  real array: (3,ndf)
!! @param[in]  chi_hat_vert  real array: (nverts,3)
!! @param[in]  panel_id      real: field giving ID of cubed sphere panel
!! @param[in]  ndf_pid       integer: number of DoFs per cell for panel_id space
!! @param[in]  undf_pid      integer: number of universal DoFs for panel_id space
!! @param[in]  map_pid       integer: DoF map for panel_id space
  subroutine assign_coordinate_sph( nlayers,       &
                                    ndf,           &
                                    nverts,        &
                                    undf,          &
                                    map,           &
                                    chi_1,         &
                                    chi_2,         &
                                    chi_3,         &
                                    column_coords, &
                                    chi_hat_node,  &
                                    chi_hat_vert,  &
                                    panel_id,      &
                                    ndf_pid,       &
                                    undf_pid,      &
                                    map_pid        )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nlayers, ndf, nverts, undf
    integer(kind=i_def), intent(in)  :: map(ndf)
    real(kind=r_def),    intent(out) :: chi_1(undf), chi_2(undf), chi_3(undf)
    real(kind=r_def),    intent(in)  :: column_coords(3,nverts,nlayers)
    real(kind=r_def),    intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)
    integer(kind=i_def), intent(in)  :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)  :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)  :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: k, df, dfk, vert
    integer(kind=i_def) :: panel
    real(kind=r_def)    :: alpha, beta, radius

    real(kind=r_def) :: interp_weight
    real(kind=r_def) :: v_x, v_y, v_z, v_a, v_b, v_r
    real(kind=r_def) :: vertex_local_coords(3,nverts)

    panel = int(panel_id(map_pid(1)))

    ! Compute the representation of the coordinate field
    do k = 0, nlayers-1
       vertex_local_coords(:,:) = column_coords(:,:,k+1)

       do df = 1, ndf
          dfk = map(df)+k
          ! Compute interpolation weights
          alpha = 0.0_r_def
          beta = 0.0_r_def
          radius = 0.0_r_def
          do vert = 1,nverts
             interp_weight = &
                   (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))
             v_x = vertex_local_coords(1,vert)
             v_y = vertex_local_coords(2,vert)
             v_z = vertex_local_coords(3,vert)

             call xyz2alphabetar(v_x,v_y,v_z,panel,v_a,v_b,v_r)
             alpha = alpha + interp_weight*v_a
             beta = beta + interp_weight*v_b
             radius = radius + interp_weight*v_r
          end do

          chi_1(dfk) = alpha
          chi_2(dfk) = beta
          chi_3(dfk) = radius - scaled_radius

        end do

     end do

   end subroutine assign_coordinate_sph

end module assign_coordinate_field_mod
