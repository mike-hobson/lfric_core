!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> Contains the diagnostic field metadata definition for the diagnostics miniapp
!-------------------------------------------------------------------------------
! The module name consists of:
!   * the science section name
!   * the field group name
!   * meta_mod
! all joined by '__' (double underscores)
module diagnostics_colours__diagnostics__meta_mod

    ! The required use statements include:
    !   * field_meta_data_type to hold the metadata for each field
    !   * enums that define allowed values for certain metadata items
    !       (eg W3, STANDARD_TIMESTEP & BILINEAR below)
    !   * any dimensions that fields live on
    use diagnostics_mod,                only: field_meta_data_type
    use constants_mod,                  only: REAL_TYPE
    use fs_continuity_mod,              only: W3
    use time_step_enum_mod,             only: STANDARD_TIMESTEP
    use interpolation_enum_mod,         only: BILINEAR
    use vertical_dimensions_mod,        only: model_height_dimension

    implicit none

    private

    !> @brief Type that holds the metadata for the diagnostics miniapp
    !> diagnostic fields
    !>
    !> The name of the type uses the same format as the module name
    !> (with _mod replaced with _type)
    !> The type contains an instance of field_meta_data_type for each
    !> field in the field group
    type, public :: diagnostics_colours__diagnostics__meta_type
        type(field_meta_data_type), public :: hex
    end type diagnostics_colours__diagnostics__meta_type

    interface diagnostics_colours__diagnostics__meta_type
        module procedure diagnostics_colours__diagnostics__meta_constructor
    end interface

contains

    !> @brief Constructor for diagnostic miniapp's diagnostics field metadata
    !>
    !> Contains the metadata type for the field group and instantiates each
    !> field with the appropriate metadata
    function diagnostics_colours__diagnostics__meta_constructor() result(self)

        implicit none

        type(diagnostics_colours__diagnostics__meta_type) :: self

        self%hex = field_meta_data_type(&
            unique_id = "diagnostics__hex", &
            long_name = "", &
            units = "1", &
            function_space = W3, &
            order = 0, &
            io_driver = "write_field_face", &
            trigger = "_checksum: true;", &
            description = "A colourful diagnostic field", &
            data_type = REAL_TYPE, &
            time_step = STANDARD_TIMESTEP, &
            recommended_interpolation = BILINEAR, &
            packing = 0, &
            vertical_dimension = model_height_dimension(), &
            standard_name = "hex")

        end function diagnostics_colours__diagnostics__meta_constructor

end module diagnostics_colours__diagnostics__meta_mod