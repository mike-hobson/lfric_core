!-------------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> Contains the non-spatial field metadata definition for the diagnostics miniapp
!-------------------------------------------------------------------------------
! The module name consists of:
!   * the science section name
!   * the field group name
!   * meta_mod
! all joined by '__' (double underscores)
module colours__non_spatial__meta_mod

    ! The required use statements include:
    !   * field_meta_data_type to hold the metadata for each field
    !   * enums that define allowed values for certain metadata items
    !       (eg W3, STANDARD_TIMESTEP & BILINEAR below)
    !   * any dimensions that fields live on
    use diagnostics_mod,                only: field_meta_data_type
    use constants_mod,                  only: REAL_TYPE, str_def, &
                                              str_short, r_def
    use fs_continuity_mod,              only: W3
    use time_step_enum_mod,             only: STANDARD_TIMESTEP
    use interpolation_enum_mod,         only: BILINEAR
    use non_spatial_dimension_mod,      only: non_spatial_dimension_type, &
                                              NUMERICAL, CATEGORICAL

    implicit none

    private

    !> @brief Type that holds the metadata for the diagnostics miniapp
    !> non-spatial fields
    !>
    !> The name of the type uses the same format as the module name
    !> (with _mod replaced with _type)
    !> The type contains an instance of field_meta_data_type for each
    !> field in the field group and a string containing the field group name,
    !> in the form <section>__<group>
    type, public :: colours__non_spatial__meta_type
        type(field_meta_data_type), public :: mutable_numbers, &
                                              mutable_categories, &
                                              immutable_both
        character(str_def) :: name = "colours__non_spatial"
    end type colours__non_spatial__meta_type

    interface colours__non_spatial__meta_type
        module procedure colours__non_spatial__meta_constructor
    end interface

contains

    !> @brief Constructor for diagnostic miniapp's non-spatial field metadata
    !>
    !> Contains the metadata type for the field group and instantiates each
    !> field with the appropriate metadata
    function colours__non_spatial__meta_constructor() result(self)

        implicit none

        type(colours__non_spatial__meta_type) :: self

        self%mutable_categories = field_meta_data_type( &
            unique_id = "colours__mutable_categories", &
            long_name = "mutable_categories", &
            units = "arbitrary units", &
            function_space = W3, &
            order = 0, &
            io_driver = "write_field_single_face", &
            trigger = "__checksum: true;", &
            description = "A diagnostic field with a mutable categorical &
                    &non-spatial dimension", &
            data_type = REAL_TYPE, &
            time_step = STANDARD_TIMESTEP, &
            recommended_interpolation = BILINEAR, &
            packing = 0, &
            non_spatial_dimension = [non_spatial_dimension_type( &
                    dimension_name = "mutable_categories", &
                    dimension_category = CATEGORICAL, &
                    help_text = "Mutable categories help text")])

        self%mutable_numbers = field_meta_data_type( &
            unique_id = "colours__mutable_numbers", &
            long_name = "mutable_numbers", &
            units = "arbitrary units", &
            function_space = W3, &
            order = 0, &
            io_driver = "write_field_single_face", &
            trigger = "__checksum: true;", &
            description = "A diagnostic field with a mutable numerical &
                    &non-spatial dimension", &
            data_type = REAL_TYPE, &
            time_step = STANDARD_TIMESTEP, &
            recommended_interpolation = BILINEAR, &
            packing = 0, &
            non_spatial_dimension = [non_spatial_dimension_type( &
                    dimension_name = "mutable_numbers", &
                    dimension_category = NUMERICAL, &
                    help_text = "Mutable numbers help text", &
                    non_spatial_units = "other units")])

        self%immutable_both = field_meta_data_type( &
            unique_id = "colours__immutable_both", &
            long_name = "immutable_both", &
            units = "arbitrary units", &
            function_space = W3, &
            order = 0, &
            io_driver = "write_field_single_face", &
            trigger = "__checksum: true;", &
            description = "A diagnostic field with immutable categorical &
                    &and numerical non-spatial dimensions", &
            data_type = REAL_TYPE, &
            time_step = STANDARD_TIMESTEP, &
            recommended_interpolation = BILINEAR, &
            packing = 0, &
            non_spatial_dimension = [non_spatial_dimension_type( &
                    dimension_name = "immutable_categories", &
                    dimension_category = CATEGORICAL, &
                    help_text = "Immutable categories help text", &
                    label_definition = [character(str_short) :: 'a','b','c']), &
                                     non_spatial_dimension_type( &
                    dimension_name = "immutable_numbers", &
                    dimension_category = NUMERICAL, &
                    help_text = "Immutable numbers help text", &
                    axis_definition = [real(r_def) :: 1,2,3], &
                    non_spatial_units = "other units")])

    end function colours__non_spatial__meta_constructor

end module colours__non_spatial__meta_mod
