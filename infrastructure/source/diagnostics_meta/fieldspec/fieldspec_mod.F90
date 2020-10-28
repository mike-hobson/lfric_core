!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the fieldspec class
!>
!> @details A fieldspec contains all relevant information to create a field


module fieldspec_mod

  use constants_mod,        only: i_def, str_def
  use linked_list_data_mod, only: linked_list_data_type

  implicit none

  private

  type, extends(linked_list_data_type), public :: fieldspec_type

    private

    !> Unique id used by the diagnostic system to identify the field
    character(str_def) :: unique_id
    !> Other information used to create a field
    integer(i_def)     :: mesh_id
    integer(i_def)     :: function_space
    integer(i_def)     :: order
    integer(i_def)     :: field_kind
    integer(i_def)     :: field_type
    integer(i_def)     :: io_driver

  contains

    !> Getter to return the unique_id
    procedure, public :: get_unique_id

    !> Getter to return the mesh_id
    procedure, public :: get_mesh_id

    !> Getter to return the function_space
    procedure, public :: get_function_space

    !> Getter to return the order
    procedure, public :: get_order

    !> Getter to return the kind
    procedure, public :: get_kind

    !> Getter to return the type
    procedure, public :: get_type

    !> Getter to return the io_driver
     procedure, public :: get_io_driver

  end type fieldspec_type

  interface fieldspec_type
    module procedure fieldspec_constructor
  end interface


contains


  !> Construct a <code>fieldspec_type</code> object.
  !>
  !> @param [in] unique_id A unique identifer for the field
  !> @param [in] mesh_id The mesh to create field with
  !> @param [in] function_space The function space to create field with
  !> @param [in] order The order of the field
  !> @param [in] field_kind The kind of the field
  !> @param [in] field_type The type of the field
  !> @return self the fieldspec object
  !>
  function fieldspec_constructor( unique_id, mesh_id,     &
                                  function_space, order,          &
                                  field_kind, field_type, &
                                  io_driver ) &
                                 result(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    character(*),               intent(in)    :: unique_id
    integer(i_def),             intent(in)    :: mesh_id
    integer(i_def),             intent(in)    :: function_space
    integer(i_def),             intent(in)    :: order
    integer(i_def),             intent(in)    :: field_kind
    integer(i_def),             intent(in)    :: field_type
    integer(i_def),             intent(in)    :: io_driver

    type(fieldspec_type), target :: self

    self%unique_id      = trim(unique_id)
    self%mesh_id        = mesh_id
    self%function_space = function_space
    self%order          = order
    self%field_kind     = field_kind
    self%field_type     = field_type
    self%io_driver      = io_driver

  end function fieldspec_constructor


  !> Getter for unique_id
  !> @param[in]  self  fieldspec_type
  !> @return unique_id

  function get_unique_id(self) result(unique_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    character(str_def) :: unique_id

    unique_id = trim(self%unique_id)

  end function get_unique_id

  !> Getter for mesh_id
  !> @param[in]  self  fieldspec_type
  !> @return mesh_id
  function get_mesh_id(self) result(mesh_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: mesh_id

    mesh_id = self%mesh_id

  end function get_mesh_id

  !> Getter for function_space
  !> @param[in]  self  fieldspec_type
  !> @return function_space
  function get_function_space(self) result(function_space)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: function_space

    function_space = self%function_space

  end function get_function_space

  !> Getter for order
  !> @param[in]  self  fieldspec_type
  !> @return order
  function get_order(self) result(order)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: order

    order = self%order

  end function get_order

  !> Getter for field_kind
  !> @param[in]  self  fieldspec_type
  !> @return field_kind
  function get_kind(self) result(field_kind)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: field_kind

    field_kind = self%field_kind

  end function get_kind

  !> Getter for field_type
  !> @param[in]  self  fieldspec_type
  !> @return field_type
  function get_type(self) result(field_type)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: field_type

    field_type = self%field_type

  end function get_type

  !> Getter for io_driver
  !> @param[in]  self  fieldspec_type
  !> @return io_driver
  function get_io_driver(self) result(io_driver)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: io_driver

    io_driver = self%io_driver

  end function get_io_driver


end module fieldspec_mod
