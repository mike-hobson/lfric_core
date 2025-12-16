!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Provides access to the state of the Conway's Game of Life
!> @details  Extends the abstract external field class to provide access to
!>           either a) set the state, or b) write the current state to the log

module conway_state_mod

  use abstract_external_field_mod, only: abstract_external_field_type
  use constants_mod,            only: i_def, r_def, l_def, str_def, imdi
  use integer_field_mod,        only: integer_field_type, &
                                      integer_field_proxy_type
  use field_collection_mod,     only: field_collection_type
  use function_space_mod,       only: function_space_type
  use log_mod,                  only: log_event,       &
                                      LOG_LEVEL_INFO,  &
                                      log_scratch_space
  use mesh_mod,                 only: mesh_type
  use reference_element_mod,    only: W, S, E, N

  implicit none

  private

type, extends(abstract_external_field_type), public :: conway_state_type
  private
  ! The cell id at which to insert the item into the field
  integer(i_def) :: cell
  ! A pointer to an inteeger field (the fieldf in the abstract is a real field)
  type(integer_field_type), pointer :: lfric_integer_field
contains
  private
  !> Initialises the object
  procedure, public :: initialise
  !> Dummy that copies data from the LFRic field, but not neede here
  procedure, public :: copy_from_lfric
  !> Copies state data into the LFRic field
  procedure, public :: copy_to_lfric
  !> Manually tidies up
  procedure, public :: clear
  !> Tidies up on destruction
  final             :: finalise
end type  conway_state_type

  contains

  !> @brief Initialises the external field used for coupling
  !> @param [in] lfric_field_ptr Pointer to an lfric field
  !> @param [in] sorting_index   Index to sort data for coupling
  !
  subroutine initialise( self, lfric_integer_field_ptr, cell )
  implicit none

  class(conway_state_type),      intent(inout) :: self
  type(integer_field_type), pointer, intent(in)    :: lfric_integer_field_ptr
  integer(i_def),                    intent(in)    :: cell

!  This external field is based on an integer field - the abstract
!  contains a real field - so there is no need to 
!  call self%abstract_external_field_initialiser(lfric_field_ptr)
  self%lfric_integer_field => lfric_integer_field_ptr

  self%cell = cell

  end subroutine initialise


  !>@brief Accesses the state and writes it out to the log
  !>@param return_code Optional return code from the copy_from procedure
  !
  subroutine copy_from_lfric(self, return_code)

  implicit none

  class(conway_state_type), intent(inout) :: self
  integer(i_def), optional,        intent(out)   :: return_code

  write(log_scratch_space, '(A)' ) &
    "copy_from_lfric for the conway_state object is not written, yet"
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine copy_from_lfric


  !> @brief copies data into the state field
  !> @param return_code The return code from the copy_to procedure
  !
  subroutine copy_to_lfric( self, return_code )
  implicit none

  class(conway_state_type), intent(inout) :: self
  integer(i_def), optional,     intent(out)   :: return_code

  type( integer_field_proxy_type )       :: current_state_proxy
  integer(i_def)                         :: cell
  integer(i_def), pointer                :: cell_dofmap(:)
  type(function_space_type), pointer     :: fs
  type(mesh_type), pointer               :: mesh

  current_state_proxy = self%lfric_integer_field%get_proxy()
  fs => self%lfric_integer_field%get_function_space()
  mesh => self%lfric_integer_field%get_mesh()

  ! Insert a glider
  cell = self%cell
  cell_dofmap => fs%get_cell_dofmap(cell)
  current_state_proxy%data(cell_dofmap(1)) = 1
  cell = mesh%get_cell_next(E, cell)
  current_state_proxy%data(cell_dofmap(1)) = 1
  cell = mesh%get_cell_next(E, cell)
  current_state_proxy%data(cell_dofmap(1)) = 1
  cell = mesh%get_cell_next(S, cell)
  current_state_proxy%data(cell_dofmap(1)) = 1
  cell = mesh%get_cell_next(S, cell)
  cell = mesh%get_cell_next(W, cell)
  current_state_proxy%data(cell_dofmap(1)) = 1
  call log_event( 'conway: glider added to state', LOG_LEVEL_INFO )

  end subroutine copy_to_lfric


  ! Finaliser/Clear
  !
  !> @brief Deallocates the memory associated with the object.
  subroutine clear(self)
  implicit none
  class(conway_state_type), intent(inout) :: self

  nullify(self%lfric_integer_field)

  end subroutine clear

  subroutine finalise(self)
  implicit none
  type(conway_state_type), intent(inout) :: self

  call self%clear()

  end subroutine finalise

end module conway_state_mod
