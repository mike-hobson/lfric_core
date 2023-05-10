!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   Defines a container object (namelist_item_type) to hold a concrete
!>          key-value pair object via their allocatable abstract type
!>          (keyvalue_type).
!> @details This container is required so that arrays of key-value pair objects
!>          can be created. Concrete objects of different data-types
!>          can be contained in a single array of the container object
!>          (namelist_item_type). The following concrete keyvalue types
!>          are supported in this module.
!>
!>          - i32_keyvalue_type,   i32_arr_keyvalue_type
!>          - i64_keyvalue_type,   i64_arr_keyvalue_type
!>          - r32_keyvalue_type,   r32_arr_keyvalue_type
!>          - r64_keyvalue_type,   r64_arr_keyvalue_type
!>          - logic_keyvalue_type, logic_arr_keyvalue_type
!>          - str_keyvalue_type,   str_arr_keyvalue_type
!>
!>          This container object allows concrete keyvalue types to be
!>          initialised and their values retrieved. However, the container
!>          object (namelist_item_type) does not expose the concrete types
!>          and so prevents data being altered after initialisation.
!>
!>          All concrete types are aliased so that the only public methods
!>          are:
!>          - key()
!>          - initialise( key, value )
!>          - value( output_variable )
!----------------------------------------------------------------------------
module namelist_item_mod

  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

  use constants_mod,  only: imdi, rmdi, cmdi, str_def
  use log_mod,        only: log_event, log_scratch_space, LOG_LEVEL_ERROR
  use keyvalue_mod,   only: keyvalue_type,                                  &
                            i32_keyvalue_type,     i64_keyvalue_type,       &
                            i32_arr_keyvalue_type, i64_arr_keyvalue_type,   &
                            r32_keyvalue_type,     r64_keyvalue_type,       &
                            r32_arr_keyvalue_type, r64_arr_keyvalue_type,   &
                            logic_keyvalue_type,   logic_arr_keyvalue_type, &
                            str_keyvalue_type,     str_arr_keyvalue_type
  implicit none

  private

  public :: namelist_item_type

  !=========================================
  ! Namelist item type to containe abstract
  !=========================================
  type :: namelist_item_type

    private

    class(keyvalue_type), allocatable :: keyvalue_pair

  contains
    procedure, private :: init_i32
    procedure, private :: init_i64
    procedure, private :: init_i32_arr
    procedure, private :: init_i64_arr
    procedure, private :: init_r32
    procedure, private :: init_r64
    procedure, private :: init_r32_arr
    procedure, private :: init_r64_arr
    procedure, private :: init_logic
    procedure, private :: init_logic_arr
    procedure, private :: init_str
    procedure, private :: init_str_arr

    procedure, private :: value_i32
    procedure, private :: value_i64
    procedure, private :: value_i32_arr
    procedure, private :: value_i64_arr
    procedure, private :: value_r32
    procedure, private :: value_r64
    procedure, private :: value_r32_arr
    procedure, private :: value_r64_arr
    procedure, private :: value_logic
    procedure, private :: value_logic_arr
    procedure, private :: value_str
    procedure, private :: value_str_arr

    procedure :: get_key

    generic   :: get_value  => value_i32,     value_i64,       &
                               value_i32_arr, value_i64_arr,   &
                               value_r32,     value_r64,       &
                               value_r32_arr, value_r64_arr,   &
                               value_logic,   value_logic_arr, &
                               value_str,     value_str_arr
    generic   :: initialise => init_i32,      init_i64,        &
                               init_i32_arr,  init_i64_arr,    &
                               init_r32,      init_r64,        &
                               init_r32_arr,  init_r64_arr,    &
                               init_logic,    init_logic_arr,  &
                               init_str,      init_str_arr

  end type namelist_item_type

contains

!=========================================
! 2.0 Initialisers
!=========================================
! Initialisers for concrete keyvalue pair types
!
! 2.1 Initialisers for SCALAR VALUES
!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 32-bit integer.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  32-bit integer to assign to the "key"
subroutine init_i32( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int32), intent(in) :: value

  allocate(i32_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (i32_keyvalue_type)
    call clay%initialise( trim(key), value )
  end select

  return
end subroutine init_i32


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 64-bit integer.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  64-bit integer to assign to the "key"
subroutine init_i64( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int64), intent(in) :: value

  allocate(i64_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (i64_keyvalue_type)
    call clay%initialise( trim(key), value )
  end select

  return
end subroutine init_i64


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 32-bit real.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  32-bit real to assign to the "key"
subroutine init_r32( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  real(real32), intent(in) :: value

  allocate(r32_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (r32_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_r32


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 64-bit real.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  64-bit real to assign to the "key"
subroutine init_r64( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  real(real64), intent(in) :: value

  allocate(r64_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (r64_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_r64


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for logical.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  Logical value to assign to the "key"
subroutine init_logic( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  logical,      intent(in) :: value

  allocate(logic_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (logic_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_logic


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for string value.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  String value to assign to the "key"
subroutine init_str( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  character(*), intent(in) :: value

  allocate(str_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (str_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_str


!=========================================
! 2.2 Initialisers for ARRAY VALUES
!=========================================

!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 32-bit integer array.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  32-bit integer array to assign to the "key"
subroutine init_i32_arr( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int32), intent(in) :: value(:)

  allocate(i32_arr_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (i32_arr_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_i32_arr


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 64-bit integer array.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  64-bit integer array to assign to the "key"
subroutine init_i64_arr( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int64), intent(in) :: value(:)

  allocate(i64_arr_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (i64_arr_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_i64_arr


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 32-bit real array.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  32-bit real array to assign to the "key"
subroutine init_r32_arr( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  real(real32), intent(in) :: value(:)

  allocate(r32_arr_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (r32_arr_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_r32_arr


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for 64-bit real array.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  64-bit real array to assign to the "key"
subroutine init_r64_arr( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  real(real64), intent(in) :: value(:)

  allocate(r64_arr_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (r64_arr_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_r64_arr


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for logical array.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  Logical array to assign to the "key"
subroutine init_logic_arr( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  logical,      intent(in) :: value(:)

  allocate(logic_arr_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (logic_arr_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_logic_arr


!===========================================================
!> @brief Initialises abstract key-value pair to concrete
!>        keyvalue object for string array.
!> @param[in] key    String variable to use as the "key".
!> @param[in] value  String array to assign to the "key"
subroutine init_str_arr( self, key, value )

  implicit none

  class(namelist_item_type), intent(inout) :: self

  character(*), intent(in) :: key
  character(*), intent(in) :: value(:)

  allocate(str_arr_keyvalue_type :: self%keyvalue_pair)

  select type( clay => self%keyvalue_pair )
  class is (str_arr_keyvalue_type)
    call clay%initialise( key, value )
  end select

  return
end subroutine init_str_arr



!=========================================
! 3.0 Getters
!=========================================
! Retrieve data associated with keys for
! for concrete keyvalue pair types
!
!===========================================================
! 3.1 Querys the key string for the object
!===========================================================
!> @brief Returns the key string value used by the object.
!> @return key    String variable to use as the "key". This is
!>                the string set on initialising the abstract
!>                object to a concrete.
function get_key( self ) result( key )

  implicit none

  class(namelist_item_type), intent(in) :: self

  character(:), allocatable :: key

  if (allocated(self%keyvalue_pair)) then
    key = self%keyvalue_pair%get_key()
  else
    allocate( key, source=trim(cmdi) )
  end if

  return
end function get_key


!===========================================================
! 3.2 Getters for SCALAR VALUES
!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  32-bit integer value held in object.
subroutine value_i32( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  integer(int32),            intent(out) :: value

  select type( clay => self%keyvalue_pair )
  class is (i32_keyvalue_type)
    value = clay%value
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected i32_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_i32


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  64-bit integer value held in object.
subroutine value_i64( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  integer(int64),            intent(out) :: value

  select type( clay => self%keyvalue_pair )
  class is (i64_keyvalue_type)
    value = clay%value
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected i64_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_i64


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  32-bit real value held in object.
subroutine value_r32( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  real(real32),              intent(out) :: value

  select type (clay => self%keyvalue_pair)
  class is (r32_keyvalue_type)
    value = clay%value
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected r32_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_r32


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  64-bit real value held in object.
subroutine value_r64( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  real(real64),              intent(out) :: value

  select type (clay => self%keyvalue_pair)
  class is (r64_keyvalue_type)
    value = clay%value
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected r64_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_r64


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  Logical value held in object.
subroutine value_logic( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  logical,                   intent(out) :: value

  select type (clay => self%keyvalue_pair)
  class is (logic_keyvalue_type)
    value = clay%value
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected logic_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_logic


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  String value held in object.
subroutine value_str( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  character(*),              intent(out) :: value

  select type( clay => self%keyvalue_pair )
  class is (str_keyvalue_type)
    value = trim(clay%value)
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected str_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_str


!===========================================================
! 3.3 Getters for ARRAY VALUES
!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  32-bit integer array value held in object.
subroutine value_i32_arr( self, value )

  implicit none

  class(namelist_item_type),   intent(in)  :: self
  integer(int32), allocatable, intent(out) :: value(:)

  select type( clay => self%keyvalue_pair )
  class is (i32_arr_keyvalue_type)
    allocate(value, source=clay%value)
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected i32_arr_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_i32_arr


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  64-bit integer array value held in object.
subroutine value_i64_arr( self, value )

  implicit none

  class(namelist_item_type),   intent(in)  :: self
  integer(int64), allocatable, intent(out) :: value(:)

  select type( clay => self%keyvalue_pair )
  class is (i64_arr_keyvalue_type)
    allocate(value, source=clay%value)
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected i64_arr_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_i64_arr


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  32-bit real array value held in object.
subroutine value_r32_arr( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  real(real32), allocatable, intent(out) :: value(:)

  select type( clay => self%keyvalue_pair )
  class is (r32_arr_keyvalue_type)
    allocate( value(size(clay%value)) )
    value(:) = clay%value(:)
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected r32_arr_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_r32_arr


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  64-bit real array value held in object.
subroutine value_r64_arr( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  real(real64), allocatable, intent(out) :: value(:)

  select type( clay => self%keyvalue_pair )
  class is (r64_arr_keyvalue_type)
    allocate( value( size(clay%value)) )
    value(:) = clay%value(:)
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected r64_arr_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_r64_arr


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  Logical array value held in object.
subroutine value_logic_arr( self, value )

  implicit none

  class(namelist_item_type), intent(in)  :: self
  logical, allocatable,      intent(out) :: value(:)

  select type( clay => self%keyvalue_pair )
  class is (logic_arr_keyvalue_type)
    allocate( value( size(clay%value)) )
    value(:) = clay%value(:)
  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected logic_arr_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_logic_arr


!===========================================================
!> @brief Extracts the value held by the concrete type.
!> @param[out] value  String array value held in object.
subroutine value_str_arr( self, value )

  implicit none

  class(namelist_item_type), intent(in) :: self

  !> @todo This was applied with #3547. This would have been
  !>       similar to the scalar string: i.e.
  !>
  !>       character(*), allocatable, intent(out) :: value(:)
  !>
  !>       However, the revision of the Intel compiler on the XC40
  !>       produced unexpected behaviour so the length has been
  !>       limited to str_def. This should be revisited when the
  !>       XC40 compilers are later than 17.0.0.098/5.
  character(str_def), allocatable, intent(out) :: value(:)

  integer :: arr_len
  integer :: i

  select type( clay => self%keyvalue_pair )
  class is (str_arr_keyvalue_type)
    arr_len = size(clay%value)
    allocate( value(arr_len) )
    do i=1, arr_len
      value(i) = trim(clay%value(i))
    end do

  class default
    write( log_scratch_space, '(A)' ) &
        'Object is not the expected str_arr_keyvalue_type.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  return
end subroutine value_str_arr

end module namelist_item_mod
