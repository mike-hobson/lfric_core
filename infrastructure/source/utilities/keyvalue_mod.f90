!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   Defines key-value pair objects.
!> @details An abstract key-value pair object is defined (keyvalue_type) along
!>          with the following concrete types.
!>
!> * keyvalue_type (abstract)
!>     - i32_keyvalue_type,       i64_keyvalue_type
!>     - r32_keyvalue_type,       r64_keyvalue_type
!>     - logic_keyvalue_type,     str_keyvalue_type
!>     - i32_arr_keyvalue_type,   i64_arr_keyvalue_type
!>     - r32_arr_keyvalue_type,   r64_arr_keyvalue_type
!>     - logic_arr_keyvalue_type, str_arr_keyvalue_type
!>
!> Concrete types can only be initialised once, however the value is directly
!> accessible without a getter/setter.
!----------------------------------------------------------------------------
module keyvalue_mod

  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

  use constants_mod, only: imdi, cmdi, str_longlong, str_def
  use log_mod,       only: log_event, log_scratch_space, &
                           LOG_LEVEL_INFO, LOG_LEVEL_ERROR

  implicit none

  private

  public :: keyvalue_type
  public :: i32_keyvalue_type,     i64_keyvalue_type
  public :: i32_arr_keyvalue_type, i64_arr_keyvalue_type
  public :: r32_keyvalue_type,     r64_keyvalue_type
  public :: r32_arr_keyvalue_type, r64_arr_keyvalue_type
  public :: logic_keyvalue_type,   logic_arr_keyvalue_type
  public :: str_keyvalue_type,     str_arr_keyvalue_type

  !=========================================
  ! Key-Value pair abstract type
  !=========================================
  type, abstract :: keyvalue_type

    private
    character(:), allocatable :: key

  contains
    procedure :: keyvalue_initialise
    procedure :: get_key => get_keyvalue_key
  end type keyvalue_type


  !===============================================
  ! Concrete Types of Abstract: keyvalue_type
  !===============================================

  !===============================================
  ! i32_keyvalue_type: 32-bit integer key-value pair
  !===============================================
  type, extends(keyvalue_type) :: i32_keyvalue_type
    integer(int32) :: value = -huge(0_int32)
  contains
    procedure :: initialise => init_i32_keyvalue
  end type i32_keyvalue_type

  !===============================================
  ! i64_keyvalue_type: 64-bit integer key-value pair
  !===============================================
  type, extends(keyvalue_type) :: i64_keyvalue_type
    integer(int64) :: value = -huge(0_int64)
  contains
    procedure :: initialise => init_i64_keyvalue
  end type i64_keyvalue_type

  !===============================================
  ! r32_keyvalue_type: 32-bit real key-value pair
  !===============================================
  type, extends(keyvalue_type) :: r32_keyvalue_type
    real(real32) :: value = -huge(0.0_real32)
  contains
    procedure :: initialise => init_r32_keyvalue
  end type r32_keyvalue_type

  !==============================================
  ! r64_keyvalue_type: 64-bit real key-value pair
  !==============================================
  type, extends(keyvalue_type) :: r64_keyvalue_type
    real(real64) :: value = -huge(0.0_real64)
  contains
    procedure :: initialise => init_r64_keyvalue
  end type r64_keyvalue_type

  !==============================================
  ! logic_keyvalue_type: Logical key-value pair
  !==============================================
  type, extends(keyvalue_type) :: logic_keyvalue_type
    logical :: value = .false.
  contains
    procedure :: initialise => init_logic_keyvalue
  end type logic_keyvalue_type

  !==============================================
  ! str_keyvalue_type: string key-value pair
  !==============================================
  type, extends(keyvalue_type) :: str_keyvalue_type
    character(str_longlong) :: value
  contains
    procedure :: initialise => init_str_keyvalue
  end type str_keyvalue_type

  !===============================================
  ! i32_keyvalue_type: 32-bit integer key-value pair
  !===============================================
  type, extends(keyvalue_type) :: i32_arr_keyvalue_type
    integer(int32), allocatable :: value(:)
  contains
    procedure :: initialise => init_i32_arr_keyvalue
  end type i32_arr_keyvalue_type

  !===============================================
  ! i64_keyvalue_type: 64-bit integer key-value pair
  !===============================================
  type, extends(keyvalue_type) :: i64_arr_keyvalue_type
    integer(int64), allocatable :: value(:)
  contains
    procedure :: initialise => init_i64_arr_keyvalue
  end type i64_arr_keyvalue_type

  !===============================================
  ! i32_keyvalue_type: 32-bit integer key-value pair
  !===============================================
  type, extends(keyvalue_type) :: r32_arr_keyvalue_type
    real(real32), allocatable :: value(:)
  contains
    procedure :: initialise => init_r32_arr_keyvalue
  end type r32_arr_keyvalue_type

  !===============================================
  ! i64_keyvalue_type: 64-bit integer key-value pair
  !===============================================
  type, extends(keyvalue_type) :: r64_arr_keyvalue_type
    real(real64), allocatable :: value(:)
  contains
    procedure :: initialise => init_r64_arr_keyvalue
  end type r64_arr_keyvalue_type

  !==============================================
  ! logic_arr_keyvalue_type: Logical key-value pair
  !==============================================
  type, extends(keyvalue_type) :: logic_arr_keyvalue_type
    logical, allocatable :: value(:)
  contains
    procedure :: initialise => init_logic_arr_keyvalue
  end type logic_arr_keyvalue_type

  !==============================================
  ! str_arr_keyvalue_type: string array key-value pair
  !==============================================
  type, extends(keyvalue_type) :: str_arr_keyvalue_type
    !> @todo This was applied with #3547. This would have been
    !>       similar to the scalar string: i.e.
    !>
    !>       character(str_longlong), allocatable :: value(:)
    !>
    !>       However, the revision of the Intel compiler on the XC40
    !>       produced unexpected behaviour so the length has been
    !>       limited to str_def. This should be revisited when the
    !>       XC40 compilers are later than 17.0.0.098/5.
    character(str_def), allocatable :: value(:)
  contains
    procedure :: initialise => init_str_arr_keyvalue
  end type str_arr_keyvalue_type

contains


!=========================================
! Abstract Initialiser
!=========================================
subroutine keyvalue_initialise( self, key )

  implicit none

  class(keyvalue_type), intent(inout) :: self
  character(*),         intent(in)    :: key

  if ( allocated(self%key) ) then
    write( log_scratch_space,'(A)' ) &
        'Type already initialised as ' // trim(self%key)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  else
    allocate( self%key, source=trim(key) )
  end if

  return
end subroutine keyvalue_initialise


! @brief Returns key string for the pair
function get_keyvalue_key( self ) result( key )

  implicit none

  class(keyvalue_type), intent(in) :: self

  character(:), allocatable :: key

  allocate( key, source=trim(self%key) )

  return
end function get_keyvalue_key


!=================================================
! Concrete initialisers
! These set the key/value of the concrete types
!=================================================

! 32-bit integer key-value pair initialiser
!==============================================================
!> @brief Initialiser for 32-bit integer key-value pair object
!> @param [in] key      String used for key
!> @param [in] value    32-bit integer value
subroutine init_i32_keyvalue( self, key, value )

  implicit none

  class(i32_keyvalue_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int32), intent(in) :: value

  call self%keyvalue_initialise( key )
  self%value = value

  return
end subroutine init_i32_keyvalue




! 64-bit integer key-value pair initialiser
!==============================================================
!> @brief Initialiser for 64-bit integer key-value pair object
!> @param [in] key      String used for key
!> @param [in] value    64-bit integer value
subroutine init_i64_keyvalue( self, key, value )

  implicit none

  class(i64_keyvalue_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int64), intent(in) :: value

  call self%keyvalue_initialise( key )
  self%value = value

  return
end subroutine init_i64_keyvalue


! 32-bit real key-value pair initialiser
!==============================================================
!> @brief Initialiser for 32-bit real key-value pair object
!> @param [in] key    String used for key
!> @param [in] value  32-bit real value
subroutine init_r32_keyvalue( self, key, value )

  implicit none

  class(r32_keyvalue_type), intent(inout) :: self

  character(*), intent(in) :: key
  real(real32), intent(in) :: value

  call self%keyvalue_initialise( key )
  self%value = value

  return
end subroutine init_r32_keyvalue


! 64-bit real key-value pair initialiser
!==============================================================
!> @brief Initialiser for 64-bit real key-value pair object
!> @param [in] key    String used for key
!> @param [in] value  64-bit real value
subroutine init_r64_keyvalue( self, key, value )

  implicit none

  class(r64_keyvalue_type), intent(inout) :: self

  character(*), intent(in) :: key
  real(real64), intent(in) :: value

  call self%keyvalue_initialise( key )
  self%value = value

  return
end subroutine init_r64_keyvalue


! Logical key-value pair initialiser
!==============================================================
!> @brief Initialiser for logical key-value pair object
!> @param [in] key    String used for key
!> @param [in] value  Logical value
subroutine init_logic_keyvalue( self, key, value )

  implicit none

  class(logic_keyvalue_type), intent(inout) :: self

  character(*), intent(in) :: key
  logical,      intent(in) :: value

  call self%keyvalue_initialise( key )
  self%value = value

  return
end subroutine init_logic_keyvalue


! Logical key-value pair initialiser
!==============================================================
!> @brief Initialiser for logical key-value pair object
!> @param [in] key    String used for key
!> @param [in] value  Logical value
subroutine init_str_keyvalue( self, key, value )

  implicit none

  class(str_keyvalue_type), intent(inout) :: self

  character(*), intent(in) :: key
  character(*), intent(in) :: value

  call self%keyvalue_initialise( key )
  self%value = trim(value)

  return
end subroutine init_str_keyvalue


! 32-bit integer key-value pair initialiser
!==============================================================
!> @brief Initialiser for 32-bit integer key-value pair object
!> @param [in] key      String used for key
!> @param [in] value    32-bit integer value
subroutine init_i32_arr_keyvalue( self, key, value )

  implicit none

  class(i32_arr_keyvalue_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int32), intent(in) :: value(:)

  call self%keyvalue_initialise( key )
  allocate( self%value, source=value )

  return
end subroutine init_i32_arr_keyvalue


! 64-bit integer key-value pair initialiser
!==============================================================
!> @brief Initialiser for 64-bit integer key-value pair object
!> @param [in] key      String used for key
!> @param [in] value    64-bit integer value
subroutine init_i64_arr_keyvalue( self, key, value )

  implicit none

  class(i64_arr_keyvalue_type), intent(inout) :: self

  character(*),   intent(in) :: key
  integer(int64), intent(in) :: value(:)

  call self%keyvalue_initialise( key )
  allocate( self%value, source=value )

  return
end subroutine init_i64_arr_keyvalue


! 32-bit integer key-value pair initialiser
!==============================================================
!> @brief Initialiser for 32-bit integer key-value pair object
!> @param [in] key      String used for key
!> @param [in] value    32-bit integer value
subroutine init_r32_arr_keyvalue( self, key, value )

  implicit none

  class(r32_arr_keyvalue_type), intent(inout) :: self

  character(*),   intent(in) :: key
  real(real32), intent(in) :: value(:)

  call self%keyvalue_initialise( key )
  allocate( self%value, source=value )

  return
end subroutine init_r32_arr_keyvalue


! 64-bit integer key-value pair initialiser
!==============================================================
!> @brief Initialiser for 64-bit integer key-value pair object
!> @param [in] key      String used for key
!> @param [in] value    64-bit integer value
subroutine init_r64_arr_keyvalue( self, key, value )

  implicit none

  class(r64_arr_keyvalue_type), intent(inout) :: self

  character(*),   intent(in) :: key
  real(real64), intent(in) :: value(:)

  call self%keyvalue_initialise( key )
  allocate( self%value, source=value )

  return
end subroutine init_r64_arr_keyvalue



subroutine init_logic_arr_keyvalue( self, key, value )

  implicit none

  class(logic_arr_keyvalue_type), intent(inout) :: self

  character(*), intent(in) :: key
  logical,      intent(in) :: value(:)

  call self%keyvalue_initialise( key )
  allocate( self%value(size(value)) )
  self%value(:) = value(:)

  return
end subroutine init_logic_arr_keyvalue

subroutine init_str_arr_keyvalue( self, key, value )

  implicit none

  class(str_arr_keyvalue_type), intent(inout) :: self

  character(*), intent(in) :: key
  character(*), intent(in) :: value(:)

  integer :: arr_len, i

  arr_len  = size(value)

  call self%keyvalue_initialise( key )
  allocate( self%value(arr_len) )
  do i=1, arr_len
    self%value(i) = trim(value(i))
  end do

  return
end subroutine init_str_arr_keyvalue

end module keyvalue_mod
