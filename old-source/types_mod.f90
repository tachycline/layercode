module Types
! double precision
  integer, parameter,public :: r8=selected_real_kind(15,300)
! single precision
  integer, parameter,public :: r4=selected_real_kind(1,1)
! integer types
  integer, parameter,public :: i1=selected_int_kind(1)
  integer, parameter,public :: i2=selected_int_kind(3)
  integer, parameter,public :: i3=selected_int_kind(5)
  integer, parameter,public :: i4=selected_int_kind(10)

end module Types
