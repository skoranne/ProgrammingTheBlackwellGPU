!!! File    : interval.f90
!!! Author  : Sandeep Koranne (C) 2025.
!!! Purpose : Modern Fortran programming
module Interval
  implicit none
  public
  type interval
     integer :: lo, hi
   contains
     procedure, pass :: valid, add, mul, write
  end type interval
  interface operator(+)
     module procedure add
  end interface operator(+)
  interface operator(*)
     module procedure mul
  end interface operator(*)  
  !! Please make sure that assignment subroutines are (o,i)
  interface assignment(=)
     module procedure interval_from_integer, integer_from_interval
  end interface assignment(=)
contains
  subroutine integer_from_interval(a,retval)
    type(interval), intent(in) :: a
    integer :: retval
    retval = (a%lo+a%hi)/2
  end subroutine integer_from_interval
  subroutine interval_from_integer(o,i)
    type(interval) :: o
    integer, intent(in) :: i
    o%lo = i; o%hi = i
  end subroutine interval_from_integer
  logical function valid(a)
    class(interval), intent(in) :: a
    valid = a%hi > a%lo !we dont want point intervals
  end function valid
  type(interval) function add(a,b)
    class(interval), intent(in) :: a, b
    add%lo = min( a%lo, b%lo )
    add%hi = max( a%hi, b%hi )
  end function add
  !Not a trivial function to write as the output can be invalid
  type(interval) function mul(a,b)
    class(interval), intent(in) :: a, b
    mul%lo = max( a%lo, b%lo )
    mul%hi = min( a%hi, b%hi )
  end function mul
  subroutine write(a)
    class(interval), intent(in) :: a
    write (*,*) '[lo: ', a%lo, ' hi: ', a%hi, ' ]'
  end subroutine write
end module Interval

program main
  use Interval
  implicit none
  type(interval) :: c
  type(interval) :: a = interval(5,10), b=interval(7,8)
  call a%write()
  call b%write()
  c = add(a,b)
  call c%write()
  c = 3
  call c%write()
  c = a+b
  call c%write()
  c = a*b
  call c%write()
end program main

 ! [lo:             5  hi:            10  ]
 ! [lo:             7  hi:             8  ]
 ! [lo:             5  hi:            10  ]
 ! [lo:             3  hi:             3  ]
 ! [lo:             5  hi:            10  ]
 ! [lo:             7  hi:             8  ]
