!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! File    : pi_agm.f90
!!! Author  : Sandeep Koranne, 2025. Based on L. N. Trefethen TDA pi-agm
!!! Purpose : Calculates Pi: 3.14159265358979323846264338327950288419716939937510
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  use iso_fortran_env
  implicit none
  !integer,parameter :: mydp = real64 !selected_real_kind(16)
  integer,parameter :: mydp = real128 !selected_real_kind(16)
  real(kind=mydp) :: x,y,p,s
  real(mydp),parameter :: PI = 4.0_mydp*atan(1.0_mydp)
  integer :: i
  integer,parameter :: N=20
  y = sqrt(sqrt(2.0_mydp))
  x = (y+1.0_mydp/y)/(2.0_mydp)
  p = 2.0_mydp + sqrt(2.0_mydp)
  print *,'                        3.14159265358979323846264338327950288419716939937510'
  do i=1,N
     write (*,'(I2,A,F20.18,A,F20.18,A,F20.18)') i,'th approximation to π: ', p, ' ',PI,' err= ',abs(p-PI)
     p = p*(1.0_mydp + x)/(1.0_mydp + y)
     s = sqrt(x)
     y = (y*s+1.0_mydp/s)/(1.0_mydp+y)
     x = (s+1.0_mydp/s)/2.0_mydp
  end do
  stop
end program main
