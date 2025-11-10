!!! File    : diagonal.f90
!!! Author  : Sandeep Koranne (C)
!!! Purpose : Example of pointer usage, and also a delightful example
program main
  use iso_fortran_env
  implicit none
  integer, parameter :: N = 9
  complex(kind=real64), pointer :: storage(:), matrix(:,:), main_diagonal(:)
  real(kind=real64), pointer :: stability(:)
  complex(kind=real64), allocatable :: v_new(:)
  integer            :: i,j
  allocate( storage(N*N) )
  allocate( v_new(N) )
  v_new = cmplx(1.0, 0.0)
  write (*,*) Cnorm2(v_new)
  matrix(1:N, 1:N) => storage
  main_diagonal(1:N) => storage(::N+1) ! this is unique to Fortran
  do i=1,N
     do j=1,N
        matrix(i,j) = cmplx(i+j,i*j)
     end do
  end do
  write (*,*) "Main diagonal:"
  write (*,*) main_diagonal
  stability => main_diagonal%im  
  write (*,*) "Re part of main diagonal: ", norm2(stability)!, norm2(v_new)
  write (*,*) stability  ! this is such an elegant formalism
  deallocate( storage )
contains
  pure real(real64) function Cnorm2(x)
    complex(kind=real64), intent(in) :: x(:)
    Cnorm2 = sqrt( dot_product(x,x) )
  end function Cnorm2
  
end program main
