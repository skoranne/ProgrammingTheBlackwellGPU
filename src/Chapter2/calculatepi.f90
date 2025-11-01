!!! File   : calculatepi.f90
!!! Author : Sandeep Koranne (C)

program main
  use omp_lib
  implicit none
  integer :: n1, n2, i, n
  integer, allocatable :: seed(:)
  real(8) :: pi, x1, x2

  n1 = 1000*1000

  ! Make space for a deterministic seed
  call random_seed(size=n)
  allocate(seed(n))

  !$omp parallel firstprivate(seed), private(x1, x2, i), shared(n1), reduction(+:n2)
  n2 = 0

  ! Deterministic seed, distinct on each thread.
  seed = [ (i, i=0, size(seed)-1) ] + omp_get_thread_num()
  call random_seed(put=seed)

  !$omp do
  do i = 1, n1
    call random_number(x1)
    call random_number(x2)
    if (x1**2 + x2**2 <= 1.0_8) n2 = n2 + 1
  end do
  !$omp end do

  !$omp end parallel
  pi = 4.0_8*dble(n2)/dble(n1)
  print *, pi
end program main
