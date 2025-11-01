!!! Example code to show List comprehensions in Fortran
module PythagoreanTriple
  type Triple
     integer :: a,b,c
  end type Triple
contains
  function ConvertToTriple(A,N) result(retval)
    logical,dimension(:), intent(in) :: A
    integer,intent(in)               :: N
    integer,dimension(:),allocatable :: pos
    integer,dimension(:),allocatable :: iota
    integer                          :: i,alloc_stat
    type(Triple),dimension(:),allocatable :: retval
    iota=[ (i, i=1,size(A)) ]
    pos = pack( iota, A )
    allocate(retval(count(A)),stat=alloc_stat)
    if( alloc_stat /= 0 ) stop "Not able to allocate memory"
    do concurrent (i=1:count(A))
       associate (a=>retval(i)%a, b=>retval(i)%b, c=>retval(i)%c)
         c = ceiling(pos(i)/real(N**2))
         a = (pos(i)-(c-1)*N**2)
         b = ceiling(a/real(N))
         a = mod(a,N)
       end associate
    end do    
  end function ConvertToTriple
  pure logical function CheckTriple(a,b,c) result(retval)
    integer,intent(in) :: a,b,c
    integer            :: s
    s = a+b+c
    retval =  a/real(s-a)+b/real(s-b)+c/real(s-c) == 4
  end function CheckTriple
  pure logical function CheckPTriple(a,b,c) result(retval)
    integer,intent(in) :: a,b,c
    retval = c**2 == a**2+b**2
  end function CheckPTriple
end module PythagoreanTriple

program GeneratePythagoreanTriple
  use iso_fortran_env
  use PythagoreanTriple
  use omp_lib
  implicit none
  integer :: a,b,c,total_found
  integer :: time_1, time_2, delta_t, countrate, countmax
  real(real64) :: secs
  !500 => 772
  !5000 => 11362 in real    0m11.083s nvfortran -mp=multicore -stdpar=multicore -O3
  !5000 => 11362 in real    0m3.801s -parallel -fopenmp -O3
  !5000 => real    1m17.624s

  !50000  Number of triples found =      163010
  !IntelLLVM 1h
  !50000 nvfortran 10m OpenMP 
  !integer, parameter :: N = 50000
  integer, parameter :: N = 50000
  logical,dimension(:),allocatable :: val
  type(Triple),dimension(:),allocatable :: answer
  integer :: alloc_stat
  call system_clock(count_max=countmax, count_rate=countrate)
  call system_clock(time_1)
  
  total_found = 0
  !do concurrent(c=1:N) reduce(+:total_found)
  !   do concurrent(b=1:N) reduce(+:total_found)
  !      do concurrent(a=1:N) reduce(+:total_found)
  !! $ qmp parallel do shared(total_found) default(none)
  !$omp target teams map(tofrom: total_found)
  !$omp distribute parallel do SIMD private(a,b,c) reduction(+: total_found)
  do c=1,N
     do b=1,N 
        do a=1,N
           !There is a HUGE runtime difference between calling this function vs inline of the code
           !if(CheckPTriple(a,b,c)) total_found = total_found+1
           if( c**2 == a**2+b**2 ) total_found = total_found+1
           !where( CheckPTriple(a,b,c) ) total_found = total_found+1
        end do
     end do
  end do
  !$omp end target teams
  call system_clock(time_2)
  delta_t = time_2-time_1
  secs = real(delta_t)/real(countrate)
  print *, secs,' seconds'
  !write(6,'(1I5,13X,1A,F10.2)')OMP_GET_MAX_THREADS(),'|',secs
  
  !val = [ ((( c**2 == a**2+b**2, a=1,N),b=1,N),c=1,N)]
  !val = [ ((( CheckTriple(a,b,c), a=1,N),b=1,N),c=1,N)]
  print *,'Number of triples found =', total_found
  !answer = ConvertToTriple( val, N )
  !print *,'--------- Pythagorean Triples upto ',N,'---------'
  !write (*,'(3A5)') 'a','b','c'
  !write (*,'(A)') '---------------'
  !write (*,'(3I5)') answer
end program GeneratePythagoreanTriple
  
! Number of triples found =           4
! --------- Pythagorean Triples upto           10 ---------
!    a    b    c
!---------------
!    4    3    5
!    3    4    5
!    8    6   10
!    6    8   10
