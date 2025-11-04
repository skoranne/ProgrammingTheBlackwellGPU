!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! File    : test_intrinsics.f90
!!! Author  : Sandeep Koranne, 2025. 
!!! Purpose : Become familiar with builtin intrinsics of modern Fortran
!!! UpperTriangular: if(row > col) B(row,col) = 0
!!! LowerTriangular: if(row < col) B(row,col) = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!include 'mkl_vsl.f90'
module TestIntrinsics
  use iso_fortran_env
  use Utilities
  implicit none
contains
  subroutine ConvertToBanded(B,m1,m2)
    real,dimension(:,:),intent(inout) :: B
    integer,intent(in) :: m1,m2
    integer :: row,col
    if( size(B,dim=1) /= size(B,dim=2) ) error stop "B should be a square matrix"
    do col=1,size(B,dim=2)
       do row=1,size(B,dim=1)
          if( (row > col+m1) .or. (col > row+m2) ) then
             B(row,col) = 0
          end if
       end do
    end do
    return
  end subroutine ConvertToBanded
  function EncodeBanded(B,m1,m2) result(retval)
    real,dimension(:,:),intent(in) :: B
    integer,intent(in) :: m1,m2
    real,dimension(:,:) :: retval(size(B,dim=1),m1+m2+1)
    integer :: row,col,i
    if( size(B,dim=1) /= size(B,dim=2) ) error stop "B should be a square matrix"    
    retval(:,m1+1)=[(B(i,i),i=1,size(B,dim=1))] !diagonals
    do col=1,m2
       retval(:,m1+1+col) = [(B(i,i+col),i=1,size(B,dim=1))]
    end do
    do col=1,m1
       retval(:,m1+1-col) = eoshift([(B(i+col,i),i=1,size(B,dim=1))],dim=1,shift=-col)
    end do
  end function EncodeBanded
  function EncodeBandedLAPACK(B,m1,m2) result(retval)
    !AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
    real,dimension(:,:),intent(in) :: B
    integer,intent(in) :: m1,m2
    real,dimension(:,:) :: retval(2*m1+m2+1,size(B,dim=2))
    integer :: row,col,i1,i2,k,M,N
    if( size(B,dim=1) /= size(B,dim=2) ) error stop "B should be a square matrix"
    N = size(B,dim=1)
    M = m1+m2+1 !total number of rows
    retval = 0.0
    do col=1,N ! this is from the LINPACK/LAPACK book pp 2.3
       i1 = max(1,col-m2)
       i2 = min(N,col+m1)
       do row=i1,i2
          k = row-col+M
          retval(k,col) = B(row,col)
       end do
    end do
  end function EncodeBandedLAPACK
  
  function BandedMultiply(Z,m1,m2,x) result(retval)
    real,dimension(:,:),intent(in) :: Z
    integer,intent(in) :: m1,m2
    real,dimension(size(Z,dim=1)),intent(in) :: x
    real,dimension(size(Z,dim=1)) :: retval
    real,dimension(:,:),allocatable :: y
    integer :: row,col,i,N
    N = size(Z,dim=1)
    if( N /= size(x) ) error stop "A and x should have same dim"
    y = spread(x,dim=2,ncopies=m1+m2+1)
    y = transpose(eoshift(y,dim=1,shift=ArithmeticSequenceInt(-m1,1,m1+m2+1)))
    retval = [( dot_product(z(i,:),y(:,i)),i=1,N)]
  end function BandedMultiply
  function BandedMultiplySerial(Z,m1,m2,x) result(retval)
    real,dimension(:,:),intent(in) :: Z
    integer,intent(in) :: m1,m2
    real,dimension(size(Z,dim=1)),intent(in) :: x
    real,dimension(size(Z,dim=1)) :: retval
    real,dimension(:,:),allocatable :: y
    integer :: row,col,i,k,N
    N = size(Z,dim=1)
    if( N /= size(x) ) error stop "A and x should have same dim"
    retval = 0.0
    do row=1,N
       k = row-m1-1
       do col=max(1,1-k),min(m1+m2+1,N-k)
          retval(row) = retval(row) + Z(row,col)*x(col+k)
       end do
    end do
  end function BandedMultiplySerial
  subroutine TestDGETRF()
    integer :: m,n,kl,ku,ldab,info,k,i,j
    integer,dimension(:),allocatable :: ipiv
    real(real64),dimension(:,:),allocatable :: ab
    read (*,*) !skip first line
    read (*,*) m,n,kl,ku
    ldab = 2*kl + ku + 1
    allocate(ab(ldab,n))
    allocate(ipiv(max(m,n)))
    k = kl + ku + 1
    read (*, *)((ab(k+i-j,j),j=max(i-kl,1),min(i+ku,n)), i=1, m)
    print *,'LAPACK Banded matrix:',shape(ab)
    do i=1,ldab
       write (*,'(9F8.4)') ab(i,:)
    end do
    print *,'Transpose trick'
    write (*,'(9F8.4)') transpose(ab)
    !     Factorize A
    call dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)

    print *,'TestDGETRF Info=',info
    if( info < 0 ) then
       print *,'Row ',info,' has problem'
       error stop "Something wrong in DGBTRF"
    else if( info >0 ) then
       print *,'Row ',info,' has problem'
       error stop "Something wrong in DGBTRF"
    end if
    print *,'LAPACK factored matrix:',shape(ab)
    write (*,'(9F8.4)') transpose(ab)
    print *,'IPIV=',ipiv
    deallocate(ab)
    deallocate(ipiv)
  end subroutine TestDGETRF
  
end module TestIntrinsics

program main
  use iso_fortran_env
  use TestIntrinsics
  use LAPACK95
  !use mkl_vsl
  implicit none
  integer,parameter :: N = 9
  integer :: i,j,k,row,col
  integer,dimension(N*N) :: A = [(i,i=1,N*N)]
  real(real64),dimension(:,:),allocatable :: B

  integer,dimension(4) :: x = [(i,i=1,4)]
  integer,dimension(:,:),allocatable :: y
  real,dimension(:,:),allocatable :: z
  real,dimension(:,:),allocatable :: b1,b4
  real,dimension(:),allocatable   :: b2,b3,b5,b6
  real,dimension(N) :: x0
  integer :: lb,ub
  integer :: info
  integer,parameter    :: seed = 32125
  integer,dimension(N) :: ipiv
  character(len=20) :: format_str
  character(len=20) :: encoded_str  
  write (format_str,100) N
  100 format('(',I1,'F8.2)')
  print *,'The format_str is:',format_str
  !call TestDGETRF()
  call random_number(x0)
  x0 = ArithmeticSequence(1.0,1.0,N)
  !x is a column vector
  print *,'x='
  write (*,'(1I4)') x
  y = spread(x,dim=1,ncopies=2)  
  print *,'y=spread(x,dim=1,ncopies=2)'
  write (*,'(2I4)') y

  y = spread(x,dim=2,ncopies=3)  
  print *,'y=spread(x,dim=2,ncopies=3)'
  write (*,'(4I4)') y
  call srand(seed)
  B = reshape(A,[N,N])
  call random_number(B)
!#define SAND
#if defined SAND
  do col=1,N
     do row=1,N
        B(row,col) = rand() !row*10+col
     end do
  end do
#endif
  
  write (*,format_str) transpose(B)+0.0
  lb=3;ub=2
  write (encoded_str,100) lb+ub+1  
  !make B into a Banded matrix
  z=B+0.0
  call ConvertToBanded(z,m1=lb,m2=ub)
  print *,'z='
  write (*,format_str) transpose(z)
  b3 = matmul( z, x0 )  !original answer
  !z=transpose(z) !for printing

  !what is the size of bandencoding of z considering there are few NZ
  b1 = EncodeBanded(z,m1=lb,m2=ub)
  b4 = EncodeBandedLAPACK(z,m1=lb,m2=ub)  
  print *,'Encoded banded matrix:'
  write (*,encoded_str) transpose(b1)
  print *,'LAPACK encoded matrix:',shape(b4)
  write (*,format_str) transpose(b4)    
  !deallocate(b1) !is this a memory leak as it is an automatic returned array
  b2 = BandedMultiply(b1,m1=lb,m2=ub,x=x0)
  !write (*,format_str) x0
  print *,'Error=',norm2(b2-b3)
  b5 = BandedMultiplySerial(b1,m1=lb,m2=ub,x=x0)
  print *,'Error Serial=',norm2(b5-b3)
  !lets factor b4 banded matrix into LU using LAPACK
  call dgbtrf(N, N, lb, ub, b4, 2*lb+ub+1, ipiv, info)
  if( info < 0 ) then
     print *,'Row ',info,' has problem'
     error stop "Something wrong in DGBTRF"
  else if( info >0 ) then
     print *,'Row ',info,' has problem'
     !error stop "Something wrong in DGBTRF"
  end if
  print *,'LAPACK factored matrix:'
  write (*,format_str) transpose(b4)    
  print *,'IPIV=',ipiv
end program main

