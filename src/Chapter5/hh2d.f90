!!! File    : hh2d.f90
!!! Author  : Sandeep Koranne (C)
!!! Purpose : Explanation of 2D Householder Reflector
!!!         : given a vector x we want to project it onto |x|e_1 or -|x|e_1
!!!         : whichever is closest. To do this, we calculate v = x-|x|e_1
!!!         : and build H = (I - vv*/|v*v|) then reflect across this by 2
!!!         : vv* is an outer-product, x and v are column vectors of dim n

module Householder2d
  use iso_fortran_env
  implicit none
contains
  function Identity(n)
    integer, intent(in) :: n
    real(kind=real64) :: Identity(n,n)
    integer :: j
    Identity = 0.0_real64
    do j = 1,n
       Identity(j,j) = 1.0_real64
    end do
  end function Identity

  !> Householder Matrix Formation Routine for First Element Reflection
  !>
  !> This subroutine computes the Householder reflection matrix H that zeros out
  !> all elements of the input vector x except the first one. The transformation
  !> is defined as H = I - βvvᵀ where v is the Householder vector and β is a scalar.
  !> The resulting matrix H reflects vector x onto a multiple of the first standard
  !> basis vector e₁, such that H·x = ±||x||·e₁.
  !>
  !> The subroutine handles numerical stability by choosing the sign of the
  !> first element to avoid cancellation errors when computing the Householder vector.
  !>
  !> @param x Input vector of length n (intent in)
  !> @param H Output Householder (intent out)
  !> @param beta Scalar coefficient used in the Householder transformation (intent out)
  !>
  !> @note This implementation includes diagnostic output for debugging purposes
  !> @note The matrix H is guaranteed to be orthogonal (HᵀH = HHᵀ = I)
  !> @note The transformation preserves the Euclidean norm of vectors

  subroutine HH2D(x,H,beta)
    real(kind=real64), intent(in)  :: x(:)
    real(kind=real64), intent(out) :: H(size(x,1),size(x,1))    
    real(kind=real64), intent(out) :: beta
    real(kind=real64) :: norm_x2, norm_x
    real(kind=real64) :: O(size(x,1),size(x,1))
    real(kind=real64) :: v(size(x,1))   
    beta = 0.0_real64
    norm_x = norm2(x)
    !v = x-|x|e1
    v = x
    if (x(1) >= 0.0_real64) then
       beta = -norm_x
    else
       beta = norm_x
    end if
    v(1) = x(1) + beta
    norm_x2 = dot_product(v,v)
    H = Identity(size(x,1))
    O = outer_prod(v)
    beta = 2.0_real64/norm_x2
    H = H - (beta)*O
  end subroutine HH2D
  !> Compute the Householder reflection matrix.
  !>
  !> This function constructs a Householder reflection matrix H such that when applied
  !> to the input vector x, it zeros out all elements from index k+1 to m, while
  !> preserving the first k-1 elements. The transformation is defined as:
  !>   H = I - βvvᵀ
  !> where v is the Householder vector and β = 2 / (vᵀv).
  !>
  !> The algorithm chooses the sign of the norm to avoid numerical cancellation
  !> and returns the identity matrix if the norm of the subvector is below machine
  !> precision (1.0e-12).
  !>
  !> @param x Input vector of length m (1-based indexing), intent(in)
  !> @param k Column index (1-based) from which to start the Householder transformation (intent(in))
  !> @return H Householder reflection matrix of size m×m (result)
  !>
  !> @pre x must be a real vector of at least k elements
  !> @pre k must be between 1 and> @post H is an orthogonal matrix (H^T * H = I)
  !> @post H * x = (x₁, x₂, ..., xₖ, 0, 0, ..., 0)ᵀ
  !>
  !> @example
  !> For x = [1.0, 2.0, 3.0, 4.0] and k = 2:
  !>   - The subvector [2.0, 3.0, 4.0] is reflected onto a multiple of e₂
  !>   - Resulting H will zero out elements 3 and 4 of x
  !>
  !> @see Householder transformation, QR decomposition, orthogonal matrices
  !>
  !> @author Sandeep Koranne
  !> @date 11/13/2025
  !> @version 1.0
  !> @since 2025
  !>
  !> @note This implementation follows standard numerical linear algebra practices.
  !>       It is used in QR decomposition algorithms to introduce zeros below the diagonal.

  function HHReflector(x, k) result(H)
    implicit none
    integer, intent(in) :: k
    real(kind=real64), intent(in) :: x(:)
    real(kind=real64) :: H(size(x,1), size(x,1))
    real(kind=real64) :: v(size(x,1)), beta, norm_x
    integer :: i, j,m
    m = size(x,1)
    ! Initialize identity matrix
    H = Identity(m)
    ! Compute Householder vector for the subvector starting at position k
    ! This creates a reflector that zeros out elements k+1 to m
    norm_x = sqrt(sum(x(k:m)**2))
    if (abs(norm_x) < 1.0e-12_real64) then
       return
    end if
    if (x(k) >= 0.0_real64) then
       beta = -norm_x
    else
       beta = norm_x
    end if
    v = 0.0_real64
    v(k) = x(k) + beta
    v(k+1:) = x(k+1:)
    beta = 2.0_real64 / sum(v(k:m)**2)
    do i = k, m
       do j = k, m
          H(i, j) = H(i, j) - beta * v(i) * v(j)
       end do
    end do
  end function HHReflector

  function outer_prod(v) result(O)
    real(kind=real64),intent(in) :: v(:)
    real(kind=real64) :: O(size(v,1),size(v,1))
    integer :: i, j
    do i = 1, size(v,1)
       do j = 1, size(v,1)
          O(i,j) = v(i) * v(j)
       end do
    end do
  end function outer_prod
end module Householder2d


program main
  use Householder2d
  use iso_fortran_env
  implicit none
  integer, parameter :: n = 4
  real(kind=real64), allocatable :: x(:)
  real(kind=real64), allocatable :: y(:,:)
  real(kind=real64) :: beta
  allocate(x(n), y(n,n))
  !x = [3.0_real64, 4.0_real64]
  x = [3.0_real64,2.0_real64,1.0_real64,3.0_real64]
  call HH2D( x, y, beta )
  write (*,*) '|x| = ', norm2(x)
  write (*,'(4F6.2)') x
  write (*,'(4(4F6.2,/))') transpose(y)  
  write (*,*) ' beta = ', beta
  write (*,'(4F8.5)') matmul(y,x)
  y = HHReflector(x,2)
  write (*,'(4F10.5)') matmul(y,x)  
  deallocate(x,y)
end program main
