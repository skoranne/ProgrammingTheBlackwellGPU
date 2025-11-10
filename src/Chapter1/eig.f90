!!! File    : eig.f90
!!! Author  : Sandeep Koranne (C)
!!! Purpose : Demonstrate matmul offloading for Rayleigh iteration

program eigenvalue_calculator
  use iso_fortran_env, only: real64
  implicit none
  integer :: i
  complex(real64), allocatable :: matrix(:,:)
  complex(real64), allocatable :: eigenvals(:), eigs_ans(:)
  complex(real64), allocatable :: gershgorin_bounds(:)
  ! Initialize a 3x3 complex matrix for testing
  !   3×3 Matrix{ComplexF64}:
  !  1.0+0.0im  0.0+1.0im  1.0+0.5im
  !  0.0+0.5im  2.0+0.0im  1.0+0.0im
  !  1.0+0.0im  0.0+1.0im  3.0+0.0im

  ! julia> ei=eigvals(A)
  ! 3-element Vector{ComplexF64}:
  !  0.7867103000148132 - 0.023986000103288857im
  !  1.6305038358200563 - 0.7339786518673036im
  !  3.5827858641651327 + 0.7579646519705924im
  ! ZGEEV Example Program Results
  ! 
  ! Eigenvalues
  !    ( 0.7867,-0.0240) ( 1.6305,-0.7340) ( 3.5828, 0.7580)

  allocate(matrix(10,10))
  ! matrix = reshape([ &
  !      (1.0d0, 0.0d0), (0.0d0, 1.0d0), (1.0d0, 0.5d0), &
  !      (0.0d0, 0.5d0), (2.0d0, 0.0d0), (1.0d0, 0.0d0), &
  !      (1.0d0, 0.0d0), (0.0d0, 1.0d0), (3.0d0, 0.0d0) &
  !      ], [3,3])

  ! Calculate eigenvalues using Rayleigh power iteration
  call generate_hermitian_matrix(matrix,10)
  eigenvals = calculate_eigenvalues(matrix)
  eigs_ans = eigs(matrix)
  ! Calculate Gershgorin bounds
  gershgorin_bounds = calculate_gershgorin_bounds(matrix)
  write (*,*) 'GEigenvalues: ', maxval(abs(eigenvals))
  write(*,*) 'Gershgorin bounds:'
  do i = 1, size(gershgorin_bounds)/2
     write(*,'(1x,2f8.4,1x,"to",1x,2f8.4)') &
          gershgorin_bounds(2*i-1)%re, gershgorin_bounds(2*i-1)%im, &
          gershgorin_bounds(2*i)%re, gershgorin_bounds(2*i)%im
  end do

  ! Clean up
  deallocate(matrix, eigenvals, eigs_ans, gershgorin_bounds)

contains
  !! Hermitian matrices result in Real Eigenvalues.
  subroutine generate_hermitian_matrix(matrix, size)
    use iso_fortran_env
    implicit none
    integer, intent(in) :: size
    complex(real64), intent(out) :: matrix(size, size)
    integer :: i, j
    real(real64) :: rand_valA, rand_valB
    matrix = (0.0d0, 0.0d0)
    ! Fill upper triangular part with random complex numbers
    do i = 1, size
       do j = i, size
          ! Generate random real and imaginary parts
          call random_number(rand_valA)
          call random_number(rand_valB)
          write (*,*) 'Rng: ', rand_valA, ' ', rand_valB
          matrix(i, j) = cmplx(rand_valA, rand_valB)
          ! For Hermitian matrix, matrix(i,j) = conjg(matrix(j,i))
          if (i /= j) then
             matrix(j, i) = conjg(matrix(i, j))
          end if
       end do
    end do
    ! Make diagonal real (optional - for better numerical properties)
    do i = 1, size
       matrix(i, i) = cmplx(real(matrix(i, i)), 0.0d0)
    end do
  end subroutine generate_hermitian_matrix

  function eigs(M) result(eigenvals)
    implicit none    
    complex(real64), intent(in) :: M(:,:)
    complex(real64), allocatable :: eigenvals(:)
    integer, parameter :: nb = 64, nin = 5, nout = 6
    integer :: i, ifail, info, lda, ldvr, lwork, n
    complex(kind=real64), allocatable :: vr(:,:), work(:)
    complex(kind=real64) :: dummy(1,1)
    real(kind=real64), allocatable :: rwork(:)
    n = size(M,1)
    allocate(vr(n,n), eigenvals(n),rwork(2*n))
    lwork = -1
    call zgeev('No left vectors', 'Vectors (right)', n, M, n, eigenvals, dummy, 1, vr, n, dummy, lwork, rwork, info)
    lwork = max((nb+1)*n, nint(real(dummy(1,1))))
    allocate (work(lwork))

    call zgeev('No left vectors', 'Vectors (right)', n, M, n, eigenvals, dummy, 1, vr, n, work, lwork, rwork, info)
    if( info == 0 ) then
       write (*,*) 'ZGEEV Eigenvalues: ', maxval(abs(eigenvals))
       write (*,*) eigenvals
    end if
  end function eigs

  ! Pure function to calculate eigenvalues using Rayleigh power iteration
  function calculate_eigenvalues(M) result(eigenvals)
    implicit none
    complex(real64), intent(in) :: M(:,:)
    complex(real64), allocatable :: eigenvals(:)

    integer :: n, i, j, iter
    integer, parameter :: max_iter = 1000
    real(real64), parameter :: tolerance = 1e-3_real64
    complex(real64), allocatable :: v(:), v_new(:), temp_vec(:)
    complex(real64) :: lambda, lambda_old, max_eigenval
    complex(real64), allocatable :: matrix(:,:)
    real(real64) :: current_norm
    logical :: converged

    n = size(M, 1)
    allocate(eigenvals(n))
    allocate(v(n))
    allocate(v_new(n))
    allocate(temp_vec(n))
    allocate(matrix(n,n))
    matrix = M
    ! Calculate all eigenvalues using power iteration
    do i = 1, 1 !deflate is not working
       ! Initialize for this eigenvalue calculation
       v = cmplx(1.0d0, 0.0d0) / sqrt(real(n, real64))
       lambda = cmplx(0.0d0, 0.0d0)
       converged = .false.

       ! Rayleigh quotient iteration for this eigenvalue
       do iter = 1, max_iter
          v_new = matmul(matrix,v) ! 3x3 x 3x1 := 3x1
          lambda_old = lambda
          lambda = dot_product(v, v_new)! / dot_product(v, v)
          ! Check convergence λ
          !write (*,'(A10,I4,1X,A3,F6.3,A2,F6.3)') 'Iteration: ',iter, 'λ ', lambda%re,'+i', lambda%im
          if (abs(lambda - lambda_old) < tolerance) then
             converged = .true.
             exit
          end if
          current_norm = Cnorm2(v_new)
          v = v_new / Cnorm2(v_new) 
       end do

       eigenvals(i) = lambda

       ! Deflate the matrix for next eigenvalue
       ! This is a simplified approach - full deflation would be more complex
       if (i < n) then
          matrix = matrix - lambda * identity_matrix(n)
       end if
    end do

    ! Clean up
    deallocate(v, v_new, temp_vec, matrix)
  end function calculate_eigenvalues

  ! Pure function to calculate Gershgorin circle bounds
  pure function calculate_gershgorin_bounds(M) result(bounds)
    use gershgorin
    implicit none
    complex(real64), intent(in)  :: M(:,:)
    complex(real64), allocatable :: bounds(:)
    type(gershgorin_disk),allocatable :: disks(:)
    integer :: i
    disks = ComputeGershgorinDisk(M)
    allocate(bounds(2*size(M,1)))
    do i = 1, size(M,1)
       bounds(2*i-1) = disks(i)%center - disks(i)%radius
       bounds(2*i) = disks(i)%center + disks(i)%radius       
    end do
    deallocate(disks) !! good idea to clear
  end function calculate_gershgorin_bounds

  ! Helper function to create identity matrix
  pure function identity_matrix(n) result(I)
    implicit none
    integer, intent(in) :: n
    complex(real64) :: I(n,n)
    integer :: j
    I = cmplx(0.0d0, 0.0d0)
    do j = 1, n
       I(j,j) = cmplx(1.0d0, 0.0d0)
    end do
  end function identity_matrix
  pure real(real64) function Cnorm2(x)
    complex(kind=real64), intent(in) :: x(:)
    Cnorm2 = sqrt( dot_product(x,x) )
  end function Cnorm2

end program eigenvalue_calculator
