!!! File   : disk.f90
!!! Author : Sandeep Koranne (C) 2025. All rights reserved
!!! Purpose: Implementation of Poisson/Dirichlet problem on disk
!!! HDF5_FC=ifx /scratch1/skoranne/INTELCOMPILERS/bin/h5fc -O3 disk.f90 -o disk.exe
!!! h5ls and h5stat show file info, and also use Julia HDF5 package
module poisson_solver
  use iso_fortran_env
  use hdf5
  implicit none
  private
  public :: SolvePoissonDirichlet
  real(real64), parameter :: pi = 4*atan(1.0_real64)    
contains
  pure elemental real(real64) function sincube(angle)
    real(real64),intent(in) :: angle
    sincube = sin(angle)**3
  end function sincube
  
  function SolvePoissonDirichlet(nr, ntheta,B) result(soln)
    integer, intent(in) :: nr, ntheta
    real(real64), intent(in) :: B
    real(real64), allocatable :: u(:)
    ! Grid setup
    real(real64), allocatable :: r(:), theta(:)
    real(real64), allocatable :: A(:, :)
    real(real64), allocatable :: g(:,:)
    real(real64), allocatable :: soln(:)
    integer :: ir, itheta
    allocate(r(1:nr))
    allocate(theta(1:ntheta))
    allocate(soln(nr*ntheta))
    allocate(u(nr*ntheta))
    ! Generate grid points
    call GenerateGrid( nr, ntheta, B, r, theta)
    ! Initialize Green's matrix
    allocate(A(ntheta*nr,ntheta)) !skoranne this is a change
    call ComputePoissonKernel( nr, ntheta, B, r, theta, A )
    ! Boundary data
    allocate(g(ntheta,1))
    g(:,1) = sin(2*theta)
    !g(:,1) = (sin(theta))**3 ! sin(3x) = 3 sin(x) - 4 sin^3(x)
    do ir=1,nr
       !print *,'soln at r=',r(ir),'ntheta=',ntheta,' ',(ir-1)*ntheta+1,':',(ir)*ntheta
       !soln((ir-1)*ntheta+1:ir*ntheta) = (r(ir)*3.0_real64/4.0_real64)*sin(theta)-(r(ir)**3)*0.25_real64*sin(3*theta)
       soln((ir-1)*ntheta+1:ir*ntheta) = r(ir)*sin(2*theta)
       !write (*,'(8F8.4)') soln((ir-1)*ntheta+1:ir*ntheta)
    end do
    !print *,'soln at boundary='
    !write (*,'(8F8.4)') sin(theta)
    !write (*,'(8F8.4)') (B*3.0_real64/4.0_real64)*sin(theta)-B**3*0.25_real64*sin(3*theta)
    !print *,'soln=',soln
    ! Solve the linear system
    call ApplyOperator(nr, ntheta, B, A, g, soln, u)
  end function SolvePoissonDirichlet

!!!==============================================================================
!!! subroutine GenerateGrid
!!!
!!! Generates structured radial and angular grid points for polar coordinate systems
!!!
!!! This subroutine creates uniform grid points for 2D polar coordinates, typically
!!! used in computational physics and numerical simulations involving cylindrical
!!! symmetry. The grid points are generated with uniform spacing in both radial
!!! and angular directions.
!!!
!!! The subroutine generates interior grid points only, excluding the boundary
!!! points (0 and B) to facilitate numerical methods that require interior points
!!! for finite difference or finite element calculations.
!!! Parameters:
!!!   nr      - Number of radial grid points (interior points only)
!!!   ntheta  - Number of angular grid points (uniform distribution around circle)
!!!   B       - Boundary parameter defining radial extent (0 to B)
!!!
!!! Output Parameters:
!!!   r       - Array of radial grid points (allocated, size nr)
!!!   theta   - Array of angular grid points (allocated, size ntheta)
!!!
!!! Grid Generation Details:
!!!   Radial Grid: Uniform spacing from dr to B (excluding boundaries)
!!!                where dr = B / (nr + 1)
!!!   Angular Grid: Uniform distribution from 0 to 2π radians
!!!                 with spacing 2π/ntheta radians
!!!
!!! Example Usage:
!!!   Call GenerateGrid(10, 20, 5.0_real64, radial_points, angular_points)
!!!   ! Generates 10 radial points from 0.5 to 5.0 and 20 angular points from 0 to 2π
!!!
  
  subroutine GenerateGrid(nr, ntheta, B, r, theta)
    integer, intent(in) :: nr, ntheta
    real(real64), intent(in) :: B
    real(real64), allocatable :: r(:), theta(:)
    ! Generate radial grid points (interior points only, excluding boundary)
    real(real64) :: dr 
    integer :: i
    dr = B / (nr + 1)
    r = [(i*dr, i=1, nr)]
    ! Generate angular grid points
    theta = [(2*pi*(i-1)*1.0_real64/ntheta, i=1, ntheta)]
    !print *,'dr=',dr,' r =',r,' theta=',theta
  end subroutine GenerateGrid

!!!==============================================================================
!!! subroutine ComputePoissonKernel
!!!
!!! Computes the Poisson kernel for a circular domain used in solving Laplace's equation
!!!
!!! This subroutine generates the Poisson kernel matrix A that relates boundary values to
!!! interior point values in a circular domain. The kernel is fundamental for harmonic
!!! analysis and boundary value problems in circular geometries.
!!!
!!! The Poisson kernel is defined as: P(r,θ,θ') = (B² - r²) / (B² + r² - 2Br cos(θ-θ'))
!!! where B is the boundary radius, r is the radial coordinate, and θ is the angular coordinate
!!! Parameters:
!!!   nr      - Number of radial grid points
!!!   ntheta  - Number of angular grid points
!!!   B       - Boundary radius of the circular domain
!!!   r(:)    - Array of radial grid points (input)
!!!   theta(:) - Array of angular grid points (input)
!!!   A(:, :) - Output Poisson kernel matrix (nr*ntheta × ntheta)
!!!
!!! Input Validation:
!!!   - Validates that output matrix A has correct dimensions (nr*ntheta × ntheta)
!!!
!!! Algorithm:
!!!   1. Pre-computes numerator = (B² - r²)/(2π) and denominator = B² + r²
!!!   2. For each boundary angular point, computes angular differences
!!!   3. Computes kernel values using the Poisson kernel formula
!!!   4. Normalizes each column to sum to 1 for numerical stability
!!!
!!! Output Matrix Structure:
!!!   - Rows correspond to interior grid points (radial × angular indexing)
!!!   - Columns correspond to boundary angular points
!!!   - A(i,j) represents kernel value for interior point i due to boundary point j
!!!
!!! Mathematical Properties:
!!!   - Kernel exhibits symmetry
!!!   - Values are positive for valid domain configurations
!!!   - Column normalization ensures probability interpretation
!!!
!!! Applications:
!!!   - Solving Dirichlet problems for Laplace's equation
!!!   - Harmonic analysis in circular domains
!!!   - Potential theory computations
!!!   - Numerical methods for electrostatics and fluid dynamics
!!!
!!! Author: Sandeep Koranne (C)
!!! Date: April 2024
!!! Version: 1.0
!!!
!!!==============================================================================
  
  ! nr is the number of radial divisions, ntheta the angular divisions
  ! nr*ntheta is the total number of unknowns
  ! function value on the boundary is known
  ! Poisson Kernel := \frac{B^2-r^2}{B^2+r^2 - 2 B r cos(phi-theta)
  ! we organize the matrix by phi values (on the boundary)
  subroutine ComputePoissonKernel(nr, ntheta, B, r, theta, A)
    integer, intent(in) :: nr, ntheta
    real(real64), intent(in) :: B ! boundary radius
    real(real64), intent(in) :: r(:), theta(:)
    real(real64), intent(out) :: A(:, :)
    integer :: i, btheta
    real(real64), allocatable :: numerator(:), denominator(:), cosdiff(:),row_sum(:)
    if( size(A,dim=1) /= nr*ntheta .or. size(A,dim=2) /= ntheta ) then
       stop "Error in A dimensions"
    end if
    numerator = (1.0_real64/(2.0_real64*pi))*(B**2 - r**2)
    denominator = B**2 + r**2
    !print *,'Numerator:'
    !write (*,'(8F8.4)') numerator
    !print *,'Denominator:'    
    !write (*,'(8F8.4)') denominator
    !print *,''
    do btheta = 1, ntheta
       cosdiff  = [(cos(2.0_real64*pi/ntheta*(btheta-1-i)),i=0,ntheta-1)]
       do i = 1, nr
          A((i-1)*ntheta+1:(i*ntheta),btheta) = numerator(i)/(denominator(i) - 2.0*B*r(i)*cosdiff)
       end do
    end do
    !for some reason the matrix A is not normalized
    !print *,'size(A)',size(A,dim=1),' ',size(A,dim=2),'row_sum=',sum(A,dim=2)
    row_sum = sum(A,dim=2) !sum along dim=2 columns

    do i=1,size(A,dim=2)
       A(:,i) = A(:,i)/row_sum
    end do
    !print *,'size(A)',size(A,dim=1),' ',size(A,dim=2),'row_sum=',sum(A,dim=2)    
  end subroutine ComputePoissonKernel

!!!==============================================================================
!!! subroutine ApplyOperator
!!!
!!! Applies the Poisson operator to boundary data and compares with exact solution
!!!
!!! This subroutine computes the numerical solution to the Poisson equation using the
!!! pre-computed Poisson kernel matrix A. It applies the operator to boundary data g,
!!! compares the result with the exact solution soln, and provides detailed output
!!! including error analysis and data export.
!!!
!!! The subroutine performs:
!!!   1. Matrix-vector multiplication: u = A * g(:,1) (applying Poisson operator)
!!!   2. Comparison with exact solution
!!!   3. Grid point generation for output
!!!   4. Detailed error analysis and printing
!!!   5. Data export to HDF5 file
!!!
!!! Parameters:
!!!   nr      - Number of radial grid points
!!!   ntheta  - Number of angular grid points
!!!   B       - Boundary radius of the circular domain
!!!   A(:, :) - Poisson kernel matrix (nr*ntheta × ntheta)
!!!   g(:,:)  - Boundary data matrix (size: nr*ntheta × ntheta)
!!!   soln(:) - Exact (size: nr*ntheta)
!!!   u(:)    - Output computed solution vector (size: nr*ntheta)
!!!
!!! Algorithm:
!!!   1. Computes u = A * g(:,1) using matrix multiplication
!!!   2. Generates radial and angular grid points using GenerateGrid
!!!   3. Performs element-wise comparison between computed and exact solutions
!!!   4. Calculates L2 norm of error: ‖u-soln‖₂
!!!   5. Exports all data to HDF5 file for further analysis
!!!
!!! Output:
!!!   - Detailed comparison of computed vs exact solution at each grid point
!!!   - L2 error norm between computed and exact solutions
!!!   - Data export to HDF5 file with all grid points and solution data
!!!
!!! Mathematical Context:
!!!   - Solves Δu = f in circular domain with boundary data g
!!!   - Uses Poisson kernel method: u = A * g
!!!   - Compares numerical solution with analytical solution
!!!
!!! Dependencies:
!!!   - GenerateGrid: For grid point generation
!!!   - PrintDataToHDF5: For data export functionality
!!!
!!! Usage:
!!! Called after Poisson kernel matrix A is computed
!!!   - Requires boundary data g and exact solution soln
!!!
!!! Author: Sandeep Koranne (C)
!!! Date: April 2024
!!! Version: 1.0
!!!
!!!==============================================================================
  
  subroutine ApplyOperator(nr, ntheta, B, A, g, soln,u)
    integer, intent(in) :: nr, ntheta
    real(real64), intent(in) :: B
    real(real64), intent(in) :: A(:, :)
    real(real64), intent(in) :: g(:,:)
    real(real64), intent(in) :: soln(:)    
    real(real64), intent(out) :: u(:)
    integer :: i,j
    real(real64), allocatable :: r(:), theta(:)
    
    !print *,'A'
    !by default Fortran writes arrays in column order to screen
    !write( * , "(*(g0.3))" ) ( (A(i,j)," ",j=1,size(A,dim=2)), new_line("A"), i=1,size(A,dim=1) )
    u=matmul(A,g(:,1))
    !print *,'g',' size(g) = ', size(g,dim=1)
    !write (*,'(8F8.4)') g
    print *,'-------------- Soln -------------------'
    !write (*,'(8F8.4)') u
    !write (*,'(8F8.4)') soln
    allocate(r(1:nr))
    allocate(theta(1:ntheta))
    ! Generate grid points for output writing
    call GenerateGrid( nr, ntheta, B, r, theta)
    do i=1,nr
       do j=1,ntheta
          associate (&
               &op_sol => u((i-1)*ntheta+j), &
               &nm_sol => soln((i-1)*ntheta+j))
            !write (*,'(7(A,F8.4))') 'r= ',r(i),' θ= ',theta(j),&
            !     ' u= ',op_sol,' B= ',B,' bdry= ',sin(2*theta(j)),' soln= ', nm_sol,' Δ=',abs(op_sol-nm_sol)
          end associate
       end do
    end do
    print *,'‖u-soln‖₂=',norm2(u-soln)
    call PrintDataToHDF5(nr,ntheta,A,B,r,theta,soln,u,g(:,1))
  end subroutine ApplyOperator

  subroutine PrintDataToHDF5(nr,ntheta,A,B,r,theta,soln,u,g)
    integer, intent(in) :: nr, ntheta
    real(real64), intent(in) :: A(:,:),B,r(:),theta(:),soln(:),u(:),g(:)
    integer :: total_unknowns
    character(len=*),parameter :: fileName = "poisson.h5"
    character(len=*),parameter :: groupName= "/PoissonGroup"
    character(len=*),parameter :: gridDsetName = "/PoissonGroup/U"
    character(len=*),parameter :: gridSizeDsetName = "GridSize"
    character(len=*),parameter :: linopDsetName = "Operator"
    character(len=*),parameter :: solnDsetName = "Soln"
    character(len=*),parameter :: bdryDsetName = "BoundaryValue"    
    character(len=*),parameter :: rDsetName = "R"
    character(len=*),parameter :: thetaDsetName = "Theta"    
    integer(HID_T):: fileid, gid !this is the order file,group,space,set
    integer,parameter :: num_datasets = 8
    !since we have many datasets, we will index them
    integer(HSIZE_T),dimension(num_datasets,2) :: DID
    integer       :: error
    integer(HID_T):: atype_id
    integer(HSIZE_T),dimension(2) :: dims
    dims(1) = nr
    dims(2) = ntheta
    total_unknowns = nr*ntheta
    if( total_unknowns /= size(soln) .or. total_unknowns /= size(u) ) stop "Dimensions mismatched"

    call h5open_f(error)
    call h5fcreate_f(fileName, H5F_ACC_TRUNC_F, fileid, error)
    call h5gcreate_f(fileid, groupName, gid, error)
    dims(1)=1
    call h5screate_simple_f(1,dims(1),DID(6,1),error)
    dims(2)=nr
    call h5acreate_f(gid,"nr",H5T_NATIVE_INTEGER,DID(6,1),DID(6,2),error)
    CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
    CALL h5tset_size_f(atype_id, 1, error)
    !call h5awrite_f(DID(6,1),atype_id,dims(2),error)
    call h5aclose_f(DID(6,2),error)
    call h5tclose_f(atype_id,error)

    dims(1)=2
    call h5screate_simple_f(1,dims,DID(7,1),error)
    call h5dcreate_f(gid, gridSizeDsetName, H5T_NATIVE_INTEGER, DID(7,1), DID(7,2), error)
    call h5dwrite_f(DID(7,2),H5T_NATIVE_INTEGER, [nr,ntheta], dims, error)
    call h5dclose_f(DID(7,2),error)    
    call h5sclose_f(DID(7,1),error)
    
    dims(1)=nr
    dims(2)=ntheta
    call h5screate_simple_f(2,dims,DID(1,1),error)
    call h5dcreate_f(fileid, gridDsetName, H5T_IEEE_F64LE, DID(1,1), DID(1,2), error)
    call h5dwrite_f(DID(1,2),H5T_IEEE_F64LE, reshape(u,[nr,ntheta]), dims, error)
    call h5dclose_f(DID(1,2),error)    
    call h5sclose_f(DID(1,1),error)

    
    call h5screate_simple_f(2,dims,DID(3,1),error)
    call h5dcreate_f(fileid, solnDsetName, H5T_IEEE_F64LE, DID(3,1), DID(3,2), error)
    call h5dwrite_f(DID(3,2),H5T_IEEE_F64LE, reshape(soln,[nr,ntheta]), dims, error)
    call h5dclose_f(DID(3,2),error)        
    call h5sclose_f(DID(3,1),error)
    
    dims(2) = nr*ntheta !watch this is transpose of A
    dims(1) = ntheta
    call h5screate_simple_f(2,dims,DID(2,1),error)
    call h5dcreate_f(gid, linopDsetName, H5T_IEEE_F64LE, DID(2,1), DID(2,2), error)
    call h5dwrite_f(DID(2,2),H5T_IEEE_F64LE,A,dims,error)
    call h5dclose_f(DID(2,2),error)        
    call h5sclose_f(DID(2,1),error)

    dims(1) = nr
    call h5screate_simple_f(1,dims,DID(4,1),error)
    call h5dcreate_f(gid, rDsetName, H5T_IEEE_F64LE, DID(4,1), DID(4,2), error)
    call h5dwrite_f(DID(4,2),H5T_IEEE_F64LE,r,dims,error)
    call h5dclose_f(DID(4,2),error)        
    call h5sclose_f(DID(4,1),error)

    dims(1) = ntheta
    call h5screate_simple_f(1,dims,DID(5,1),error)
    call h5dcreate_f(gid, thetaDsetName, H5T_IEEE_F64LE, DID(5,1), DID(5,2), error)
    call h5dwrite_f(DID(5,2),H5T_IEEE_F64LE,theta,dims,error)
    call h5dclose_f(DID(5,2),error)        
    call h5sclose_f(DID(5,1),error)

    call h5screate_simple_f(1,dims,DID(8,1),error)
    call h5dcreate_f(gid, bdryDsetName, H5T_IEEE_F64LE, DID(8,1), DID(8,2), error)
    call h5dwrite_f(DID(8,2),H5T_IEEE_F64LE, g, dims, error)
    call h5dclose_f(DID(8,2),error)        
    call h5sclose_f(DID(8,1),error)
    
    call h5gclose_f(gid,error)
    call h5fclose_f(fileid,error)
    call h5close_f(error)
  end subroutine PrintDataToHDF5
  
end module poisson_solver

!!!==============================================================================
!!! program main
!!!
!!! Main program for solving Poisson's equation using Dirichlet boundary conditions
!!!
!!! This program serves as the entry point for computing the solution to the Poisson
!!! equation in a circular domain using the Poisson kernel method. It handles command
!!! line arguments for grid parameters and orchestrates the solution computation. The program:
!!!   1. Parses command line arguments for grid dimensions
!!!   2. Initializes grid parameters (nr, ntheta)
!!!   3. Calls the Poisson solver routine to compute the solution
!!!   4. Provides output indicating successful computation
!!!
!!! Command Line Arguments:
!!!   - First argument: Number of radial grid points (nr)
!!!   - Second argument: Number of angular grid points (ntheta)
!!!   - If no arguments provided, defaults to nr=4, ntheta=8
!!!
!!! Input Parameters:
!!!   - nr: Number of radial grid points (default: 4)
!!!   - ntheta: Number of angular grid points (default: 8)
!!!
!!! Output:
!!!   - Prints grid dimensions used
!!!   - Displays success message upon completion
!!!
!!! Dependencies:
!!!   - poisson_solver module: Contains SolvePoissonDirichlet subroutine
!!!   - iso_fortran_env: For standard environment features
!!!
!!! Usage Examples:
!!!   ./program                    ! Uses default grid sizes
!!!   ./program 32 64              ! Uses 32×64 grid
!!!   ./program 16 32              ! Uses 16× grid
!!!
!!! Mathematical Context:
!!!   - Solves Δu = f in circular domain with Dirichlet boundary conditions
!!!   - Uses Poisson kernel method for numerical solution
!!!   - Implements harmonic analysis in polar coordinates
!!!
!!! Author: Sandeep Koranne (C) 2025
!!! Date: April 2024
!!! Version: 1.0
!!!
!!!==============================================================================

program main
  use poisson_solver
  use iso_fortran_env
  implicit none
  integer :: num_args,stat
  !integer :: nr = 128, ntheta = 128
  integer :: nr=4, ntheta =8
  integer,dimension(10) :: args
  integer :: i
  real(real64), allocatable :: soln(:)
  character(len=255) :: arg_str
  args = 0
  num_args = command_argument_count()
  if( num_args > 0 ) then
     do i=1,min(num_args,10)
        call get_command_argument( i, arg_str )
        read( arg_str, *, iostat=stat ) args(i)
        if( stat /= 0 ) then
           print *,'Arg conversion failed: ', i,' ', trim(arg_str)
           args(i) = -1
        end if
     end do
  end if
  print *,'nargc = ', num_args
  if( args(1) > 0 ) nr = args(1)
  if( args(1) > 0 ) ntheta = args(2)  
  print *,'nr=',nr,'ntheta=',ntheta
  soln = SolvePoissonDirichlet(nr, ntheta,1.0_real64)
  !print *,'α'
  !Δ
  print *, "Solution to Poisson problem Δ computed successfully."
end program main
