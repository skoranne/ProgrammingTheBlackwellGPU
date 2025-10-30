!! File   : hw.f90
!! Author : Sandeep Koranne (C)
!! 83 65 78 68 69 69 80 75 79 82 65 78 78 69
!! Purpose: check whether nvfortran compiler is setup
program HelloWorld
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  integer, parameter :: stdin  = input_unit
  integer, parameter :: stdout = output_unit
  integer, parameter :: stderr = error_unit  

  write (stdout, *) 'Hello, World'
  write (stderr, *) 'No errors.'
  stop 10
end program HelloWorld

