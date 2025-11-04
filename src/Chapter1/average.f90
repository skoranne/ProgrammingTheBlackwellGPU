!* *qwen3-coder:latest* >> PROMPT: generate a Fortran program that reads 3 integers and then computes and prints their sum and average. Use modern Fortran idioms like iso_fortran_env and also prompt the user for Input.


!** [qwen3-coder:latest: RESPONSE]



!*** FINISHED#+begin_src fortran
program sum_and_average
    use iso_fortran_env, only: input_unit, output_unit, real64
    implicit none
    
    integer, parameter :: dp = real64
    integer :: a, b, c
    real(dp) :: sum_val, average
    
    ! Prompt user for input
    write(output_unit, '(A)') 'Enter three integers:'
    write(output_unit, '(A)', advance='no') 'First integer: '
    read(input_unit, *) a
    write(output_unit, '(A)', advance='no') 'Second integer: '
    read(input_unit, *) b
    write(output_unit, '(A)', advance='no') 'Third integer: '
    read(input_unit, *) c
    
    ! Compute sum and average
    sum_val = real(a, dp) + real(b, dp) + real(c, dp)
    average = sum_val / 3.0_dp
    
    ! Print results
    write(output_unit, '(A, I0)') 'Sum: ', int(sum_val)
    write(output_unit, '(A, F0.2)') 'Average: ', average
    
end program sum_and_average
!#+end_src
