program test_check_kywd_compatible_pass2
 use mr_keyword, only: ist, check_kywd_compatible
 use mr_keyword, only: X2C, casscf, casscf_prog, hf_prog, skiphf
 use test_utils, only: reset_mr_keyword_defaults
 implicit none

 call reset_mr_keyword_defaults()

 X2C = .true.
 casscf = .true.
 casscf_prog = 'orca'
 hf_prog = 'gaussian'
 skiphf = .false.
 ist = 0
 call check_kywd_compatible()

 write(6,'(A)') 'PASS: check_kywd_compatible pass2'
end program test_check_kywd_compatible_pass2
