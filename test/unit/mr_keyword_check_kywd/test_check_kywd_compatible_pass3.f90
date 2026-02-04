program test_check_kywd_compatible_pass3
 use mr_keyword, only: ist, check_kywd_compatible
 use mr_keyword, only: casci, casci_prog
 use test_utils, only: reset_mr_keyword_defaults
 implicit none

 call reset_mr_keyword_defaults()

 casci = .true.
 casci_prog = 'orca'
 ist = 0
 call check_kywd_compatible()

 write(6,'(A)') 'PASS: check_kywd_compatible pass3'
end program test_check_kywd_compatible_pass3
