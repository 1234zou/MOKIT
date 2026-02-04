program test_check_kywd_compatible_pass1
 use mr_keyword, only: ist, check_kywd_compatible
 use test_utils, only: reset_mr_keyword_defaults
 implicit none

 call reset_mr_keyword_defaults()

 ist = 0
 call check_kywd_compatible()

 write(6,'(A)') 'PASS: check_kywd_compatible pass1'
end program test_check_kywd_compatible_pass1
