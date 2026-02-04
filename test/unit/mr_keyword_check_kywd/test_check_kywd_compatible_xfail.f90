program test_check_kywd_compatible_xfail
 use mr_keyword, only: ist, new_mult, check_kywd_compatible
 use test_utils, only: reset_mr_keyword_defaults
 implicit none

 call reset_mr_keyword_defaults()

 ! Case: invalid NewMult with ist != 5 should fail
 ist = 1
 new_mult = 2
 call check_kywd_compatible()

 write(6,'(A)') 'FAIL: check_kywd_compatible did not stop on invalid NewMult'
 stop 2
end program test_check_kywd_compatible_xfail
