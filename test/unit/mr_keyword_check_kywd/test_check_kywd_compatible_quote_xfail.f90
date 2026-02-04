program test_check_kywd_compatible_quote_xfail
 use mr_keyword, only: ist, check_kywd_compatible
 use mr_keyword, only: mcpdft, mcpdft_prog
 use test_utils, only: reset_mr_keyword_defaults
 implicit none

 call reset_mr_keyword_defaults()

 mcpdft = .true.
 mcpdft_prog = "'pyscf'"
 ist = 0
 call check_kywd_compatible()

 write(6,'(A)') 'FAIL: check_kywd_compatible did not stop on single quotes'
 stop 2
end program test_check_kywd_compatible_quote_xfail
