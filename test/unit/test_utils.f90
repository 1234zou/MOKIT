module test_utils
 implicit none
contains
 subroutine assert_int_equal(name, got, expected)
  character(len=*), intent(in) :: name
  integer, intent(in) :: got, expected
  if(got /= expected) then
   call fail(trim(name)//' expected='//trim(int_to_str(expected))// &
             ' got='//trim(int_to_str(got)))
  end if
 end subroutine assert_int_equal

 subroutine assert_real_equal(name, got, expected)
  character(len=*), intent(in) :: name
  real(kind=8), intent(in) :: got, expected
  if(got /= expected) then
   call fail(trim(name)//' expected='//trim(real_to_str(expected))// &
             ' got='//trim(real_to_str(got)))
  end if
 end subroutine assert_real_equal

 subroutine assert_real_close(name, got, expected, atol)
  character(len=*), intent(in) :: name
  real(kind=8), intent(in) :: got, expected, atol
  if(abs(got-expected) > atol) then
   call fail(trim(name)//' expected='//trim(real_to_str(expected))// &
             ' got='//trim(real_to_str(got)))
  end if
 end subroutine assert_real_close

 subroutine assert_path_suffix(name, value, suffix)
  character(len=*), intent(in) :: name, value, suffix
  integer :: i
  i = len_trim(value) - len_trim(suffix) + 1
  if(i < 1) call fail(trim(name)//' shorter than suffix')
  if(value(i:) /= suffix) then
   call fail(trim(name)//' expected suffix='//trim(suffix)// &
             ' got='//trim(value))
  end if
 end subroutine assert_path_suffix

 subroutine assert_string_equal(name, got, expected)
  character(len=*), intent(in) :: name, got, expected
  if(trim(got) /= trim(expected)) then
   call fail(trim(name)//' expected='//trim(expected)//' got='//trim(got))
  end if
 end subroutine assert_string_equal

 function int_to_str(i) result(out)
  integer, intent(in) :: i
  character(len=32) :: out
  write(out,'(I0)') i
 end function int_to_str

 function real_to_str(x) result(out)
  real(kind=8), intent(in) :: x
  character(len=64) :: out
  write(out,'(ES24.16E3)') x
 end function real_to_str

 subroutine fail(msg)
  character(len=*), intent(in) :: msg
  write(6,'(A)') 'FAIL: '//trim(msg)
  stop 1
 end subroutine fail

 subroutine reset_mr_keyword_defaults()
  use mr_keyword, only: ist, new_mult
  use mr_keyword, only: readrhf, readuhf, readno, on_thres, uno_thres, iroot, nstate
  use mr_keyword, only: nmr, icss, DKH2, X2C, RI, F12, CIonly
  use mr_keyword, only: casci, casscf, dmrgci, dmrgscf, gvb, c_gvb_conv, c_fcgvb
  use mr_keyword, only: gvb_prog, localm, casci_prog, casscf_prog, dmrgci_prog, dmrgscf_prog
  use mr_keyword, only: mcpdft, mcpdft_prog
  use mr_keyword, only: caspt2, hf_prog, skiphf
  ist = 0
  new_mult = 0
  readrhf = .false.
  readuhf = .false.
  readno = .false.
  on_thres = 0.2d0
  uno_thres = 1d-5
  iroot = 0
  nstate = 0
  nmr = .false.
  icss = .false.
  DKH2 = .false.
  X2C = .false.
  RI = .false.
  F12 = .false.
  CIonly = .false.
  casci = .false.
  casscf = .false.
  dmrgci = .false.
  dmrgscf = .false.
  gvb = .false.
  c_gvb_conv = .false.
  c_fcgvb = .false.
  gvb_prog = 'gamess'
  localm = 'pm'
  casci_prog = 'pyscf'
  casscf_prog = 'pyscf'
  dmrgci_prog = 'pyscf'
  dmrgscf_prog = 'pyscf'
  mcpdft = .false.
  mcpdft_prog = 'openmolcas'
  caspt2 = .false.
  hf_prog = 'gaussian'
  skiphf = .false.
 end subroutine reset_mr_keyword_defaults
end module test_utils
