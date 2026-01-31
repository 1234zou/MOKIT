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
end module test_utils
