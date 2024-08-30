! written by jxzou at 20230728: a wrapper of mkl2fch and fch2amo for ORCA->Amesp

program mkl2amo
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, narg, irel, SYSTEM, RENAME
 character(len=240), allocatable :: str_arg(:)
 character(len=240) :: mklname, fchname, inpname0, inpname, orbfile0, orbfile

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,A)') ' ERROR in program mkl2amo: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2amo a.mkl'
  write(6,'(A)')   ' Example 2: mkl2amo a.mkl b.aip'
  write(6,'(A)')   ' Example 3: mkl2amo a.mkl -sfx2c'
  write(6,'(A,/)') ' Example 4: mkl2amo a.mkl b.aip -sfx2c'
  stop
 end if

 allocate(str_arg(narg))
 do i = 1, narg, 1
  call getarg(i, str_arg(i))
 end do ! for i

 mklname = str_arg(1)
 call require_file_exist(mklname)
 call find_specified_suffix(mklname, '.mkl', i)
 inpname = mklname(1:i-1)//'.aip'
 k = 2; irel = -1

 if(narg > 1) then
  i = LEN_TRIM(str_arg(2))
  if(str_arg(2)(i-3:i) == '.aip') then
   inpname = str_arg(2)
   k = 3
  end if
  if(k <= narg) then
   if(TRIM(str_arg(k)) == '-sfx2c') then
    irel = -3
   else
    write(6,'(/,A)') 'ERROR in subroutine mkl2fch: wrong command line argument!'
    write(6,'(A)')   'Example: mkl2amo a.mkl -sfx2c'
    stop
   end if
  end if
 end if

 deallocate(str_arg)
 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)

 i = SYSTEM('fch2amo '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2amo: failed to call utility fch2amo.'
  write(6,'(A)') 'mkl2amo is a wrapper of mkl2fch and fch2amo, so fch2amo is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.aip'
 orbfile0 = fchname(1:i-1)//'.amo'

 call find_specified_suffix(inpname, '.aip', i)
 orbfile = inpname(1:i-1)//'.amo'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
 i = RENAME(TRIM(orbfile0), TRIM(orbfile))
 write(6,'(/,A)') 'Conversion done. You need to use `a2m` for .amo -> .mo'
end program mkl2amo

