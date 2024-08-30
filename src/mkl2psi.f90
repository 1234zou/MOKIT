! written by jxzou at 20230728: a wrapper of mkl2fch and fch2psi for ORCA ->
! PSI4

program mkl2psi
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, narg, irel, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname
 character(len=240) :: a_file0, b_file0, a_file, b_file
 character(len=240), allocatable :: str_arg(:)
 logical :: alive

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,A)') ' ERROR in program mkl2psi: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2psi a.mkl'
  write(6,'(A)')   ' Example 2: mkl2psi a.mkl b.inp'
  write(6,'(A)')   ' Example 3: mkl2psi a.mkl -dkh2'
  write(6,'(A)')   ' Example 4: mkl2psi a.mkl -sfx2c'
  write(6,'(A,/)') ' Example 5: mkl2psi a.mkl b.inp -sfx2c'
  stop
 end if

 allocate(str_arg(narg))
 do i = 1, narg, 1
  call getarg(i, str_arg(i))
 end do ! for i

 mklname = str_arg(1)
 call require_file_exist(mklname)
 call find_specified_suffix(mklname, '.mkl', i)
 inpname = mklname(1:i-1)//'.inp'
 k = 2; irel = -1

 if(narg > 1) then
  i = LEN_TRIM(str_arg(2))
  if(str_arg(2)(i-3:i) == '.inp') then
   inpname = str_arg(2)
   k = 3
  end if
  if(k <= narg) then
   select case(TRIM(str_arg(k)))
   case('-sfx2c')
    irel = -3
   case('-dkh2')
    irel = 2
   case default
    write(6,'(/,A)') 'ERROR in subroutine mkl2psi: wrong command line argument!'
    write(6,'(A)')   'Example: mkl2psi a.mkl -sfx2c'
    stop
   end select
  end if
 end if

 deallocate(str_arg)
 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)

 i = SYSTEM('fch2psi '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2psi: failed to call utility fch2psi.'
  write(6,'(A)') 'mkl2psi is a wrapper of mkl2fch and fch2psi, so fch2psi is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'_psi.inp'
 a_file0 = fchname(1:i-1)//'.A'
 b_file0 = fchname(1:i-1)//'.B'

 call find_specified_suffix(inpname, '.inp', i)
 a_file = inpname(1:i-1)//'.A'
 b_file = inpname(1:i-1)//'.B'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
 i = RENAME(TRIM(a_file0), TRIM(a_file))
 inquire(file=TRIM(b_file0), exist=alive)
 if(alive) i = RENAME(TRIM(b_file0), TRIM(b_file))

 call modify_orbname_in_psi4_inp(inpname, a_file)
end program mkl2psi

! modify the orbital filename in PSI4 .inp file
subroutine modify_orbname_in_psi4_inp(inpname, a_file)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1, b_file
 character(len=240), intent(in) :: inpname, a_file

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'

 call find_specified_suffix(a_file, '.A', i)
 b_file = a_file(1:i-1)//'.B'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:9) == 'scfenergy') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  if(buf(1:10) == 'scf_wfn.Ca') buf="scf_wfn.Ca().load('"//TRIM(a_file)//"')"
  if(buf(1:10) == 'scf_wfn.Cb') buf="scf_wfn.Cb().load('"//TRIM(b_file)//"')"
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_orbname_in_psi4_inp

