! written by jxzou at 20260116: ab initio PBC calculation helper
! Currently only implement: scanning various thresholds for a target CP2K calculation

program main
 implicit none
 character(len=240) :: inpname

 inpname = ''
 call getarg(1, inpname)
 call scan_eps_default_in_cp2k_inp(inpname)
end program main

subroutine set_eps_default_in_cp2k_inp(inpname, eps_default)
 implicit none
 integer :: i, k, fid, fid1, RENAME
 real(kind=8), intent(in) :: eps_default
!f2py intent(in) :: eps_default
 character(len=49), parameter :: error_warn = 'ERROR in subroutine set_eps_defa&
                                              &ult_in_cp2k_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(INDEX(buf,'&QS') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&QS" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = INDEX(buf,'EPS_DEFAULT')
  if(k > 0) exit
  if(INDEX(buf,'&END QS') > 0) then
   i = -1; exit
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"EPS_DEFAULT" cannot be found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A,ES7.1)') buf(1:k+10)//' ', eps_default

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_eps_default_in_cp2k_inp

subroutine set_eps_scf_in_cp2k_inp(inpname, eps_scf)
 implicit none
 integer :: i, k, fid, fid1, RENAME
 real(kind=8), intent(in) :: eps_scf
!f2py intent(in) :: eps_scf
 character(len=45), parameter :: error_warn = 'ERROR in subroutine set_eps_scf_&
                                              &in_cp2k_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(INDEX(buf,'&SCF') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&SCF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = INDEX(buf,'EPS_SCF')
  if(k > 0) exit
  if(INDEX(buf,'&END SCF') > 0) then
   i = -1; exit
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"EPS_SCF" cannot be found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A,ES7.1)') buf(1:k+6)//' ', eps_scf

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_eps_scf_in_cp2k_inp

subroutine set_mgrid_cutoff_in_cp2k_inp(inpname, cutoff)
 implicit none
 integer :: i, k, fid, fid1, RENAME
 integer, intent(in) :: cutoff
!f2py intent(in) :: cutoff
 character(len=50), parameter :: error_warn = 'ERROR in subroutine set_mgrid_cu&
                                              &toff_in_cp2k_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(INDEX(buf,'&MGRID') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&MGRID" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = INDEX(buf,'CUTOFF')
  if(k>0 .and. INDEX(buf,'REL_CUTOFF')==0) exit
  if(INDEX(buf,'&END MGRID') > 0) then
   i = -1; exit
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"CUTOFF" cannot be found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A,I0)') buf(1:k+5)//' ', cutoff

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_mgrid_cutoff_in_cp2k_inp

subroutine set_mgrid_rel_cutoff_in_cp2k_inp(inpname, rel_cutoff)
 implicit none
 integer :: i, k, fid, fid1, RENAME
 integer, intent(in) :: rel_cutoff
!f2py intent(in) :: rel_cutoff
 character(len=54), parameter :: error_warn = 'ERROR in subroutine set_mgrid_re&
                                              &l_cutoff_in_cp2k_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(INDEX(buf,'&MGRID') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&MGRID" cannot be found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = INDEX(buf,'REL_CUTOFF')
  if(k > 0) exit
  if(INDEX(buf,'&END MGRID') > 0) then
   i = -1; exit
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"REL_CUTOFF" cannot be found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A,I0)') buf(1:k+9)//' ', rel_cutoff

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_mgrid_rel_cutoff_in_cp2k_inp

subroutine increase_int_in_cp2k_inpname(inpname, new_inp)
 implicit none
 integer :: i, j, k, m
 character(len=240), intent(in) :: inpname
 character(len=240), intent(out) :: new_inp

 call find_specified_suffix(inpname, '.in', i)
 k = INDEX(inpname, '_', back=.true.)
 if(k == 0) then
  new_inp = inpname(1:i-1)//'_1.inp'
 else
  read(inpname(k+1:i-1),fmt=*,iostat=j) m
  if(j == 0) then
   write(new_inp,'(A,I0,A)') inpname(1:k),m+1,'.inp'
  else
   new_inp = inpname(1:i-1)//'_1.inp'
  end if
 end if
end subroutine increase_int_in_cp2k_inpname

subroutine scan_eps_default_in_cp2k_inp(inpname)
 implicit none
 integer :: i, k
 integer, parameter :: nstep = 4
 integer, parameter :: nproc = 64
 real(kind=8) :: eps_default, e0, e1
 real(kind=8), parameter :: eps_default0 = 1d-9
 real(kind=8), parameter :: e_diff = 1d-6
 character(len=240) :: old_inp, new_inp, outname
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname
 logical :: converged

 converged = .false.; e0 = 1d0
 eps_default = eps_default0; old_inp = inpname

 do i = 1, nstep, 1
  call increase_int_in_cp2k_inpname(old_inp, new_inp)
  call find_specified_suffix(new_inp, '.in', k)
  outname = new_inp(1:k-1)//'.out'
  call sys_copy_file(TRIM(inpname), TRIM(new_inp), .false.)
  eps_default = eps_default*0.1d0
  call set_eps_default_in_cp2k_inp(new_inp, eps_default)
  call submit_cp2k_job(new_inp, nproc)
  call read_scf_e_from_cp2k_out(outname, e1)
  if(DABS(e0-e1) < 1d-6) then
   converged = .true.
   exit
  end if
  old_inp = new_inp
  e0 = e1
 end do ! for i

 if(.not. converged) then
  write(6,'(/,A)') 'ERROR in subroutine scan_eps_default_in_cp2k_inp: scanning &
                   &EPS_DEFAULT does'
  write(6,'(A)') 'not lead to a converged result. This is not likely to happen.&
                 & Please check your file '
  write(6,'(A)') TRIM(inpname)
  stop
 end if
 write(6,'(A,ES7.1)') 'eps_default=', eps_default
end subroutine scan_eps_default_in_cp2k_inp

