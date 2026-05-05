! written by jxzou at 20260116: ab initio PBC calculation helper
! Currently only implement: scanning various thresholds/cutoffs for a target CP2K calculation

! The order for checking:
! 1) EPD_DEFAULT, EPS_SCF
! 2) CUTOFF, REL_CUTOFF
! 3) NGRIDS
! 4) USE_FINER_GRID

!program main
! implicit none
! integer :: conv_param
! character(len=240) :: inpname
!
! inpname = ''
! call getarg(1, inpname)
! call scan_cp2k_conv_param(inpname, 'EPS_DEFAULT', 1e-6, 64, conv_param)
!end program main

! set/activate WFN_RESTART_FILE_NAME and SCF_GUESS in &FORCE_EVAL/&DFT
subroutine set_wfn_restart_file_in_cp2k_inp(inpname, wfn_file)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 character(len=54), parameter :: error_warn = 'ERROR in subroutine set_wfn_rest&
                                              &art_file_in_cp2k_inp: '
 character(len=240) :: buf0, buf, inpname1, wfn_name
 character(len=240), intent(in) :: inpname, wfn_file
!f2py intent(in) :: inpname, wfn_file
 logical :: found

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 k = LEN_TRIM(wfn_file)
 ! if wfn_file is an empty string, use CP2K default xxx-RESTART.wfn
 if(k == 0) then
  wfn_name = inpname(1:i-1)//'-RESTART.wfn'
 else
  wfn_name = wfn_file
 end if

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:4) == '&DFT') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&DFT" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 found = .false.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:21) == 'WFN_RESTART_FILE_NAME') then
    found = .true.; exit
   end if
   if(buf(1:3) == '&QS') exit
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(found) then
  write(buf,'(A,I0,A)') '(',j-k,'X,A)'
  write(fid1,fmt=TRIM(buf)) 'WFN_RESTART_FILE_NAME '//TRIM(wfn_name)
 else
  write(buf,'(A,I0,A)') '(',j-k,'X,A)'
  write(fid1,fmt=TRIM(buf)) 'WFN_RESTART_FILE_NAME '//TRIM(wfn_name)
  write(fid1,fmt=TRIM(buf)) '&QS'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:4) == '&SCF') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&SCF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:7) == 'EPS_SCF') then
    write(fid1,'(A)') TRIM(buf0)
    exit
   end if
   if(buf(1:9) == 'SCF_GUESS') cycle
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 write(buf,'(A,I0,A)') '(',j-k,'X,A)'
 write(fid1,fmt=TRIM(buf)) 'SCF_GUESS RESTART'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:8) == '&END SCF') then
    write(fid1,'(A)') TRIM(buf0)
    exit
   end if
   if(buf(1:9) == 'SCF_GUESS') cycle
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&END SCF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_wfn_restart_file_in_cp2k_inp

! set/modify &FORCE_EVAL/&DFT/&QS/EPS_DEFAULT in a specified CP2K input file
subroutine set_eps_default_in_cp2k_inp(inpname, int_eps)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: int_eps ! 12 means 1e-12
!f2py intent(in) :: int_eps
 character(len=49), parameter :: error_warn = 'ERROR in subroutine set_eps_defa&
                                              &ult_in_cp2k_inp: '
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 if(int_eps<9 .or. int_eps>14) then
  write(6,'(/,A,I0)') error_warn//'invalid int_eps=', int_eps
  write(6,'(A)') 'The allowed range is [9,14], which means [1e-9,1e-14].'
  stop
 end if

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:3) == '&QS') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&QS" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:11) == 'EPS_DEFAULT') exit
   if(buf(1:7) == '&END QS') then
    i = -1; exit
   end if
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"EPS_DEFAULT" cannot be found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(buf,'(A,I0,A)') '(',j-k,'X,A,I0)'
 write(fid1,fmt=TRIM(buf)) 'EPS_DEFAULT 1e-', int_eps

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_eps_default_in_cp2k_inp

! set/modify &FORCE_EVAL/&DFT/&SCF/EPS_SCF in a specified CP2K input file
subroutine set_eps_scf_in_cp2k_inp(inpname, int_eps)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: int_eps ! 8 means 1e-8
!f2py intent(in) :: int_eps
 character(len=45), parameter :: error_warn = 'ERROR in subroutine set_eps_scf_&
                                              &in_cp2k_inp: '
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 if(int_eps<5 .or. int_eps>9) then
  write(6,'(/,A,I0)') error_warn//'invalid int_eps=', int_eps
  write(6,'(A)') 'The allowed range is [5,9], which means [1e-5,1e-9].'
  stop
 end if

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:4) == '&SCF') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&SCF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:7) == 'EPS_SCF') exit
   if(buf(1:8) == '&END SCF') then
    i = -1; exit
   end if
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"EPS_SCF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(buf,'(A,I0,A)') '(',j-k,'X,A,I0)'
 write(fid1,fmt=TRIM(buf)) 'EPS_SCF 1e-', int_eps

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_eps_scf_in_cp2k_inp

! set/modify NGRIDS in &FORCE_EVAL/&DFT/&MGRID
subroutine set_mgrid_ngrids_in_cp2k_inp(inpname, ngrids)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: ngrids
!f2py intent(in) :: ngrids
 character(len=50), parameter :: error_warn = 'ERROR in subroutine set_mgrid_ng&
                                              &rids_in_cp2k_inp: '
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname
 logical :: found

 if(ngrids<3 .or. ngrids>8) then
  write(6,'(/,A,I0)') error_warn//'invalid NGRIDS=', ngrids
  write(6,'(A)') 'The allowed range is [3,8].'
  stop
 end if

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:6) == '&MGRID') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&MGRID" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 found = .false.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:6) == 'NGRIDS') then
    found = .true.; exit
   end if
   if(buf(1:10) == '&END MGRID') exit
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(found) then
  write(buf,'(A,I0,A)') '(',j-k,'X,A,I0)'
  write(fid1,fmt=TRIM(buf)) 'NGRIDS ', ngrids
 else
  write(buf,'(A,I0,A)') '(',j-k+1,'X,A,I0)'
  write(fid1,fmt=TRIM(buf)) 'NGRIDS ', ngrids
  write(buf,'(A,I0,A)') '(',j-k,'X,A)'
  write(fid1,fmt=TRIM(buf)) '&END MGRID'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_mgrid_ngrids_in_cp2k_inp

! set/modify CUTOFF in &FORCE_EVAL/&DFT/&MGRID
subroutine set_mgrid_cutoff_in_cp2k_inp(inpname, cutoff)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: cutoff
!f2py intent(in) :: cutoff
 character(len=50), parameter :: error_warn = 'ERROR in subroutine set_mgrid_cu&
                                              &toff_in_cp2k_inp: '
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 if(cutoff<100 .or. cutoff>3000) then
  write(6,'(/,A,I0)') error_warn//'invalid CUTOFF=', cutoff
  write(6,'(A)') 'The allowed range is [100,3000].'
  stop
 end if

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:6) == '&MGRID') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&MGRID" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:6) == 'CUTOFF') exit
   if(buf(1:10) == '&END MGRID') then
    i = -1; exit
   end if
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"CUTOFF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(buf,'(A,I0,A)') '(',j-k,'X,A,I0)'
 write(fid1,fmt=TRIM(buf)) 'CUTOFF ', cutoff

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_mgrid_cutoff_in_cp2k_inp

! set/modify REL_CUTOFF in &FORCE_EVAL/&DFT/&MGRID
subroutine set_mgrid_rel_cutoff_in_cp2k_inp(inpname, rel_cutoff)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: rel_cutoff
!f2py intent(in) :: rel_cutoff
 character(len=54), parameter :: error_warn = 'ERROR in subroutine set_mgrid_re&
                                              &l_cutoff_in_cp2k_inp: '
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 if(rel_cutoff<10 .or. rel_cutoff>200) then
  write(6,'(/,A,I0)') error_warn//'invalid REL_CUTOFF=', rel_cutoff
  write(6,'(A)') 'The allowed range is [10,200].'
  stop
 end if

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   if(buf(1:6) == '&MGRID') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&MGRID" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf0
  if(i /= 0) exit
  j = LEN_TRIM(buf0)
  if(j > 0) then
   buf = ADJUSTL(buf0)
   k = LEN_TRIM(buf)
   if(buf(1:10) == 'REL_CUTOFF') exit
   if(buf(1:10) == '&END MGRID') then
    i = -1; exit
   end if
  end if
  write(fid1,'(A)') TRIM(buf0)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"REL_CUTOFF" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(buf,'(A,I0,A)') '(',j-k,'X,A,I0)'
 write(fid1,fmt=TRIM(buf)) 'REL_CUTOFF ', rel_cutoff

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_mgrid_rel_cutoff_in_cp2k_inp

! set finer grid in &FORCE_EVAL/&DFT/&XC/&XC_GRID in a specified CP2K .inp file
subroutine set_finer_grid_in_cp2k_inp(inpname, finer)
 implicit none
 integer :: i, j, k, n_indent, fid, fid1, RENAME
 character(len=1) :: str1
 character(len=48), parameter :: error_warn = 'ERROR in subroutine set_finer_gr&
                                              &id_in_cp2k_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname
 logical :: found1, found2, alive1, alive2
 logical, intent(in) :: finer
!f2py intent(in) :: finer

 str1 = 'F'
 if(finer) str1 = 'T'

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  k = LEN_TRIM(buf)
  if(k > 0) then
   buf = ADJUSTL(buf)
   k = LEN_TRIM(buf)
   if(buf(1:3) == '&XC') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&XC" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 found1 = .false.; found2 = .false.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  j = LEN_TRIM(buf)
  if(j > 0) then
   buf = ADJUSTL(buf)
   k = LEN_TRIM(buf)
   if(buf(1:8) == '&XC_GRID') found1 = .true.
   if(buf(1:14) == 'USE_FINER_GRID') found2 = .true.
   alive1 = (k == 7)
   alive2 = (k>7 .and. (buf(8:8)==' ' .or. buf(8:8)=='#'))
   if((buf(1:7)=='&END XC') .and. (alive1 .or. alive2)) exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&END XC" cannot be found in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 if(found1) then ! &XC_GRID exists
  if(found2) then ! &XC_GRID/USE_FINER_GRID both exist
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    j = LEN_TRIM(buf)
    if(j > 0) then
     buf = ADJUSTL(buf)
     if(buf(1:14) == 'USE_FINER_GRID') exit
    end if
   end do ! for while
   do while(.true.)
    BACKSPACE(fid1)
    BACKSPACE(fid1)
    read(fid1,'(A)') buf
    j = LEN_TRIM(buf)
    if(j > 0) then
     buf = ADJUSTL(buf)
     k = LEN_TRIM(buf)
     n_indent = j - k
     if(buf(1:14) == 'USE_FINER_GRID') exit
    end if
   end do ! for while
   BACKSPACE(fid1)
   write(buf,'(A,I0,A)') '(',n_indent,'X,A)'
   write(fid1,fmt=TRIM(buf)) 'USE_FINER_GRID '//str1
  else            ! &XC_GRID exists, but USE_FINER_GRID does not exist
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    j = LEN_TRIM(buf)
    if(j > 0) then
     buf = ADJUSTL(buf)
     k = LEN_TRIM(buf)
     if(k==12 .and. buf(1:12)=='&END XC_GRID') exit
     alive1 = (k == 3)
     alive2 = (k>3 .and. (buf(4:4)==' ' .or. buf(4:4)=='#'))
     if(buf(1:3)=='&XC' .and. (alive1 .or. alive2)) then
      write(6,'(/,A)') error_warn//'"&END XC_GRID" not found in'
      write(6,'(A)') 'file '//TRIM(inpname)
      close(fid)
      close(fid1,status='delete')
      stop
     end if
    end if
   end do ! for while
   do while(.true.)
    BACKSPACE(fid1)
    BACKSPACE(fid1)
    read(fid1,'(A)') buf
    j = LEN_TRIM(buf)
    if(j > 0) then
     buf = ADJUSTL(buf)
     k = LEN_TRIM(buf)
     n_indent = j - k
     if(buf(1:12) == '&END XC_GRID') exit
    end if
   end do ! for while
   BACKSPACE(fid1)
   write(buf,'(2(A,I0),A)') '(',n_indent+1,'X,A,/,',n_indent,'X,A)'
   write(fid1,fmt=TRIM(buf)) 'USE_FINER_GRID '//str1, '&END XC_GRID'
  end if
 else            ! &XC_GRID does not exist
  if(found2) then ! &XC_GRID does not exist, but USE_FINER_GRID exists
   write(6,'(/,A)') error_warn//'"USE_FINER_GRID" exists but'
   write(6,'(A)') '"&XC_GRID" does not exist, which seems impossible. Please ch&
                  &eck the file'
   write(6,'(A)') TRIM(inpname)
   close(fid)
   close(fid1,status='delete')
   stop
  else            ! neither of &XC_GRID/USE_FINER_GRID exists
   BACKSPACE(fid1)
   n_indent = j - k
   write(buf,'(2(A,I0),A)') '(',n_indent+1,'X,A,/,',n_indent+2,'X,A)'
   write(fid1,fmt=TRIM(buf)) '&XC_GRID', 'USE_FINER_GRID '//str1
   write(buf,'(2(A,I0),A)') '(',n_indent+1,'X,A,/,',n_indent,'X,A)'
   write(fid1,fmt=TRIM(buf)) '&END XC_GRID', '&END XC'
  end if
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine set_finer_grid_in_cp2k_inp

! increase/update the integer appears in the CP2K input filename
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

! scan convergence parameters in a specified CP2K input file (.inp)
! (EPS_DEFAULT, EPS_SCF, CUTOFF, REL_CUTOFF, NGRIDS, USE_FINER_GRID)
subroutine scan_cp2k_conv_param(inpname, param_str, e_diff, nproc, conv_param)
 implicit none
 integer :: i, k, icase
 integer, parameter :: param0(6) = [10,5,400,40,4,0]
 integer, parameter :: stepsize(6) = [1,1,200,20,1,1]
 integer, parameter :: max_nstep(6) = [4,4,15,4,3,2]
 integer, intent(in) :: nproc
!f2py intent(in) :: nproc
 integer, intent(out) :: conv_param
!f2py intent(out) :: conv_param
 real(kind=8) :: e0, e1
 real(kind=8), intent(in) :: e_diff
!f2py intent(in) :: e_diff
 character(len=42), parameter :: error_warn='ERROR in subroutine scan_cp2k_conv&
                                            &_param: '
 character(len=240) :: old_inp, new_inp, outname, wfn_file
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname
 character(len=*), intent(in) :: param_str
!f2py intent(in) :: param_str
 logical :: converged

 converged = .false.; e0 = 0d0; old_inp = inpname
 call find_specified_suffix(inpname, '.in', i)
 wfn_file = inpname(1:i-1)//'-RESTART.wfn'

 icase = 0
 select case(TRIM(param_str))
 case('EPS_DEFAULT')
  icase = 1
 case('EPS_SCF')
  icase = 2
 case('CUTOFF')
  icase = 3
 case('REL_CUTOFF')
  icase = 4
 case('NGRIDS')
  icase = 5
 case('FINER_GRID','USE_FINER_GRID')
  icase = 6
 case default
  write(6,'(/,A)') error_warn//' invalid param_str="'//TRIM(param_str)//'"'
  write(6,'(A)') 'Currently supported keywords: "EPS_DEFAULT", "EPS_SCF", "CUTO&
                 &FF", "REL_CUTOFF"'
  write(6,'(A)') '"NGRIDS", "USE_FINER_GRID".'
  stop
 end select

 do i = 1, max_nstep(icase), 1
  call increase_int_in_cp2k_inpname(old_inp, new_inp)
  call find_specified_suffix(new_inp, '.in', k)
  outname = new_inp(1:k-1)//'.out'
  call sys_copy_file(TRIM(old_inp), TRIM(new_inp), .false.)
  if(i == 1) call set_wfn_restart_file_in_cp2k_inp(new_inp, wfn_file)
  conv_param = param0(icase) + (i-1)*stepsize(icase)
  select case(icase)
  case(1)
   call set_eps_default_in_cp2k_inp(new_inp, conv_param)
  case(2)
   call set_eps_scf_in_cp2k_inp(new_inp, conv_param)
  case(3)
   call set_mgrid_cutoff_in_cp2k_inp(new_inp, conv_param)
  case(4)
   call set_mgrid_rel_cutoff_in_cp2k_inp(new_inp, conv_param)
  case(5)
   call set_mgrid_ngrids_in_cp2k_inp(new_inp, conv_param)
  case(6)
   call set_finer_grid_in_cp2k_inp(new_inp, (conv_param==1))
  end select
  call submit_cp2k_job(new_inp, nproc)
  call read_scf_e_from_cp2k_out(outname, e1)
  if(i>1 .and. DABS(e0-e1)<e_diff) then
   converged = .true.
   exit
  end if
  old_inp = new_inp
  e0 = e1
 end do ! for i

 if(.not. converged) then
  write(6,'(/,A)') error_warn//'scanning "'//TRIM(param_str)//'"'
  write(6,'(A)') 'does not lead to a converged result. Please check files '&
               //TRIM(inpname)
  write(6,'(A)') 'and '//TRIM(outname)
  stop
 end if

 ! use the previous convergence parameter
 conv_param = param0(icase) + (i-2)*stepsize(icase)

 select case(icase)
 case(1,2)
  write(6,'(A,I0)') 'Recommended '//TRIM(param_str)//' 1e-', conv_param
 case(3,4,5)
  write(6,'(A,I0)') 'Recommended '//TRIM(param_str)//' ', conv_param
 case(6)
  write(6,'(A,L1)') 'Recommended USE_FINER_GRID ', (conv_param==1)
 end select
end subroutine scan_cp2k_conv_param

