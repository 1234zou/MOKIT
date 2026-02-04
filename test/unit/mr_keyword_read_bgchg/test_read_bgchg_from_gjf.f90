program test_read_bgchg_from_gjf
 use mol, only: natom, nuc, coor, nbgchg, bgcharge, ptchg_e, nuc_pt_e
 use mr_keyword, only: chgname, rigid_scan, relaxed_scan
 use phys_cons, only: Bohr_const
 use test_utils, only: assert_int_equal, assert_real_equal, assert_real_close, &
  assert_path_suffix, assert_string_equal, fail
 implicit none
 interface
  subroutine read_bgchg_from_gjf(gjfname, no_coor)
   character(len=240), intent(in) :: gjfname
   logical, intent(in) :: no_coor
  end subroutine read_bgchg_from_gjf
 end interface
 real(kind=8), parameter :: tol = 1d-10
 real(kind=8) :: expected_ptchg, expected_nuc_pt
 real(kind=8) :: r12, r1, r2
 logical :: exists

 call test_case_simple('mr_keyword_read_bgchg/read_bgchg_from_gjf.gjf')
 call test_case_second('mr_keyword_read_bgchg/read_bgchg_from_gjf2.gjf')

 write(6,'(A)') 'PASS: read_bgchg_from_gjf'
contains
 subroutine reset_state()
  natom = 1
  if(allocated(nuc)) deallocate(nuc)
  if(allocated(coor)) deallocate(coor)
  if(allocated(bgcharge)) deallocate(bgcharge)
  allocate(nuc(natom))
  allocate(coor(3,natom))
  nuc(1) = 1
  coor(:,1) = 0d0
  nbgchg = 0
  ptchg_e = 0d0
  nuc_pt_e = 0d0
 end subroutine reset_state

 subroutine test_case_simple(fname)
  character(len=*), intent(in) :: fname
  character(len=240) :: gjf
  call reset_state()
  rigid_scan = .false.
  relaxed_scan = .false.

  gjf = ' '
  gjf = fname
  call read_bgchg_from_gjf(gjf, .false.)
  call assert_string_equal('chgname', TRIM(chgname), 'mr_keyword_read_bgchg/read_bgchg_from_gjf.chg')

  call assert_int_equal('nbgchg', nbgchg, 2)
  call assert_real_equal('bgcharge(1,1)', bgcharge(1,1), 0d0)
  call assert_real_equal('bgcharge(2,1)', bgcharge(2,1), 0d0)
  call assert_real_equal('bgcharge(3,1)', bgcharge(3,1), 1d0)
  call assert_real_equal('bgcharge(4,1)', bgcharge(4,1), 0.5d0)
  call assert_real_equal('bgcharge(1,2)', bgcharge(1,2), 0d0)
  call assert_real_equal('bgcharge(2,2)', bgcharge(2,2), 0d0)
  call assert_real_equal('bgcharge(3,2)', bgcharge(3,2), 2d0)
  call assert_real_equal('bgcharge(4,2)', bgcharge(4,2), -0.5d0)

  r12 = 1d0 / Bohr_const
  expected_ptchg = (0.5d0 * (-0.5d0)) / r12
  call assert_real_close('ptchg_e', ptchg_e, expected_ptchg, tol)

  r1 = 1d0 / Bohr_const
  r2 = 2d0 / Bohr_const
  expected_nuc_pt = 0.5d0 / r1 + (-0.5d0) / r2
  call assert_real_close('nuc_pt_e', nuc_pt_e, expected_nuc_pt, tol)

  call assert_path_suffix('chgname', chgname, '.chg')
  inquire(file=trim(chgname), exist=exists)
  if(.not. exists) call fail('chg file not created: '//trim(chgname))
 end subroutine test_case_simple

 subroutine test_case_second(fname)
  character(len=*), intent(in) :: fname
  character(len=240) :: gjf
  call reset_state()
  rigid_scan = .false.
  relaxed_scan = .false.

  gjf = ' '
  gjf = fname
  call read_bgchg_from_gjf(gjf, .false.)
  call assert_string_equal('chgname', TRIM(chgname), 'mr_keyword_read_bgchg/read_bgchg_from_gjf2.chg')

  call assert_int_equal('nbgchg', nbgchg, 2)
  call assert_real_equal('bgcharge(1,1)', bgcharge(1,1), 0d0)
  call assert_real_equal('bgcharge(2,1)', bgcharge(2,1), -8.4035d0)
  call assert_real_equal('bgcharge(3,1)', bgcharge(3,1), 0d0)
  call assert_real_equal('bgcharge(4,1)', bgcharge(4,1), -1.125d0)
  call assert_real_equal('bgcharge(1,2)', bgcharge(1,2), 0d0)
  call assert_real_equal('bgcharge(2,2)', bgcharge(2,2), 8.4035d0)
  call assert_real_equal('bgcharge(3,2)', bgcharge(3,2), 0d0)
  call assert_real_equal('bgcharge(4,2)', bgcharge(4,2), -1.125d0)

  r12 = (2d0 * 8.4035d0) / Bohr_const
  expected_ptchg = ((-1.125d0) * (-1.125d0)) / r12
  call assert_real_close('ptchg_e', ptchg_e, expected_ptchg, tol)

  r1 = 8.4035d0 / Bohr_const
  expected_nuc_pt = (-1.125d0) / r1 + (-1.125d0) / r1
  call assert_real_close('nuc_pt_e', nuc_pt_e, expected_nuc_pt, tol)

  call assert_path_suffix('chgname', chgname, '.chg')
  inquire(file=trim(chgname), exist=exists)
  if(.not. exists) call fail('chg file not created: '//trim(chgname))
 end subroutine test_case_second

end program test_read_bgchg_from_gjf
