! written by jxzou at 20231109: partition a large molecule into monomers
!  and generate initial MOs as well as amplitudes from fragement calculations

!program frag_qc
! implicit none
! integer :: i, wfn_type
! integer, parameter :: n = 7
! integer :: wfn_type0(n)
! character(len=240) :: fchname0(n), fchname
! character(len=24) :: data_string
! logical :: pos(n)
!
! wfn_type0 = [1,3,3,3,3,3,3]
! pos = [.true.,.true.,.false.,.true.,.false.,.true.,.false.]
! wfn_type = 1
! fchname0(1) = 'Cu_Imidazole-2_old.fch'
! do i = 2, n, 1
!  fchname0(i) = 'Cu_atom_6-31Gdp.fch'
! end do ! for i
! fchname = 'Cu_Imidazole-2_ds_dm.fch'
!
! call fdate(data_string)
! write(6,'(A)') TRIM(data_string)
!
!! call direct_sum_frag_fock_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
!! call direct_sum_frag_mo_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
! call direct_sum_frag_dm_in_fch(n, fchname0, fchname)
!
! call fdate(data_string)
! write(6,'(/,A)') TRIM(data_string)
!end program frag_qc

program frag_qc
 implicit none
 integer :: wfn_type
 integer, parameter :: n = 2
 integer :: wfn_type0(n)
 character(len=240) :: fchname0(n), fchname
 character(len=24) :: data_string
 logical :: pos(n)

 wfn_type0 = [1,1]
 pos = [.true.,.true.]
 wfn_type = 1
 fchname0 = ['Cu_Imidazole-1_o.fch','Cu_Imidazole-2_UPBE_UwB97X_UwB97M-V_r.fch']
 fchname = 'Cu_Imidazole-3_only_r.fch'

 call fdate(data_string)
 write(6,'(A)') TRIM(data_string)

! call direct_sum_frag_fock_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
 call direct_sum_frag_mo_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
! call direct_sum_frag_dm_in_fch(n, fchname0, fchname)

 call fdate(data_string)
 write(6,'(/,A)') TRIM(data_string)
end program frag_qc

