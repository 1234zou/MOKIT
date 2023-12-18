! written by jxzou at 20231109: partition a large molecule into monomers
!  and generate initial MOs as well as amplitudes from fragement calculations

program frag_qc
 implicit none
 integer :: wfn_type
 integer, parameter :: n = 3
 integer :: wfn_type0(n)
 character(len=240) :: fchname0(3), fchname
 character(len=24) :: data_string
 logical :: pos(n)

 wfn_type0 = 1; wfn_type=1; pos = .true.
 fchname0 = ['DMC_trimer1_rhf.fch','DMC_trimer2_rhf.fch','DMC_trimer3_rhf.fch']
 fchname = 'DMC_trimer_rhf_conver6.fch'

 call fdate(data_string)
 write(6,'(A)') TRIM(data_string)

! call direct_sum_frag_fock_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
 call direct_sum_frag_mo_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)

 call fdate(data_string)
 write(6,'(/,A)') TRIM(data_string)
end program frag_qc
