---------------------------------------------------
Example 1: transfer RHF MOs from Gaussian to GAMESS
---------------------------------------------------
Run 'fch2inp 01-water_cc-pVTZ_6D10F.fch' in Shell, a file 01-water_cc-pVTZ_6D10F.inp would be generated.

(Usually the keywords in .inp file are already appropriate, but do check if advanced options are desired)

---------------------------------------------------
Example 2: transfer UHF MOs from Gaussian to GAMESS
---------------------------------------------------
Run 'fch2inp 02-bis_6D10F.fchk -uhf' in Shell, a file 02-bis_6D10F.inp would be generated.

----------------------------------------------------------------
Example 2.1: transfer UHF MOs from Gaussian to MOLCAS/OpenMolcas
----------------------------------------------------------------
Run 'fch2inporb 02-bis_6D10F.fchk' in Shell, two files(.input and .INPORB) would be generated.

-------------------------------------------------------
Example 3: transfer paired LMOs from Gaussian to GAMESS
-------------------------------------------------------
Run 'fch2inp 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.fch -gvb 3' in Shell, the .inp would be generated.
You can run
'gms 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.inp 01 8 >& 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.gms'
to submit this job. This is a GVB(3) job of GAMESS with paired LMOs as initial guess.

-----------------------------------------------------
Example 3.1: transfer GVB NOs from GAMESS to Gaussian
-----------------------------------------------------
After Example 3 finished, you can run the following 2 steps
(1) cp 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.fch 03-ben_cc-pVDZ_6D10F_gvb3.fch
(2) dat2fch 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.dat 03-ben_cc-pVDZ_6D10F_gvb3.fch -gvb 3
The GVB(3) orbitals are in 03-ben_cc-pVDZ_6D10F_gvb3.fch file, which can be visualizied by GaussView.

-------------------------------------------------
Example 4: transfer RHF MOs from Gaussian to ORCA
-------------------------------------------------
Run 'fch2mkl 04-PhCOOH_def2TZVP_rhf.fchk' in Shell, two files(.mkl and .inp) would be generated.
Then run 'orca_2mkl 04-PhCOOH_def2TZVP_rhf -gbw' in Shell, you can obtain 04-PhCOOH_def2TZVP_rhf.gbw file.
Finally run 'orca 04-PhCOOH_def2TZVP_rhf.inp >& 04-PhCOOH_def2TZVP_rhf.out', you will find that the RHF energy
can be exactly reproduced using the utility fch2mkl.

--------------------------------------------------
Example 5: transfer GVB NOs from Gaussian to PySCF
--------------------------------------------------
In Example 3.1, you already got the file 03-ben_cc-pVDZ_6D10F_gvb3.fch, now run
'python 05-ben_cc-pVDZ_6D10F_gvb3_2CAS66.py >& 05-ben_cc-pVDZ_6D10F_gvb3_2CAS66.out'
The GVB NOs will be read and used as the initial guess of CAS(6,6). After CASSCF finished,
CASSCF NOs will be exported.

(Note: the basis sets in .py file are obtained by running
'fch2inp 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.fch && bas_gms2py 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.inp'
in Shell. Do not use the built-in basis sets of PySCF.)

--------------------------------------------------
Example 6: transfer UHF MOs from Gaussian to PySCF
--------------------------------------------------
Run 'python 06-bis_6D10F.py' in Shell. The UHF canonical MOs in 02-bis_6D10F.fchk will be read in.
One will find that the UHF energy immediately converges after little change.

