*[中文版说明请点击README_zh.md](README_zh.md)*
# Molecular Orbital KIT (MOKIT)
MOKIT offers various utilities and modules to transfer MOs among various quantum
chemistry software packages. Besides, the automr program in MOKIT can set up and
run common multi-reference calculations in a block-box way.

A list of important utilities along with their functions are shown below
![MOKIT utilities with their functions](doc/orbital_transfer.png)

With MOKIT, one can perform multi-reference calculations in a quite simple way,
and utilize the best modules of each program. E.g.

  UHF(UNO) -> CASSCF -> CASPT2  
  Gaussian&emsp;&emsp;PySCF&emsp;&emsp;OpenMolcas  
or  
  UHF(UNO) -> GVB   -> CASSCF -> NEVPT2  
  Gaussian&emsp;&emsp;GAMESS&emsp;&emsp;PySCF&emsp;&emsp;PySCF  
or   
  RHF      -> GVB   -> CASSCF -> ic-MRCISD+Q  
  Gaussian&emsp;&emsp;GAMESS&emsp;&emsp;PySCF&emsp;&emsp;OpenMolcas

Negligible energy loss(usually<1e-6 a.u., for the same wave function method in two
programs) are ensured during transferring MOs, since the basis order of angular
momentum up to H(i.e. l=5) are considered.

Pre-built `Windows* OS` executables of 20 utilities are provided in [Releases](https://gitlab.com/jxzou/mokit/-/releases).

Note that although MOKIT aims to make the multi-reference calculations block-box,
the users are still required to have practical experiences of quantum chemistry
computations (e.g. familiar with routine DFT calculations in Gaussian). You are
encouraged to learn how to use Gaussian if you are a fresh hand.

Feb 20, 2021

Installation
------------

* Prerequisites
    - Fortran compiler: ifort(>=2015) or gfortran(>=4.8.5 and <8.0)
    - Intel MKL (installing Intel compilers recommended)
    - f2py (installing Anaconda Python3 recommended)

* Compile all modules

        cd src
        make all

* Compile individual utility or module  
  E.g.

        make fch2inp

* After 'make all', you need to set environment variables `MOKIT_ROOT`, `PATH` and `PYTHONPATH`.  
  E.g. if MOKIT is installed in /home/$USER/software/mokit, the following should be set in ~/.bashrc:

        export MOKIT_ROOT=/home/$USER/software/mokit
        export PATH=$MOKIT_ROOT/bin:$PATH
        export PYTHONPATH=$MOKIT_ROOT/lib:$PYTHONPATH
        export ORCA=/opt/orca-4.2.1/orca
        export GMS=/home/$USER/software/gamess/rungms

  (Remember to modify the `ORCA` and `GMS` paths to suit your local environment)

* The original GAMESS code can only deal with GVB <=12 pairs. But nowadays we
  can do hundreds of pairs. To go beyond 12 pairs, please read Section 4.4.9 in
  [manual](doc/).

Quick Start
-----------
* Each utility is self-explanatory. For example, run `fch2inp` in Shell,
  you will find

   ERROR in subroutine fch2inp: wrong command line arguments!  
   Example 1 (R(O)HF, CAS): fch2inp a.fch  
   Example 2 (UHF)        : fch2inp a.fch -uhf  
   Example 3 (GVB)        : fch2inp a.fch -gvb [npair]  
   Example 4 (ROGVB)      : fch2inp a.fch -gvb [npair] -open [nopen]

* For usages of modules in lib/, see [examples/utilities/readme.txt](examples/utilities/readme.txt)

* The input syntax of the automr program is like Gaussian gjf. For example, the input
  file '00-h2o_cc-pVDZ_1.5.gjf' of the water molecule at d(O-H) = 1.5 A is shown below

        %mem=4GB
        %nprocshared=4
        #p CASSCF/cc-pVDZ
        
        mokit{}
        
        0 1
        O      -0.23497692    0.90193619   -0.068688
        H       1.26502308    0.90193619   -0.068688
        H      -0.73568721    2.31589843   -0.068688

  Run

        automr 00-h2o_cc-pVDZ_1.5.gjf >& 00-h2o_cc-pVDZ_1.5.out

  in Shell. The automr program will successively perform HF, GVB, and CASSCF
  computations. See [examples/](examples/) for more examples.

QC Packages can be called by AutoMR
----------
* [Gaussian](http://gaussian.com/)
* [PySCF](https://github.com/pyscf/pyscf)
* [GAMESS](https://www.msg.chem.iastate.edu/gamess/index.html)
* [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)
* [Molpro](https://www.molpro.net/)
* [ORCA](https://orcaforum.kofo.mpg.de)
* [BDF](http://182.92.69.169:7226/Introduction)
* [PSI4](https://github.com/psi4/psi4/)

Some Tips
---------
* To avoid unnecessary errors, you must specify keywords 'nosymm int=nobasistransform'
  in Gaussian .gjf file, if you want to provide a .fch(k) file to AutoMR.

Bug Report
----------
* If you find any bug frequently occurs, please contact the author jxzou via njumath[at]sina.cn,
  with your .fch(k) file and output files attached. The author may not answer or update code
  frequently since he is being postponed due to his PhD program.

* You can also open an issue on the [Issues](https://gitlab.com/jxzou/mokit/-/issues) page.

TODO
----
* MOs trasferring among BAGEL, NWCHEM, Dalton, etc.

* Support more multireference methods like ic-MRCC

* Develop/Implement black-box strategies of excited state calculations

Citation
--------
* If you use MOKIT in your work, please cite the gitlab page of MOKIT as

   Jingxiang Zou, Molecular Orbital Kit (MOKIT), https://gitlab.com/jxzou/mokit (accessed month day, year)

* If you use MOKIT to peform calculations involving GVB, please also cite the following two papers

   DOI: 10.1021/acs.jctc.8b00854; DOI: 10.1021/acs.jpca.0c05216.

* Please read the [manual](doc/) for details of citation. Your proper citation
  would be a great encouragement to developer.

Disclaimer
----------
Copyright (c) 2021 jxzou

All rights reserved.

Redistribution and use in source and binary forms are permitted provided that the above copyright notice and this paragraph are duplicated in all such forms and that any documentation, advertising materials, and other materials related to such distribution and use acknowledge that the software was developed by the author. The name of the authors may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

