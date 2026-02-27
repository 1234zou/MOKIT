*[中文说明](README_zh.md)*
*&emsp;&emsp;[Download the program](https://gitlab.com/jxzou/mokit/-/archive/master/mokit-master.zip)*
*&emsp;&emsp;[GitLab repo](https://gitlab.com/jxzou/MOKIT)*
*&emsp;&emsp;[GitHub mirror](https://github.com/1234zou/MOKIT)*
*&emsp;&emsp;[NJU-Git mirror](https://git.nju.edu.cn/jxzou/mokit)*
*&emsp;&emsp;[Documentation](https://doc.mokit.xyz/)*

# Molecular Orbital KIT (MOKIT)
MOKIT offers various utilities and modules to transform MOs among various quantum chemistry software packages. Besides, the `automr` program in MOKIT can set up and run multi-reference quantum chemistry calculations in a black-box way.

A list of important utilities along with their functions are shown below  
![MOKIT utilities with their functions](doc/orbital_transfer.png)

With MOKIT, one can perform multi-reference calculations in a quite simple way, and utilize the best modules of each program. E.g.

  UHF(UNO) -> CASSCF -> CASPT2  
  Gaussian&emsp;&emsp;PySCF&emsp;&emsp;OpenMolcas  
or  
  UHF(UNO) -> GVB   -> CASSCF -> NEVPT2  
  Gaussian&emsp;&emsp;GAMESS&emsp;&emsp;PySCF&emsp;&emsp;PySCF  
or   
  RHF      -> GVB   -> CASSCF -> ic-MRCISD+Q  
  Gaussian&emsp;GAMESS&emsp;PySCF&emsp;OpenMolcas

Negligible energy loss(usually<1e-6 a.u., for the same wave function method in two programs) are ensured during transferring MOs, since the basis order of angular momentum up to H(i.e. *l*=5) are considered.

Note that although MOKIT aims to make the multi-reference calculations black-box, the users are still required to have practical experiences of quantum chemistry computations, e.g. familiar with routine DFT calculations in [Gaussian](https://gaussian.com). You are encouraged to learn how to use Gaussian if you are a fresh hand.

Feb 27, 2026

Dependencies
------------
Dependencies on quantum chemistry packages are different for each executable or module. Here the minimum requirements for binary executables `automr`, `frag_guess_wfn` and Python modules `py2xxx` are listed:
1. automr: GAMESS, PySCF
2. frag_guess_wfn: Gaussian
3. Most of the utilities do not depend on quantum chemistry packages except that the modules `py2gau`, `py2orca`, `py2molpro`, etc, work with PySCF installed.

> Note: automr will not call the original GAMESS executable `gamess.00.x`, but call `gamess.01.x`. To generate it, please read Section [4.4.10](https://doc.mokit.xyz/chap4-4.html#4410-gvb_prog) in MOKIT manual and use the provided Shell script to automatically modify GAMESS code and re-compile.

Installation
------------
You can choose one of the four options shown below to install MOKIT on Linux or MacOS, and they are for full functionalities. If you only want the utility `frag_guess_wfn` or some other utility like `fch2mkl`, see [here](https://doc.mokit.xyz/chap2-2.html#223-only-want-frag_guess_wfn) for an even easier way to install.
Pre-built `Windows OS` executables for a few utilities are provided in [Releases](https://gitlab.com/jxzou/mokit/-/releases). But they are outdated compared with master branch, and there's no way to use full functionality on Windows.

### Option 1: Install from conda (for Linux and MacOS)
This is the easiest way, but network is required to auto-download the requirements. Creating a new environment before installing is highly recommended, to avoid changing your base environment. You can create the environment and install in one go like
```
conda create -n mokit-py311 python=3.11 mokit -c mokit -c conda-forge
conda activate mokit-py311
```
For Linux x86-64, you can use any version of Python 3.9-3.12. But for MacOS arm64, only 3.11 is available. MacOS users can also [use the homebrew-toolchains](https://doc.mokit.xyz/chap2-2.html#option-2-use-homebrew-toolchains-for-macos-only) to install MOKIT.

For more details about conda channels and how to update/uninstall MOKIT using conda, please read [here](https://doc.mokit.xyz/chap2-2.html#option-1-install-from-conda-for-linux-and-macos). You need to activate the environment `mokit-py311` before using MOKIT, and you can run `conda deactivate` to exit the environment when you do not use it. Please read [more details](https://doc.mokit.xyz/chap2-4.html) to install and use MOKIT on a Cluster(集群). If you have no access to network, but still don't want to compile the source code, you can try Option 2 below.

### Option 2: Use Pre-compiled MOKIT

Pre-built `Linux` executables can be downloaded [here](https://doc.mokit.xyz/chap2-2.html#222-pre-built-linux-executables-and-libraries).

* Prerequisites: Python3 environment and NumPy.

* A detailed guide for choosing the version of pre-built artifacts and resolving dependencies can be found [here](https://doc.mokit.xyz/chap2-2.html#222-pre-built-linux-executables-and-libraries)

* After downloading the pre-built artifacts, you need to set the following environment variables (assuming MOKIT is put in `$HOME/software/mokit`) in your `~/.bashrc`:

```bash
export MOKIT_ROOT=$HOME/software/mokit
export PATH=$MOKIT_ROOT/bin:$PATH
export PYTHONPATH=$MOKIT_ROOT:$PYTHONPATH
export LD_LIBRARY_PATH=$MOKIT_ROOT/mokit/lib:$LD_LIBRARY_PATH
export GMS=$HOME/software/gamess/rungms
```
The `LD_LIBRARY_PATH` is needed since the OpenBLAS dynamic library is put there. Remember to modify the `GMS` path to suit your local environment. Note that you need to exit the terminal and re-login, in order to activate newly written environment variables.

### Option 3: Build from Source
The link to latest version of MOKIT source code can be found [here](https://doc.mokit.xyz/chap2-3.html).

* Prerequisites
    - Fortran compiler: `ifort`(>=2017) or `gfortran`(>=4.8.5)
    - Intel MKL(recommended) or [OpenBLAS](https://github.com/xianyi/OpenBLAS)
    - f2py (installing Anaconda Python3 recommended)

* Compile all modules
```
cd src
make all
cd ..
pip install -e . --prefix=.
```

* After `make all`, you need to set environment variables `MOKIT_ROOT`, `PATH` and `PYTHONPATH`. E.g. if MOKIT is installed in `$HOME/software/mokit`, the following should be set in `~/.bashrc`:

```bash
export MOKIT_ROOT=$HOME/software/mokit
export PATH=$MOKIT_ROOT/bin:$PATH
export PYTHONPATH=$MOKIT_ROOT:$PYTHONPATH
export GMS=$HOME/software/gamess/rungms
```
Remember to modify the `GMS` path to suit your local environment. Note that you need to exit the terminal and re-login, in order to activate newly written environment variables.

Quick Start
-----------
* Each utility is self-explanatory. For example, run `fch2inp` in Shell, you will find
```
 ERROR in program fch2inp: wrong command line arguments!
 Example 1 (R(O)HF/UHF/CAS): fch2inp h2o.fch
 Example 2 (SF-CIS)        : fch2inp high_spin.fch -sfcis
 Example 3 (SF-TDDFT)      : fch2inp high_spin.fch -sf
 Example 4 (MRSF-CIS)      : fch2inp triplet.fch -mrsfcis
 Example 5 (MRSF-TDDFT)    : fch2inp triplet.fch -mrsf
 Example 6 (GVB)           : fch2inp h2o.fch -gvb [Npair]
```
You can search a utility and read their documentations [here](https://doc.mokit.xyz/chap4-5.html).

* For usages of modules in mokit/lib/, see [examples/utilities/readme.txt](examples/utilities/readme.txt)

* The input syntax of the `automr` program is just Gaussian input file (.gjf). For example, the input file `00-h2o_cc-pVDZ_1.5.gjf` of the water molecule at d(O-H) = 1.5 A is shown below
```
%mem=4GB
%nprocshared=4
#p CASSCF/cc-pVDZ

mokit{}

0 1
O      -0.23497692    0.90193619   -0.068688
H       1.26502308    0.90193619   -0.068688
H      -0.73568721    2.31589843   -0.068688
```

Run the following commands
```
automr 00-h2o_cc-pVDZ_1.5.gjf >00-h2o_cc-pVDZ_1.5.out 2>&1
```

The `automr` program will successively perform HF, GVB, and CASSCF computations. The active space will be automatically determined as (4,4) during computations. See [Quick Start at documentation](https://doc.mokit.xyz/chap3_quick.html) and [User Guide](https://doc.mokit.xyz/chap4_guide.html) for more information. See [examples/](examples/) for more examples.


MOKIT is able to transform MOs among these QC Packages
----------
* [Amesp](https://amesp.xyz)
* [BDF](https://bdf-manual.readthedocs.io/en/latest/index.html)
* [CFOUR](https://cfour.uni-mainz.de/cfour/index.php?n=Main.HomePage)
* [CP2K](https://github.com/cp2k/cp2k)
* [Dalton](https://gitlab.com/dalton/dalton)
* [GAMESS](https://www.msg.chem.iastate.edu/gamess/index.html)
* [Gaussian](http://gaussian.com)
* [Molpro](https://www.molpro.net)
* [MRCC](https://mrcc.hu/index.php)
* [NWChem](https://github.com/nwchemgit/nwchem)
* [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)
* [OpenQP](https://github.com/Open-Quantum-Platform/openqp)
* [ORCA](https://orcaforum.kofo.mpg.de)
* [PSI4](https://github.com/psi4/psi4)
* [PySCF](https://github.com/pyscf/pyscf)
* [Q-Chem](https://www.q-chem.com)
* [QM4D](https://www.qm4d.org)
* [REST](https://gitee.com/restgroup)
* [Turbomole](https://www.turbomole.org)


Some Tips
---------
* To avoid unnecessary errors, it is strongly recommended to specify keywords `nosymm int=nobasistransform` in Gaussian input file (.gjf), if you want to provide a .fch(k) file to `automr`.


Bug Report
----------
* If you find any bug frequently occurs or have any suggestions, you can open an issue on the [Issues](https://gitlab.com/jxzou/mokit/-/issues) page.

* You can also contact the developer jxzou via E-mail njumath[at]sina.cn, with related files (gjf, fch, out, etc) attached.


TODO
----
* MOs transferring among BAGEL, COLUMBUS, etc.

* Develop/Implement robust black-box strategies of excited state calculations

Citation
--------
* Currently we have not published the paper of MOKIT program. If you use (any module or utility of) MOKIT in your work, please cite MOKIT as

   Jingxiang Zou, Molecular Orbital Kit (MOKIT), https://gitlab.com/jxzou/mokit (accessed month day, year)

* If you use MOKIT to perform calculations involving GVB, citing the following two papers would be appreciated

   DOI: [10.1021/acs.jctc.8b00854](https://www.doi.org/10.1021/acs.jctc.8b00854); DOI: [10.1021/acs.jpca.0c05216](https://www.doi.org/10.1021/acs.jpca.0c05216).

* If you use MOKIT in your work, please cite MOKIT in the main body of your paper. Citing MOKIT only in Supporting Information of your paper is insufficient. EndNote citation files can be found [here](https://gitlab.com/jxzou/mokit/-/tree/master/doc?ref_type=heads). More details and examples of citation can be found in [manual](https://doc.mokit.xyz/chap1-2.html). 您的规范引用是对开发者的极大鼓励。您可以使用MOKIT为其他人做计算（包括代算）甚至是构建AI智能体，但务必提醒他/她在发表文章时恰当地引用MOKIT和计算中用到的量子化学软件。

* Click [here](https://doc.mokit.xyz/citing.html) to see published papers which cited MOKIT.

