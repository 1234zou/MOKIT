*[中文说明](README_zh.md)*
*&emsp;&emsp;[Download the program](https://gitlab.com/jxzou/mokit/-/archive/master/mokit-master.zip)*
*&emsp;&emsp;[GitLab repo](https://gitlab.com/jxzou/MOKIT)*
*&emsp;&emsp;[GitHub mirror](https://github.com/1234zou/MOKIT)*
*&emsp;&emsp;[NJU-Git mirror](https://git.nju.edu.cn/jxzou/mokit)*
*&emsp;&emsp;[Documentation](https://jeanwsr.gitlab.io/mokit-doc-mdbook/)*

# Molecular Orbital KIT (MOKIT)
MOKIT offers various utilities and modules to transfer MOs among various quantum chemistry software packages. Besides, the `automr` program in MOKIT can set up and run common multi-reference calculations in a black-box way.

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

Note that although MOKIT aims to make the multi-reference calculations black-box, the users are still required to have practical experiences of quantum chemistry computations (e.g. familiar with routine DFT calculations in Gaussian). You are encouraged to learn how to use Gaussian if you are a fresh hand.

May 22, 2025

Dependencies
------------
Dependencies on quantum chemistry packages are different for each executable or module. Here the minimum requirements for binary executables `automr`, `frag_guess_wfn` and Python modules `py2xxx` are listed:
1. `automr`: GAMESS, PySCF
2. `frag_guess_wfn`: Gaussian
3. Most of the utilities do not depend on quantum chemistry packages except that the modules `py2gau`, `py2orca`, `py2molpro`, etc, work with PySCF installed.

Note that the original GAMESS code can only deal with GVB <=12 pairs. But nowadays we can do hundreds of pairs. To go beyond 12 pairs, please read Section 4.4.10 in [manual](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap4-4.html#4410-gvb_prog).

Installation
------------
You can choose one of the four options shown below to install MOKIT on Linux or MacOS, and they are for full functionalities. If you only want the utility `frag_guess_wfn` or other binary utilities, like `fch2mkl`, see [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap2-2.html#223-only-want-frag_guess_wfn) for an even easier way to install.
Pre-built `Windows OS` executables for a few utilities are provided in [Releases](https://gitlab.com/jxzou/mokit/-/releases). But they are outdated compared with master branch, and there's no way to use full functionality on Windows.

### Option 1: Install from conda (for Linux only)
This is the easiest way, but network is required to auto-download the requirements (like Intel MKL). And, creating a new environment before installing is highly recommended, to avoid changing your base environment.
```
conda create -n mokit-py39 python=3.9 # 3.8~3.11 are available
conda activate mokit-py39
conda install mokit -c mokit
```

You need to keep the environment `mokit-py39` activated when using MOKIT. If you want to install MOKIT with conda-forge channel, please read [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap2-2.html#use-mokit-with-conda-forge-channel). If you have no access to network, but still don't want to compile MOKIT manually, you can try option 3.

### Option 2: Use homebrew-toolchains (for MacOS only)
* Prerequisites: 
    - You need to install [homebrew](https://brew.sh) on your mac. See [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap2-2.html#optional-2-use-homebrew-toolchains-for-macos-only) for more tips.
    - You need to install conda via brew and install numpy in base env. via pip, as follows

```
brew install --cask miniconda
conda init bash # (or zsh) 
conda activate base
pip install numpy
```

Then 
`brew install ansatzx/homebrew-mokit/mokit`

Or `brew tap ansatzx/homebrew-mokit` and then `brew install mokit`.

Finally, follow caveats guides, add the following in your zsh(bash/fish etc.) profile.
```zsh
export MOKIT_ROOT="$(brew --prefix)/Cellar/mokit/master"
export PATH=$MOKIT_ROOT/bin:$PATH
export PYTHONPATH=$MOKIT_ROOT:$PYTHONPATH
export LD_LIBRARY_PATH=$MOKIT_ROOT/mokit/lib:$LD_LIBRARY_PATH
```

### Option 3: Use Pre-compiled MOKIT

Pre-built `Linux` executables can be downloaded [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap2-2.html#222-pre-built-linux-executables-and-libraries).

* Prerequisites: Python3 environment and NumPy.

* A detailed guide for choosing the version of pre-built artifacts and resolving dependencies can be found [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap2-2.html#222-pre-built-linux-executables-and-libraries)

* After downloading the pre-built artifacts, you need to set the following environment variables (assuming MOKIT is put in `$HOME/software/mokit`) in your `~/.bashrc`:

```bash
export MOKIT_ROOT=$HOME/software/mokit
export PATH=$MOKIT_ROOT/bin:$PATH
export PYTHONPATH=$MOKIT_ROOT:$PYTHONPATH
export LD_LIBRARY_PATH=$MOKIT_ROOT/mokit/lib:$LD_LIBRARY_PATH
export GMS=$HOME/software/gamess/rungms
```
The `LD_LIBRARY_PATH` is needed since the OpenBLAS dynamic library is put there. Remember to modify the `GMS` path to suit your local environment. 
Note that you need to exit the terminal and re-login, in order to activate newly written environment variables.

### Option 4: Build from Source
The link to latest version of MOKIT source code can be found [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap2-3.html).

* Prerequisites
    - Fortran compiler: `ifort`(>=2017) or `gfortran`(>=4.8.5)
    - Intel MKL(recommended) or [OpenBLAS](https://github.com/xianyi/OpenBLAS)
    - f2py (installing Anaconda Python3 recommended)

* Compile all modules
```
cd src
make all
```

* After `make all`, you need to set environment variables `MOKIT_ROOT`, `PATH` and `PYTHONPATH`. E.g. if MOKIT is installed in `$HOME/software/mokit`, the following should be set in `~/.bashrc`:

```bash
export MOKIT_ROOT=$HOME/software/mokit
export PATH=$MOKIT_ROOT/bin:$PATH
export PYTHONPATH=$MOKIT_ROOT:$PYTHONPATH
export GMS=$HOME/software/gamess/rungms
```

Remember to modify the `GMS` path to suit your local environment. 
Note that you need to exit the terminal and re-login, in order to activate newly written environment variables.

Quick Start
-----------
* Each utility is self-explanatory. For example, run `fch2inp` in Shell, you will find
```
 ERROR in subroutine fch2inp: wrong command line arguments!  
 Example 1 (R(O)HF, UHF, CAS): fch2inp a.fch  
 Example 2 (GVB)             : fch2inp a.fch -gvb [npair]  
 Example 3 (ROGVB)           : fch2inp a.fch -gvb [npair] -open [nopen]
```
You can search a utility and read their documentations [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap4-5.html).

* For usages of modules in mokit/lib/, see [examples/utilities/readme.txt](examples/utilities/readme.txt)

* The input syntax of the `automr` program is like Gaussian gjf. For example, the input file `00-h2o_cc-pVDZ_1.5.gjf` of the water molecule at d(O-H) = 1.5 A is shown below
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

Run
```
automr 00-h2o_cc-pVDZ_1.5.gjf >& 00-h2o_cc-pVDZ_1.5.out
```
in Shell. The `automr` program will successively perform HF, GVB, and CASSCF computations. The active space will be automatically determined as (4,4) during computations. 
See [Quick Start at documentation](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap3_quick.html) and [User Guide](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap4_guide.html) for more information.
See [examples/](examples/) for more examples.

QC Packages can be called
----------
* [Gaussian](http://gaussian.com)
* [PySCF](https://github.com/pyscf/pyscf)
* [GAMESS](https://www.msg.chem.iastate.edu/gamess/index.html)
* [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)
* [Molpro](https://www.molpro.net)
* [ORCA](https://orcaforum.kofo.mpg.de)
* [BDF](http://182.92.69.169:7226/Introduction)
* [PSI4](https://github.com/psi4/psi4)
* [Dalton](https://gitlab.com/dalton/dalton)
* [Q-Chem](https://www.q-chem.com)

Some Tips
---------
* To avoid unnecessary errors, you must specify keywords `nosymm int=nobasistransform` in Gaussian .gjf file, if you want to provide a .fch(k) file to `automr`.

Bug Report
----------
* If you find any bug frequently occurs or have any suggestions, you can open an
  issue on the [Issues](https://gitlab.com/jxzou/mokit/-/issues) page.

* You can also contact the developer jxzou via E-mail njumath[at]sina.cn, with
  related files (gjf, fch, out, etc) attached.

TODO
----
* MOs transferring among NWCHEM, BAGEL, COLUMBUS, etc.

* Develop/Implement robust black-box strategies of excited state calculations

Citation
--------
* Currently we have not published the paper of MOKIT program. If you use (any module or utility of) MOKIT in your work, please cite MOKIT as

   Jingxiang Zou, Molecular Orbital Kit (MOKIT), https://gitlab.com/jxzou/mokit (accessed month day, year)

* If you use MOKIT to perform calculations involving GVB, citing the following two papers would be appreciated

   DOI: [10.1021/acs.jctc.8b00854](https://www.doi.org/10.1021/acs.jctc.8b00854); DOI: [10.1021/acs.jpca.0c05216](https://www.doi.org/10.1021/acs.jpca.0c05216).

* If you use MOKIT in your work, please cite MOKIT in the main body of your paper. Citing MOKIT only in Supporting Information of your paper is insufficient. EndNote citation files can be found [here](https://gitlab.com/jxzou/mokit/-/tree/master/doc?ref_type=heads). More details and examples of citation can be found in [manual](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap1-2.html). 您的规范引用是对开发者的极大鼓励。您可以使用MOKIT为其他人做计算（包括代算），但务必提醒他/她在发表文章时恰当地引用MOKIT和计算中用到的量子化学软件。

* Click [here](https://jeanwsr.gitlab.io/mokit-doc-mdbook/citing.html) to see published papers which cited MOKIT.

