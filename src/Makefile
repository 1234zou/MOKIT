# written by jxzou at 20200322
# updated by jxzou at 20200515: delete 'install', add $(BIN)

BIN = ../bin

# --- gfortran ---
#F90    = gfortran
#FFLAGS = -O2 -cpp
#MKL_FLAGS = -lmkl_rt
#F2PY   = f2py
#F2_FLAGS =
#F2_MKL_FLAGS = --link-lapack_opt

# --- ifort ----
F90    = ifort
FFLAGS = -O2 -fpp
MKL_FLAGS = -mkl
F2PY   = f2py
F2_FLAGS = --fcompiler=intelem --compiler=intelem
F2_MKL_FLAGS = --link-lapack_opt

OBJ_add_bgcharge_to_inp = file_op.o string_manipulate.o read_gms_inp.o add_bgcharge_to_inp.o
OBJ_align_orbitals = rwwfn.o align_orbitals.o
OBJ_automr = file_op.o string_manipulate.o mr_keyword.o util_wrapper.o rwwfn.o \
             read_fch.o rwgeom.o do_hf.o do_gvb.o do_cas.o do_mrpt2.o do_mrpt3.o \
             do_mcpdft.o do_mrci.o do_mrcc.o call_qc_calc_int.o automr.o
OBJ_auto_pair = auto_pair.f90
OBJ_bas_fch2py = file_op.o string_manipulate.o util_wrapper.o rwwfn.o bas_fch2py.o
OBJ_bas_gau2molcas = file_op.o string_manipulate.o read_gms_inp.o read_fch.o bas_gau2molcas.o
OBJ_bas_gms2bdf = file_op.o string_manipulate.o read_gms_inp.o bas_gms2bdf.o
OBJ_bas_gms2dal = file_op.o string_manipulate.o read_gms_inp.o bas_gms2dal.o
OBJ_bas_gms2molcas = file_op.o string_manipulate.o read_gms_inp.o bas_gms2molcas.o
OBJ_bas_gms2molpro = file_op.o string_manipulate.o read_gms_inp.o bas_gms2molpro.o
OBJ_bas_gms2psi = file_op.o string_manipulate.o read_gms_inp.o bas_gms2psi.o
OBJ_bas_gms2py = file_op.o string_manipulate.o read_gms_inp.o bas_gms2py.o
OBJ_bdf2fch = file_op.o string_manipulate.o rwwfn.o read_fch.o bdf2fch.o
OBJ_bdf2mkl = file_op.o bdf2mkl.o
OBJ_dal2fch = file_op.o string_manipulate.o rwwfn.o dal2fch.o
OBJ_dat2fch = file_op.o string_manipulate.o read_fch.o rwwfn.o dat2fch.o
OBJ_extract_noon2fch = file_op.o extract_noon2fch.o
OBJ_fch2bdf = file_op.o string_manipulate.o rwwfn.o util_wrapper.o read_fch.o \
              rwgeom.o read_gms_inp.o fch2bdf.o
OBJ_fch2com = file_op.o string_manipulate.o rwwfn.o util_wrapper.o read_fch.o fch2com.o
OBJ_fch2dal = file_op.o string_manipulate.o rwwfn.o util_wrapper.o read_fch.o fch2dal.o
OBJ_fch2mkl = file_op.o string_manipulate.o read_fch.o fch2mkl.o
OBJ_fch2inp = file_op.o string_manipulate.o read_fch.o fch2inp.o
OBJ_fch2inporb = file_op.o string_manipulate.o rwwfn.o util_wrapper.o read_fch.o fch2inporb.o
OBJ_fch2psi = file_op.o string_manipulate.o rwwfn.o util_wrapper.o read_fch.o fch2psi.o
OBJ_fch2qm4d = file_op.o string_manipulate.o read_fch.o rwgeom.o fch2qm4d.o
OBJ_fch_u2r = file_op.o fch_u2r.o
OBJ_fch_mo_copy = rwwfn.o fch_mo_copy.o
OBJ_frag_guess_wfn = file_op.o string_manipulate.o read_fch.o util_wrapper.o \
                     mr_keyword.o do_hf.o rwwfn.o rwgeom.o read_gms_inp.o \
                     call_qc_calc_int.o frag_guess_wfn.o
OBJ_chk2py = chk2py.f90
OBJ_fch2py = fch2py.f90
OBJ_find_mc_pair = file_op.o find_mc_pair.o
OBJ_gvb_correct_pairs = string_manipulate.o rwwfn.o gvb_correct_pairs.o
OBJ_gvb_sort_pairs = file_op.o gvb_sort_pairs.o
OBJ_mkl2fch = file_op.o string_manipulate.o read_fch.o read_mkl.o rwwfn.o mkl2fch.o
OBJ_mkl2gjf = file_op.o string_manipulate.o read_fch.o read_mkl.o mkl2gjf.o
OBJ_orb2fch = file_op.o string_manipulate.o rwwfn.o read_fch.o orb2fch.o
OBJ_replace_xyz_in_inp = file_op.o string_manipulate.o read_fch.o rwgeom.o \
                         replace_xyz_in_inp.o
OBJ_xml2fch = file_op.o string_manipulate.o rwwfn.o read_fch.o xml2fch.o
OBJ_py2fch = py2fch.f90
OBJ_pyuno  = pyuno.f90
OBJ_construct_vir = construct_vir.f90
OBJ_ortho  = ortho.f90
OBJ_assoc_rot = assoc_rot.f90
OBJ_assoc_loc = assoc_loc.f90
OBJ_pop = pop.f90
OBJ_lo = lo.f90
OBJ_mo_svd = mo_svd.f90
OBJ_mo_svd_format = string_manipulate.o rwwfn.o mo_svd_format.o
OBJ_solve_ON_matrix = string_manipulate.o rwwfn.o solve_ON_matrix.o
OBJ_rwgeom = string_manipulate.f90 rwgeom.f90
OBJ_rwwfn = rwwfn.f90
OBJ_excited = excited.f90

.PHONY: clean distclean add_bgcharge_to_inp align_orbitals automr bas_fch2py \
        bas_gau2molcas bas_gms2bdf bas_gms2dal bas_gms2molcas bas_gms2molpro \
        bas_gms2psi bas_gms2py bdf2fch bdf2mkl dal2fch dat2fch extract_noon2fch \
        fch2bdf fch2com fch2dal fch2mkl fch2inp fch2inporb fch2psi fch2qm4d \
        fch_u2r fch_mo_copy find_mc_pair gvb_correct_pairs gvb_sort_pairs mkl2fch \
        mkl2gjf orb2fch replace_xyz_in_inp xml2fch mo_svd mo_svd_format solve_ON_matrix \
        rwgeom rwwfn excited

%.o: %.f90
	$(F90) -c $< -o $@ $(FFLAGS)

help:
	@echo ;\
	echo " make [target], where the [target] could be" ;\
	echo ;\
	echo " help      : print the current information" ;\
	echo " automr    : generate automr, automated multireference calculations" ;\
	echo " fch2bdf   : generate fch2bdf    (Gau->BDF)" ;\
	echo " fch2com   : generate fch2com    (Gau->Molpro)" ;\
	echo " fch2dal   : generate fch2dal    (Gau->Dalton)" ;\
	echo " fch2inp   : generate fch2inp    (Gau->GAMESS)" ;\
	echo " fch2inporb: generate fch2inporb (Gau->OpenMolcas)" ;\
	echo " fch2mkl   : generate fch2mkl    (Gau->ORCA)" ;\
	echo " fch2psi   : generate fch2psi    (Gau->PSI4)" ;\
	echo " fch2py    : generate fch2py     (Gau->PySCF)" ;\
	echo " fch2qm4d  : generate fch2qm4d   (Gau->QM4D)" ;\
	echo " bdf2fch   : generate bdf2fch    (BDF->Gau)" ;\
	echo " dat2fch   : generate dat2fch    (GAMESS->Gau)" ;\
	echo " mkl2gjf   : generate mkl2gjf    (ORCA->Gau)" ;\
	echo " mkl2fch   : generate mkl2fch    (ORCA->Gau)" ;\
	echo " orb2fch   : generate orb2fch    (OpenMolcas->Gau)" ;\
	echo " xml2fch   : generate xml2fch    (Molpro->Gau)" ;\
	echo " xxx       : generate xxx or xxx.so (xxx is any individual module)" ;\
	echo " all       : all executables and *.so files" ;\
	echo " clean     : delete *.mod *.o" ;\
	echo " distclean : delete *.mod *.o and clean ../lib ../bin" ;\
	echo

add_bgcharge_to_inp: $(OBJ_add_bgcharge_to_inp)
	$(F90) $(OBJ_add_bgcharge_to_inp) -o $(BIN)/add_bgcharge_to_inp $(FFLAGS)

align_orbitals: $(OBJ_align_orbitals)
	$(F90) $(OBJ_align_orbitals) -o $(BIN)/align_orbitals $(FFLAGS)

automr: $(OBJ_automr)
	$(F90) $(OBJ_automr) -o $(BIN)/automr $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "AutoMR requires almost all modules, so run 'make all' if other modules haven't been compiled" ;\
	echo

bas_fch2py: $(OBJ_bas_fch2py)
	$(F90) $(OBJ_bas_fch2py) -o $(BIN)/bas_fch2py $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "bas_fch2py requires utilities fch2inp and bas_gms2py, please also compile them" ;\
	echo

bas_gau2molcas: $(OBJ_bas_gau2molcas)
	$(F90) $(OBJ_bas_gau2molcas) -o $(BIN)/bas_gau2molcas $(FFLAGS)

bas_gms2bdf: $(OBJ_bas_gms2bdf)
	$(F90) $(OBJ_bas_gms2bdf) -o $(BIN)/bas_gms2bdf $(FFLAGS)

bas_gms2dal: $(OBJ_bas_gms2dal)
	$(F90) $(OBJ_bas_gms2dal) -o $(BIN)/bas_gms2dal $(FFLAGS)

bas_gms2molcas: $(OBJ_bas_gms2molcas)
	$(F90) $(OBJ_bas_gms2molcas) -o $(BIN)/bas_gms2molcas $(FFLAGS)

bas_gms2molpro: $(OBJ_bas_gms2molpro)
	$(F90) $(OBJ_bas_gms2molpro) -o $(BIN)/bas_gms2molpro $(FFLAGS)

bas_gms2psi: $(OBJ_bas_gms2psi)
	$(F90) $(OBJ_bas_gms2psi) -o $(BIN)/bas_gms2psi $(FFLAGS)

bas_gms2py: $(OBJ_bas_gms2py)
	$(F90) $(OBJ_bas_gms2py) -o $(BIN)/bas_gms2py $(FFLAGS)

bdf2fch: $(OBJ_bdf2fch)
	$(F90) $(OBJ_bdf2fch) $(MKL_FLAGS) -o $(BIN)/bdf2fch $(FFLAGS)

bdf2mkl: $(OBJ_bdf2mkl)
	$(F90) $(OBJ_bdf2mkl) $(MKL_FLAGS) -o $(BIN)/bdf2mkl $(FFLAGS)
	@echo ;\
	echo "bdf2mkl requires utilities bdf2fch and fch2mkl, please also compile them" ;\
	echo

dal2fch: $(OBJ_dal2fch)
	$(F90) $(OBJ_dal2fch) -o $(BIN)/dal2fch $(FFLAGS) $(MKL_FLAGS)

dat2fch: $(OBJ_dat2fch)
	$(F90) $(OBJ_dat2fch) $(MKL_FLAGS) -o $(BIN)/dat2fch $(FFLAGS) $(MKL_FLAGS)

extract_noon2fch: $(OBJ_extract_noon2fch)
	$(F90) $(OBJ_extract_noon2fch) -o $(BIN)/extract_noon2fch

fch2bdf: $(OBJ_fch2bdf)
	$(F90) $(OBJ_fch2bdf) -o $(BIN)/fch2bdf $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "fch2bdf requires utilites fch2inp and bas_gms2bdf, please also compile them" ;\
	echo

fch2com: $(OBJ_fch2com)
	$(F90) $(OBJ_fch2com) -o $(BIN)/fch2com $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "fch2com requires utilites fch2inp and bas_gms2molpro, please also compile them" ;\
	echo

fch2dal: $(OBJ_fch2dal)
	$(F90) $(OBJ_fch2dal) -o $(BIN)/fch2dal $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "fch2dal requires utilites fch2inp and bas_gms2dal, please also compile them" ;\
	echo

fch2mkl: $(OBJ_fch2mkl)
	$(F90) $(OBJ_fch2mkl) -o $(BIN)/fch2mkl $(FFLAGS)

fch2inp: $(OBJ_fch2inp)
	$(F90) $(OBJ_fch2inp) -o $(BIN)/fch2inp $(FFLAGS)

fch2inporb: $(OBJ_fch2inporb)
	$(F90) $(OBJ_fch2inporb) -o $(BIN)/fch2inporb $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "fch2inporb requires utilities fch2inp and bas_gms2molcas, please also compile them" ;\
	echo

fch2psi: $(OBJ_fch2psi)
	$(F90) $(OBJ_fch2psi) -o $(BIN)/fch2psi $(FFLAGS) $(MKL_FLAGS)
	@echo ;\
	echo "fch2psi requires utilities fch2inp and bas_gms2psi, please also compile them" ;\
	echo

fch2qm4d: $(OBJ_fch2qm4d)
	$(F90) $(OBJ_fch2qm4d) -o $(BIN)/fch2qm4d $(FFLAGS) $(MKL_FLAGS)

fch_u2r: $(OBJ_fch_u2r)
	$(F90) $(OBJ_fch_u2r) -o $(BIN)/fch_u2r $(FFLAGS)

fch_mo_copy: $(OBJ_fch_mo_copy)
	$(F90) $(OBJ_fch_mo_copy) -o $(BIN)/fch_mo_copy $(FFLAGS) $(MKL_FLAGS)

chk2py: $(OBJ_chk2py)
	$(F2PY) -c $(OBJ_chk2py) -m chk2py $(F2_FLAGS)
	@mv chk2py*.so ../lib/

fch2py: $(OBJ_fch2py)
	$(F2PY) -c $(OBJ_fch2py) -m fch2py $(F2_FLAGS)
	@mv fch2py*.so ../lib/

find_mc_pair: $(OBJ_find_mc_pair)
	$(F90) $(OBJ_find_mc_pair) -o $(BIN)/find_mc_pair $(FFLAGS)

frag_guess_wfn: $(OBJ_frag_guess_wfn)
	$(F90) $(OBJ_frag_guess_wfn) -o $(BIN)/frag_guess_wfn $(FFLAGS) $(MKL_FLAGS)

orb2fch: $(OBJ_orb2fch)
	$(F90) $(OBJ_orb2fch) -o $(BIN)/orb2fch $(FFLAGS) $(MKL_FLAGS)

py2fch: $(OBJ_py2fch)
	$(F2PY) -c $(OBJ_py2fch) -m py2fch $(F2_FLAGS)
	@mv py2fch*.so ../lib/

gvb_correct_pairs: $(OBJ_gvb_correct_pairs)
	$(F90) $(OBJ_gvb_correct_pairs) -o $(BIN)/gvb_correct_pairs $(FFLAGS) $(MKL_FLAGS)

gvb_sort_pairs: $(OBJ_gvb_sort_pairs)
	$(F90) $(OBJ_gvb_sort_pairs) -o $(BIN)/gvb_sort_pairs $(FFLAGS)

mkl2fch: $(OBJ_mkl2fch)
	$(F90) $(OBJ_mkl2fch) -o $(BIN)/mkl2fch $(FFLAGS) $(MKL_FLAGS)

mkl2gjf: $(OBJ_mkl2gjf)
	$(F90) $(OBJ_mkl2gjf) -o $(BIN)/mkl2gjf $(FFLAGS)

replace_xyz_in_inp: $(OBJ_replace_xyz_in_inp)
	$(F90) $(OBJ_replace_xyz_in_inp) -o $(BIN)/replace_xyz_in_inp $(FFLAGS)

xml2fch: $(OBJ_xml2fch)
	$(F90) $(OBJ_xml2fch) -o $(BIN)/xml2fch $(FFLAGS) $(MKL_FLAGS)

pyuno: $(OBJ_pyuno)
	$(F2PY) -c $(OBJ_pyuno) -m uno $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv uno*.so ../lib/

construct_vir: $(OBJ_construct_vir)
	$(F2PY) -c $(OBJ_construct_vir) -m construct_vir $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv construct_vir*.so ../lib/

ortho: $(OBJ_ortho)
	$(F2PY) -c $(OBJ_ortho) -m ortho $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv ortho*.so ../lib/

assoc_rot: $(OBJ_assoc_rot)
	$(F2PY) -c $(OBJ_assoc_rot) -m assoc_rot $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv assoc_rot*.so ../lib/

assoc_loc: $(OBJ_assoc_loc)
	$(F2PY) -c $(OBJ_assoc_loc) -m assoc_loc $(F2_FLAGS)
	@mv assoc_loc*.so ../lib/

pop: $(OBJ_pop)
	$(F2PY) -c $(OBJ_pop) -m pop $(F2_FLAGS)
	@mv pop*.so ../lib/

lo: $(OBJ_lo)
	$(F2PY) -c $(OBJ_lo) -m lo $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv lo*.so ../lib/

auto_pair: $(OBJ_auto_pair)
	$(F2PY) -c $(OBJ_auto_pair) -m auto_pair $(F2_FLAGS)
	@mv auto_pair*.so ../lib/

mo_svd: $(OBJ_mo_svd)
	$(F2PY) -c $(OBJ_mo_svd) -m mo_svd $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv mo_svd*.so ../lib/

mo_svd_format: $(OBJ_mo_svd_format)
	$(F90) $(OBJ_mo_svd_format) $(MKL_FLAGS) -o $(BIN)/mo_svd

solve_ON_matrix: $(OBJ_solve_ON_matrix)
	$(F90) $(OBJ_solve_ON_matrix) $(MKL_FLAGS) -o $(BIN)/solve_ON_matrix

rwgeom: $(OBJ_rwgeom)
	$(F2PY) -c $(OBJ_rwgeom) -m rwgeom $(F2_FLAGS)
	@mv rwgeom*.so ../lib/

rwwfn: $(OBJ_rwwfn)
	$(F2PY) -c $(OBJ_rwwfn) -m rwwfn $(F2_FLAGS)
	@mv rwwfn*.so ../lib/

excited: $(OBJ_excited)
	$(F2PY) -c $(OBJ_excited) -m excited $(F2_FLAGS) $(F2_MKL_FLAGS)
	@mv excited*.so ../lib/

all: add_bgcharge_to_inp automr auto_pair bas_fch2py bas_gau2molcas bas_gms2bdf \
     bas_gms2dal bas_gms2molcas bas_gms2molpro bas_gms2psi bas_gms2py bdf2fch \
     bdf2mkl dal2fch dat2fch extract_noon2fch fch2bdf fch2com fch2dal fch2mkl \
     fch2inp fch2inporb fch2psi fch2qm4d fch_u2r fch_mo_copy frag_guess_wfn \
     chk2py fch2py find_mc_pair gvb_correct_pairs gvb_sort_pairs mkl2fch mkl2gjf \
     orb2fch replace_xyz_in_inp xml2fch py2fch pyuno construct_vir ortho assoc_rot \
     assoc_loc pop lo mo_svd mo_svd_format solve_ON_matrix rwgeom rwwfn excited

clean:
	rm -f *.o *.mod

distclean:
	rm -rf *.o *.mod *.so ../bin/* ../lib/*.so ../lib/__pycache__

