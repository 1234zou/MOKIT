# written by jxzou at 20200322
# updated by jxzou at 20200515: delete 'install', add $(BIN)

BIN = ../bin

# --- gfortran ---
#F90    = gfortran
#FFLAGS = -O2
#F2PY   = f2py
#F2_FLAGS =
#F2_MKL_FLAGS = --link-lapack_opt

# --- ifort ----
F90    = ifort
FFLAGS = -O2 -traceback
MKL_FLAGS = -mkl
F2PY   = f2py
F2_FLAGS = --fcompiler=intelem --compiler=intelem
F2_MKL_FLAGS = --link-lapack_opt

OBJ_align_orbitals = rwwfn.o align_orbitals.o
OBJ_automr = mr_keyword.o rwwfn.o rwgeom.o automr.o
OBJ_auto_pair = auto_pair.f90
OBJ_bas_fch2py = bas_fch2py.o rwwfn.o
OBJ_bas_gms2molcas = bas_gms2molcas.o
OBJ_bas_gms2py = bas_gms2py.o
OBJ_dat2fch = rwwfn.o dat2fch.o
OBJ_extract_noon2fch = extract_noon2fch.o
OBJ_fch2mkl = read_fch.o fch2mkl.o
OBJ_fch2inp = read_fch.o fch2inp.o
OBJ_fch2inporb = fch2inporb.o rwwfn.o
OBJ_fch_u2r = fch_u2r.o
OBJ_fch_mo_copy = rwwfn.o fch_mo_copy.o
OBJ_chk2py = chk2py.f90
OBJ_fch2py = fch2py.f90
OBJ_find_mc_pair = find_mc_pair.o
OBJ_gvb_correct_pairs = rwwfn.o gvb_correct_pairs.o
OBJ_gvb_sort_pairs = gvb_sort_pairs.o
OBJ_mkl2fch = rwwfn.o read_fch.o mkl2fch.o
OBJ_py2fch = py2fch.f90
OBJ_pyuno  = pyuno.f90
OBJ_construct_vir = construct_vir.f90
OBJ_ortho  = ortho.f90
OBJ_assoc_rot = assoc_rot.f90
OBJ_assoc_loc = assoc_loc.f90
OBJ_pop = pop.f90
OBJ_lo = lo.f90
OBJ_mo_svd = mo_svd.f90
OBJ_mo_svd_format = rwwfn.o mo_svd_format.o
OBJ_solve_ON_matrix = rwwfn.o solve_ON_matrix.o

.PHONY: clean distclean align_orbitals automr bas_fch2py bas_gms2molcas bas_gms2py dat2fch \
        fch2mkl fch2inp fch2inporb fch_u2r fch_mo_copy find_mc_pair gvb_correct_pairs \
        gvb_sort_pairs mkl2fch mo_svd mo_svd_format solve_ON_matrix extract_noon2fch

%.o: %.f90
	$(F90) -c $< -o $@ $(FFLAGS)

help:
	@echo ;\
	echo " make [target], where the [target] could be" ;\
	echo ;\
	echo " help      : print the current information" ;\
	echo " automr    : automatically do multireference calculations" ;\
	echo " fch2inp   : generate fch2inp" ;\
	echo " fch2inporb: generate fch2inporb" ;\
	echo " xx        : generate xx or xx.so (xx is any individual module)" ;\
	echo " all       : all executables and *.so files" ;\
	echo " clean     : delete *.mod *.o" ;\
	echo " distclean : delete *.mod *.o and clean ../lib ../bin" ;\
	echo

align_orbitals: $(OBJ_align_orbitals)
	$(F90) $(OBJ_align_orbitals) -o $(BIN)/align_orbitals $(FFLAGS)

automr: $(OBJ_automr)
	$(F90) $(OBJ_automr) -o $(BIN)/automr $(FFLAGS)
	@echo ;\
	echo "AutoMR requires almost all modules, so run 'make all' if other modules haven't been compiled" ;\
	echo

bas_fch2py: $(OBJ_bas_fch2py)
	$(F90) $(OBJ_bas_fch2py) -o $(BIN)/bas_fch2py $(FFLAGS)
	@echo ;\
	echo "bas_fch2py requires modules fch2inp and bas_gms2py, so compile them if they haven't been compiled" ;\
	echo

bas_gms2molcas: $(OBJ_bas_gms2molcas)
	$(F90) $(OBJ_bas_gms2molcas) -o $(BIN)/bas_gms2molcas $(FFLAGS)

bas_gms2py: $(OBJ_bas_gms2py)
	$(F90) $(OBJ_bas_gms2py) -o $(BIN)/bas_gms2py $(FFLAGS)

dat2fch: $(OBJ_dat2fch)
	$(F90) $(OBJ_dat2fch) -o $(BIN)/dat2fch $(FFLAGS)

fch2mkl: $(OBJ_fch2mkl)
	$(F90) $(OBJ_fch2mkl) -o $(BIN)/fch2mkl $(FFLAGS)

fch2inp: $(OBJ_fch2inp)
	$(F90) $(OBJ_fch2inp) -o $(BIN)/fch2inp $(FFLAGS)

fch2inporb: $(OBJ_fch2inporb)
	$(F90) $(OBJ_fch2inporb) -o $(BIN)/fch2inporb $(FFLAGS)
	@echo ;\
	echo "fch2inporb requires moduled fch2inp and bas_gms2molcas, so compile them if they haven't been compiled" ;\
	echo

fch_u2r: $(OBJ_fch_u2r)
	$(F90) $(OBJ_fch_u2r) -o $(BIN)/fch_u2r $(FFLAGS)

fch_mo_copy: $(OBJ_fch_mo_copy)
	$(F90) $(OBJ_fch_mo_copy) -o $(BIN)/fch_mo_copy $(FFLAGS)

chk2py: $(OBJ_chk2py)
	$(F2PY) -c $(OBJ_chk2py) -m chk2py $(F2_FLAGS)
	@mv chk2py*.so ../lib/

fch2py: $(OBJ_fch2py)
	$(F2PY) -c $(OBJ_fch2py) -m fch2py $(F2_FLAGS)
	@mv fch2py*.so ../lib/

find_mc_pair: $(OBJ_find_mc_pair)
	$(F90) $(OBJ_find_mc_pair) -o $(BIN)/find_mc_pair $(FFLAGS)

py2fch: $(OBJ_py2fch)
	$(F2PY) -c $(OBJ_py2fch) -m py2fch $(F2_FLAGS)
	@mv py2fch*.so ../lib/

gvb_correct_pairs: $(OBJ_gvb_correct_pairs)
	$(F90) $(OBJ_gvb_correct_pairs) -o $(BIN)/gvb_correct_pairs $(FFLAGS)

gvb_sort_pairs: $(OBJ_gvb_sort_pairs)
	$(F90) $(OBJ_gvb_sort_pairs) -o $(BIN)/gvb_sort_pairs $(FFLAGS)

mkl2fch: $(OBJ_mkl2fch)
	$(F90) $(OBJ_mkl2fch) -o $(BIN)/mkl2fch $(FFLAGS)

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

extract_noon2fch: $(OBJ_extract_noon2fch)
	$(F90) $(OBJ_extract_noon2fch) -o $(BIN)/extract_noon2fch

all: automr auto_pair bas_fch2py bas_gms2molcas bas_gms2py dat2fch fch2mkl \
     fch2inp fch2inporb fch_u2r fch_mo_copy chk2py fch2py find_mc_pair gvb_correct_pairs \
     gvb_sort_pairs mkl2fch py2fch pyuno construct_vir ortho assoc_rot assoc_loc pop lo \
     mo_svd mo_svd_format solve_ON_matrix extract_noon2fch

clean:
	rm -f *.o *.mod

distclean:
	rm -f *.o *.mod *.so ../bin/* ../lib/*.so

