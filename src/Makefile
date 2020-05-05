# written by jxzou at 20200322

# --- gfortran ---
#F90    = gfortran
#FFLAGS = -O2
#F2PY   = f2py
#F2_FLAGS =
#F2_MKL_FLAGS = --link-lapack_opt

# --- ifort ----
F90    = ifort
FFLAGS = -O2
MKL_FLAGS = -mkl

# --- F2PY ---
F2PY   = f2py
F2_FLAGS = --fcompiler=intelem --compiler=intelem
F2_MKL_FLAGS = --link-lapack_opt

OBJ_bas_gms2molcas = bas_gms2molcas.o
OBJ_bas_gms2py = bas_gms2py.o
OBJ_dat2fch = dat2fch.o
OBJ_extract_noon2fch = extract_noon2fch.o
OBJ_fch2mkl = read_fch.o fch2mkl.o
OBJ_fch2inp = read_fch.o fch2inp.o
OBJ_fch2inporb = fch2inporb.o rwwfn.o
OBJ_fch_u2r = fch_u2r.o
OBJ_fch_mo_copy = fch_mo_copy.o
OBJ_chk2py = chk2py.f90
OBJ_fch2py = fch2py.f90
OBJ_find_mc_pair = find_mc_pair.o
OBJ_gvb_sort_pairs = gvb_sort_pairs.o
OBJ_py2fch = py2fch.f90
OBJ_pyuno  = pyuno.f90
OBJ_construct_vir = construct_vir.f90
OBJ_ortho  = ortho.f90
OBJ_assoc_rot = assoc_rot.f90
OBJ_assoc_loc = assoc_loc.f90
OBJ_pop = pop.f90
OBJ_lo = lo.f90
OBJ_auto_pair = auto_pair.f90
OBJ_mo_svd = mo_svd.f90
OBJ_solve_ON_matrix = solve_ON_matrix.o rwwfn.o

.PHONY: clean distclean bas_gms2molcas bas_gms2py dat2fch fch2mkl fch2inp fch2inporb fch_u2r fch_mo_copy find_mc_pair gvb_sort_pairs solve_ON_matrix extract_noon2fch

%.o: %.f90
	$(F90) -c $< -o $@ $(FFLAGS)

help:
	@echo ;\
	echo " make [target], where the [target] could be" ;\
	echo ;\
	echo " help      : print the current information" ;\
	echo " fch2inp   : generate fch2inp" ;\
	echo " fch2inporb: generate fch2inporb" ;\
	echo " xx        : generate xx or xx.so (xx is any individual module)" ;\
	echo " all       : all executables and *.so files" ;\
	echo " install   : move all executables and .so files to ../bin and ../lib" ;\
	echo " clean     : delete *.mod *.o" ;\
	echo " distclean : delete *.mod *.o and clean ../lib ../bin" ;\
	echo

bas_gms2molcas: $(OBJ_bas_gms2molcas)
	$(F90) $(OBJ_bas_gms2molcas) -o bas_gms2molcas $(FFLAGS)

bas_gms2py: $(OBJ_bas_gms2py)
	$(F90) $(OBJ_bas_gms2py) -o bas_gms2py $(FFLAGS)

dat2fch: $(OBJ_dat2fch)
	$(F90) $(OBJ_dat2fch) -o dat2fch $(FFLAGS)

fch2mkl: $(OBJ_fch2mkl)
	$(F90) $(OBJ_fch2mkl) -o fch2mkl $(FFLAGS)

fch2inp: $(OBJ_fch2inp)
	$(F90) $(OBJ_fch2inp) -o fch2inp $(FFLAGS)

fch2inporb: $(OBJ_fch2inporb)
	$(F90) $(OBJ_fch2inporb) -o fch2inporb $(FFLAGS)

fch_u2r: $(OBJ_fch_u2r)
	$(F90) $(OBJ_fch_u2r) -o fch_u2r $(FFLAGS)

fch_mo_copy: $(OBJ_fch_mo_copy)
	$(F90) $(OBJ_fch_mo_copy) -o fch_mo_copy $(FFLAGS)

chk2py: $(OBJ_chk2py)
	$(F2PY) -c $(OBJ_chk2py) -m chk2py $(F2_FLAGS)

fch2py: $(OBJ_fch2py)
	$(F2PY) -c $(OBJ_fch2py) -m fch2py $(F2_FLAGS)

find_mc_pair: $(OBJ_find_mc_pair)
	$(F90) $(OBJ_find_mc_pair) -o find_mc_pair $(FFLAGS)

py2fch: $(OBJ_py2fch)
	$(F2PY) -c $(OBJ_py2fch) -m py2fch $(F2_FLAGS)

gvb_sort_pairs: $(OBJ_gvb_sort_pairs)
	$(F90) $(OBJ_gvb_sort_pairs) -o gvb_sort_pairs $(FFLAGS)

pyuno: $(OBJ_pyuno)
	$(F2PY) -c $(OBJ_pyuno) -m uno $(F2_FLAGS) $(F2_MKL_FLAGS)

construct_vir: $(OBJ_construct_vir)
	$(F2PY) -c $(OBJ_construct_vir) -m construct_vir $(F2_FLAGS) $(F2_MKL_FLAGS)

ortho: $(OBJ_ortho)
	$(F2PY) -c $(OBJ_ortho) -m ortho $(F2_FLAGS) $(F2_MKL_FLAGS)

assoc_rot: $(OBJ_assoc_rot)
	$(F2PY) -c $(OBJ_assoc_rot) -m assoc_rot $(F2_FLAGS) $(F2_MKL_FLAGS)

assoc_loc: $(OBJ_assoc_loc)
	$(F2PY) -c $(OBJ_assoc_loc) -m assoc_loc $(F2_FLAGS)

pop: $(OBJ_pop)
	$(F2PY) -c $(OBJ_pop) -m pop $(F2_FLAGS)

lo: $(OBJ_lo)
	$(F2PY) -c $(OBJ_lo) -m lo $(F2_FLAGS) $(F2_MKL_FLAGS)

auto_pair: $(OBJ_auto_pair)
	$(F2PY) -c $(OBJ_auto_pair) -m auto_pair $(F2_FLAGS)

mo_svd: $(OBJ_mo_svd)
	$(F2PY) -c $(OBJ_mo_svd) -m mo_svd $(F2_FLAGS)

solve_ON_matrix: $(OBJ_solve_ON_matrix)
	$(F90) $(OBJ_solve_ON_matrix) $(MKL_FLAGS) -o solve_ON_matrix

extract_noon2fch: $(OBJ_extract_noon2fch)
	$(F90) $(OBJ_extract_noon2fch) -o extract_noon2fch

all: bas_gms2molcas bas_gms2py dat2fch fch2mkl fch2inp fch2inporb fch_u2r fch_mo_copy chk2py \
     fch2py find_mc_pair gvb_sort_pairs py2fch pyuno construct_vir ortho assoc_rot assoc_loc \
     pop lo auto_pair mo_svd solve_ON_matrix extract_noon2fch

install:
	mv bas_gms2molcas bas_gms2py dat2fch fch2mkl fch2inp fch2inporb fch_u2r fch_mo_copy find_mc_pair ../bin;
	mv gvb_sort_pairs solve_ON_matrix extract_noon2fch ../bin;
	mv chk2py*.so fch2py*.so py2fch*.so uno*.so construct_vir*.so ortho*.so assoc_rot*.so assoc_loc*.so ../lib;
	mv pop*.so lo*.so auto_pair*.so mo_svd*.so ../lib

clean:
	rm -f *.o *.mod

distclean:
	rm -f *.o *.mod *.so ../bin/* ../lib/*

