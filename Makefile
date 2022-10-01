#Makefile for seismicgradiometry

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

FC = ifort
FFLAGS = -assume byterecl -mcmodel=medium -O3 -xHOST -no-prec-div -ipo
DEFS = -DDOUBLE -DMKL
LIBS = -liomp5 -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lnetcdff -lnetcdf -lhdf5
LIBDIR = -L${MKLROOT}/lib/intel64 -L/usr/local/netcdff-4.5.2/lib -L/usr/local/netcdf-4.7.4/lib -L/usr/local/hdf5-1.12/lib
INCDIR = -I. -I${MKLROOT}/include/intel64/lp64 -I/usr/local/netcdff-4.5.2/include




calc_bpf_coef = calc_bpf_coef.f90
calc_bpf_order = calc_bpf_order.f90
calc_lpf_coef = calc_lpf_coef.f90
calc_lpf_order = calc_lpf_order.f90
constants = constants.F90
deconvolution = deconvolution.F90
greatcircle = greatcircle.f90
lonlat_xy_conv = lonlat_xy_conv.F90
nrtype = nrtype.F90
read_sacfile = read_sacfile.F90
sac_decimation = sac_decimatino.F90
sac_deconvolve = sac_deconvlve.F90
sac_integrate = sac_integrate.F90
seismicgradiometry = seismicgradiometry.F90
sort = sort.F90
tandem = tandem.f90
grdfile_io = grdfile_io.F90

##Module
mod_constants = $(constants:.F90=.mod)
mod_greatcircle = $(greatcircle:.f90=.mod)
mod_lonlat_xy_conv = $(lonlat_xy_conv:.F90=.mod)
mod_nrtype = $(nrtype:.F90=.mod)
mod_read_sacfile = $(read_sacfile:.F90=.mod)
mod_grdfile_io = $(grdfile_io:.F90=.mod)

##Object
o_calc_bpf_coef = $(calc_bpf_coef:.f90=.o)
o_calc_bpf_order = $(calc_bpf_order:.f90=.o)
o_calc_lpf_coef = $(calc_lpf_coef:.f90=.o)
o_calc_lpf_order = $(calc_lpf_order:.f90=.o)
o_tandem = $(tandem:.f90=.o)
o_constants = $(constants:.F90=.o)
o_deconvolution = $(deconvolution:.F90=.o)
o_grdfile_io = $(grdfile_io:.F90=.o)
o_greatcircle = $(greatcircle:.f90=.o)
o_lonlat_xy_conv = $(lonlat_xy_conv:.F90=.o)
o_nrtype = $(nrtype:.F90=.o)
o_read_sacfile = $(read_sacfile:.F90=.o)
o_sac_decimation = $(sac_decimation:.F90=.o)
o_sac_deconvolve = $(sac_deconvolve:.F90=.o)
o_sac_integreate = $(sac_integrate:.F90=.o)
o_seismicgradiometry = $(seismicgradiometry:.F90=.o)
o_sort = $(sort:.F90=.o)
o_tandem = $(tandem:.f90=.o)

##Module dependency
$(mod_constants): $(constants) $(o_constants)
$(mod_greatcircle): $(greatcircle) $(o_greatcircle)
$(mod_nrtype): $(nrtype) $(o_nrtype)
$(mod_lonlat_xy_conv): $(lonlat_xy_conv) $(o_lonlat_xy_conv)
$(mod_read_sacfile): $(read_sacfile) $(o_read_sacfile)
$(mod_grdfile_io): $(grdfile_io) $(o_grdfile_io)

##Object dependency
$(o_calc_bpf_coef): $(calc_bpf_coef) $(mod_nrtype) $(mod_constants)
$(o_calc_bpf_order): $(calc_bpf_order) $(mod_nrtype) $(mod_constants)
$(o_calc_lpf_coef): $(calc_lpf_coef) $(mod_nrtype) $(mod_constants)
$(o_calc_lpf_order): $(calc_lpf_order) $(mod_nrtype) $(mod_constants)
$(o_tandem): $(tandem) $(mod_nrtype)
$(o_constants): $(constants) $(mod_nrtype)
$(o_deconvolution): $(deconvolution) $(mod_nrtype) $(mod_constants)
$(o_greatcircle): $(greatcircle) $(mod_nrtype) $(mod_constants)
$(o_lonlat_xy_conv): $(lonlat_xy_conv) $(mod_nrtype) $(mod_constants)
$(o_nrtype): $(nrtype)
$(o_read_sacfile): $(read_sacfile) $(mod_nrtype)
$(o_sac_decimation): $(sac_decimation) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile)
$(o_sac_deconvolve): $(sac_deconvolve) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile)
$(o_sac_integrate): $(sac_integrate) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile)
$(o_grdfile_io): $(grdfile_io) $(mod_nrtype)
	$(FC) $(DEFS) $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $< -c -o $@
$(o_seismicgradiometry): $(seismicgradiometry) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_grdfile_io) \
	$(mod_lonlat_xy_conv) $(mod_sort)
	$(FC) $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS) $< -c -o $@

.F90.o:
	$(FC) $(DEFS) $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $< -c -o $@
.f90.o:
	$(FC) $(DEFS) $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $< -c -o $@
.F90.mod:
	@:
.f90.mod:
	@:

seismicgradiometry: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_lonlat_xy_conv) \
	$(o_grdfile_io) $(o_read_sacfile) $(o_sort) $(o_seismicgradiometry)
	$(FC) $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS) $^ -o $@

clean:
	rm -f *.o *.mod

