# Makefile for seismicgradiometry
# Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
# Released under the MIT license.
# see https://opensource.org/licenses/MIT


.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

FC = ifort
FFLAGS = -assume byterecl -mcmodel=medium -O3 -xHOST -no-prec-div -ipo -qmkl
DEFS = -DDOUBLE -DMKL -DELLIPSE
LIBS = -lnetcdff -lnetcdf -lmkl_lapack95_lp64 -lhdf5
LIBDIR = -L${NETCDF_FORTRAN_LIB} -L${NETCDF_LIB} -L${HDF5_LIB}
INCDIR = -I. -I${NETCDF_FORTRAN_INC}


calc_bpf_coef = calc_bpf_coef.f90
calc_bpf_order = calc_bpf_order.f90
calc_lpf_coef = calc_lpf_coef.f90
calc_lpf_order = calc_lpf_order.f90
constants = constants.F90
deconvolution = deconvolution.F90
gradiometry_parameters = gradiometry_parameters.F90
greatcircle = greatcircle.f90
lonlat_xy_conv = lonlat_xy_conv.F90
nrtype = nrtype.F90
read_sacfile = read_sacfile.F90
sac_decimation = sac_decimation.F90
sac_deconvolve = sac_deconvolve.F90
sac_integrate = sac_integrate.F90
seismicgradiometry = seismicgradiometry.F90
sort = sort.F90
tandem = tandem.f90
grdfile_io = grdfile_io.F90
line_fit = line_fit.f90
typedef = typedef.F90
calc_kernelmatrix = calc_kernelmatrix.F90
geompack2 = geompack2.f90
geometry = geometry.f90

##Module
mod_nrtype = $(nrtype:.F90=.mod)
mod_constants = $(constants:.F90=.mod)
mod_gradiometry_parameters = $(gradiometry_parameters:.F90=.mod)
mod_greatcircle = $(greatcircle:.f90=.mod)
mod_lonlat_xy_conv = $(lonlat_xy_conv:.F90=.mod)
mod_read_sacfile = $(read_sacfile:.F90=.mod)
mod_grdfile_io = $(grdfile_io:.F90=.mod)
mod_deconvolution = $(deconvolution:.F90=.mod)
mod_typedef = $(typedef:.F90=.mod)
mod_calc_kernelmatrix = $(calc_kernelmatrix:.F90=.mod)

##Object
o_calc_bpf_coef = $(calc_bpf_coef:.f90=.o)
o_calc_bpf_order = $(calc_bpf_order:.f90=.o)
o_calc_lpf_coef = $(calc_lpf_coef:.f90=.o)
o_calc_lpf_order = $(calc_lpf_order:.f90=.o)
o_tandem = $(tandem:.f90=.o)
o_constants = $(constants:.F90=.o)
o_gradiometry_parameters = $(gradiometry_parameters:.F90=.o)
o_deconvolution = $(deconvolution:.F90=.o)
o_grdfile_io = $(grdfile_io:.F90=.o)
o_greatcircle = $(greatcircle:.f90=.o)
o_lonlat_xy_conv = $(lonlat_xy_conv:.F90=.o)
o_nrtype = $(nrtype:.F90=.o)
o_read_sacfile = $(read_sacfile:.F90=.o)
o_sac_decimation = $(sac_decimation:.F90=.o)
o_sac_deconvolve = $(sac_deconvolve:.F90=.o)
o_sac_integrate = $(sac_integrate:.F90=.o)
o_seismicgradiometry = $(seismicgradiometry:.F90=.o)
o_sort = $(sort:.F90=.o)
o_tandem = $(tandem:.f90=.o)
o_deconvolution = $(deconvolution:.F90=.o)
o_line_fit = $(line_fit:.f90=.o)
o_typedef = $(typedef:.F90=.o)
o_calc_kernelmatrix = $(calc_kernelmatrix:.F90=.o)
o_geompack2 = $(geompack2:.f90=.o)
o_geometry = $(geometry:.f90=.o)

##Module dependency
$(mod_constants): $(constants) $(o_constants)
$(mod_gradiometry_parameters): $(gradiometry_parameters) $(o_gradiometry_parameters)
$(mod_greatcircle): $(greatcircle) $(o_greatcircle)
$(mod_nrtype): $(nrtype) $(o_nrtype)
$(mod_lonlat_xy_conv): $(lonlat_xy_conv) $(o_lonlat_xy_conv)
$(mod_read_sacfile): $(read_sacfile) $(o_read_sacfile)
$(mod_grdfile_io): $(grdfile_io) $(o_grdfile_io)
$(mod_deconvolution): $(deconvolution) $(o_deconvolution)
$(mod_typedef): $(typedef) $(o_typedef)
$(mod_calc_kernelmatrix): $(calc_kernelmatrix) $(o_calc_kernelmatrix)

##Object dependency
$(o_calc_bpf_coef): $(calc_bpf_coef) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_calc_bpf_order): $(calc_bpf_order) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_calc_lpf_coef): $(calc_lpf_coef) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_calc_lpf_order): $(calc_lpf_order) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_tandem): $(tandem) $(mod_nrtype) $(o_nrtype)
$(o_constants): $(constants) $(mod_nrtype) $(o_nrtype)
$(o_deconvolution): $(deconvolution) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_gradiometry_parameters): $(gradiometry_parameters) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_greatcircle): $(greatcircle) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_lonlat_xy_conv): $(lonlat_xy_conv) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_nrtype): $(nrtype)
$(o_read_sacfile): $(read_sacfile) $(mod_nrtype) $(o_nrtype)
$(o_sac_decimation): $(sac_decimation) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(o_nrtype) $(o_constants) \
	$(o_read_sacfile)
$(o_sac_deconvolve): $(sac_deconvolve) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_deconvolution) \
	$(o_nrtype) $(o_constants) $(o_read_sacfile) $(o_deconvolution)
$(o_sac_integrate): $(sac_integrate) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(o_nrtype) $(o_constants) \
	$(o_read_sacfile)
$(o_grdfile_io): $(grdfile_io) $(mod_nrtype) $(o_nrtype)
$(o_typedef): $(typedef) $(mod_nrtype) $(o_nrtype)
$(o_geometry): $(geometry) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_geompack2): $(geompack2) $(mod_nrtype) $(mod_constants) $(o_nrtype) $(o_constants)
$(o_calc_kernelmatrix): $(calc_kernelmatrix) $(mod_nrtype) $(mod_constants) $(mod_typedef) $(mod_greatcircle) $(mod_sort) \
	$(o_gradiometry_parameters) $(o_nrtype) $(o_constants) $(o_typedef) $(o_gradiometry) $(o_sort)
$(o_seismicgradiometry): $(seismicgradiometry) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_grdfile_io) \
	$(mod_lonlat_xy_conv) $(mod_typedef) $(mod_gradiometry_parameters) $(mod_calc_kernelmatrix) \
	$(o_gradiometry_parameters) $(o_calc_kernelmatrix) $(o_nrtype) $(o_constants) $(o_read_sacfile) $(o_grdfile_io) \
	$(o_lonlat_xy_conv) $(o_typedef)

.F90.o:
	$(FC) $< -c -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)
.f90.o:
	$(FC) $< -c -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)
.F90.mod:
	@:
.f90.mod:
	@:

seismicgradiometry: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_lonlat_xy_conv) \
	$(o_grdfile_io) $(o_read_sacfile) $(o_sort) $(o_greatcircle) $(o_seismicgradiometry) $(o_gradiometry_parameters) \
        $(o_typedef) $(o_calc_kernelmatrix) $(o_geompack2) $(o_geometry)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

sac_decimation: $(o_nrtype) $(o_constants) $(o_calc_lpf_order) $(o_calc_lpf_coef) $(o_tandem) $(o_read_sacfile) \
	$(o_sac_decimation)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

sac_deconvolve: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_read_sacfile) \
	$(o_deconvolution) $(o_sac_deconvolve)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

sac_integrate: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_read_sacfile) \
	$(o_line_fit) $(o_sac_integrate)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)


clean:
	rm -f *.o *.mod

