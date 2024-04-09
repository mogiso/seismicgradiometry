# Makefile for seismicgradiometry
# Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
# Released under the MIT license.
# see https://opensource.org/licenses/MIT


.SUFFIXES:
.SUFFIXES: .f .F90 .f90 .o

ifeq ($(arch), ifx)
  FC = ifx
  FFLAGS = -assume byterecl -mcmodel=medium -O3 -xHOST -no-prec-div -ipo -qmkl -warn all -traceback
  DEFS = -DDOUBLE -DMKL -DPARTICLEVELOCITY
  NETCDF_FORTRAN_LIB = /usr/local/intel/lib
  NETCDF_FORTRAN_INC = /usr/local/intel/include
  NETCDF_LIB = /usr/local/intel/lib
  NETCDF_INC = /usr/local/intel/include
  LIBS = -lnetcdff -lnetcdf -lmkl_lapack95_lp64
  LIBDIR = -L${NETCDF_FORTRAN_LIB} -L${NETCDF_LIB}
  INCDIR = -I. -I${NETCDF_FORTRAN_INC}
endif

ifeq ($(arch),ifx-debug)
  FC = ifx
  FFLAGS = -assume byterecl -mcmodel=medium -O0 -qmkl -warn all -check all -fpe0 -CB -traceback -g
  DEFS = -DDOUBLE -DMKL
  NETCDF_FORTRAN_LIB = /usr/local/intel/lib
  NETCDF_FORTRAN_INC = /usr/local/intel/include
  NETCDF_LIB = /usr/local/intel/lib
  NETCDF_INC = /usr/local/intel/include
  LIBS = -lnetcdff -lnetcdf -lmkl_lapack95_lp64
  LIBDIR = -L${NETCDF_FORTRAN_LIB} -L${NETCDF_LIB}
  INCDIR = -I. -I${NETCDF_FORTRAN_INC} -I${NETCDF_INC}
endif

ifeq ($(arch),gfortran-debug)
  FC = gfortran
  FFLAGS = -O0 -Wall -fbounds-check -fbacktrace -Wuninitialized -g
  DEFS = -DDOUBLE -DSORATENA
  NETCDF_FORTRAN_LIB = /usr/lib/x86_64-linux-gnu
  NETCDF_FORTRAN_INC = /usr/include
  NETCDF_LIB = ${NETCDF_FORTRAN_LIB}
  NETCDF_INC = ${NETCDF_FORTRAN_INC}
  LIBS = -lnetcdff -lnetcdf -llapack95 -llapack -lblas
  LIBDIR = -L${NETCDF_FORTRAN_LIB} -L/usr/local/lib
  INCDIR = -I. -I${NETCDF_FORTRAN_INC} -I/usr/local/include
endif

ifeq ($(arch),gfortran)
  FC = gfortran
  FFLAGS = -O2 -frecursive
  DEFS = -DDOUBLE -DPARTICLEVELOCITY
  NETCDF_FORTRAN_LIB = /usr/lib/x86_64-linux-gnu
  NETCDF_FORTRAN_INC = /usr/include
  NETCDF_LIB = ${NETCDF_FORTRAN_LIB}
  NETCDF_INC = ${NETCDF_FORTRAN_INC}
  LIBS = -lnetcdff -lnetcdf -llapack95 -llapack -lblas
  LIBDIR = -L${NETCDF_FORTRAN_LIB} -L/usr/local/lib
  INCDIR = -I. -I${NETCDF_FORTRAN_INC} -I/usr/local/include
endif

aeluma = AutomatedEventLocationUsingaMeshofArrays.F90
calc_minmax_waveform_grd = calc_minmax_waveform_grd.F90
calc_bpf_coef = calc_bpf_coef.F90
calc_bpf_order = calc_bpf_order.F90
calc_kernelmatrix = calc_kernelmatrix.F90
calc_lpf_coef = calc_lpf_coef.F90
calc_lpf_order = calc_lpf_order.F90
calc_froude_number_grd = calc_froude_number_grd.F90
correlation = correlation.F90
constants = constants.F90
cosine_taper = cosine_taper.F90
deconvolution = deconvolution.F90
fftsg = fftsg.f
geometry = geometry.f90
geompack2 = geompack2.f90
gradiometry_parameters = gradiometry_parameters.F90
greatcircle = greatcircle.f90
grdfile_io = grdfile_io.F90
itoa = itoa.F90
line_fit = line_fit.f90
lonlat_xy_conv = lonlat_xy_conv.F90
make_waveform_from_grd = make_waveform_from_grd.F90
nrtype = nrtype.F90
read_sacfile = read_sacfile.F90
sac_decimation = sac_decimation.F90
sac_deconvolve = sac_deconvolve.F90
sac_integrate = sac_integrate.F90
seismicgradiometry = seismicgradiometry.F90
seismicgradiometry_reducingvelocity2 = seismicgradiometry_reducingvelocity2.F90
sort = sort.F90
tandem = tandem.F90
typedef = typedef.F90

##Module
mod_calc_kernelmatrix = $(calc_kernelmatrix:.F90=.mod)
mod_constants = $(constants:.F90=.mod)
mod_correlation = $(correlation:.F90=.mod)
mod_cosine_taper = taper.mod
mod_deconvolution = $(deconvolution:.F90=.mod)
mod_gradiometry_parameters = $(gradiometry_parameters:.F90=.mod)
mod_greatcircle = $(greatcircle:.f90=.mod)
mod_grdfile_io = $(grdfile_io:.F90=.mod)
mod_itoa = $(itoa:.F90=.mod)
mod_lonlat_xy_conv = $(lonlat_xy_conv:.F90=.mod)
mod_nrtype = $(nrtype:.F90=.mod)
mod_read_sacfile = $(read_sacfile:.F90=.mod)
mod_tandem = $(tandem:.F90=.mod)
mod_typedef = $(typedef:.F90=.mod)

##Object
o_aeluma = $(aeluma:.F90=.o)
o_calc_bpf_coef = $(calc_bpf_coef:.F90=.o)
o_calc_bpf_order = $(calc_bpf_order:.F90=.o)
o_calc_lpf_coef = $(calc_lpf_coef:.F90=.o)
o_calc_lpf_order = $(calc_lpf_order:.F90=.o)
o_calc_minmax_waveform_grd = $(calc_minmax_waveform_grd:.F90=.o)
o_calc_froude_number_grd = $(calc_froude_number_grd:.F90=.o)
o_calc_kernelmatrix = $(calc_kernelmatrix:.F90=.o)
o_constants = $(constants:.F90=.o)
o_correlation = $(correlation:.F90=.o)
o_cosine_taper = $(cosine_taper:.F90=.o)
o_deconvolution = $(deconvolution:.F90=.o)
o_fftsg = $(fftsg:.f=.o)
o_gradiometry_parameters = $(gradiometry_parameters:.F90=.o)
o_grdfile_io = $(grdfile_io:.F90=.o)
o_greatcircle = $(greatcircle:.f90=.o)
o_itoa = $(itoa:.F90=.o)
o_lonlat_xy_conv = $(lonlat_xy_conv:.F90=.o)
o_make_waveform_from_grd = $(make_waveform_from_grd:.F90=.o)
o_nrtype = $(nrtype:.F90=.o)
o_read_sacfile = $(read_sacfile:.F90=.o)
o_sac_decimation = $(sac_decimation:.F90=.o)
o_sac_deconvolve = $(sac_deconvolve:.F90=.o)
o_sac_integrate = $(sac_integrate:.F90=.o)
o_seismicgradiometry = $(seismicgradiometry:.F90=.o)
o_seismicgradiometry_reducingvelocity2 = $(seismicgradiometry_reducingvelocity2:.F90=.o)
o_sort = $(sort:.F90=.o)
o_tandem = $(tandem:.F90=.o)
o_deconvolution = $(deconvolution:.F90=.o)
o_line_fit = $(line_fit:.f90=.o)
o_typedef = $(typedef:.F90=.o)
o_geompack2 = $(geompack2:.f90=.o)
o_geometry = $(geometry:.f90=.o)

##Module dependency
$(mod_constants): $(constants) $(o_constants)
$(mod_correlation): $(correlation) $(o_correlation)
$(mod_cosine_taper): $(cosine_taper) $(o_cosine_taper)
$(mod_gradiometry_parameters): $(gradiometry_parameters) $(o_gradiometry_parameters)
$(mod_greatcircle): $(greatcircle) $(o_greatcircle)
$(mod_itoa): $(itoa) $(o_itoa)
$(mod_nrtype): $(nrtype) $(o_nrtype)
$(mod_lonlat_xy_conv): $(lonlat_xy_conv) $(o_lonlat_xy_conv)
$(mod_read_sacfile): $(read_sacfile) $(o_read_sacfile)
$(mod_grdfile_io): $(grdfile_io) $(o_grdfile_io)
$(mod_deconvolution): $(deconvolution) $(o_deconvolution)
$(mod_tandem): $(tandem) $(o_tandem)
$(mod_typedef): $(typedef) $(o_typedef)
$(mod_calc_kernelmatrix): $(calc_kernelmatrix) $(o_calc_kernelmatrix)

##Object dependency
$(o_calc_bpf_coef): $(calc_bpf_coef) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_calc_bpf_order): $(calc_bpf_order) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_calc_kernelmatrix): $(calc_kernelmatrix) $(mod_nrtype) $(mod_constants) $(mod_typedef) $(mod_greatcircle) $(mod_sort) \
        $(mod_gradiometry_parameters) \
	$(o_constants) $(o_gradiometry_parameters) $(o_sort)
$(o_calc_lpf_coef): $(calc_lpf_coef) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_calc_lpf_order): $(calc_lpf_order) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_calc_froude_number_grd): $(calc_froude_number_grd) $(mod_nrtype) $(mod_constants) $(mod_typedef) $(mod_calc_kernelmatrix) \
        $(mod_grdfile_io)
$(o_constants): $(constants) $(mod_nrtype)
$(o_correlation): $(correlation) $(mod_nrtype)
$(o_cosine_taper): $(cosine_taper) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_deconvolution): $(deconvolution) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_fftsg): $(fftsg)
$(o_geometry): $(geometry) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_geompack2): $(geompack2) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_gradiometry_parameters): $(gradiometry_parameters) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_grdfile_io): $(grdfile_io) $(mod_nrtype)
$(o_greatcircle): $(greatcircle) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_itoa): $(itoa)
$(o_lonlat_xy_conv): $(lonlat_xy_conv) $(mod_nrtype) $(mod_constants) $(o_constants)
$(o_make_waveform_from_grd): $(make_waveform_from_grd) $(mod_nrtype) $(mod_grdfile_io) $(mod_gradiometry_parameters)
$(o_nrtype): $(nrtype)
$(o_read_sacfile): $(read_sacfile) $(mod_nrtype)
$(o_sac_decimation): $(sac_decimation) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_tandem)
$(o_sac_deconvolve): $(sac_deconvolve) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_deconvolution)
$(o_sac_integrate): $(sac_integrate) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_tandem)
$(o_seismicgradiometry): $(seismicgradiometry) $(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_grdfile_io) \
	$(mod_lonlat_xy_conv) $(mod_typedef) $(mod_gradiometry_parameters) $(mod_calc_kernelmatrix) $(mod_tandem) \
        $(o_gradiometry_parameters) $(o_constants)
$(o_aeluma): $(aeluma) \
	$(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_grdfile_io) $(mod_cosine_taper) $(mod_tandem) \
	$(mod_lonlat_xy_conv) $(mod_typedef) $(mod_gradiometry_parameters) $(mod_calc_kernelmatrix) $(mod_cosine_taper) \
        $(o_gradiometry_parameters) $(o_constants)
$(o_seismicgradiometry_reducingvelocity2): $(seismicgradiometry_reducingvelocity2) \
	$(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_grdfile_io) $(mod_tandem) $(mod_itoa) \
	$(mod_lonlat_xy_conv) $(mod_typedef) $(mod_gradiometry_parameters) $(mod_calc_kernelmatrix) $(o_gradiometry_parameters)
$(o_calc_minmax_waveform_grd): $(calc_minmax_waveform_grd) \
	$(mod_nrtype) $(mod_constants) $(mod_read_sacfile) $(mod_grdfile_io) $(mod_tandem) \
	$(mod_lonlat_xy_conv) $(mod_typedef) $(mod_gradiometry_parameters) $(mod_calc_kernelmatrix)
$(o_tandem): $(tandem) $(mod_nrtype)
$(o_typedef): $(typedef) $(mod_nrtype)

seismicgradiometry: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_lonlat_xy_conv) \
	$(o_grdfile_io) $(o_read_sacfile) $(o_sort) $(o_greatcircle) $(o_seismicgradiometry) $(o_gradiometry_parameters) \
        $(o_typedef) $(o_calc_kernelmatrix) $(o_geompack2) $(o_geometry)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

aeluma: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) \
	$(o_lonlat_xy_conv) $(o_itoa) $(o_correlation) \
	$(o_grdfile_io) $(o_read_sacfile) $(o_sort) $(o_greatcircle) $(o_gradiometry_parameters) \
	$(o_cosine_taper) $(o_typedef) $(o_calc_kernelmatrix) $(o_geompack2) $(o_geometry) $(o_fftsg) \
	$(o_aeluma)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

seismicgradiometry_reducingvelocity2: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) \
	$(o_lonlat_xy_conv) $(o_itoa) \
	$(o_grdfile_io) $(o_read_sacfile) $(o_sort) $(o_greatcircle) $(o_seismicgradiometry_reducingvelocity2) \
	$(o_gradiometry_parameters) $(o_typedef) $(o_calc_kernelmatrix) $(o_geompack2) $(o_geometry)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

calc_minmax_waveform_grd: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) \
	$(o_lonlat_xy_conv) \
	$(o_grdfile_io) $(o_read_sacfile) $(o_sort) $(o_greatcircle) $(o_calc_minmax_waveform_grd) \
	$(o_gradiometry_parameters) $(o_typedef) $(o_calc_kernelmatrix) $(o_geompack2) $(o_geometry)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

sac_decimation: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_read_sacfile) \
	$(o_sac_decimation)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

sac_deconvolve: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_read_sacfile) \
	$(o_deconvolution) $(o_sac_deconvolve)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

sac_integrate: $(o_nrtype) $(o_constants) $(o_calc_bpf_order) $(o_calc_bpf_coef) $(o_tandem) $(o_read_sacfile) \
	$(o_line_fit) $(o_sac_integrate)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

make_waveform_from_grd: $(o_nrtype) $(o_grdfile_io) $(o_gradiometry_parameters) $(o_make_waveform_from_grd)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

calc_froude_number_grd: $(o_nrtype) $(o_constants) $(o_gradiometry_parameters) $(o_lonlat_xy_conv) $(o_typedef) \
        $(o_calc_kernelmatrix) $(o_geompack2) $(o_geometry) $(o_grdfile_io) $(o_calc_froude_number_grd)
	$(FC) $^ -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

.F90.mod:
	@:
.f90.mod:
	@:

.F90.o:
	$(FC) $< -c -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)
.f90.o:
	$(FC) $< -c -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)
.f.o:
	$(FC) $< -c -o $@ $(FFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) $(DEFS)

clean:
	rm -f *.o *.mod *__genmod.f90

