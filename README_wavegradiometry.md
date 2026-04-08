# mogiso/seismicgradiometry: wave gradiometry

## Description
This program package contains the code for wave gradiometry, one of the array analyses methods.

Language: Fortran 90
Input waveform format: SAC (seismic analysis code), binary

When excute the file, firstly it reads waveform files (sac-binary formatted). Then, it constructs array network from the geometry of stations
read from waveform files.

Require: [NetCDF-Fortran](https://www.unidata.ucar.edu/software/netcdf/docs-fortran/index.html "NetCDF-Fortran"),
and [LAPACK95](http://www.netlib.org/lapack95/ "LAPACK95").
NetCDF-Fortran depends on [NetCDF](https://www.unidata.ucar.edu/software/netcdf/),
and LAPACK95 depends on [LAPACK](http://www.netlib.org/lapack/). [Intel Math Kernel Library](https://software.intel.com/mkl) can be used instead of LAPACK95.

I used this program package in Ogiso and Tsushima (2022; 2025). Please refer these papers for the detail of the method.

[Ogiso, M. and H. Tsushima, Ocean-wave gradiometry: Visualizing and extracting propagation features of the 15 January 2022 tsunami wavefield with dense ocean-bottom pressure gauge arrays, Seismo. Res. Lett., 94, 626-636, 2022.](https://doi.org/10.1785/0220240358)

[Ogiso, M. and H. Tsushima, Propagation characteristics of tsunamis in an atmosphere-ocean coupled system revealed by wave gradiometry: The 2022 Hunga Tonga–Hunga Ha’apai eruption case, Seismo. Res. Lett., 96, 744-757, 2025.](https://doi.org/10.1785/0220240358)

## Compile
  $ make -f Makefile seismicgradiometry_reducingvelocity2

## Parameters
Analyses parameters are listed in gradiometry_parameters.F90
