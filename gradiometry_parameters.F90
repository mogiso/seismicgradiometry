!! Parameters of wave gradiometry
!! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

module gradiometry_parameters
  use nrtype, only : fp
  use constants, only : deg2rad
  implicit none
  private

  !!S-net OBP long-period (20-60min.)
  !real(kind = fp), public, parameter :: order = 1.0e-2_fp !!Pa -> hpa
  !real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  !real(kind = fp), public, parameter :: x_start = -350.0_fp, y_start = -600.0_fp, &
  !&                                     x_end = 350.0_fp, y_end = 600.0_fp
  !real(kind = fp), public, parameter :: center_lon = 142.5_fp, center_lat = 38.25_fp 
  !real(kind = fp), public, parameter :: dgrid_x = 20.0_fp, dgrid_y = 20.0_fp
  !real(kind = fp), public, parameter :: cutoff_dist = 80.0_fp
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (20.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (10.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp
  !integer,         public, parameter :: ntime_slowness = 61, ntime_slowness2 = (ntime_slowness - 1) / 2
  !integer,         public, parameter :: nsta_grid_max = 40, nsta_grid_min = 5  !!For S-net/DONET OBPG array
  !integer,         public, parameter :: ntime = 630
  !integer,         public, parameter :: ntime = 1024 !!testdata

  !!DONET OBP long-period (20-60min.)
  !real(kind = fp), public, parameter :: order = 1.0e-2_fp  !!Pa -> hpa
  !real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  !real(kind = fp), public, parameter :: x_start = -150.0_fp, y_start = -100.0_fp, &
  !&                                     x_end = 150.0_fp, y_end = 100.0_fp
  !real(kind = fp), public, parameter :: center_lon = 135.75_fp, center_lat = 33.2_fp 
  !real(kind = fp), public, parameter :: dgrid_x = 10.0_fp, dgrid_y = 10.0_fp
  !real(kind = fp), public, parameter :: cutoff_dist = 80.0_fp
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (20.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (10.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp
  !integer,         public, parameter :: ntime_slowness = 61, ntime_slowness2 = (ntime_slowness - 1) / 2
  !integer,         public, parameter :: nsta_grid_max = 40, nsta_grid_min = 5  !!For S-net/DONET OBPG array
  !integer,         public, parameter :: ntime = 630
  !integer,         public, parameter :: ntime = 1024 !!testdata

  !!DONET OBP short-period (6-20min.)
  !real(kind = fp), public, parameter :: order = 1.0e-2_fp  !!Pa -> hpa
  !real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  !real(kind = fp), public, parameter :: x_start = -150.0_fp, y_start = -100.0_fp, &
  !&                                     x_end = 150.0_fp, y_end = 100.0_fp
  !real(kind = fp), public, parameter :: center_lon = 135.75_fp, center_lat = 33.2_fp
  !real(kind = fp), public, parameter :: dgrid_x = 10.0_fp, dgrid_y = 10.0_fp
  !real(kind = fp), public, parameter :: cutoff_dist = 30.0_fp
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (20.0_fp * 60.0_fp), fh = 1.0_fp / (6.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (3.0_fp * 60.0_fp), ap = 0.5_fp, as = 10.0_fp
  !integer,         public, parameter :: ntime_slowness = 61, ntime_slowness2 = (ntime_slowness - 1) / 2
  !integer,         public, parameter :: nsta_grid_max = 40, nsta_grid_min = 5  !!For S-net/DONET OBPG array
  !integer,         public, parameter :: ntime = 630
  !integer,        public, parameter :: ntime = 1024 !!testdata

  !!SK-net long-period ground motion (5-10s)
  real(kind = fp), public, parameter :: order = 1.0e-6_fp  !!nm -> mm
  real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  real(kind = fp), public, parameter :: x_start = -250.0_fp, y_start = -250.0_fp, &
  &                                     x_end = 205.0_fp, y_end = 225.0_fp
  real(kind = fp), public, parameter :: center_lon = 139.0_fp, center_lat = 36.0_fp
  real(kind = fp), public, parameter :: dgrid_x = 2.5_fp, dgrid_y = 2.5_fp
  real(kind = fp), public, parameter :: cutoff_dist = 10.0_fp
  integer,         public, parameter :: nsta_grid_max = 10, nsta_grid_min = 3
  integer,         public, parameter :: ntime_slowness = 60, ntime_slowness2 = (ntime_slowness - 1) / 2
  real(kind = fp), public, parameter :: fl = 1.0_fp / 10.0_fp, fh = 1.0_fp / 5.0_fp, fs = fh * 2.0_fp, &
  &                                     ap = 0.5_fp, as = 10.0_fp
  integer,         public, parameter :: ntime = 27000
#ifdef ELLIPSE
  !!sigma_x: normal to propagation direction, sigma_y: parallel to propagation direction
  real(kind = fp), public, parameter :: sigma_x = cutoff_dist * 2.0_fp, sigma_y = cutoff_dist
  real(kind = fp), public, parameter :: evlon = 137.8910_fp, evlat = 36.6928_fp
#endif


  !!common constants
  integer,         public, parameter :: ngrid_x = int((x_end - x_start) / dgrid_x) + 1
  integer,         public, parameter :: ngrid_y = int((y_end - y_start) / dgrid_y) + 1
  !!threshold of interstation distance (vertix distance) for delaunay triangulation
  real(kind = fp), public, parameter :: interstationdistance_min = 0.05_fp
  real(kind = fp), public, parameter :: cos_taper_ratio = 0.1_fp
  real(kind = fp), public, parameter :: xcorr_min = 0.7_fp

  integer,         public, parameter :: ntimestep = 100
  integer,         public, parameter :: naddstation_array = 1
  integer,         public, parameter :: niteration_max = 5
  integer,         public, parameter :: ngradient = 1000, ngradient2 = ngradient / 2, ngradient4 = ngradient / 4

  !!for AELUMA method
  integer,         public, parameter :: ntime_fft = 1024
 
end module gradiometry_parameters
  
