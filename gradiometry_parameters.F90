!! Parameters of wave gradiometry
!! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

module gradiometry_parameters
  use nrtype, only : fp
  use constants, only : deg2rad, pi
  implicit none
  private

  real(kind = fp), public, parameter :: eps = 1.0e-10_fp

  
#ifdef SORATENA
  !!S-net region
  real(kind = fp), public, parameter :: x_start = -350.0_fp, y_start = -600.0_fp, &
  &                                     x_end   = 350.0_fp,  y_end   = 600.0_fp
  real(kind = fp), public, parameter :: center_lon = 142.5_fp, center_lat = 38.25_fp 
  real(kind = fp), public, parameter :: dgrid_x = 20.0_fp, dgrid_y = 20.0_fp
  !!For 15s-sampled data 
  integer,         public, parameter :: ntime = 2400
  integer,         public, parameter :: ntimestep = 4 
  integer,         public, parameter :: ngradient2 = 240
  integer,         public, parameter :: ntime_slowness = 241
  !!data order, array configuration
  real(kind = fp), public, parameter :: order = 1.0_fp !!already hpa
  integer,         public, parameter :: nsta_grid_min = 15, nsta_grid_max = 23  !!For S-net/DONET OBPG array
  integer,         public, parameter :: naddstation_array = 20
  real(kind = fp), public, parameter :: cutoff_dist = 80.0_fp
  real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  !!bandpass filter
  real(kind = fp), public, parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (10.0_fp * 60.0_fp),  &
  &                                     fs = 1.0_fp / (5.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (40.0_fp * 60.0_fp), fh = 1.0_fp / (10.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (5.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (180.0_fp * 60.0_fp), fh = 1.0_fp / (4.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (3.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp

#elif defined (SKNET)

  !!SK-net region
  !!SK-net long-period ground motion (5-10s)
  real(kind = fp), public, parameter :: eps = 1.0e-6_fp
  real(kind = fp), public, parameter :: x_start = -250.0_fp, y_start = -250.0_fp, &
  &                                     x_end   =  205.0_fp, y_end   =  225.0_fp
  real(kind = fp), public, parameter :: center_lon = 139.0_fp, center_lat = 36.0_fp
  real(kind = fp), public, parameter :: dgrid_x = 2.5_fp, dgrid_y = 2.5_fp
  !!For 0.01s-sampled data  (SK-net)
  integer,         public, parameter :: ntime = 18000
  integer,         public, parameter :: ntimestep = 100
  integer,         public, parameter :: ngradient2 = 1500
  integer,         public, parameter :: ntime_slowness = 60
  !!data order, array configuration
  real(kind = fp), public, parameter :: order = 1.0e-6_fp !!nm -> mm
  integer,         public, parameter :: nsta_grid_min = 10, nsta_grid_max = 4
  integer,         public, parameter :: naddstation_array = 1
  real(kind = fp), public, parameter :: cutoff_dist = 20.0_fp
  real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  !!Bandpass filter
  real(kind = fp), public, parameter :: fl = 1.0_fp / 10.0_fp, fh = 1.0_fp / 5.0_fp, fs = fh * 2.0_fp, &
  &                                     ap = 0.5_fp, as = 10.0_fp

#else

  !!S-net region
  real(kind = fp), public, parameter :: x_start = -350.0_fp, y_start = -600.0_fp, &
  &                                     x_end   = 350.0_fp,  y_end   = 600.0_fp
  real(kind = fp), public, parameter :: center_lon = 142.5_fp, center_lat = 38.25_fp 
  real(kind = fp), public, parameter :: dgrid_x = 20.0_fp, dgrid_y = 20.0_fp
  !!DONET region
  !real(kind = fp), public, parameter :: x_start = -150.0_fp, y_start = -100.0_fp, &
  !&                                     x_end   = 150.0_fp,  y_end   = 100.0_fp
  !real(kind = fp), public, parameter :: center_lon = 135.75_fp, center_lat = 33.2_fp 
  !real(kind = fp), public, parameter :: dgrid_x = 10.0_fp, dgrid_y = 10.0_fp
  !!data order, array configuration
  !!For 15s-sampled data 
  integer,         public, parameter :: ntime = 2400
  integer,         public, parameter :: ntimestep = 4 
  integer,         public, parameter :: ngradient2 = 240
  integer,         public, parameter :: ntime_slowness = 241
  !!For 6s-sampled data
  !integer,         public, parameter :: ntime = 6000
  !integer,         public, parameter :: ntimestep = 10
  !integer,         public, parameter :: ngradient2 = 600
  !integer,         public, parameter :: ntime_slowness = 601, ntime_slowness2 = (ntime_slowness - 1) / 2
  !real(kind = fp), public, parameter :: order = 1.0e-2_fp !!Pa -> hPa
  real(kind = fp), public, parameter :: order = 1.0_fp 
  integer,         public, parameter :: nsta_grid_min = 4, nsta_grid_max = 4
  integer,         public, parameter :: naddstation_array = 1
  real(kind = fp), public, parameter :: cutoff_dist = 80.0_fp !!S-net
  !!2023Philippines
  !real(kind = fp), public, parameter :: order = 1.0_fp !!2023Philippines: already hPa
  !integer,         public, parameter :: nsta_grid_min = 4, nsta_grid_max = 6
  !integer,         public, parameter :: naddstation_array = 3
  !real(kind = fp), public, parameter :: cutoff_dist = 60.0_fp  !!DONET
  real(kind = fp), public, parameter :: az_diff_max = 150.0_fp * deg2rad
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (20.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (10.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp
  real(kind = fp), public, parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (10.0_fp * 60.0_fp), &
  &                                     fs = 1.0_fp / (5.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp
  !real(kind = fp), public, parameter :: fl = 1.0_fp / (180.0_fp * 60.0_fp), fh = 1.0_fp / (30.0_fp * 60.0_fp), &
  !&                                     fs = 1.0_fp / (5.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp

#endif


  real(kind = fp), public, parameter :: error_omega = 2.0_fp * pi * 1.0_fp / ((1.0_fp / fl + 1.0_fp / fh) * 0.5_fp)

#ifdef ELLIPSE
  !!sigma_x: normal to propagation direction, sigma_y: parallel to propagation direction
  real(kind = fp), public, parameter :: sigma_x = cutoff_dist * 2.0_fp, sigma_y = cutoff_dist
  real(kind = fp), public, parameter :: evlon = 137.8910_fp, evlat = 36.6928_fp
#endif

  !!threshold of interstation distance (vertix distance) for delaunay triangulation
  real(kind = fp), public, parameter :: interstationdistance_min = 0.05_fp

  !!common constants
  integer,         public, parameter :: ngrid_x = int((x_end - x_start) / dgrid_x) + 1
  integer,         public, parameter :: ngrid_y = int((y_end - y_start) / dgrid_y) + 1
  integer,         public, parameter :: niteration_max = 15

#ifdef GREEN_CORRECTION
  real(kind = fp), public, parameter :: depth_ref = 0.1_fp  !!in km
#endif

  real(kind = fp), public, parameter :: grav_acc = 980.655_fp !!in cm/s^2

end module gradiometry_parameters
  
