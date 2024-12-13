! Parameters of AELUMA method
!! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

module aeluma_parameters
  use nrtype, only : fp
  implicit none
  private

  real(kind = fp), public, parameter :: eps = 1.0e-10_fp

  
  !!F-net aeluma method, read waveforms from shmdump
  !!temporal origin: land-area centroid
  real(kind = fp), public, parameter :: center_lon_aeluma = 137.6336_fp, center_lat_aeluma = 37.4581_fp
  integer,         public, parameter :: nadd_station_aeluma = 1
  
  !!Bandpass filter
  real(kind = fp), public, parameter :: fl = 1.0_fp / 50.0_fp, fh = 1.0_fp / 2.0_fp, fs = fh * 2.0_fp, &
  &                                     ap = 0.5_fp, as = 10.0_fp

  integer,         public, parameter :: ntime = 2400
  integer,         public, parameter :: ntimestep = 4   

  real(kind = fp), public, parameter :: interstationdistance_min = 0.05_fp
  real(kind = fp), public, parameter :: cutoff_dist = 100.0_fp
  integer,         public, parameter :: ntime_fft = 2048, ntime_fft2 = ntime_fft / 2
  integer,         public, parameter :: nsec_buf = 100
  real(kind = fp), public, parameter :: xcorr_min = 0.6_fp
  real(kind = fp), public, parameter :: cos_taper_ratio = 0.1_fp
  real(kind = fp), public, parameter :: lagtime_max = 20.0_fp * 60.0_fp, lagtime_min = -lagtime_max
 
end module aeluma_parameters
  
