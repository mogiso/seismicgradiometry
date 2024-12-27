! Parameters of AELUMA method
!! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

module aeluma_parameters
  use nrtype, only : fp
  implicit none
  private

  !!F-net aeluma method, read waveforms from shmdump
  !!temporal origin: land-area centroid
  real(kind = fp), public, parameter :: center_lon_aeluma = 137.6336_fp, center_lat_aeluma = 37.4581_fp
  integer,         public, parameter :: nadd_station_aeluma = 0
  
  !!Bandpass filter
  real(kind = fp), public, parameter :: fl = 1.0_fp / 50.0_fp, fh = 1.0_fp / 20.0_fp, fs = fh * 2.0_fp, &
  &                                     ap = 0.5_fp, as = 10.0_fp

  real(kind = fp), public, parameter :: interstationdistance_min = 5.0_fp
  real(kind = fp), public, parameter :: cutoff_dist = 100.0_fp
  integer,         public, parameter :: nsampling_int = 2
  integer,         public, parameter :: sampling_int(1 : nsampling_int) = [100, 20]
  real(kind = fp), public, parameter :: sampling_sec(1 : nsampling_int) &
  &                                     = 1.0_fp / real(sampling_int(1 : nsampling_int), kind = fp)
  integer,         public, parameter :: sampling_int_use = 5
  integer,         public, parameter :: nsec_buf  = 360
  integer,         public, parameter :: ntime_fft = 2048, ntime_fft2 = ntime_fft / 2
  integer,         public, parameter :: waveform_buf_index_max = nsec_buf * sampling_int_use
  real(kind = fp), public, parameter :: xcorr_min = 0.2_fp
  real(kind = fp), public, parameter :: cos_taper_ratio = 0.1_fp
  real(kind = fp), public, parameter :: lagtime_max = 200.0_fp, lagtime_min = -lagtime_max
  real(kind = fp), public, parameter :: order = 1.0e+6_fp
  integer,         public, parameter :: nwinch = 65536

  integer,         public, parameter :: ntime = 6000, ntimestep = 10
end module aeluma_parameters
  
