! Parameters of AELUMA method
!! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

module aeluma_parameters
  use nrtype
  use constants
  implicit none
  private

  !!F-net aeluma method, read waveforms from shmdump
  !!temporal origin: land-area centroid
  real(kind = fp),    public, parameter :: center_lon_aeluma = 137.6336_fp, center_lat_aeluma = 37.4581_fp
  integer,            public, parameter :: nadd_station_aeluma = 0
  
  !!Bandpass filter
  real(kind = fp),    public, parameter :: fl = 1.0_fp / 50.0_fp, fh = 1.0_fp / 20.0_fp, fs = fh * 2.0_fp, &
  &                                        ap = 0.5_fp, as = 10.0_fp

  real(kind = fp),    public, parameter :: interstationdistance_min = 5.0_fp
  real(kind = fp),    public, parameter :: cutoff_dist = 100.0_fp
  real(kind = fp),    public, parameter :: dtimestep = 1.0_fp
  integer,            public, parameter :: nsampling_int = 2
  integer,            public, parameter :: sampling_int(1 : nsampling_int) = [100, 20]
  real(kind = fp),    public, parameter :: sampling_sec(1 : nsampling_int) &
  &                                        = 1.0_fp / real(sampling_int(1 : nsampling_int), kind = fp)
  integer,            public, parameter :: sampling_int_use = 5
  integer,            public, parameter :: nsec_buf  = 360
  integer,            public, parameter :: ntime_fft = 2048, ntime_fft2 = ntime_fft / 2
  integer,            public, parameter :: waveform_buf_index_max = nsec_buf * sampling_int_use
  real(kind = fp),    public, parameter :: xcorr_min = 0.02_fp
  real(kind = fp),    public, parameter :: lagtime_ratio_threshold = 0.7_fp
  real(kind = fp),    public, parameter :: cos_taper_ratio = 0.1_fp
  real(kind = fp),    public, parameter :: lagtime_max = 200.0_fp, lagtime_min = -lagtime_max
  real(kind = fp),    public, parameter :: order = 1.0e+6_fp
  integer,            public, parameter :: nwinch = 65536

  integer,            public, parameter :: iwin_wave = 0  !!AELUMA_shmudmp.F90
  integer,            public, parameter :: iwin_map = 0, iwin_legend = 1, iwin_eplist = 2  !!plot_map_vector.F90
  real(kind = sp),    public, parameter :: plot_dy_eplist = 5.0_sp
  real(kind = sp),    public, parameter :: plot_x_eplist = 1.0_sp, plot_y_eplist = 27.0_sp
  !!For plot_map_vector.F90
  real(kind = sp),    public, parameter :: width = 300.0_sp, height = 300.0_sp, scale = 1.0_sp
  real(kind = sp),    public, parameter :: vector_len_ref = 0.5, vector_width = 1.5, vector_head1 = 2.5, vector_head2 = 4.0
  real(kind = fp),    public, parameter :: lon_w = 120.0_fp, lon_e = 149.0_fp, &
  &                                        lat_s = 22.5_fp,  lat_n = 48.0_fp, center_lon = 135.0_fp
  real(kind = sp),    public, parameter :: likelihood_legend_normalize   = 2.0e+3
  character(len = 5), public, parameter :: likelihood_legend_normalize_c = "x5e-4"

  !!Color palette
  integer,            public, parameter :: color_likelihood(1 : 3, 1 : 10) = reshape([255, 255, 204, &  !!scm/lajolla -T0/8/1
  &                                                                                   253, 245, 175, &
  &                                                                                   247, 216, 117, &
  &                                                                                   237, 174,  86, &
  &                                                                                   229, 136,  81, &
  &                                                                                   208,  96,  76, &
  &                                                                                   155,  68,  63, &
  &                                                                                    99,  51,  40, &
  &                                                                                    48,  34,  16, &
  &                                                                                    26,  26,   1], [3, 10])
  integer,            public, parameter :: color_correlation(1 : 3, 1 : 10) = reshape([230, 230, 240, & !!scm/acton -T0/8/1 -I
  &                                                                                    223, 213, 228, &
  &                                                                                    214, 181, 206, &
  &                                                                                    212, 153, 187, &
  &                                                                                    209, 124, 166, &
  &                                                                                    177, 103, 149, &
  &                                                                                    140,  98, 142, &
  &                                                                                     99,  77, 121, &
  &                                                                                     62,  48,  91, &
  &                                                                                     46,  33,  77], [3, 10])

  !!location estimation
  integer,            public, parameter :: nparticle = 1000, niter = 3, iteration_count_max = 1000, nepicenter = 4
  real(kind = fp),    public, parameter :: sigma_azdiff = 15.0_fp * deg2rad, sigma_azdiff2 = sigma_azdiff ** 2
  real(kind = fp),    public, parameter :: sameaz_num = 10.0_fp * deg2rad, sameaz_num2 = 20.0_fp ** 2, azweight_coef = 0.7_fp
  real(kind = fp),    public, parameter :: sigma_particle = 0.3_fp
  real(kind = fp),    public, parameter :: sigma_dist = log(100.0_fp), ot_coef = 0.7_fp, sigma_otdiff = 60.0_fp, &
  &                                        sigma_dist2 = sigma_dist ** 2, sigma_otdiff2 = sigma_otdiff ** 2
  real(kind = fp),    public, parameter :: correlation_threshold = 0.4_fp
  real(kind = fp),    public, parameter :: ref_appvelocity = 0.3_fp
  real(kind = fp),    public, parameter :: max_slowness = 0.4_fp

  integer,            public, parameter :: narray_use_min = 9
  integer,            public, parameter :: epicenter_acceptcount_threshold = 90
  integer,            public, parameter :: epicenter_renew_threshold = 150
  real(kind = fp),    public, parameter :: min_likelihood_eqobs = 0.5_fp / (pi * sigma_azdiff * sigma_otdiff) &
  &                                                             * exp(-0.5_fp * (2.5_fp ** 2))

  integer,            public, parameter :: nthread_min = 1

  integer,            public, parameter :: ntime = 6000, ntimestep = 10
end module aeluma_parameters
  
