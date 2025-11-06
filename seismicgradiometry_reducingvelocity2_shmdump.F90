!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

program seismicgradiometry_reducingvelocity2_shmdump
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use typedef
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy, xy2bl
  use gradiometry_parameters
  use calc_kernelmatrix
  use itoa
  use tandem

  implicit none

  integer :: nstation, nch, nsample, timeindex_diff_min, i, j, ii, jj, kk, ncount, n, ios, ngrad, ndecimate
  real(kind = fp)              :: dt = 1.0_fp / real(sampling_int_use, kind = fp)
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)
  type(location)               :: location_sta(1 : nwinch)
  real(kind = fp)              :: dx_east, dy_north, denominator, uu, uut, utut, max_innerproduct, &
  &                               innerproduct_tmp, relativeerror, &
  &                               uxu(1 : 2), uxut(1 : 2), numerator_slowness(1 : 2), numerator_ampterm(1 : 2), &
  &                               obs_vector(1 : nsta_grid_max), &
  &                               sigma_slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               sigma_ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), & 
  &                               waveform_real(1 : maxval(sampling_int)), &
  &                               waveform_est_tmp(1 : 4, 1 : ngradient2), &
  &                               waveform_est_tmp2(1 : ngradient2, 1 : 4), &
  &                               waveform_est_plot(1 : ngrid_x, 1 : ngrid_y), &
  &                               slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               slowness_correction(1 : 2, 1 : ngrid_x, 1 : ngrid_y), slowness_cor_prev(1 : 2), &
  &                               kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                               error_matrix(1 : 3, 1 : ngrid_x, 1 : ngrid_y), &
  &                               waveform_buf(1 : waveformbuf_index_max, 1 : nwinch)
  integer                      :: grid_stationwinch(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                               nsta_count(1 : ngrid_x, 1 : ngrid_y), timeindex_diff(1 : nsta_grid_max), &
  &                               waveform_tmp(1 : maxval(sampling_int))
  logical                      :: calc_grad, &
  &                               grid_enough_sta(1 : ngrid_x, 1 : ngrid_y), is_usewinch(1 : nwinch)
  integer,         allocatable :: station_winch(:)
  character(len = 129)         :: outfile
  character(len = 14)          :: time_char
  character(len = 2)           :: yr(1 : nsec_buf), mo(1 : nsec_buf), dy(1 : nsec_buf), &
  &                               hh(1 : nsec_buf), mm(1 : nsec_buf), ss(1 : nsec_buf)
#ifdef PARTICLEVELOCITY
  real(kind = fp)              :: particlevelocity(1 : 2, 1 : ngrid_x, 1 : ngrid_y)
  integer                      :: ncount1
  character(len = 129)         :: outfile_particle
#endif
  character(len = 255)         :: chtbl, chtbl_line

  integer                      :: filter_m(1 : nsampling_int), filter_n(1 : nsampling_int)
  real(kind = fp)              :: c(1 : nsampling_int), gn(1 : nsampling_int), station_sensitivity(1 : nwinch)
  real(kind = fp), allocatable :: h(:, :), uv(:, :)             !!For time-domain recursive filter 
  real(kind = fp)              :: ad_v_min, sensor_sens, naturalfreq, damp, stlat_tmp, stlon_tmp, stelev_tmp, &
  &                               ptime_cor, stime_cor
  integer                      :: winch_tmp, transdelay, mon_order, recflag, sensor_amp, ad_bit
  character(len = 4)           :: winch_char, comp_tmp, sensor_unit
  character(len = 8)           :: stname_tmp, stname(1 : nwinch)


  !nsta = command_argument_count()

  !allocate(waveform_obs(1 : ntime, 1 : nsta))
  !allocate(location_sta(1 : nsta), sacfile(1 : nsta), begin(1 : nsta))

  !do i = 1, nsta
  !  call get_command_argument(i, value = sacfile(i))
  !enddo
  !!read sac-formatted waveforms
  !do i = 1, nsta
  !  call read_sachdr(sacfile(i), begin = begin(i), delta = dt, npts = npts_tmp, &
  !  &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
  !  call read_sacdata(sacfile(i), ntime, waveform_obs(:, i))
  !  waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
#ifdef GREEN_CORRECTION
  !  if(location_sta(i)%depth .ne. -12.345_fp) &
  !  &  waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * (location_sta(i)%depth / depth_ref) ** 0.25_fp
#endif
  !enddo

  call get_command_argument(1, value = chtbl)

  !!read win-formatted channel table
  is_usewinch(1 : nwinch) = .false.
  open(unit = 10, file = trim(chtbl))
  nstation = 0
  do
    read(10, '(a255)', iostat = ios) chtbl_line
    if(ios .ne. 0) exit
    if(chtbl_line(1 : 1) .eq. "#") cycle
    nstation = nstation + 1
  enddo
  rewind(10)
  allocate(station_winch(1 : nstation))
  i = 1
  do
    read(10, '(a256)', iostat = ios) chtbl_line
    if(ios .ne. 0) exit
    if(chtbl_line(1 : 1) .eq. "#") cycle
    do i = 1, 255
      if(chtbl_line(j : j) .eq. "/") chtbl_line(j : j) = "-"
    enddo
    read(chtbl_line, *) winch_char, recflag, transdelay, stname_tmp, comp_tmp, mon_order, &
    &                   ad_bit, sensor_sens, sensor_unit, naturalfreq, damp, sensor_amp, ad_v_min, &
    &                   stlat_tmp, stlon_tmp, stelev_tmp, ptime_cor, stime_cor
    read(winch_char, '(z4)') station_winch(i)
    is_usewinch(station_winch(i)) = .true.
    stname(station_winch(i)) = trim(stname_tmp)
    station_sensitivity(station_winch(i)) = ad_v_min / (sensor_sens * 10.0_fp ** (sensor_amp / 20))
    location_sta(station_winch(i))%lon = stlon_tmp
    location_sta(station_winch(i))%lat = stlat_tmp
    i = i + 1
  enddo
  close(10)

  !!calculate filter parameter
  do i = 1, nsampling_int
    call calc_bpf_order(fl, fh, fs, ap, as, sampling_sec(i), filter_m(i), filter_n(i), c(i))
  enddo
  allocate(h(1 : 4 * maxval(filter_m), 1 : nsampling_int))
  do i = 1, nsampling_int
    call calc_bpf_coef(fl, fh, sampling_sec(i), filter_m(i), filter_n(i), h(:, i), c(i), gn(i))
  enddo

  !!set grid location
  do j = 1, ngrid_y
    do i = 1, ngrid_x
      location_grid(i, j)%x_east  = x_start + dgrid_x * real(i - 1, kind = fp)
      location_grid(i, j)%y_north = y_start + dgrid_y * real(j - 1, kind = fp)
      call xy2bl(location_grid(i, j)%y_north * 1000.0_fp, location_grid(i, j)%x_east * 1000.0_fp, &
      &          center_lon, center_lat, location_grid(i, j)%lon, location_grid(i, j)%lat)
    enddo
  enddo
  !!convert station longitude/latitude to x_east/y_north
  open(unit = 12, file = "station_location.txt")
  do i = 1, nstation
    call bl2xy(location_sta(station_winch(i))%lon,     location_sta(station_winch(i))%lat,    &
    &          center_lon,                             center_lat,                            &
    &          location_sta(station_winch(i))%y_north, location_sta(station_winch(i))%x_east)
    location_sta(station_winch(i))%y_north = location_sta(station_winch(i))%y_north / 1000.0_fp
    location_sta(station_winch(i))%x_east  = location_sta(station_winch(i))%x_east  / 1000.0_fp
    write(12, '(5(e15.7, 1x))') location_sta(station_winch(i))%x_east, location_sta(station_winch(i))%y_north, &
    &                           location_sta(station_winch(i))%lon,    location_sta(station_winch(i))%lat, &
    &                           location_sta(station_winch(i))%depth
  enddo
  close(12)

  !!make kernel matrix for each grid
  !call calc_kernelmatrix_circle(location_grid, location_sta, grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)
  call calc_kernelmatrix_delaunay2_shmdump(location_grid, station_winch, location_sta, naddstation_array, grid_enough_sta, &
  &                                        nsta_count, grid_stationwinch, kernel_matrix, error_matrix = error_matrix)

#ifdef PARTICLEVELOCITY
  particlevelocity(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
#endif

  !!calculate amplitude and its spatial derivatives at each grid
  !open(unit = 30, file = "log")
  !do kk = 1, int(ntime / ntimestep)
  waveform_buf(1 : waveformbuf_index_max, 1 : nwinch) = 0.0_fp
  time_loop: do

    calc_sampling_loop: do ii = 1, sampling_int_calc
      yr(1 : nsec_buf - 1) = yr(2 : nsec_buf)
      mo(1 : nsec_buf - 1) = mo(2 : nsec_buf)
      dy(1 : nsec_buf - 1) = dy(2 : nsec_buf)
      hh(1 : nsec_buf - 1) = hh(2 : nsec_buf)
      mm(1 : nsec_buf - 1) = mm(2 : nsec_buf)
      ss(1 : nsec_buf - 1) = ss(2 : nsec_buf)

      read(*, *, iostat = ios) yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf), nch
      if(ios .ne. 0) error stop
      write(0, '(a, 6(a2, 1x), a, i0)') "Reading ", &
      &                                 yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf), &
      &                                 " nch = ", nch
      if(.not. allocated(uv)) then
        allocate(uv(1 : 4 * maxval(filter_m), 1 : nwinch))
        uv(1 : 4 * maxval(filter_m), 1 : nwinch) = 0.0_fp
      endif

      waveform_buf(1 : waveformbuf_index_max - sampling_int_use, 1 : nwinch) &
      & = waveform_buf(sampling_int_use + 1 : waveformbuf_index_max, 1 : nwinch)
      !!read waveforms from stdin, filtering, decimation
      do j = 1, nch
        read(*, *) winch_char, nsample, (waveform_tmp(i), i = 1, nsample)
        read(winch_char, '(z4)') winch_tmp
        if(.not. is_usewinch(winch_tmp)) cycle
        do i = 1, nsampling_int
          if(nsample .eq. sampling_int(i)) exit
        enddo
        waveform_real(1 : nsample) = real(waveform_tmp(1 : nsample), kind = fp) * station_sensitivity(winch_tmp)
        if(i .le. nsampling_int) then
          call tandem3(waveform_real(1 : nsample), h(:, i), gn(i), 1, past_uv = uv(:, winch_tmp))
        endif
        ndecimate = sampling_int(i) / sampling_int_use
        do i = 1, nsample / ndecimate
          waveform_buf(waveformbuf_index_max - sampling_int_use + i, winch_tmp) = waveform_real(ndecimate * (i - 1) + 1)
        enddo
      enddo
    enddo calc_sampling_loop

  
    !timeindex = ntimestep * (kk - 1) + 1
    !write(0, '(a, i0, a)') "Time index = ", kk, " Calculate amplitudes and their gradients at each grid"
    !if(timeindex .lt. 1 .or. timeindex - ntimestep .gt. ntime) cycle

    !call int_to_char(kk, 4, ctimeindex)
    write(time_char, '(a, 6(i2.2))') "20", yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf)
    outfile = "slowness_gradiometry_" // trim(time_char) // ".dat"
    open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 10, status = "replace")
    ncount = 1
#ifdef PARTICLEVELOCITY
    outfile_particle = "particlevelocity_gradiometry_" // trim(time_char) // ".dat"
    open(unit   = 11, file = trim(outfile_particle), form = "unformatted", &
    &    access = "direct", recl = 4 * 4, status = "replace")
    ncount1 = 1
#endif

    !!initial step: estimate slowness vector without reducing velocity
    !!First, estimate spatial gradients
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        slowness(1 : 2, ii, jj)            = 0.0_fp
        slowness_correction(1 : 2, ii, jj) = 0.0_fp
        sigma_slowness(1 : 2, ii, jj)      = 0.0_fp
        ampterm(1 : 2, ii, jj)             = 0.0_fp
        sigma_ampterm(1 : 2, ii, jj)       = 0.0_fp
        waveform_est_plot(ii, jj)          = 0.0_fp
        !!calculate slowness vector using conventional array analysis
        if(grid_enough_sta(ii, jj) .eqv. .false.) cycle

        !!calculate spatial gradients of wavefield iteratively
        n = 0
        do
          calc_grad = .true.
          n = n + 1
          if(n .gt. niteration_max) exit
          !print *, "iterative reducing, count = ", n
          waveform_est_tmp(1 : 4, 1 : ngradient2) = 0.0_fp
          waveform_est_tmp2(1 : ngradient2, 1 : 4) = 0.0_fp

          do i = 1, nsta_count(ii, jj)
            dx_east  = location_grid(ii, jj)%x_east  - location_sta(grid_stationwinch(i, ii, jj))%x_east
            dy_north = location_grid(ii, jj)%y_north - location_sta(grid_stationwinch(i, ii, jj))%y_north
            timeindex_diff(i) = int((dx_east * slowness(1, ii, jj) + dy_north * slowness(2, ii, jj)) / dt)
          enddo
          timeindex_diff_min = minval(timeindex_diff)

          ngrad = 0
          do j = 1, ngradient2
            obs_vector(1 : nsta_grid_max) = 0.0_fp
            do i = 1, nsta_count(ii, jj)
              if(waveformbuf_index_max - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min .lt. 1 .or. &
                 waveformbuf_index_max - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min .gt. waveformbuf_index_max) then
                calc_grad = .false.
                exit
              endif
              obs_vector(i) &
              &  = waveform_buf(waveformbuf_index_max - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min, &
              &                 grid_stationwinch(i, ii, jj))
            enddo
            if(calc_grad .eqv. .false.) exit

            waveform_est_tmp(1 : 3, j) &
            &  = matmul(kernel_matrix(1 : 3, 1 : nsta_count(ii, jj), ii, jj), obs_vector(1 : nsta_count(ii, jj)))
            ngrad = ngrad + 1
          enddo

#ifdef PARTICLEVELOCITY
          if(n .eq. 1) then
            !!unit: [hPa(cm)/s]
            do j = ngradient2 - ntimestep + 1, ngradient2
              particlevelocity(1 : 2, ii, jj) = particlevelocity(1 : 2, ii, jj) &
              &                               - waveform_est_tmp(2 : 3, j) * 1.0e-5_fp & !!hPa(cm)/km -> hPa(cm)/cm
              &                               * grav_acc * dt
            enddo
            write(11, rec = ncount1) real(x_start + dgrid_x * real(ii - 1, kind = fp), kind = sp), &
            &                        real(y_start + dgrid_y * real(jj - 1, kind = fp), kind = sp), &
            &                        real(particlevelocity(1, ii, jj), kind = sp), &
            &                        real(particlevelocity(2, ii, jj), kind = sp)
            ncount1 = ncount1 + 1
          endif
#endif
          !!calculate slowness and app. geom. spreading terms at each grid
          if(ngrad .le. 2) cycle  !!if the number of data is small, do not calculate gradiometry coefficients

          !!Time derivative
          do j = 1, ngrad
            if(j .gt. 1 .and. j .lt. ngrad) then
              waveform_est_tmp(4, j) = (waveform_est_tmp(1, j + 1) - waveform_est_tmp(1, j - 1)) / dt * 0.5_fp
            elseif(j .eq. 1) then
              waveform_est_tmp(4, j) = (waveform_est_tmp(1, j + 1) - waveform_est_tmp(1, j))     / dt
            elseif(j .eq. ngrad) then
              waveform_est_tmp(4, j) = (waveform_est_tmp(1, j)     - waveform_est_tmp(1, j - 1)) / dt
            endif
            do i = 1, 4
              waveform_est_tmp2(j, i) = waveform_est_tmp(i, j)
            enddo
          enddo
          if(calc_grad .eqv. .false.) exit
          !if(n .eq. 1) waveform_est_plot(ii, jj) = waveform_est_tmp(1, ngradient2)
          waveform_est_plot(ii, jj) = waveform_est_tmp(1, ngrad)


          !!estimate slowness term
          slowness_cor_prev(1 : 2) = slowness_correction(1 : 2, ii, jj)
          ampterm(1 : 2, ii, jj) = 0.0_fp
          uu   = dot_product(waveform_est_tmp2(1 : ngrad, 1), waveform_est_tmp2(1 : ngrad, 1))
          uut  = dot_product(waveform_est_tmp2(1 : ngrad, 1), waveform_est_tmp2(1 : ngrad, 4))
          utut = dot_product(waveform_est_tmp2(1 : ngrad, 4), waveform_est_tmp2(1 : ngrad, 4))
          denominator = uu * utut - uut * uut
          do i = 1, 2
            uxu(i)  = dot_product(waveform_est_tmp2(1 : ngrad, i + 1), waveform_est_tmp2(1 : ngrad, 1))
            uxut(i) = dot_product(waveform_est_tmp2(1 : ngrad, i + 1), waveform_est_tmp2(1 : ngrad, 4))
            numerator_slowness(i) = uu   * uxut(i) - uut * uxu(i)
            numerator_ampterm(i)  = utut * uxu(i)  - uut * uxut(i)

            slowness_correction(i, ii, jj) = -numerator_slowness(i) / denominator
            ampterm(i, ii, jj)             =  numerator_ampterm(i)  / denominator
          enddo
          slowness(1 : 2, ii, jj) = slowness(1 : 2, ii, jj) + slowness_correction(1 : 2, ii, jj)


          !if(slowness_correction(1, ii, jj) * slowness_correction(1, ii, jj) &
          !&  + slowness_correction(2, ii, jj) * slowness_correction(2, ii, jj) .lt. eps) exit
          if(  (slowness_correction(1, ii, jj) - slowness_cor_prev(1)) ** 2 &
          &  + (slowness_correction(2, ii, jj) - slowness_cor_prev(2)) ** 2 .lt. eps) exit
        enddo

        if(calc_grad .eqv. .false.) cycle

#ifdef PARTICLEVELOCITY
        do j = 1, ngrad
          particlevelocity(1 : 2, ii, jj) = particlevelocity(1 : 2, ii, jj) &
          &                                  + waveform_est_tmp(2 : 3, j) * 1e-2_fp / 1e+3_fp &   !!cm/km
          &                                  * grav_acc * 60.0_fp / real(ntimestep, kind = fp)
        enddo
        write(11, rec = ncount1) real(x_start + dgrid_x * real(ii - 1, kind = fp), kind = sp), &
        &                        real(y_start + dgrid_y * real(jj - 1, kind = fp), kind = sp), &
        &                        real(particlevelocity(1, ii, jj), kind = sp), &
        &                        real(particlevelocity(2, ii, jj), kind = sp)
        ncount1 = ncount1 + 1
#endif
        !!error estimation
        max_innerproduct = 0.0_fp
        do i = 1, nsta_count(ii, jj)
          dx_east  = location_grid(ii, jj)%x_east  - location_sta(grid_stationwinch(i, ii, jj))%x_east
          dy_north = location_grid(ii, jj)%y_north - location_sta(grid_stationwinch(i, ii, jj))%y_north
          innerproduct_tmp = (slowness(1, ii, jj) * dx_east + slowness(2, ii, jj) * dy_north) ** 2
          if(innerproduct_tmp .ge. max_innerproduct) then
            !print *, dx_east, dy_north, innerproduct_tmp
            max_innerproduct = innerproduct_tmp
          endif
        enddo
        relativeerror = 0.5_fp * error_omega ** 2 * max_innerproduct
        !print *, ii, jj, error_omega, max_innerproduct, relativeerror

        do j = 1, ngrad
          do i = 1, 2
            sigma_slowness(i, ii, jj) = sigma_slowness(i, ii, jj) &
                                      !!dB/du
            &                         + ((2.0_fp * waveform_est_tmp(1, j)     * uxut(i) &
            &                                    - waveform_est_tmp(4, j)     * uxu (i) &
            &                                    - waveform_est_tmp(i + 1, j) * uut) / denominator &
            &                           - 2.0_fp * numerator_slowness(i) / (denominator ** 2) &
            &                                    * (waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut)) ** 2 &
            &                         * (waveform_est_tmp(1, j) * relativeerror) ** 2 * error_matrix(1, ii, jj) &
                                      !!dB/dut
            &                         + ((waveform_est_tmp(i + 1, j) * uu - waveform_est_tmp(1, j) * uxu(i)) / denominator &
            &                           - 2.0_fp * (numerator_slowness(i) / denominator ** 2) &
            &                                    * (waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut)) ** 2 &
            &                         * (waveform_est_tmp(4, j) * 2.0_fp / (dt ** 2) * relativeerror ) ** 2 &
            &                         * error_matrix(1, ii, jj) &
                                      !!dB/dux
            &                         + ((waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut) / denominator) ** 2 &
            &                         * (waveform_est_tmp(i + 1, j) * relativeerror) ** 2 * error_matrix(i + 1, ii, jj)
            sigma_ampterm(i, ii, jj) = sigma_ampterm(i, ii, jj) &
                                     !!dA/du
            &                        + ((waveform_est_tmp(i + 1, j) * utut - waveform_est_tmp(4, j) * uxut(i)) / denominator &
            &                           - 2.0_fp * (numerator_ampterm(i) / denominator ** 2) &
            &                                    * (waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut)) ** 2 &
            &                        * (waveform_est_tmp(1, j) * relativeerror) ** 2 * error_matrix(1, ii, jj) &
                                     !!dA/dut
            &                        + ((2.0_fp * waveform_est_tmp(4, j)     * uxu (i) &
            &                                   - waveform_est_tmp(1, j)     * uxut(i) &
            &                                   - waveform_est_tmp(i + 1, j) * uut) / denominator &
            &                          - 2.0_fp * numerator_ampterm(i) / (denominator ** 2) &
            &                                   * (waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut)) ** 2 &
            &                        * (waveform_est_tmp(4, j) * 2.0_fp / (dt ** 2) * relativeerror)  ** 2 &
            &                        * error_matrix(1, ii, jj) &
                                     !!dA/dux
            &                        + ((waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut) / denominator) ** 2 &
            &                        * (waveform_est_tmp(i + 1, j) * relativeerror) ** 2 * error_matrix(i + 1, ii, jj)
          enddo
        enddo

        !print *, "grid index = ", ii, jj, (n - 1)
        !print *, "gradiometry slowness nocor", slowness_correction(1, ii, jj), slowness_correction(2, ii, jj)
        !print *, "gradiometry slowness nocor", slowness(1, ii, jj), slowness(2, ii, jj)
        !app_velocity = 1.0_fp / sqrt(slowness(1, ii, jj) ** 2 + slowness(2, ii, jj) ** 2)
        !backazimuth = atan2(slowness(2, ii, jj), slowness(1, ii, jj)) * rad2deg
        !print *, "gradiometry slowness cor", app_velocity, backazimuth
        !app_velocity = ampterm(1, ii, jj) * sin(backazimuth * deg2rad) + ampterm(2, ii, jj) * cos(backazimuth * deg2rad)
        !print *, "gradiometry ampterm", ampterm(1, ii, jj), ampterm(2, ii, jj), app_velocity

        write(10, rec = ncount) real(x_start + dgrid_x * real(ii - 1, kind = fp), kind = sp), &
        &                       real(y_start + dgrid_y * real(jj - 1, kind = fp), kind = sp), &
        &                       real(slowness(1, ii, jj), kind = sp),       real(slowness(2, ii, jj), kind = sp), & 
        &                       real(sigma_slowness(1, ii, jj), kind = sp), real(sigma_slowness(2, ii, jj), kind = sp), &
        &                       real(ampterm(1, ii, jj), kind = sp),        real(ampterm(2, ii, jj), kind = sp), &
        &                       real(sigma_ampterm(1, ii, jj), kind = sp),  real(sigma_ampterm(2, ii, jj), kind = sp)
        ncount = ncount + 1

      enddo
    enddo
    close(10)
#ifdef PARTICLEVELOCITY
    close(11)
#endif

    outfile = "amplitude_gradiometry_" // trim(time_char) // ".grd"
    outfile = trim(outfile)
    call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est_plot, outfile, &
    &                        nanval = 0.0_fp)

  enddo time_loop
  close(30)


  stop
end program seismicgradiometry_reducingvelocity2_shmdump

