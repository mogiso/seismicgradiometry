!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

program seismicgradiometry_reducingvelocity
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use typedef
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy, xy2bl
  use gradiometry_parameters
  use calc_kernelmatrix
  use correlation
  use itoa
  use taper

  implicit none

  real(kind = fp), parameter :: eps = 1.0e-5_fp

  integer :: nsta, npts_tmp, timeindex, timeindex_diff, i, j, k, ii, jj, kk, icount, m, n, ngrad_mean
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)
  type(location),  allocatable :: location_sta(:)
  real(kind = fp), allocatable :: begin(:)
  real(kind = fp), allocatable :: waveform_obs(:, :),       & !!(1 : ntime, 1 : nsta)
  &                               waveform_est(:, :, :, :), & !!(4(u, dudx, dudy, dudt), ngrid_x, ngrid_y, ntime/ntimestep)
  &                               slowness_est_matrix(:, :, :, :), &
  &                               lagtime(:), &
  &                               maxval_xcorr(:)
  real(kind = fp), allocatable :: h(:), uv(:, :)             !!For time-domain recursive filter 
  real(kind = fp)              :: dt, c, gn, dx_east, dy_north, app_velocity, backazimuth, &
  &                               obs_vector(1 : nsta_grid_max), &
  &                               xcorr(1 : ntime_fft), &
  &                               minval_xcorr(1 : ngrid_x, 1 : ngrid_y), &
  &                               sigma_slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               sigma_ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), & 
  &                               taper_window(1 : ntime_fft), &
  &                               waveform_fft(1 : ntime_fft, 1 : 2), &
  &                               waveform_est_tmp(1 : 3, -ntime_fft / 4 : ntime_fft / 4), &
  &                               waveform_est_tmp2(1 : ntime_slowness, 1 : 4), &
  &                               waveform_est_plot(1 : ngrid_x, 1 : ngrid_y), &
  &                               slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               slowness_xcorr(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y)
  integer                      :: grid_stationindex(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                               nsta_count(1 : ngrid_x, 1 : ngrid_y), &
  &                               max_xcorr(1), &
  &                               nsta_correlation(1 : ngrid_x, 1 : ngrid_y)
  logical                      :: calc_grad_mean, &
  &                               grid_enough_sta(1 : ngrid_x, 1 : ngrid_y)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  character(len = 4) :: ctimeindex


  nsta = command_argument_count()

  allocate(waveform_obs(1 : ntime, 1 : nsta), waveform_est(1 : 4, 1 : ngrid_x, 1 : ngrid_y, 1 : int(ntime / ntimestep)))
  allocate(location_sta(1 : nsta), sacfile(1 : nsta), begin(1 : nsta))
  waveform_est(1 : 4, 1 : ngrid_x, 1 : ngrid_y, 1 : int(ntime / ntimestep)) = 0.0_fp

  do i = 1, nsta
    call get_command_argument(i, value = sacfile(i))
  enddo

  !!read sac-formatted waveforms
  do i = 1, nsta
    call read_sachdr(sacfile(i), begin = begin(i), delta = dt, npts = npts_tmp, &
    &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
    call read_sacdata(sacfile(i), ntime, waveform_obs(:, i))
    !waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
  enddo
  call cosine_taper(cos_taper_ratio, ntime_fft, taper_window)

  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(1 : 4 * m), uv(1 : 4 * m, 1 : nsta))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  !do i = 1, nsta
  !  call tandem2(waveform_obs(:, i), waveform_obs(:, i), ntime, h, m, 1, gn, uv)
  !enddo
  deallocate(h, uv)

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
  do i = 1, nsta
    call bl2xy(location_sta(i)%lon, location_sta(i)%lat, center_lon, center_lat, &
    &          location_sta(i)%y_north, location_sta(i)%x_east)
    location_sta(i)%y_north = location_sta(i)%y_north / 1000.0_fp
    location_sta(i)%x_east  = location_sta(i)%x_east / 1000.0_fp
    write(12, '(5(e15.7, 1x))') location_sta(i)%x_east, location_sta(i)%y_north, location_sta(i)%lon, location_sta(i)%lat, &
    &                           location_sta(i)%depth
  enddo
  close(12)

  !!make kernel matrix for each grid
  !call calc_kernelmatrix_circle(location_grid, location_sta, grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)

  call calc_kernelmatrix_delaunay2(location_grid, location_sta, naddstation_array, grid_enough_sta, &
  &                                nsta_count, grid_stationindex, kernel_matrix, &
  &                                nsta_correlation = nsta_correlation, slowness_est_matrix = slowness_est_matrix)

  !!calculate amplitude and its spatial derivatives at each grid
  !do kk = 1, int(ntime / ntimestep)
  do kk = 800, 1000
    timeindex = ntimestep * (kk - 1) + 1
    if(timeindex - ntime_fft + 1 .lt. 1) cycle
    if(kk - ntime_slowness + 1 .lt. 1) cycle
    write(0, '(a, i0, a)') "Time index = ", kk, " Calculate amplitudes and their gradients at each grid"
    waveform_est(1 : 4, 1 : ngrid_x, 1 : ngrid_y, kk) = 0.0_fp
    sigma_slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_sp
    sigma_ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_sp

    if(timeindex .lt. 1 .or. timeindex - ntimestep .gt. ntime) cycle

    call int_to_char(kk, 4, ctimeindex)
    outfile = "slowness_gradiometry_" // trim(ctimeindex) // ".dat"
    open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 10, status = "replace")
    icount = 1

    !!initial step: estimate slowness vector without reducing velocity
    !!First, estimate spatial gradients
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        slowness(1 : 2, ii, jj) = 0.0_fp
        ampterm(1 : 2, ii, jj) = 0.0_fp
        !!calculate slowness vector using conventional array analysis
        if(grid_enough_sta(ii, jj) .eqv. .false.) cycle
        minval_xcorr(ii, jj) = 1.0_fp
        !correlation
        allocate(lagtime(1 : nsta_correlation(ii, jj)), maxval_xcorr(1 : nsta_correlation(ii, jj)))
        k = 1
        do j = 1, nsta_count(ii, jj) - 1
          waveform_fft(1 : ntime_fft, 1) &
          &  = waveform_obs(timeindex - ntime_fft + 1 : timeindex, grid_stationindex(j, ii, jj)) * taper_window(1 : ntime_fft)
          do i = j + 1, nsta_count(ii, jj)
            waveform_fft(1 : ntime_fft, 2) &
            &  = waveform_obs(timeindex - ntime_fft + 1 : timeindex, grid_stationindex(i, ii, jj)) * taper_window(1 : ntime_fft)

            call correlation_fft(waveform_fft, ntime_fft, xcorr)
            max_xcorr = maxloc(xcorr)
            if(xcorr(max_xcorr(1)) .le. minval_xcorr(ii, jj)) minval_xcorr(ii, jj) = xcorr(max_xcorr(1))
            lagtime(k) = real(max_xcorr(1) - ntime_fft2, kind = fp) * dt
            !print '(3(i0, 1x), 2(1x, f8.5))', k, grid_stationindex(i, ii, jj), grid_stationindex(j, ii, jj), &
            !&                                 xcorr(max_xcorr(1)), lagtime(k)
            k = k + 1
          enddo
        enddo
        if(minval_xcorr(ii, jj) .le. xcorr_min) then
          deallocate(lagtime, maxval_xcorr)
          cycle
        endif

        slowness_xcorr(1 : 2, ii, jj) &
        &  = matmul(slowness_est_matrix(1 : 2, 1 : nsta_correlation(ii, jj), ii, jj), lagtime(1 : nsta_correlation(ii, jj)))
        app_velocity = 1.0_fp / sqrt(slowness_xcorr(1, ii, jj) ** 2 + slowness_xcorr(2, ii, jj) ** 2)
        backazimuth = atan2(slowness_xcorr(2, ii, jj), slowness_xcorr(1, ii, jj)) * rad2deg
        print *, ii, jj, app_velocity, backazimuth
        deallocate(lagtime, maxval_xcorr)

        !!calculate spatial gradients of wavefield
        !set up observation vector
        calc_grad_mean = .true.
        ngrad_mean = 0
        do j = -ntime_fft4, ntime_fft4, 1
          obs_vector(1 : nsta_grid_max) = 0.0_fp
          do i = 1, nsta_count(ii, jj)
            dx_east  = location_grid(ii, jj)%x_east  - location_sta(grid_stationindex(i, ii, jj))%x_east
            dy_north = location_grid(ii, jj)%y_north - location_sta(grid_stationindex(i, ii, jj))%y_north
            timeindex_diff = int((dx_east * slowness_xcorr(1, ii, jj) + dy_north * slowness_xcorr(2, ii, jj)) * dt + 0.5_fp)
            if(timeindex - ntime_fft2 + j + timeindex_diff .lt. 1 .or. &
            &  timeindex - ntime_fft2 + j + timeindex_diff .gt. timeindex) then
              calc_grad_mean = .false.
              exit
            endif
            obs_vector(i) = waveform_obs(timeindex - ntime_fft2 + j + timeindex_diff, grid_stationindex(i, ii, jj))
            !print *, kk, timeindex - ntime_fft2 + j + timeindex_diff, grid_stationindex(i, ii, jj), obs_vector(i)
          enddo

          if(calc_grad_mean .eqv. .false.) exit

          waveform_est_tmp(1 : 3, j) &
          & = matmul(kernel_matrix(1 : 3, 1 : nsta_count(ii, jj), ii, jj), obs_vector(1 : nsta_count(ii, jj)))
          ngrad_mean = ngrad_mean + 1
        enddo

        waveform_est(1, ii, jj, kk) = waveform_est_tmp(1, 0)
        waveform_est_plot(ii, jj) = waveform_est_tmp(1, 0)

        do i = 1, ngrad_mean
          waveform_est(2 : 3, ii, jj, kk) = waveform_est(2 : 3, ii, jj, kk) + waveform_est_tmp(2 : 3, -ntime_fft4 + (i - 1))
        enddo
        waveform_est(2 : 3, ii, jj, kk) = waveform_est(2 : 3, ii, jj, kk) / real(ngrad_mean, kind = fp)
        !!Time derivative: backward differentiation
        if(kk - 1 .gt. 1) then
          waveform_est(4, ii, jj, kk) &
          &  = (waveform_est(1, ii, jj, kk) - waveform_est(1, ii, jj, kk - 1)) / (dt * real(ntimestep, kind = fp))
        endif

        !!calculate slowness and app. geom. spreading terms at each grid
        do j = 1, 4
          do i = 1, ntime_slowness
            waveform_est_tmp2(ntime_slowness - (i - 1), j) = waveform_est(j, ii, jj, kk - (i - 1))
          enddo
        enddo

        !!estimate slowness term
        do i = 1, 2
          if( (dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 1))   &
          & *  dot_product(waveform_est_tmp2(:, 4),     waveform_est_tmp2(:, 4)))  &
          & - (dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))   &
          & *  dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))) .lt. eps) exit
          slowness(i, ii, jj) = &
          &  + ((dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 1))   &
          &  *   dot_product(waveform_est_tmp2(:, i + 1), waveform_est_tmp2(:, 4)))  &
          &  -  (dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))   &
          &  *   dot_product(waveform_est_tmp2(:, 2),     waveform_est_tmp2(:, 1)))) &
          &  / ((dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 1))   &
          &  *   dot_product(waveform_est_tmp2(:, 4),     waveform_est_tmp2(:, 4)))  &
          &  -  (dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))   &
          &  *   dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))))
          ampterm(i, ii, jj) &
          &  = ((dot_product(waveform_est_tmp2(:, 4),     waveform_est_tmp2(:, 4))   &
          &  *   dot_product(waveform_est_tmp2(:, i + 1), waveform_est_tmp2(:, 1)))  &
          &  -  (dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))   &
          &  *   dot_product(waveform_est_tmp2(:, i + 1), waveform_est_tmp2(:, 4)))) &
          &  / ((dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 1))   &
          &  *   dot_product(waveform_est_tmp2(:, 4),     waveform_est_tmp2(:, 4)))  &
          &  -  (dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))   &
          &  *   dot_product(waveform_est_tmp2(:, 1),     waveform_est_tmp2(:, 4))))
        enddo
        print *, "gradiometry slowness nocor", slowness(1, ii, jj), slowness(2, ii, jj)
        slowness(1 : 2, ii, jj) = -slowness(1 : 2, ii, jj) + slowness_xcorr(1 : 2, ii, jj)
        app_velocity = 1.0_fp / sqrt(slowness(1, ii, jj) ** 2 + slowness(2, ii, jj) ** 2)
        backazimuth = atan2(slowness(2, ii, jj), slowness(1, ii, jj)) * rad2deg
        print *, "gradiometry slowness cor", app_velocity, backazimuth
        app_velocity = ampterm(1, ii, jj) * sin(backazimuth * deg2rad) + ampterm(2, ii, jj) * cos(backazimuth * deg2rad)
        print *, "gradiometry ampterm", ampterm(1, ii, jj), ampterm(2, ii, jj), app_velocity
        

        write(10, rec = icount) real(x_start + dgrid_x * real(ii - 1, kind = fp), kind = sp), &
        &                       real(y_start + dgrid_y * real(jj - 1, kind = fp), kind = sp), &
        &                       real(-slowness(1, ii, jj), kind = sp),      real(-slowness(2, ii, jj), kind = sp), &
        &                       real(sigma_slowness(1, ii, jj), kind = sp), real(sigma_slowness(2, ii, jj), kind = sp), & 
        &                       real(ampterm(1, ii, jj), kind = sp),        real(ampterm(2, ii, jj), kind = sp), &
        &                       real(sigma_ampterm(1, ii, jj), kind = sp),  real(sigma_ampterm(2, ii, jj), kind = sp)
        icount = icount + 1

      enddo
    enddo
    close(10)

    outfile = "amplitude_gradiometry_" // trim(ctimeindex) // ".grd"
    outfile = trim(outfile)
    call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est_plot, outfile)

  enddo


  stop
end program seismicgradiometry_reducingvelocity

