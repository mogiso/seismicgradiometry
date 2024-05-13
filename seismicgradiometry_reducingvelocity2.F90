!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

program seismicgradiometry_reducingvelocity2
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

  integer :: nsta, npts_tmp, timeindex, timeindex_diff_min, i, j, ii, jj, kk, ncount, m, n, ngrad
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)
  type(location),  allocatable :: location_sta(:)
  real(kind = fp), allocatable :: begin(:), waveform_obs(:, :)
  real(kind = fp), allocatable :: h(:), uv(:, :)             !!For time-domain recursive filter 
  real(kind = fp)              :: dt, c, gn, dx_east, dy_north, denominator, uu, uut, utut, min_dx, min_dy, relativeerror, &
  &                               uxu(1 : 2), uxut(1 : 2), numerator_slowness(1 : 2), numerator_ampterm(1 : 2), &
  &                               obs_vector(1 : nsta_grid_max), &
  &                               sigma_slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               sigma_ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), & 
  &                               waveform_est_tmp(1 : 4, 1 : ngradient2), &
  &                               waveform_est_tmp2(1 : ngradient2, 1 : 4), &
  &                               waveform_est_plot(1 : ngrid_x, 1 : ngrid_y), &
  &                               slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               slowness_correction(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                               kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                               error_matrix(1 : 3, 1 : ngrid_x, 1 : ngrid_y)
  integer                      :: grid_stationindex(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                               nsta_count(1 : ngrid_x, 1 : ngrid_y), timeindex_diff(1 : nsta_grid_max)
  logical                      :: calc_grad, &
  &                               grid_enough_sta(1 : ngrid_x, 1 : ngrid_y)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  character(len = 4) :: ctimeindex
#ifdef PARTICLEVELOCITY
  real(kind = fp)              :: particlevelocity(1 : 2, 1 : ngrid_x, 1 : ngrid_y)
  integer                      :: ncount1
  character(len = 129)         :: outfile_particle
#endif

  nsta = command_argument_count()

  allocate(waveform_obs(1 : ntime, 1 : nsta))
  allocate(location_sta(1 : nsta), sacfile(1 : nsta), begin(1 : nsta))

  do i = 1, nsta
    call get_command_argument(i, value = sacfile(i))
  enddo

  !!read sac-formatted waveforms
  do i = 1, nsta
    call read_sachdr(sacfile(i), begin = begin(i), delta = dt, npts = npts_tmp, &
    &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
    call read_sacdata(sacfile(i), ntime, waveform_obs(:, i))
#ifdef ENVELOPE
    !waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order * order
#else
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
#endif
#ifdef GREEN_CORRECTION
    if(location_sta(i)%depth .ne. -12.345_fp) &
    &  waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * (location_sta(i)%depth / depth_ref) ** 0.25_fp
#endif
  enddo

#ifdef ENVELOPE
#else
  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(1 : 4 * m), uv(1 : 4 * m, 1 : nsta))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  do i = 1, nsta
    call tandem3(waveform_obs(:, i), h, gn, 1, past_uv = uv(:, i))
  enddo
  !!inverse filtering
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  do i = 1, nsta
    call tandem3(waveform_obs(:, i), h, gn, -1, past_uv = uv(:, i))
  enddo
  deallocate(h, uv)
#endif

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
  &                                nsta_count, grid_stationindex, kernel_matrix, error_matrix = error_matrix)

#ifdef PARTICLEVELOCITY
  particlevelocity(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
#endif

  !!calculate amplitude and its spatial derivatives at each grid
  !open(unit = 30, file = "log")
  do kk = 1, int(ntime / ntimestep)
    timeindex = ntimestep * (kk - 1) + 1
    write(0, '(a, i0, a)') "Time index = ", kk, " Calculate amplitudes and their gradients at each grid"
    sigma_slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
    sigma_ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp

    if(timeindex .lt. 1 .or. timeindex - ntimestep .gt. ntime) cycle

    call int_to_char(kk, 4, ctimeindex)
    outfile = "slowness_gradiometry_" // trim(ctimeindex) // ".dat"
    open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 10, status = "replace")
    ncount = 1
#ifdef PARTICLEVELOCITY
    outfile_particle = "particlevelocity_gradiometry_" // trim(ctimeindex) // ".dat"
    open(unit   = 11, file = trim(outfile_particle), form = "unformatted", &
    &    access = "direct", recl = 4 * 4, status = "replace")
    ncount1 = 1
#endif

    !!initial step: estimate slowness vector without reducing velocity
    !!First, estimate spatial gradients
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        slowness(1 : 2, ii, jj) = 0.0_fp
        slowness_correction(1 : 2, ii, jj) = 0.0_fp
        ampterm(1 : 2, ii, jj) = 0.0_fp
        waveform_est_plot(ii, jj) = 0.0_fp
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
            dx_east  = location_grid(ii, jj)%x_east  - location_sta(grid_stationindex(i, ii, jj))%x_east
            dy_north = location_grid(ii, jj)%y_north - location_sta(grid_stationindex(i, ii, jj))%y_north
            timeindex_diff(i) = int((dx_east * slowness(1, ii, jj) + dy_north * slowness(2, ii, jj)) / dt)
            if(i .eq. 1) then
              min_dx = abs(dx_east)
              min_dy = abs(dy_north)
            elseif(i .gt. 1) then
              if(abs(dx_east) .le. min_dx) min_dx = abs(dx_east)
              if(abs(dy_north) .le. min_dy) min_dy = abs(dy_north)
            endif
          enddo
          timeindex_diff_min = minval(timeindex_diff)

          ngrad = 0
          do j = 1, ngradient2
            obs_vector(1 : nsta_grid_max) = 0.0_fp
            do i = 1, nsta_count(ii, jj)
              if(timeindex - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min .lt. 1 .or. &
                 timeindex - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min .gt. ntime) then
                calc_grad = .false.
                exit
              endif
              obs_vector(i) &
              &  = waveform_obs(timeindex - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min, &
              &                 grid_stationindex(i, ii, jj))
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
          slowness_correction(1 : 2, ii, jj) = 0.0_fp
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
            if(denominator .lt. eps) exit

            slowness_correction(i, ii, jj) = -numerator_slowness(i) / denominator
            ampterm(i, ii, jj) = numerator_ampterm(i) / denominator
          enddo
          slowness(1 : 2, ii, jj) = slowness(1 : 2, ii, jj) + slowness_correction(1 : 2, ii, jj)
          !!error estimation
          relativeerror = (pi * fh) ** 2 &
          &             * ((slowness(1, ii, jj) * min_dx) ** 2 + (slowness(2, ii, jj) * min_dy) ** 2)
          do j = 1, ngrad
            do i = 1, 2
              sigma_slowness(i, ii, jj) = sigma_slowness(i, ii, jj) &
                                        !!dB/du
              &                         + ((2.0_fp * waveform_est_tmp(1, j)     * uxut(i) &
              &                                    - waveform_est_tmp(4, j)     * uxu (i) &
              &                                    - waveform_est_tmp(i + 1, j) * uut) / denominator &
              &                           - 2.0_fp * (numerator_slowness(i) / denominator ** 2) &
              &                                    * (waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut)) ** 2 &
              &                         * (waveform_est_tmp(1, j) * relativeerror * error_matrix(1, ii, jj)) ** 2 &
                                        !!dB/dut
              &                         + ((waveform_est_tmp(i + 1, j) * uu - waveform_est_tmp(1, j) * uxu(i)) &
              &                           / denominator &
              &                           - 2.0_fp * (numerator_slowness(i) / denominator ** 2) &
              &                                    * (waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut)) ** 2 &
              &                         * (waveform_est_tmp(4, j) &
              &                         * (2.0_fp * relativeerror * error_matrix(1, ii, jj) / dt)) ** 2 &
                                        !!dB/dux
              &                         + ((waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut) &
              &                           / denominator) ** 2 &
              &                         * (waveform_est_tmp(i + 1, j) * (relativeerror * error_matrix(i + 1, ii, jj))) ** 2
              sigma_ampterm(i, ii, jj) = sigma_ampterm(i, ii, jj) &
                                       !!dA/du
              &                        + ((waveform_est_tmp(i + 1, j) * utut - waveform_est_tmp(4, j) * uxut(i)) &
              &                           / denominator &
              &                           - 2.0_fp * (numerator_ampterm(i) / denominator ** 2) &
              &                                    * (waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut)) ** 2 &
              &                        * (waveform_est_tmp(1, j) * relativeerror * error_matrix(1, ii, jj)) ** 2 &
                                       !!dA/dut
              &                        + ((2.0_fp * waveform_est_tmp(4, j)     * uxu (i) &
              &                                   - waveform_est_tmp(1, j)     * uxut(i) &
              &                                   - waveform_est_tmp(i + 1, j) * uut) / denominator &
              &                          - 2.0_fp * (numerator_ampterm(i) / denominator ** 2) &
              &                                   * (waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut)) ** 2 &
              &                        * (waveform_est_tmp(4, j) &
              &                        * (2.0_fp * relativeerror * error_matrix(1, ii, jj) / dt)) ** 2 &
                                       !!dA/dux
              &                        + ((waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut) &
              &                          / denominator) ** 2 &
              &                        * (waveform_est_tmp(i + 1, j) * (relativeerror * error_matrix(i + 1, ii, jj))) ** 2
            enddo
          enddo

          if(slowness_correction(1, ii, jj) * slowness_correction(1, ii, jj) &
          &  + slowness_correction(2, ii, jj) * slowness_correction(2, ii, jj) .lt. eps) exit
        enddo
        if(calc_grad .eqv. .false.) cycle


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

    outfile = "amplitude_gradiometry_" // trim(ctimeindex) // ".grd"
    outfile = trim(outfile)
    call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est_plot, outfile, &
    &                        nanval = 0.0_fp)

  enddo
  close(30)


  stop
end program seismicgradiometry_reducingvelocity2

