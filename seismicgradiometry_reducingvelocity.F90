!! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
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

  implicit none

  real(kind = fp), parameter :: eps = 1.0e-5_fp

  integer :: nsta, npts_tmp, info, i, j, ii, jj, kk, i3, j3, m, n, icount, jcount
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)
  type(location), allocatable :: location_sta(:)
  real(kind = fp), allocatable :: begin(:)
  real(kind = fp), allocatable :: waveform_obs(:, :)         !!(1 : ntime, 1 : nsta)
  real(kind = fp), allocatable :: waveform_est(:, :, :, :)   !!(1 : 3 (u, dudx, dudy), 1 : ngrid_x, 1 : ngrid_y, 1 : ntime)
  real(kind = fp), allocatable :: waveform_est_dt(:, :, :) !!(1 : ngrid_x, 1 : ngrid_y, 1 : ntime)
  real(kind = fp), allocatable :: h(:)
  real(kind = fp) :: weight(1 : nsta_grid_max, 1 : nsta_grid_max), &
  &                  kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), obs_vector(1 : nsta_grid_max), &
  &                  slowness_x(1 : ngrid_x, 1 : ngrid_y), slowness_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  sigma_slowness_x(1 : ngrid_x, 1 : ngrid_y), sigma_slowness_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  waveform_est_tmp(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_est_dt_tmp(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_est_dx(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_est_dy(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_errest_tmp(1 : ntime_slowness - 1), &
  &                  waveform_errest_dt_tmp(1 : ntime_slowness - 1), &
  &                  waveform_errest_dx(1 : ntime_slowness - 1), &
  &                  waveform_errest_dy(1 : ntime_slowness - 1), &
  &                  param_x_err_tmp(1 : 2, -ntime_slowness2 : ntime_slowness2), &
  &                  param_y_err_tmp(1 : 2, -ntime_slowness2 : ntime_slowness2), &
  &                  param_x_err(1 : 2), param_y_err(1 : 2), &
  &                  ampterm_x(1 : ngrid_x, 1 : ngrid_y), ampterm_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  sigma_ampterm_x(1 : ngrid_x, 1 : ngrid_y), sigma_ampterm_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  waveform_est_max(1 : ngrid_x, 1 : ngrid_y), waveform_est_dt_max(1 : ngrid_x, 1 : ngrid_y), &
  &                  dt, c, gn
  integer :: grid_stationindex(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), nsta_count(1 : ngrid_x, 1 : ngrid_y)
  logical :: grid_enough_sta(1 : ngrid_x, 1 : ngrid_y)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  character(len = 4) :: time_index


  nsta = command_argument_count()

  allocate(waveform_obs_buf(1 : ntime, 1 : nsta), &
  &        waveform_est_dt(1 : ngrid_x, 1 : ngrid_y, 1 : ntime_slowness), &
  &        waveform_est(1 : 3, 1 : ngrid_x, 1 : ngrid_y, 1 : ntime_slowness))
  allocate(location_sta(1 : nsta), sacfile(1 : nsta), begin(1 : nsta))

  do i = 1, nsta
    call get_command_argument(i, value = sacfile(i))
  enddo

  !!read sac-formatted waveforms
  do i = 1, nsta
    call read_sachdr(sacfile(i), begin = begin(i), delta = dt, npts = npts_tmp, &
    &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
    call read_sacdata(sacfile(i), ntime, waveform_obs(:, i))
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
  enddo

  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(1 : 4 * m), uv(1 : 4 * m, 1 : nsta))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp

  !!set grid location
  do j = 1, ngrid_y
    do i = 1, ngrid_x
      location_grid(i, j)%x_east = x_start + dgrid_x * real(i - 1, kind = fp)
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
    location_sta(i)%x_east = location_sta(i)%x_east / 1000.0_fp
    write(12, '(5(e15.7, 1x))') location_sta(i)%x_east, location_sta(i)%y_north, location_sta(i)%lon, location_sta(i)%lat, &
    &                           location_sta(i)%depth
  enddo
  close(12)

  !!make kernel matrix for each grid
  !call calc_kernelmatrix_circle(location_grid, location_sta, grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)
  call calc_kernelmatrix_delaunay(location_grid, location_sta, grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)


  !!calculate amplitude and its spatial derivatives at each grid
  do kk = 1, int(ntime / ntimestep)
    time_index = ntimestep * (kk - 1) + 1
    if(time_index - ntimestep + 1 .lt. 1) cycle

    !!initial step: estimate slowness vector without reducing velocity
    !!First, estimate spatial gradients
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        !!calculate slowness vector using conventional array analysis
        if(grid_enough_sta(ii, jj) .eqv. .false.) cycle
        !correlation
        nsta_correlation = 0
        do j = 1, nsta_count(ii, jj) - 1
          do i = j + 1, nsta_count(ii, jj)
            nsta_correlation = nsta_correlation + 1
          enddo
        enddo
        allocate(lagtime(1 : nsta_correlation))
        nlagtime = 1
        do j = 1, nsta_count(ii, jj) - 1
          waveform_fft(1 : ntime_fft, 1) &
          &  = waveform_obs(time_index - ntime_fft + 1 : time_index, grid_stationindex(j, ii, jj))
          do i = j + 1, nsta_count(ii, jj)
            waveform_fft(1 : ntime_fft, 2) &
            &  = waveform_obs(time_index - ntime_fft + 1 : time_index, grid_stationindex(i, ii, jj))

            call correlation_fft(waveform_fft, ntime_fft, xcorr)
            max_xcorr = maxloc(xcorr) - int(ntime_fft / 2)
            lagtime(nlagtime) = real(max_xcorr(1), kind = fp) * dt
            nlagtime = nlagtime + 1
          enddo
        enddo
        slowness_xcorr(1 : 2, ii, jj) &
        &  = matmul(slowness_est_matrix(1 : 2, 1 : nsta_correlation, ii, jj), lagtime(1 : nsta_correlation)
        deallocate(lagtime)

        !!calculate spatial gradients of wavefield
        !set up observation vector
        obs_vector(1 : nsta_grid_max) = 0.0_fp
        do j = 1, nsta_count(ii, jj)
          dx_east = location_grid(ii, jj)%x_east - location_sta(grid_stationindex(j, ii, jj))%x_east
          dy_north = location_grid(ii, jj)%y_north - location_sta(grid_stationindex(j, ii, jj))%y_north
          time_diff = int((dx_east * slowness_xcorr(1, ii, jj) + dy_north * slowness_xcorr(2, ii, jj)) * dt + 0.5_fp)

          obs_vector(j) = waveform_obs(time_index - ntimestep + j, grid_stationindex(i, ii, jj))
        enddo





          waveform_est_tmp(1 : 3, j, ii, jj, kk) &
          & = matmul(kernel_matrix(1 : 3, 1 : nsta_grid_max, ii, jj), obs_vector(1 : nsta_grid_max))
        enddo
        do j = 1, ntimestep
          waveform_est(1 : 3, ii, jj, kk) = waveform_est(1 : 3, ii, jj, kk) + waveform_est_tmp(1 : 3, j, ii, jj, kk)
        enddo
        waveform_est(1 : 3, ii, jj, kk) = waveform_est(1 : 3, ii, jj, kk) / real(ntimestep, kind = fp)
      enddo
    enddo
    !!Then, 
      


  do j = 1, ntime
    write(time_index, '(i4)') j
    write(0, '(a)') "calculate amplitude distribution for time " // time_index
    do i = 1, 4
      if(time_index(i : i) .eq. " ") time_index(i : i) = "0"
    enddo
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        if(grid_enough_sta(ii, jj) .eqv. .false.) cycle
        obs_vector(1 : nsta_grid_max) = 0.0_fp
        do i = 1, nsta_count(ii, jj)
          obs_vector(i) = waveform_obs(j, grid_stationindex(i, ii, jj))
        enddo
        waveform_est(1 : 3, ii, jj, j) &
        &  = matmul(kernel_matrix(1 : 3, 1 : nsta_grid_max, ii, jj), obs_vector(1 : nsta_grid_max))
      enddo
    enddo
    outfile = "amplitude_gradiometry_" // trim(time_index) // ".grd"
    call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est(1, :, :, j), outfile)
    !outfile = "amp_dx_gradiometry_" // trim(time_index) // ".grd"
    !call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est(2, :, :, j), outfile)
    !outfile = "amp_dy_gradiometry_" // trim(time_index) // ".grd"
    !call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est(3, :, :, j), outfile)
  enddo

  !!calculate time derivatives of estimated wavefield
  waveform_est_max(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  waveform_est_dt_max(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  do i = 1, ntime
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        if(i .eq. 1) then
          waveform_est_dt(ii, jj, i) = (waveform_est(1, ii, jj, i + 1) - waveform_est(1, ii, jj, i)) / dt
        elseif(i .eq. ntime) then
          waveform_est_dt(ii, jj, i) = (waveform_est(1, ii, jj, i) - waveform_est(1, ii, jj, i - 1)) / dt
        else
          waveform_est_dt(ii, jj, i) = (waveform_est(1, ii, jj, i + 1) - waveform_est(1, ii, jj, i - 1)) / dt * 0.5_fp
        endif
        if(abs(waveform_est_dt(ii, jj, i)) .ge. waveform_est_dt_max(ii, jj)) then
          waveform_est_dt_max(ii, jj) = abs(waveform_est_dt(ii, jj, i))
        endif
        if(abs(waveform_est(1, ii, jj, i)) .ge. waveform_est_max(ii, jj)) then
          waveform_est_max(ii, jj) = abs(waveform_est(1, ii, jj, i))
        endif
      enddo
    enddo
  enddo

  !!calculate slownesses at each grid
  do j = 1, ntime
    write(time_index, '(i4)') j
    do i = 1, 4
      if(time_index(i : i) .eq. " ") time_index(i : i) = "0"
    enddo
    write(0, '(a)') "calculate slowness distribution for time " // time_index
    outfile = "slowness_gradiometry_" // trim(time_index) // ".dat"
    open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 10, status = "replace")
    if(j - ntime_slowness2 .lt. 1 .or. j + ntime_slowness2 .gt. ntime) then
      close(10)
      cycle
    endif
    icount = 1
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        sigma_slowness_x(ii, jj) = 0.0_fp
        sigma_slowness_y(ii, jj) = 0.0_fp
        sigma_ampterm_x(ii, jj) = 0.0_fp
        sigma_ampterm_y(ii, jj) = 0.0_fp
        if(grid_enough_sta(ii, jj) .eqv. .false.) cycle
        do i = -ntime_slowness2, ntime_slowness2
          waveform_est_tmp(i) = waveform_est(1, ii, jj, i + j)
          waveform_est_dt_tmp(i) = waveform_est_dt(ii, jj, i + j)
          waveform_est_dx(i) = waveform_est(2, ii, jj, i + j)
          waveform_est_dy(i) = waveform_est(3, ii, jj, i + j)
        enddo
        if(( (dot_product(waveform_est_tmp, waveform_est_tmp)    * dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp)) &
        &  - (dot_product(waveform_est_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_tmp, waveform_est_dt_tmp))) &
        &  / (waveform_est_max(ii, jj) ** 2 * waveform_est_dt_max(ii, jj) ** 2) &
        &  .gt. eps) then

          !!estimate slowness term
          slowness_x(ii, jj) &
          &  = ((dot_product(waveform_est_tmp, waveform_est_tmp)    * dot_product(waveform_est_dx,     waveform_est_dt_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_dx,     waveform_est_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp)    * dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_tmp,    waveform_est_dt_tmp)))
          slowness_y(ii, jj) &
          &  = ((dot_product(waveform_est_tmp, waveform_est_tmp)    * dot_product(waveform_est_dy,     waveform_est_dt_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_dy,     waveform_est_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp)    * dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_tmp,    waveform_est_dt_tmp)))

          !!estimate amplitude term
          ampterm_x(ii, jj) &
          &  = ((dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_dx,     waveform_est_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp)    * dot_product(waveform_est_dx,     waveform_est_dt_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp)       * dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp)    * dot_product(waveform_est_tmp,    waveform_est_dt_tmp)))
          ampterm_y(ii, jj) &
          &  = ((dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp) * dot_product(waveform_est_dy,     waveform_est_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp)    * dot_product(waveform_est_dy,     waveform_est_dt_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp)       * dot_product(waveform_est_dt_tmp, waveform_est_dt_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_dt_tmp)    * dot_product(waveform_est_tmp,    waveform_est_dt_tmp)))

          !!error estimation using jackknife method
          param_x_err(1 : 2) = 0.0_fp
          param_y_err(1 : 2) = 0.0_fp
          do j3 = -ntime_slowness2, ntime_slowness2
            jcount = 1
            do i3 = -ntime_slowness2, ntime_slowness2
              !print *, i3, j3, jcount
              if(i3 .ne. j3) then
                waveform_errest_tmp(jcount) = waveform_est_tmp(i3)
                waveform_errest_dt_tmp(jcount) = waveform_est_dt_tmp(i3)
                waveform_errest_dx(jcount) = waveform_est_dx(i3)
                waveform_errest_dy(jcount) = waveform_est_dy(i3)
                jcount = jcount + 1
              endif
            enddo
            !!slowness
            param_x_err_tmp(1, j3) &
            &  = ((dot_product(waveform_errest_tmp,    waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dx,     waveform_errest_dt_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp)  &
            &  *   dot_product(waveform_errest_dx,     waveform_errest_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,    waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dt_tmp, waveform_errest_dt_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp)))
            param_y_err_tmp(1, j3) &
            &  = ((dot_product(waveform_errest_tmp,    waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dy,     waveform_errest_dt_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_dy,     waveform_errest_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,    waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dt_tmp, waveform_errest_dt_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp)))
            !!amplitude
            param_x_err_tmp(2, j3) &
            &  = ((dot_product(waveform_errest_dt_tmp, waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_dx,     waveform_errest_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp)  &
            &  *   dot_product(waveform_errest_dx,     waveform_errest_dt_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,    waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dt_tmp, waveform_errest_dt_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp)))
            param_y_err_tmp(2, j3) &
            &  = ((dot_product(waveform_errest_dt_tmp, waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_dy,     waveform_errest_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_dy,     waveform_errest_dt_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,    waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dt_tmp, waveform_errest_dt_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp) &
            &  *   dot_product(waveform_errest_tmp,    waveform_errest_dt_tmp)))
            param_x_err(1 : 2) = param_x_err(1 : 2) + param_x_err_tmp(1 : 2, j3) 
            param_y_err(1 : 2) = param_y_err(1 : 2) + param_y_err_tmp(1 : 2, j3) 
          enddo
          !!calculate mean
          param_x_err(1 : 2) = param_x_err(1 : 2) / real(ntime_slowness, kind = fp)
          param_y_err(1 : 2) = param_y_err(1 : 2) / real(ntime_slowness, kind = fp)
          !!calculate jackknife variance
          do j3 = -ntime_slowness2, ntime_slowness2
            sigma_slowness_x(ii, jj) = sigma_slowness_x(ii, jj) + (param_x_err(1) - param_x_err_tmp(1, j3)) ** 2
            sigma_slowness_y(ii, jj) = sigma_slowness_y(ii, jj) + (param_y_err(1) - param_y_err_tmp(1, j3)) ** 2
            sigma_ampterm_x(ii, jj) = sigma_ampterm_x(ii, jj) + (param_x_err(2) - param_x_err_tmp(2, j3)) ** 2
            sigma_ampterm_y(ii, jj) = sigma_ampterm_y(ii, jj) + (param_y_err(2) - param_y_err_tmp(2, j3)) ** 2
          enddo
          sigma_slowness_x(ii, jj) = sqrt(sigma_slowness_x(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
          sigma_slowness_y(ii, jj) = sqrt(sigma_slowness_y(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
          sigma_ampterm_x(ii, jj) = sqrt(sigma_ampterm_x(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
          sigma_ampterm_y(ii, jj) = sqrt(sigma_ampterm_y(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
        else
          slowness_x(ii, jj) = 0.0_fp
          slowness_y(ii, jj) = 0.0_fp
          ampterm_x(ii, jj) = 0.0_fp
          ampterm_y(ii, jj) = 0.0_fp
        endif
        write(10, rec = icount) real(x_start + dgrid_x * real(ii - 1, kind = fp), kind = sp), &
        &                       real(y_start + dgrid_y * real(jj - 1, kind = fp), kind = sp), &
        &                       real(-slowness_x(ii, jj), kind = sp), real(-slowness_y(ii, jj), kind = sp), &
        &                       real(sigma_slowness_x(ii, jj), kind = sp), real(sigma_slowness_y(ii, jj), kind = sp), & 
        &                       real(ampterm_x(ii, jj), kind = sp), real(ampterm_y(ii, jj), kind = sp), &
        &                       real(sigma_ampterm_x(ii, jj), kind = sp), real(sigma_ampterm_y(ii, jj), kind = sp)
        icount = icount + 1
      enddo
    enddo
    close(10)
  enddo 


  stop
end program seismicgradiometry_reducingvelocity

