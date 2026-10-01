!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

program seismicgradiometry_reducingvelocity2_lonlat
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use typedef
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy, xy2bl
  use itoa
  use tandem
  use m_readini
  use calc_kernelmatrix_lonlat

  implicit none

  real(kind = fp) :: dgrid_lon, dgrid_lat, lon_w, lon_e, lat_s, lat_n, center_lon, center_lat, cutoff_dist, &
  &                  dt, fl, fh, fs, ap, as, order, eps, direction
  integer         :: ngradient2, ntimestep, nsta_grid_min, nsta_grid_max, naddstation_array, niteration_max

  type(location), allocatable :: location_grid(:, :), &
  &                              location_sta(:),     &
  &                              location_gridlocal(:), &
  &                              location_stalocal(:, :)

  real(kind = fp), allocatable :: h(:), uv(:, :), &
  &                               waveform_obs(:, :), &
  &                               begin(:), &
  &                               obsvector(:), &
  &                               slowness(:, :), &
  &                               slowness_correction(:, :), &
  &                               sigma_slowness(:, :), &
  &                               ampterm(:, :), &
  &                               sigma_ampterm(:, :), & 
  &                               waveform_est_tmp(:, :), &
  &                               waveform_est_tmp2(:, :), &
  &                               waveform_est_plot(:, :), &
  &                               geospread(:, :), &
  &                               radpattern(:, :), &
  &                               kernel_matrix(:, :, :, :), &
  &                               kernel_matrix_local(:, :, :), &
  &                               error_matrix(:, :, :), &
  &                               error_matrix_local(:, :)
  integer, allocatable         :: grid_stationindex(:, :, :), &
  &                               grid_stationindex_local(:, :), &
  &                               nsta_count(:, :), &
  &                               applicable_gridindex_lon(:), &
  &                               applicable_gridindex_lat(:), &
  &                               nstalocal_count(:), &
  &                               timeindex_diff(:)
  logical, allocatable         :: grid_enough_sta(:, :)
  character(len = 129), allocatable :: sacfile(:)

  integer :: nsta, npts_tmp, timeindex, timeindex_diff_min, i, j, k, ii, jj, ncount, m, n, ngrad, unitnum, unitnum2, index_tmp, &
  &          ngrid_calc, ngrid_lon, ngrid_lat, ntime, info
  real(kind = fp)              :: c, gn, dx_east, dy_north, denominator, uu, uut, utut, max_innerproduct, &
  &                               innerproduct_tmp, relativeerror, error_omega, slowness_cor_prev(1 : 2), &
  &                               uxu(1 : 2), uxut(1 : 2), numerator_slowness(1 : 2), numerator_ampterm(1 : 2)
  logical                      :: calc_grad
  character(len = 129) :: outfile, paramfile
  character(len = 4) :: ctimeindex

  call get_command_argument(1, value = paramfile)


  !!All parameters must be presented in the file
  call readini__strict_mode(.true.) 
  !!read parameter file
  open(newunit = unitnum, file = trim(paramfile))
  call readini(unitnum, "lon_w", lon_w)
  call readini(unitnum, "lon_e", lon_e)
  call readini(unitnum, "lat_s", lat_s)
  call readini(unitnum, "lat_n", lat_n)
  call readini(unitnum, "center_lon", center_lon)
  call readini(unitnum, "center_lat", center_lat)
  call readini(unitnum, "dgrid_lon", dgrid_lon)
  call readini(unitnum, "dgrid_lat", dgrid_lat)
  call readini(unitnum, "nsta_grid_min", nsta_grid_min)
  call readini(unitnum, "nsta_grid_max", nsta_grid_max)
  call readini(unitnum, "cutoff_dist", cutoff_dist)
  call readini(unitnum, "naddstation_array", naddstation_array)
  call readini(unitnum, "dt", dt)
  call readini(unitnum, "fl", fl); fl = 1.0_fp / (fl * 60.0_fp)   !!minute -> second, period -> frequency
  call readini(unitnum, "fh", fh); fh = 1.0_fp / (fh * 60.0_fp)   !!minute -> second, period -> frequency
  call readini(unitnum, "fs", fs); fs = 1.0_fp / (fs * 60.0_fp)   !!minute -> second, period -> frequency
  call readini(unitnum, "ap", ap)
  call readini(unitnum, "as", as)
  call readini(unitnum, "ngradient2", ngradient2)
  call readini(unitnum, "order", order)
  call readini(unitnum, "ntimestep", ntimestep)
  call readini(unitnum, "niteration_max", niteration_max)
  call readini(unitnum, "eps", eps)
  close(unitnum)
  error_omega = 2.0_fp * pi * 1.0_fp / ((1.0_fp / fl + 1.0_fp / fh) * 0.5_fp)

  write(0, '(a, 3(e15.7, 1x))') "fl, fh, fs (Hz) = ", fl, fh, fs
  write(0, '(a, 3(e15.7, 1x))') "fl, fh, fs (s) = ", 1.0_fp / fl, 1.0_fp / fh, 1.0_fp / fs
  write(0, '(a, 4(e15.7, 1x))') "lon_w, lon_e, lat_s, lat_n = ", lon_w, lon_e, lat_s, lat_n
  write(0, '(a, 2(f4.2, 1x))')  "dgrid_lon, dgrid_lat = ", dgrid_lon, dgrid_lat

  nsta = command_argument_count() - 1
  allocate(location_sta(1 : nsta), sacfile(1 : nsta), begin(1 : nsta))
  do i = 1, nsta
    call get_command_argument(i + 1, value = sacfile(i))
  enddo

  !!read sac-formatted waveforms
  ntime = huge(0)
  do i = 1, nsta
    call read_sachdr(sacfile(i), begin = begin(i), npts = npts_tmp, &
    &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
    if(npts_tmp .le. ntime) ntime = npts_tmp
  enddo
  allocate(waveform_obs(1 : ntime, 1 : nsta))
  do i = 1, nsta
    call read_sacdata(sacfile(i), ntime, waveform_obs(:, i))
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
  enddo

  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(1 : 4 * m), uv(1 : 4 * m, 1 : nsta))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  do i = 1, nsta
    call tandem3(waveform_obs(:, i), h, gn, 1, past_uv = uv(:, i))
  enddo
  !uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  !do i = 1, nsta
  !  call tandem3(waveform_obs(:, i), h, gn, -1, past_uv = uv(:, i))
  !enddo
  deallocate(h, uv)

  ngrid_lon = int((lon_e - lon_w) / dgrid_lon + 0.5_fp) + 1
  ngrid_lat = int((lat_n - lat_s) / dgrid_lat + 0.5_fp) + 1
  !!set grid location for constructing arrays
  allocate(location_grid(1 : ngrid_lon, 1 : ngrid_lat), &
  &        nsta_count(1 : ngrid_lon, 1 : ngrid_lat), &
  &        grid_stationindex(1 : 3 + naddstation_array, 1 : ngrid_lon, 1 : ngrid_lat), &
  &        grid_enough_sta(1 : ngrid_lon, 1 : ngrid_lat), &
  &        kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_lon, 1 : ngrid_lat), &
  &        error_matrix(1 : 3, 1 : ngrid_lon, 1 : ngrid_lat))
  do j = 1, ngrid_lat
    do i = 1, ngrid_lon
      location_grid(i, j)%lon = lon_w + dgrid_lon * real(i - 1, kind = fp)
      location_grid(i, j)%lat = lat_s + dgrid_lat * real(j - 1, kind = fp)
      call bl2xy(location_grid(i, j)%lon,     location_grid(i, j)%lat, center_lon, center_lat, &
      &          location_grid(i, j)%y_north, location_grid(i, j)%x_east)
      !!convert meters to kilometers
      location_grid(i, j)%x_east  = location_grid(i, j)%x_east  * 1.0e-3_fp
      location_grid(i, j)%y_north = location_grid(i, j)%y_north * 1.0e-3_fp
    enddo
  enddo

  !!convert station longitude/latitude to x_east/y_north
  open(newunit = unitnum, file = "station_location.txt")
  do i = 1, nsta
    call bl2xy(location_sta(i)%lon, location_sta(i)%lat, center_lon, center_lat, &
    &          location_sta(i)%y_north, location_sta(i)%x_east)
    !! convert meters to kilometers
    location_sta(i)%y_north = location_sta(i)%y_north * 1.0e-3_fp
    location_sta(i)%x_east  = location_sta(i)%x_east  * 1.0e-3_fp
    write(unitnum, '(3(e15.7, 1x))') location_sta(i)%lon, location_sta(i)%lat, location_sta(i)%depth
  enddo
  close(unitnum)

  !!make kernel matrix for each grid
  call calc_kernelmatrix_delaunay2(paramfile, location_grid, location_sta, grid_enough_sta, &
  &                                nsta_count, grid_stationindex, kernel_matrix = kernel_matrix, error_matrix = error_matrix)

  !!count grids where the wave gradiometry is applicable
  ngrid_calc = 0
  do j = 1, ngrid_lat
    do i = 1, ngrid_lon
      if(grid_enough_sta(i, j) .eqv. .true.) ngrid_calc = ngrid_calc + 1
    enddo
  enddo
  write(0, '(a, i0)') "The number of grids where the analysis is applicable = ", ngrid_calc

  !!substitute grid indices which the analysis is applicable to the array
  !!recalculate array geometry (global coordinate -> grid-local coordinate)
  allocate(applicable_gridindex_lon(1 : ngrid_calc), &
  &        applicable_gridindex_lat(1 : ngrid_calc))
  allocate(location_gridlocal(1 : ngrid_calc),                            &
  &        location_stalocal(1 : nsta_grid_max, 1 : ngrid_calc),          &
  &        nstalocal_count(1 : ngrid_calc),                               &
  &        grid_stationindex_local(1 : nsta_grid_max, 1 : ngrid_calc),    &
  &        kernel_matrix_local(1 : 3, 1 : nsta_grid_max, 1 : ngrid_calc), &
  &        error_matrix_local(1 : 3, 1 : ngrid_calc))

  index_tmp = 1
  do jj = 1, ngrid_lat
    do ii = 1, ngrid_lon
      if(grid_enough_sta(ii, jj) .eqv. .true.) then
        applicable_gridindex_lon(index_tmp) = ii
        applicable_gridindex_lat(index_tmp) = jj

        nstalocal_count(index_tmp) = nsta_count(ii, jj)
        grid_stationindex_local(1 : nstalocal_count(index_tmp), index_tmp) = grid_stationindex(1 : nsta_count(ii, jj), ii, jj)

        location_gridlocal(index_tmp)%lon = location_grid(ii, jj)%lon
        location_gridlocal(index_tmp)%lat = location_grid(ii, jj)%lat
        call bl2xy(location_gridlocal(index_tmp)%lon,     location_gridlocal(index_tmp)%lat,    &
        &          location_grid(ii, jj)%lon,             location_grid(ii, jj)%lat,            &
        &          location_gridlocal(index_tmp)%y_north, location_gridlocal(index_tmp)%x_east)
        !!convert meters to kilometers (but must be zero)
        location_gridlocal(index_tmp)%x_east  = location_gridlocal(index_tmp)%x_east  * 1.0e-3_fp
        location_gridlocal(index_tmp)%y_north = location_gridlocal(index_tmp)%y_north * 1.0e-3_fp

        do i = 1, nstalocal_count(index_tmp)
          location_stalocal(i, index_tmp)%lon = location_sta(grid_stationindex(i, ii, jj))%lon
          location_stalocal(i, index_tmp)%lat = location_sta(grid_stationindex(i, ii, jj))%lat
          call bl2xy(location_stalocal(i, index_tmp)%lon,     location_stalocal(i, index_tmp)%lat,    &
          &          location_gridlocal(index_tmp)%lon,       location_gridlocal(index_tmp)%lat,      &
          &          location_stalocal(i, index_tmp)%y_north, location_stalocal(i, index_tmp)%x_east)
          !!convert meters to kilometers
          location_stalocal(i, index_tmp)%x_east  = location_stalocal(i, index_tmp)%x_east  * 1.0e-3_fp
          location_stalocal(i, index_tmp)%y_north = location_stalocal(i, index_tmp)%y_north * 1.0e-3_fp
        enddo
        call calculate_kernelmatrix(nstalocal_count(index_tmp), &
        &                           location_stalocal(1 : nstalocal_count(index_tmp), index_tmp), &
        &                           location_gridlocal(index_tmp), &
        &                           cutoff_dist, info, &
        &                           kernel_matrix_local(1 : 3, 1 : nstalocal_count(index_tmp), index_tmp), &
        &                           error_matrix_local(1 : 3, index_tmp))
        index_tmp = index_tmp + 1
      endif
    enddo
  enddo
  deallocate(grid_enough_sta, nsta_count, grid_stationindex, location_grid, location_sta, kernel_matrix, error_matrix)

  open(newunit = unitnum, file = "station_grid.txt")
  do j = 1, ngrid_calc
    write(unitnum, '(4(e15.7, 1x), i0)') location_gridlocal(j)%lon,    location_gridlocal(j)%lat, &
    &                                    location_gridlocal(j)%x_east, location_gridlocal(j)%y_north, nstalocal_count(j)
    do i = 1, nstalocal_count(j)
      write(unitnum, '(4(e15.7, 1x))') location_stalocal(i, j)%lon,    location_stalocal(i, j)%lat, &
      &                                location_stalocal(i, j)%x_east, location_stalocal(i, j)%y_north
    enddo
    write(unitnum, '(a)') ">"
  enddo
  close(unitnum)

  !!calculate amplitude and its spatial derivatives at each grid
  allocate(obsvector(1 : nsta_grid_max),    &
  &        slowness(1 : 2, 1 : ngrid_calc), &
  &        slowness_correction(1 : 2, 1 : ngrid_calc), &
  &        sigma_slowness(1 : 2, 1 : ngrid_calc), &
  &        ampterm(1 : 2, 1 : ngrid_calc), &
  &        sigma_ampterm(1 : 2, 1 : ngrid_calc), &
  &        waveform_est_tmp(1 : 4, 1 : ngradient2), &
  &        waveform_est_tmp2(1 : ngradient2, 1 : 4), &
  &        waveform_est_plot(1 : ngrid_lon, 1 : ngrid_lat), &
  &        timeindex_diff(1 : nsta_grid_max), &
  &        geospread(1 : ngrid_lon, 1 : ngrid_lat), &
  &        radpattern(1 : ngrid_lon, 1 : ngrid_lat))

  do ii = 1, int(ntime / ntimestep)
    timeindex = ntimestep * (ii - 1) + 1
    write(0, '(a, i0, a)') "Time index = ", ii, " Calculate amplitudes and their gradients at each grid"

    if(timeindex .lt. 1 .or. timeindex - ntimestep .gt. ntime) cycle

    call int_to_char(ii, 4, ctimeindex)
    outfile = "slowness_gradiometry_" // trim(ctimeindex) // ".dat"
    open(newunit = unitnum2, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 10, status = "replace")
    ncount = 1

    !!Estimate slowness vector with reducing velocity
    slowness(1 : 2, 1 : ngrid_calc)                 = 0.0_fp
    slowness_correction(1 : 2, 1 : ngrid_calc)      = 0.0_fp
    sigma_slowness(1 : 2, 1 : ngrid_calc)           = 0.0_fp
    ampterm(1 : 2, 1 : ngrid_calc)                  = 0.0_fp
    sigma_ampterm(1 : 2, 1 : ngrid_calc)            = 0.0_fp
    waveform_est_plot(1 : ngrid_lon, 1 : ngrid_lat) = 0.0_fp
    geospread(1 : ngrid_lon, 1 : ngrid_lat)         = 0.0_fp
    radpattern(1 : ngrid_lon, 1 : ngrid_lat)        = 0.0_fp

    do k = 1, ngrid_calc
      !!calculate spatial gradients of wavefield iteratively
      n = 0
      do
        calc_grad = .true.
        n = n + 1
        if(n .gt. niteration_max) exit
        !print *, "iterative reducing, count = ", n
        waveform_est_tmp(1 : 4, 1 : ngradient2) = 0.0_fp
        waveform_est_tmp2(1 : ngradient2, 1 : 4) = 0.0_fp

        do i = 1, nstalocal_count(k)
          dx_east  = location_gridlocal(k)%x_east  - location_stalocal(i, k)%x_east
          dy_north = location_gridlocal(k)%y_north - location_stalocal(i, k)%y_north
          timeindex_diff(i) = int((dx_east * slowness(1, k) + dy_north * slowness(2, k)) / dt)
        enddo
        timeindex_diff_min = minval(timeindex_diff)

        ngrad = 0
        do j = 1, ngradient2
          obsvector(1 : nsta_grid_max) = 0.0_fp
          do i = 1, nstalocal_count(k)
            if(timeindex - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min .lt. 1 .or. &
               timeindex - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min .gt. ntime) then
              calc_grad = .false.
              exit
            endif
            obsvector(i) &
            &  = waveform_obs(timeindex - ngradient2 + j - timeindex_diff(i) + timeindex_diff_min, grid_stationindex_local(i, k))
          enddo
          if(calc_grad .eqv. .false.) exit

          waveform_est_tmp(1 : 3, j) &
          &  = matmul(kernel_matrix_local(1 : 3, 1 : nstalocal_count(k), k), obsvector(1 : nstalocal_count(k)))
          ngrad = ngrad + 1
        enddo

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
        if(n .eq. 1) waveform_est_plot(applicable_gridindex_lon(k), applicable_gridindex_lat(k)) = waveform_est_tmp(1, ngrad)

        !!estimate slowness term
        slowness_cor_prev(1 : 2) = slowness_correction(1 : 2, k)
        ampterm(1 : 2, k) = 0.0_fp
        uu   = dot_product(waveform_est_tmp2(1 : ngrad, 1), waveform_est_tmp2(1 : ngrad, 1))
        uut  = dot_product(waveform_est_tmp2(1 : ngrad, 1), waveform_est_tmp2(1 : ngrad, 4))
        utut = dot_product(waveform_est_tmp2(1 : ngrad, 4), waveform_est_tmp2(1 : ngrad, 4))
        denominator = uu * utut - uut * uut
        do i = 1, 2
          uxu(i)  = dot_product(waveform_est_tmp2(1 : ngrad, i + 1), waveform_est_tmp2(1 : ngrad, 1))
          uxut(i) = dot_product(waveform_est_tmp2(1 : ngrad, i + 1), waveform_est_tmp2(1 : ngrad, 4))
          numerator_slowness(i) = uu   * uxut(i) - uut * uxu(i)
          numerator_ampterm(i)  = utut * uxu(i)  - uut * uxut(i)

          slowness_correction(i, k) = -numerator_slowness(i) / denominator
          ampterm(i, k)             =  numerator_ampterm(i)  / denominator
        enddo
        slowness(1 : 2, k) = slowness(1 : 2, k) + slowness_correction(1 : 2, k)


        if ((slowness_correction(1, k) - slowness_cor_prev(1)) ** 2 &
        & + (slowness_correction(2, k) - slowness_cor_prev(2)) ** 2 .lt. eps) exit
      enddo

      if(calc_grad .eqv. .false.) cycle

      !!error estimation
      max_innerproduct = 0.0_fp
      do i = 1, nstalocal_count(k)
        dx_east  = location_gridlocal(k)%x_east  - location_stalocal(i, grid_stationindex_local(i, k))%x_east
        dy_north = location_gridlocal(k)%y_north - location_stalocal(i, grid_stationindex_local(i, k))%y_north
        innerproduct_tmp = (slowness(1, k) * dx_east + slowness(2, k) * dy_north) ** 2
        if(innerproduct_tmp .ge. max_innerproduct) then
          !print *, dx_east, dy_north, innerproduct_tmp
          max_innerproduct = innerproduct_tmp
        endif
      enddo
      relativeerror = 0.5_fp * error_omega ** 2 * max_innerproduct
      !print *, ii, jj, error_omega, max_innerproduct, relativeerror

      do j = 1, ngrad
        do i = 1, 2
          sigma_slowness(i, k) = sigma_slowness(i, k) &
                                    !!dB/du
          &                         + ((2.0_fp * waveform_est_tmp(1, j)     * uxut(i) &
          &                                    - waveform_est_tmp(4, j)     * uxu (i) &
          &                                    - waveform_est_tmp(i + 1, j) * uut) / denominator &
          &                           - 2.0_fp * numerator_slowness(i) / (denominator ** 2) &
          &                                    * (waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut)) ** 2 &
          &                         * (waveform_est_tmp(1, j) * relativeerror) ** 2 * error_matrix_local(1, k) &
                                    !!dB/dut
          &                         + ((waveform_est_tmp(i + 1, j) * uu - waveform_est_tmp(1, j) * uxu(i)) / denominator &
          &                           - 2.0_fp * (numerator_slowness(i) / denominator ** 2) &
          &                                    * (waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut)) ** 2 &
          &                         * (waveform_est_tmp(4, j) * 2.0_fp / (dt ** 2) * relativeerror ) ** 2 &
          &                         * error_matrix_local(1, k) &
                                    !!dB/dux
          &                         + ((waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut) / denominator) ** 2 &
          &                         * (waveform_est_tmp(i + 1, j) * relativeerror) ** 2 * error_matrix_local(i + 1, k)

          sigma_ampterm(i, k) = sigma_ampterm(i, k) &
                                   !!dA/du
          &                        + ((waveform_est_tmp(i + 1, j) * utut - waveform_est_tmp(4, j) * uxut(i)) / denominator &
          &                           - 2.0_fp * (numerator_ampterm(i) / denominator ** 2) &
          &                                    * (waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut)) ** 2 &
          &                        * (waveform_est_tmp(1, j) * relativeerror) ** 2 * error_matrix_local(1, k) &
                                   !!dA/dut
          &                        + ((2.0_fp * waveform_est_tmp(4, j)     * uxu (i) &
          &                                   - waveform_est_tmp(1, j)     * uxut(i) &
          &                                   - waveform_est_tmp(i + 1, j) * uut) / denominator &
          &                          - 2.0_fp * numerator_ampterm(i) / (denominator ** 2) &
          &                                   * (waveform_est_tmp(4, j) * uu - waveform_est_tmp(1, j) * uut)) ** 2 &
          &                        * (waveform_est_tmp(4, j) * 2.0_fp / (dt ** 2) * relativeerror)  ** 2 &
          &                        * error_matrix_local(1, k) &
                                   !!dA/dux
          &                        + ((waveform_est_tmp(1, j) * utut - waveform_est_tmp(4, j) * uut) / denominator) ** 2 &
          &                        * (waveform_est_tmp(i + 1, j) * relativeerror) ** 2 * error_matrix_local(i + 1, k)
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

      write(unitnum2, rec = ncount) real(lon_w + dgrid_lon * real(applicable_gridindex_lon(k) - 1, kind = fp), kind = sp), &
      &                             real(lat_s + dgrid_lat * real(applicable_gridindex_lat(k) - 1, kind = fp), kind = sp), &
      &                             real(slowness(1, k),       kind = sp), real(slowness(2, k),       kind = sp), & 
      &                             real(sigma_slowness(1, k), kind = sp), real(sigma_slowness(2, k), kind = sp), &
      &                             real(ampterm(1, k),        kind = sp), real(ampterm(2, k),        kind = sp), &
      &                             real(sigma_ampterm(1, k),  kind = sp), real(sigma_ampterm(2, k),  kind = sp)
      ncount = ncount + 1

      direction = atan2(slowness(1, k), slowness(2, k))
      geospread(applicable_gridindex_lon(k), applicable_gridindex_lat(k)) &
      &  = (ampterm(1, k) * sin(direction) + ampterm(2, k) * cos(direction)) * 1.0e+2_fp
      radpattern(applicable_gridindex_lon(k), applicable_gridindex_lat(k)) &
      &  = (ampterm(1, k) * cos(direction) - ampterm(2, k) * sin(direction)) * 1.0e+2_fp
    enddo
    close(unitnum2)

    outfile = "amplitude_gradiometry_" // trim(ctimeindex) // ".grd"
    outfile = trim(outfile)
    call write_grdfile_fp_2d(lon_w, lat_s, dgrid_lon, dgrid_lat, ngrid_lon, ngrid_lat, waveform_est_plot, outfile, &
    &                        nanval = 0.0_fp, xlabel = "Longitude", ylabel = "Latitude", zlabel = "Pressure")
    outfile = "geospread_gradiometry_" // trim(ctimeindex) // ".grd"
    outfile = trim(outfile)
    call write_grdfile_fp_2d(lon_w, lat_s, dgrid_lon, dgrid_lat, ngrid_lon, ngrid_lat, geospread, outfile, &
    &                        nanval = 0.0_fp, xlabel = "Longitude", ylabel = "Latitude", zlabel = "Geo_spread")
    outfile = "radpattern_gradiometry_" // trim(ctimeindex) // ".grd"
    outfile = trim(outfile)
    call write_grdfile_fp_2d(lon_w, lat_s, dgrid_lon, dgrid_lat, ngrid_lon, ngrid_lat, radpattern, outfile, &
    &                        nanval = 0.0_fp, xlabel = "Longitude", ylabel = "Latitude", zlabel = "Rad_pattern")

  enddo

  stop
end program seismicgradiometry_reducingvelocity2_lonlat


