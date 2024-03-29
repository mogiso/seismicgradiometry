! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

program seismicgradiometry
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy, xy2bl
  use sort, only : bubblesort
  use greatcircle, only : greatcircle_dist
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif

  implicit none

  type location
    real(kind = fp) :: lon, lat, x_east, y_north, depth
  end type location

  real(kind = fp), parameter :: eps = 1.0e-5_fp

  !!For S-net/DONET OBP
  real(kind = fp), parameter :: order = 1.0e-2_fp                             !!Pa -> hpa
  real(kind = fp), parameter :: az_diff_max = 150.0_fp * deg2rad

  real(kind = fp), parameter :: x_start = -350.0_fp, y_start = -600.0_fp, &
  &                             x_end = 350.0_fp, y_end = 600.0_fp
  real(kind = fp), parameter :: center_lon = 142.5_fp, center_lat = 38.25_fp   !!S-net
  real(kind = fp), parameter :: dgrid_x = 20.0_fp, dgrid_y = 20.0_fp          !!S-net test
  real(kind = fp), parameter :: cutoff_dist = 80.0_fp                         !!S-net test
  real(kind = fp), parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (20.0_fp * 60.0_fp), &
  &                             fs = 1.0_fp / (10.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp  !!S-net test

  !real(kind = fp), parameter :: x_start = -150.0_fp, y_start = -100.0_fp, &
  !&                             x_end = 150.0_fp, y_end = 100.0_fp
  !real(kind = fp), parameter :: center_lon = 135.75_fp, center_lat = 33.2_fp   !!DONET test
  !real(kind = fp), parameter :: dgrid_x = 10.0_fp, dgrid_y = 10.0_fp          !!DONET test
  !real(kind = fp), parameter :: cutoff_dist = 80.0_fp                         !!DONET 20-60min
  !real(kind = fp), parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (20.0_fp * 60.0_fp), &
  !&                             fs = 1.0_fp / (10.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp  !!DONET long-period test
  !real(kind = fp), parameter :: cutoff_dist = 30.0_fp                         !!DONET test 6-20min
  !real(kind = fp), parameter :: fl = 1.0_fp / (20.0_fp * 60.0_fp), fh = 1.0_fp / (6.0_fp * 60.0_fp), &
  !&                             fs = 1.0_fp / (3.0_fp * 60.0_fp), ap = 0.5_fp, as = 10.0_fp  !!DONET short-period test

  integer, parameter :: ntime_slowness = 61, ntime_slowness2 = (ntime_slowness - 1) / 2
  integer, parameter :: nsta_grid_max = 40, nsta_grid_min = 5  !!For S-net/DONET OBPG array
  integer, parameter :: ntime = 630
  !integer, parameter :: ntime = 1024 !!testdata

  !!For SK-net Long-period motion
  !real(kind = fp), parameter :: az_diff_max = 150.0_fp * deg2rad
  !real(kind = fp), parameter :: x_start = -125.0_fp, y_start = -165.0_fp, &
  !&                             x_end = 205.0_fp, y_end = 165.0_fp
  !real(kind = fp), parameter :: center_lon = 139.0_fp, center_lat = 36.0_fp
  !real(kind = fp), parameter :: dgrid_x = 2.0_fp, dgrid_y = 2.0_fp
  !real(kind = fp), parameter :: cutoff_dist = 5.0_fp
  !real(kind = fp), parameter :: order = 1.0e-6_fp  !!nm -> mm
  !integer, parameter :: nsta_grid_max = 10, nsta_grid_min = 3
  !integer, parameter :: ntime_slowness = 21, ntime_slowness2 = (ntime_slowness - 1) / 2
  !real(kind = fp), parameter :: fl = 1.0_fp / 10.0_fp, fh = 1.0_fp / 5.0_fp, fs = 1.0_fp / 2.0_fp, &
  !&                             ap = 0.5_fp, as = 5.0_fp
  !integer, parameter :: ntime = 2100

#ifdef ELLIPSE
  real(kind = fp), parameter :: sigma_x = cutoff_dist, sigma_y = cutoff_dist * 2.0_fp
                                !!sigma_x: parallel to propagation direction, sigma_y: normal to propagation direction
  real(kind = fp), parameter :: evlon = 137.8910_fp, evlat = 36.6928_fp
  real(kind = fp)            :: propagation_direction
#endif


  integer, parameter :: ngrid_x = int((x_end - x_start) / real(dgrid_x, kind = fp)) + 1
  integer, parameter :: ngrid_y = int((y_end - y_start) / real(dgrid_y, kind = fp)) + 1
  integer :: nsta, npts_tmp, info, i, j, ii, jj, kk, i3, j3, m, n, icount, jcount
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)
  type(location), allocatable :: location_sta(:)
  real(kind = fp), allocatable :: begin(:)
  real(kind = fp), allocatable :: waveform_obs(:, :)         !!(1 : ntime, 1 : nsta)
  real(kind = fp), allocatable :: waveform_est(:, :, :, :)   !!(1 : 3 (u, dudx, dudy), 1 : ngrid_x, 1 : ngrid_y, 1 : ntime)
  real(kind = fp), allocatable :: waveform_est_diff(:, :, :) !!(1 : ngrid_x, 1 : ngrid_y, 1 : ntime)
  real(kind = fp), allocatable :: h(:)
  real(kind = fp), allocatable :: azimuth_order(:)
  real(kind = fp) :: dist_grid_sta(1 : nsta_grid_max), g(1 : nsta_grid_max, 1 : 3), dist_tmp, &
  &                  azimuth_grid_sta(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), az_diff_tmp, &
  &                  dist_x_tmp, dist_y_tmp, weight(1 : nsta_grid_max, 1 : nsta_grid_max), g_tmp(1 : 3, 1 : 3), &
  &                  kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), obs_vector(1 : nsta_grid_max), &
  &                  slowness_x(1 : ngrid_x, 1 : ngrid_y), slowness_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  sigma_slowness_x(1 : ngrid_x, 1 : ngrid_y), sigma_slowness_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  waveform_est_tmp(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_est_diff_tmp(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_est_dx(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_est_dy(-ntime_slowness2 : ntime_slowness2), &
  &                  waveform_errest_tmp(1 : ntime_slowness - 1), &
  &                  waveform_errest_diff_tmp(1 : ntime_slowness - 1), &
  &                  waveform_errest_dx(1 : ntime_slowness - 1), &
  &                  waveform_errest_dy(1 : ntime_slowness - 1), &
  &                  param_x_err_tmp(1 : 2, -ntime_slowness2 : ntime_slowness2), &
  &                  param_y_err_tmp(1 : 2, -ntime_slowness2 : ntime_slowness2), &
  &                  param_x_err(1 : 2), param_y_err(1 : 2), &
  &                  amp_term_x(1 : ngrid_x, 1 : ngrid_y), amp_term_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  sigma_amp_term_x(1 : ngrid_x, 1 : ngrid_y), sigma_amp_term_y(1 : ngrid_x, 1 : ngrid_y), &
  &                  waveform_est_max(1 : ngrid_x, 1 : ngrid_y), waveform_est_diff_max(1 : ngrid_x, 1 : ngrid_y), &
  &                  dt, c, gn
  integer :: grid_stationindex(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), nsta_count(1 : ngrid_x, 1 : ngrid_y), ipiv(3)
  logical :: grid_enough_sta(1 : ngrid_x, 1 : ngrid_y)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  character(len = 4) :: time_index


  nsta = command_argument_count()

  allocate(waveform_obs(1 : ntime, 1 : nsta), waveform_est_diff(1 : ngrid_x, 1 : ngrid_y, 1 : ntime), &
  &        waveform_est(1 : 3, 1 : ngrid_x, 1 : ngrid_y, 1 : ntime))
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
  allocate(h(4 * m))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  !apply filter
  do i = 1, nsta
    call tandem1(waveform_obs(:, i), waveform_obs(:, i), ntime, h, m, 1)
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * gn
    call tandem1(waveform_obs(:, i), waveform_obs(:, i), ntime, h, m, -1)
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * gn
  enddo

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
  do kk = 1, ngrid_y
    do jj = 1, ngrid_x
      nsta_count(jj, kk) = 0
      grid_enough_sta(jj, kk) = .false.
      grid_stationindex(1 : nsta_grid_max, jj, kk) = 0
      azimuth_grid_sta(1 : nsta_grid_max, jj, kk) = 1.0e+15

#ifdef ELLIPSE
      dist_grid_sta(1 : nsta_grid_max) = 0.0_fp
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      propagation_direction = propagation_direction - pi
      !print *, location_grid(jj, kk)%lon, location_grid(jj, kk)%lat, propagation_direction * rad2deg
      do ii = 1, nsta
        dist_x_tmp =   (location_sta(ii)%y_north - location_grid(jj, kk)%y_north) * cos(propagation_direction) &
        &            + (location_sta(ii)%x_east  - location_grid(jj, kk)%x_east)  * sin(propagation_direction)
        dist_y_tmp = - (location_sta(ii)%y_north - location_grid(jj, kk)%y_north) * sin(propagation_direction) &
        &            + (location_sta(ii)%x_east  - location_grid(jj, kk)%x_east)  * cos(propagation_direction)
        dist_tmp = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
        if(dist_tmp .lt. exp(-2.0_fp)) cycle

        nsta_count(jj, kk) = nsta_count(jj, kk) + 1
        do j = 1, nsta_grid_max
          if(dist_tmp .ge. dist_grid_sta(j)) then
            do i = nsta_grid_max, j + 1, -1
              dist_grid_sta(i) = dist_grid_sta(i - 1)
              grid_stationindex(i, jj, kk) = grid_stationindex(i - 1, jj, kk)
              azimuth_grid_sta(i, jj, kk) = azimuth_grid_sta(i - 1, jj, kk)
            enddo
            dist_grid_sta(j) = dist_tmp
            grid_stationindex(j, jj, kk) = ii
            azimuth_grid_sta(j, jj, kk) = atan2(location_sta(ii)%x_east  - location_grid(jj, kk)%x_east, &
            &                                   location_sta(ii)%y_north - location_grid(jj, kk)%y_north)
            exit
          endif
        enddo
      enddo
      if(nsta_count(jj, kk) .ge. nsta_grid_min) grid_enough_sta(jj, kk) = .true.
      if(nsta_count(jj, kk) .gt. nsta_grid_max) nsta_count(jj, kk) = nsta_grid_max
      !if(grid_enough_sta(jj, kk) .eqv. .true.) print *, jj, kk, nsta_count(jj, kk)
 
#else

      !!count usable stations for each grid
      dist_grid_sta(1 : nsta_grid_max) = 1.0e+15
      do ii = 1, nsta
        dist_tmp = sqrt((location_sta(ii)%x_east  - location_grid(jj, kk)%x_east)  ** 2 &
        &             + (location_sta(ii)%y_north - location_grid(jj, kk)%y_north) ** 2)
        if(dist_tmp .gt. cutoff_dist) cycle
        nsta_count(jj, kk) = nsta_count(jj, kk) + 1
        do j = 1, nsta_grid_max
          if(dist_tmp .le. dist_grid_sta(j)) then
            do i = nsta_grid_max, j + 1, -1
              dist_grid_sta(i) = dist_grid_sta(i - 1)
              grid_stationindex(i, jj, kk) = grid_stationindex(i - 1, jj, kk)
              azimuth_grid_sta(i, jj, kk) = azimuth_grid_sta(i - 1, jj, kk)
            enddo
            dist_grid_sta(j) = dist_tmp
            grid_stationindex(j, jj, kk) = ii
            !azimuth_grid_sta(j, jj, kk) = atan2(location_sta(ii)%x_east  - location_grid(jj, kk)%x_east, &
            !&                                   location_sta(ii)%y_north - location_grid(jj, kk)%y_north)
            call greatcircle_dist(location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
            &                     location_sta(ii)%lat, location_sta(ii)%lon, &
            &                     azimuth = azimuth_grid_sta(j, jj, kk))
            exit
          endif
        enddo
      enddo
      if(nsta_count(jj, kk) .ge. nsta_grid_min) grid_enough_sta(jj, kk) = .true.
      if(nsta_count(jj, kk) .gt. nsta_grid_max) nsta_count(jj, kk) = nsta_grid_max

#endif

      allocate(azimuth_order(1 : nsta_count(jj, kk)))
      azimuth_order(1 : nsta_count(jj, kk)) = azimuth_grid_sta(1 : nsta_count(jj, kk), jj, kk)
      call bubblesort(azimuth_order)
      az_diff_tmp = 2.0_fp * pi - (azimuth_order(1) - azimuth_order(nsta_count(jj, kk)))
      do i = 2, nsta_count(jj, kk)
        if(azimuth_order(i - 1) - azimuth_order(i) .gt. az_diff_tmp) then
          az_diff_tmp = azimuth_order(i - 1) - azimuth_order(i)
        endif
      enddo
      if(az_diff_tmp .gt. az_diff_max) grid_enough_sta(jj, kk) = .false.
      deallocate(azimuth_order)

      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle

      g(1 : nsta_grid_max, 1 : 3) = 0.0_fp
      weight(1 : nsta_grid_max, 1 : nsta_grid_max) = 0.0_fp
      do i = 1, nsta_count(jj, kk)
        g(i, 1) = 1.0_fp
        g(i, 2) = location_sta(grid_stationindex(i, jj, kk))%x_east - location_grid(jj, kk)%x_east
        g(i, 3) = location_sta(grid_stationindex(i, jj, kk))%y_north - location_grid(jj, kk)%y_north
#ifdef ELLIPSE
        dist_x_tmp =   (location_sta(grid_stationindex(i, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * cos(propagation_direction) &
        &            + (location_sta(grid_stationindex(i, jj, kk))%x_east  &
        &            -  location_grid(jj, kk)%x_east)  * sin(propagation_direction)
        dist_y_tmp = - (location_sta(grid_stationindex(i, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * sin(propagation_direction) &
        &            + (location_sta(grid_stationindex(i, jj, kk))%x_east  &
        &            -  location_grid(jj, kk)%x_east)  * cos(propagation_direction)

        weight(i, i) = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
#else
        weight(i, i) = exp(-(g(i, 2) ** 2 + g(i, 3) ** 2) / (cutoff_dist ** 2))
#endif
      enddo
      g_tmp = matmul(transpose(g), matmul(weight, g))
#ifdef MKL
      call getrf(g_tmp, ipiv = ipiv, info = info)
      !write(0, '(a, i0)') "GETRF info = ", info
      call getri(g_tmp, ipiv, info = info)
      !write(0, '(a, i0)') "GETRI info = ", info
#else
      call LA_GETRF(g_tmp, ipiv, info = info)
      !write(0, '(a, i0)') "LA_GETRF info = ", info
      call LA_GETRI(g_tmp, ipiv, info = info)
      !write(0, '(a, i0)') "LA_GETRI info = ", info
#endif
      if(info .ne. 0) then
        grid_enough_sta(jj, kk) = .false.
        cycle
      endif

      kernel_matrix(1 : 3, 1 : nsta_grid_max, jj, kk) = matmul(g_tmp, matmul(transpose(g), weight))

      !if(grid_enough_sta(jj, kk) .eqv. .true.) print *, jj, kk, nsta_count(jj, kk)
    enddo
  enddo

  !!calculate amplitude and its spatial derivatives at each grid
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
    outfile = "amp_dx_gradiometry_" // trim(time_index) // ".grd"
    call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est(2, :, :, j), outfile)
    outfile = "amp_dy_gradiometry_" // trim(time_index) // ".grd"
    call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, waveform_est(3, :, :, j), outfile)
  enddo

  !!calculate time derivatives of estimated wavefield
  waveform_est_max(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  waveform_est_diff_max(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  do i = 1, ntime
    do jj = 1, ngrid_y
      do ii = 1, ngrid_x
        if(i .eq. 1) then
          waveform_est_diff(ii, jj, i) = (waveform_est(1, ii, jj, i + 1) - waveform_est(1, ii, jj, i)) / dt
        elseif(i .eq. ntime) then
          waveform_est_diff(ii, jj, i) = (waveform_est(1, ii, jj, i) - waveform_est(1, ii, jj, i - 1)) / dt
        else
          waveform_est_diff(ii, jj, i) = (waveform_est(1, ii, jj, i + 1) - waveform_est(1, ii, jj, i - 1)) / dt * 0.5_fp
        endif
        if(abs(waveform_est_diff(ii, jj, i)) .ge. waveform_est_diff_max(ii, jj)) then
          waveform_est_diff_max(ii, jj) = abs(waveform_est_diff(ii, jj, i))
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
        sigma_amp_term_x(ii, jj) = 0.0_fp
        sigma_amp_term_y(ii, jj) = 0.0_fp
        if(grid_enough_sta(ii, jj) .eqv. .false.) cycle
        do i = - ntime_slowness2, ntime_slowness2
          waveform_est_tmp(i) = waveform_est(1, ii, jj, i + j)
          waveform_est_diff_tmp(i) = waveform_est_diff(ii, jj, i + j)
          waveform_est_dx(i) = waveform_est(2, ii, jj, i + j)
          waveform_est_dy(i) = waveform_est(3, ii, jj, i + j)
        enddo
        if(((dot_product(waveform_est_tmp, waveform_est_tmp) * dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp)) &
        &  - (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_tmp, waveform_est_diff_tmp))) &
        &  / (waveform_est_max(ii, jj) ** 2 * waveform_est_diff_max(ii, jj) ** 2) &
        &  .gt. eps .and. &
        &  ((dot_product(waveform_est_tmp, waveform_est_tmp) * dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp)) &
        &  - (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_tmp, waveform_est_diff_tmp))) &
        &  / (waveform_est_max(ii, jj) ** 2 * waveform_est_diff_max(ii, jj) ** 2) &
        &  .gt. eps) then

          !!estimate slowness term
          slowness_x(ii, jj) &
          &  = ((dot_product(waveform_est_tmp, waveform_est_tmp)      * dot_product(waveform_est_dx, waveform_est_diff_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_dx, waveform_est_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp) * dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_tmp, waveform_est_diff_tmp)))
          slowness_y(ii, jj) &
          &  = ((dot_product(waveform_est_tmp, waveform_est_tmp)      * dot_product(waveform_est_dy, waveform_est_diff_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_dy, waveform_est_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp) * dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_tmp, waveform_est_diff_tmp)))

          !!estimate amplitude term
          amp_term_x(ii, jj) &
          &  = ((dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_dx, waveform_est_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_dx, waveform_est_diff_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp) * dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_tmp, waveform_est_diff_tmp)))
          amp_term_y(ii, jj) &
          &  = ((dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_dy, waveform_est_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_dy, waveform_est_diff_tmp))) &
          &  / ((dot_product(waveform_est_tmp, waveform_est_tmp) * dot_product(waveform_est_diff_tmp, waveform_est_diff_tmp)) &
          &  -  (dot_product(waveform_est_tmp, waveform_est_diff_tmp) * dot_product(waveform_est_tmp, waveform_est_diff_tmp)))

          !!error estimation using jackknife method
          param_x_err(1 : 2) = 0.0_fp
          param_y_err(1 : 2) = 0.0_fp
          do j3 = -ntime_slowness2, ntime_slowness2
            jcount = 1
            do i3 = -ntime_slowness2, ntime_slowness2
              !print *, i3, j3, jcount
              if(i3 .ne. j3) then
                waveform_errest_tmp(jcount) = waveform_est_tmp(i3)
                waveform_errest_diff_tmp(jcount) = waveform_est_diff_tmp(i3)
                waveform_errest_dx(jcount) = waveform_est_dx(i3)
                waveform_errest_dy(jcount) = waveform_est_dy(i3)
                jcount = jcount + 1
              endif
            enddo
            !!slowness
            param_x_err_tmp(1, j3) &
            &  = ((dot_product(waveform_errest_tmp,      waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dx,       waveform_errest_diff_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp)  &
            &  *   dot_product(waveform_errest_dx,       waveform_errest_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,      waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_diff_tmp, waveform_errest_diff_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp)))
            param_y_err_tmp(1, j3) &
            &  = ((dot_product(waveform_errest_tmp,      waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_dy,       waveform_errest_diff_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_dy,       waveform_errest_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,      waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_diff_tmp, waveform_errest_diff_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp)))
            !!amplitude
            param_x_err_tmp(2, j3) &
            &  = ((dot_product(waveform_errest_diff_tmp, waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_dx,       waveform_errest_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp)  &
            &  *   dot_product(waveform_errest_dx,       waveform_errest_diff_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,      waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_diff_tmp, waveform_errest_diff_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp)))
            param_y_err_tmp(2, j3) &
            &  = ((dot_product(waveform_errest_diff_tmp, waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_dy,       waveform_errest_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_dy,       waveform_errest_diff_tmp))) &
            &  / ((dot_product(waveform_errest_tmp,      waveform_errest_tmp) &
            &  *   dot_product(waveform_errest_diff_tmp, waveform_errest_diff_tmp)) &
            &  -  (dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp) &
            &  *   dot_product(waveform_errest_tmp,      waveform_errest_diff_tmp)))
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
            sigma_amp_term_x(ii, jj) = sigma_amp_term_x(ii, jj) + (param_x_err(2) - param_x_err_tmp(2, j3)) ** 2
            sigma_amp_term_y(ii, jj) = sigma_amp_term_y(ii, jj) + (param_y_err(2) - param_y_err_tmp(2, j3)) ** 2
          enddo
          sigma_slowness_x(ii, jj) = sqrt(sigma_slowness_x(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
          sigma_slowness_y(ii, jj) = sqrt(sigma_slowness_y(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
          sigma_amp_term_x(ii, jj) = sqrt(sigma_amp_term_x(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
          sigma_amp_term_y(ii, jj) = sqrt(sigma_amp_term_y(ii, jj) &
          &                             * real(ntime_slowness - 1, kind = fp) / real(ntime_slowness, kind = fp))
        else
          slowness_x(ii, jj) = 0.0_fp
          slowness_y(ii, jj) = 0.0_fp
          amp_term_x(ii, jj) = 0.0_fp
          amp_term_y(ii, jj) = 0.0_fp
        endif
        write(10, rec = icount) real(x_start + dgrid_x * real(ii - 1, kind = fp), kind = sp), &
        &                       real(y_start + dgrid_y * real(jj - 1, kind = fp), kind = sp), &
        &                       real(-slowness_x(ii, jj), kind = sp), real(-slowness_y(ii, jj), kind = sp), &
        &                       real(sigma_slowness_x(ii, jj), kind = sp), real(sigma_slowness_y(ii, jj), kind = sp), & 
        &                       real(amp_term_x(ii, jj), kind = sp), real(amp_term_y(ii, jj), kind = sp), &
        &                       real(sigma_amp_term_x(ii, jj), kind = sp), real(sigma_amp_term_y(ii, jj), kind = sp)
        icount = icount + 1
      enddo
    enddo
    close(10)
  enddo 


  stop
end program seismicgradiometry

