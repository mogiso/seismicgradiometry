!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

program AutomatedEventLocationUsingaMeshofArrays
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use typedef
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy, xy2bl
  use gradiometry_parameters
  use calc_kernelmatrix
  use correlation
  use taper
  use tandem, only : tandem3
  use itoa

  implicit none

  integer :: nsta, npts_tmp, timeindex, i, j, k, ii, jj, m, n, ntriangle, nmesh_flag, nmesh_flag_tmp, ncount, npair_tmp
  type(location),  allocatable :: location_sta(:), triangle_center(:)
  real(kind = fp), allocatable :: begin(:)
  real(kind = fp), allocatable :: waveform_obs(:, :),       & !!(1 : ntime, 1 : nsta)
  &                               slowness_est_matrix(:, :, :), &
  &                               slowness(:, :), &
  &                               lagtime(:), &
  &                               minval_xcorr(:)
  real(kind = fp), allocatable :: h(:), uv(:, :)             !!For time-domain recursive filter 
  integer,         allocatable :: triangle_stationindex(:, :), nsta_count(:), tnbr(:, :)
  logical,         allocatable :: xcorr_flag(:)
  real(kind = fp)              :: dt, c, gn, app_velocity, backazimuth, &
  &                               xcorr(1 : ntime_fft), &
  &                               taper_window(1 : ntime_fft), &
  &                               waveform_fft(1 : ntime_fft, 1 : 2)
  integer                      :: max_xcorr(1)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  character(len = 4)   :: ctimeindex


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
    !waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
  enddo
  call cosine_taper(cos_taper_ratio, ntime_fft, taper_window)

  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(1 : 4 * m), uv(1 : 4 * m, 1 : nsta))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  do i = 1, nsta
    call tandem3(waveform_obs(:, i), h, gn, 1, uv(:, i))
  enddo
  deallocate(h, uv)

  !!convert station longitude/latitude to x_east/y_north
  open(unit = 12, file = "station_location.txt")
  do i = 1, nsta
    call bl2xy(location_sta(i)%lon, location_sta(i)%lat, center_lon, center_lat, &
    &          location_sta(i)%y_north, location_sta(i)%x_east)
    location_sta(i)%y_north = location_sta(i)%y_north / 1000.0_fp
    location_sta(i)%x_east  = location_sta(i)%x_east  / 1000.0_fp
    write(12, '(5(e15.7, 1x))') location_sta(i)%x_east, location_sta(i)%y_north, &
    &                           location_sta(i)%lon, location_sta(i)%lat, location_sta(i)%depth
  enddo
  close(12)

  !!make kernel matrix for each triangle
  call calc_slowness_est_matrix_delaunay(location_sta, naddstation_array, ntriangle, &
  &                                      triangle_center, slowness_est_matrix,       &
  &                                      triangle_stationindex, nsta_count, tnbr)

  !!calculate amplitude and its spatial derivatives at each grid
  !open(unit = 30, file = "log")
  allocate(minval_xcorr(1 : ntriangle), slowness(1 : 2, 1 : ntriangle))
  do jj = 1, int(ntime / ntimestep)
    timeindex = ntimestep * (jj - 1) + 1
    if(timeindex - ntime_fft + 1 .lt. 1) cycle
    write(0, '(a, i0, a)') "Time index = ", jj, " Calculate cross-correlation and slowness on each array"

    if(timeindex .lt. 1 .or. timeindex - ntimestep .gt. ntime) cycle

    call int_to_char(jj, 4, ctimeindex)
    outfile = "slowness_aeluma_" // trim(ctimeindex) // ".dat"
    open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 7, status = "replace")
    ncount = 1
   
    !!initial step: estimate slowness vector without reducing velocity
    !!First, estimate spatial gradients
    allocate(xcorr_flag(1 : ntriangle))
    do ii = 1, ntriangle
      minval_xcorr(ii) = 1.0_fp
      slowness(1 : 2, ii) = 0.0_fp
      !!calculate slowness vector using conventional array analysis
      if(nsta_count(ii) .eq. 0) cycle
      !correlation
      npair_tmp = 1
      do j = 1, nsta_count(ii) - 1
        npair_tmp = npair_tmp * j
      enddo 
      allocate(lagtime(1 : npair_tmp))
      k = 1
      do j = 1, nsta_count(ii) - 1
        waveform_fft(1 : ntime_fft, 1) &
        &  = waveform_obs(timeindex - ntime_fft + 1 : timeindex, triangle_stationindex(j, ii)) * taper_window(1 : ntime_fft)
        do i = j + 1, nsta_count(ii)
          waveform_fft(1 : ntime_fft, 2) &
          &  = waveform_obs(timeindex - ntime_fft + 1 : timeindex, triangle_stationindex(i, ii)) * taper_window(1 : ntime_fft)

          call correlation_fft(waveform_fft, ntime_fft, xcorr)
          max_xcorr = maxloc(xcorr)
          if(xcorr(max_xcorr(1)) .le. minval_xcorr(ii)) minval_xcorr(ii) = xcorr(max_xcorr(1))
          lagtime(k) = real(max_xcorr(1) - ntime_fft2, kind = fp) * dt
          !print '(3(i0, 1x), 2(1x, f8.5))', k, grid_stationindex(i, ii, jj), grid_stationindex(j, ii, jj), &
          !&                                 xcorr(max_xcorr(1)), lagtime(k)
          k = k + 1
        enddo
      enddo
      if(maxval(lagtime) .gt. lagtime_max) then
        xcorr_flag(ii) = .false.
        deallocate(lagtime)
        cycle
      endif
      if(minval(lagtime) .lt. lagtime_min) then
        xcorr_flag(ii) = .false.
        deallocate(lagtime)
        cycle
      endif
      if(minval_xcorr(ii) .le. xcorr_min) then
        xcorr_flag(ii) = .false.
        deallocate(lagtime)
        cycle
      endif

      xcorr_flag(ii) = .true.
      slowness(1 : 2, ii) &
      &  = matmul(slowness_est_matrix(1 : 2, 1 : npair_tmp, ii), lagtime(1 : npair_tmp))
      write(10, rec = ncount) real(triangle_center(ii)%x_east,  kind = sp), &
      &                       real(triangle_center(ii)%y_north, kind = sp), &
      &                       real(triangle_center(ii)%lon,     kind = sp), &
      &                       real(triangle_center(ii)%lat,     kind = sp), &
      &                       real(slowness(1, ii),             kind = sp), &
      &                       real(slowness(2, ii),             kind = sp), &
      &                       real(minval_xcorr(ii),            kind = sp)
      ncount = ncount + 1
      deallocate(lagtime)
    enddo

    !do j = 1, ntriangle
    !  nmesh_flag = 0
    !  if(xcorr_flag(j) .eqv. .true.) then
    !    nmesh_flag_tmp = 0
    !    do i = 1, 3
    !      if(xcorr_flag(tnbr(i, j)) .eqv. .true.) nmesh_flag_tmp = nmesh_flag_tmp + 1
    !    enddo
    !    if(nmesh_flag_tmp .gt. 1) nmesh_flag = nmesh_flag + 1
    !  endif
    !enddo
          
    deallocate(xcorr_flag)
    close(10)
  enddo


  stop
end program AutomatedEventLocationUsingaMeshofArrays

