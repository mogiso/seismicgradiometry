program AELUMA_shmdump
  use nrtype, only : fp, sp
  use tandem
  use typedef
  use calc_kernelmatrix, only : calc_slowness_est_matrix_delaunay_shmdump
  use aeluma_parameters
  use lonlat_xy_conv
  use taper
  use correlation
  use itoa

  implicit none
  integer,         parameter :: iwin = 0
  real(kind = sp), parameter :: window_width = 350.0_sp, window_height = 300.0_sp, scale = 1.0_sp, &
  &                             width = 320.0_sp, height = 270.0_sp, plotscale = 3.0_sp
  real(kind = sp)            :: plot_x0, plot_y0, plot_x1, plot_y1, plot_yref, dheight, dwidth
  real(kind = fp)            :: maxamp

  character(len = 4) :: winch_char, comp_tmp
  integer            :: i, j, ii, jj, ios, nsample, nch, nstation, winch_tmp, ndecimate, ntriangle, npair_tmp, &
  &                     recflag, transdelay, mon_order, ad_bit, sensor_amp, narray_success
  integer            :: yr(1 : nsec_buf), mo(1 : nsec_buf), dy(1 : nsec_buf), &
  &                     hh(1 : nsec_buf), mm(1 : nsec_buf), ss(1 : nsec_buf), &
  &                     waveform_tmp (1 : maxval(sampling_int)), max_xcorr(1)
  integer, allocatable :: station_winch(:), winch_index(:)
  real(kind = fp)    :: waveform_real(1 : maxval(sampling_int)), waveform_fft(1 : ntime_fft, 1 : 2), &
  &                     taper_window(1 : sampling_int_use * nsec_buf), xcorr(-ntime_fft2 + 1 : ntime_fft2), &
  &                     station_sensitivity(1 : nwinch)
  logical            :: is_usewinch(1 : nwinch)
  character(len = 6) :: stname(1 : nwinch), stname_tmp
  character(len = 255) :: chtbl, chtbl_line, sensor_unit
  real(kind = fp)    :: ad_v_min, sensor_sens, naturalfreq, damp, &
  &                     stlat_tmp, stlon_tmp, stelev_tmp, ptime_cor, stime_cor 
  real(kind = fp), allocatable :: slowness_matrix(:, :, :), slowness(:, :), lagtime(:), minval_xcorr(:), waveform_buf(:, :)
  integer,         allocatable :: triangle_stationwinch(:, :), nsta_count(:), tnbr(:, :)
  logical,         allocatable :: xcorr_flag(:)

  type(location)     :: location_sta(1 : nwinch)
  type(location), allocatable :: triangle_center(:)
  !!band-pass filter
  integer,         parameter :: filter_mode = 1
  real(kind = fp), allocatable :: h(:, :), uv(:, :)
  real(kind = fp)              :: gn(1 : nsampling_int), c(1 : nsampling_int)
  integer                      :: m(1 : nsampling_int), n(1 : nsampling_int)

  call getarg(1, chtbl)

  !!filter design
  do i = 1, nsampling_int
    call calc_bpf_order(fl, fh, fs, ap, as, sampling_sec(i), m(i), n(i), c(i))
  enddo
  allocate(h(1 : 4 * maxval(m), 1 : nsampling_int))
  do i = 1, nsampling_int
    call calc_bpf_coef(fl, fh, sampling_sec(i), m(i), n(i), h(:, i), c(i), gn(i))
  enddo
  !!cosine taper
  call cosine_taper(cos_taper_ratio, nsec_buf * sampling_int_use, taper_window)

  !!read win-formatted channel table
  is_usewinch(1 : nwinch) = .false.
  open(unit = 10, file = trim(chtbl))
  nstation = 0
  do 
    read(10, '(a256)', iostat = ios) chtbl_line
    if(ios .ne. 0) exit
    if(chtbl_line(1 : 1) .eq. "#") cycle
    nstation = nstation + 1
  enddo
  rewind(10)
  allocate(station_winch(1 : nstation))
  i = 1
  do
    read(10, '(a255)', iostat = ios) chtbl_line
    if(ios .ne. 0) exit
    if(chtbl_line(1 : 1) .eq. "#") cycle
    do j = 1, 255
      if(chtbl_line(j : j) .eq. "/") chtbl_line(j : j) = "-"
    enddo
    read(chtbl_line, *) winch_char, recflag, transdelay, stname_tmp, comp_tmp, mon_order, &
    &                   ad_bit, sensor_sens, sensor_unit, naturalfreq, damp, sensor_amp, ad_v_min, &
    &                   stlat_tmp, stlon_tmp, stelev_tmp, ptime_cor, stime_cor
    read(winch_char, '(z4)') station_winch(i)
    is_usewinch(station_winch(i)) = .true.
    stname(station_winch(i)) = trim(stname_tmp)
    station_sensitivity(station_winch(i)) = ad_v_min / (sensor_sens * 10.0_fp ** (sensor_amp / 20))
    location_sta(station_winch(i))%lat = stlat_tmp
    location_sta(station_winch(i))%lon = stlon_tmp
    write(0, '(i0, a, 1x, z4, 1x, a)') i, " used station winch and name = ", station_winch(i), trim(stname(station_winch(i)))
    i = i + 1
  enddo
  close(10)
  !!convert (lon, lat) to (x, y)
  do i = 1, nstation
    call bl2xy(location_sta(station_winch(i))%lon,     location_sta(station_winch(i))%lat,    &
    &          center_lon_aeluma,                      center_lat_aeluma,                     &
    &          location_sta(station_winch(i))%y_north, location_sta(station_winch(i))%x_east)
    location_sta(station_winch(i))%y_north = location_sta(station_winch(i))%y_north * 1.0e-3_fp
    location_sta(station_winch(i))%x_east  = location_sta(station_winch(i))%x_east  * 1.0e-3_fp
  enddo
  !!calculate triads
  call calc_slowness_est_matrix_delaunay_shmdump(location_sta, station_winch, nadd_station_aeluma, &
  &                                              ntriangle, triangle_center, slowness_matrix, triangle_stationwinch, &
  &                                              nsta_count, tnbr)
  allocate(minval_xcorr(1 : ntriangle), xcorr_flag(1 : ntriangle), slowness(1 : 2, 1 : ntriangle))

  !!open waveform canvas
  call pc_plotinit(iwin, "Waveform monitor", 0.0, 0.0, window_width, window_height, scale)
  call pc_setbkcolor(iwin, 255, 255, 255)
  call pc_setcolor(iwin, 0, 0, 0)
  dheight = (2.0_sp * height - window_height) / real(nstation + 1, kind = sp)
  dwidth  = (2.0_sp * width - window_width) / real(waveform_buf_index_max, kind = sp)

  !!read waveforms from standard input, then conduct analysis
  time_loop: do
    !!read date from stdin
    read(*, *, iostat = ios) yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf), nch
    if(ios .ne. 0) error stop
    write(0, '(a, 6(i2.2, 1x))') "Reading ", &
    &                             yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf)
    write(0, '(a, i0)') "nch = ", nch

    if(.not. allocated(waveform_buf)) then
      allocate(waveform_buf(1 : waveform_buf_index_max, 1 : nwinch))
      waveform_buf(1 : waveform_buf_index_max, 1 : nwinch) = 0.0_fp
    endif
    if(.not. allocated(winch_index)) then
      allocate(winch_index(1 : nwinch))
      winch_index(1 : nwinch) = 0
    endif
    if(.not. allocated(uv)) then
      allocate(uv(1 : 4 * maxval(m), 1 : nwinch))
      uv(1 : 4 * maxval(m), 1 : nwinch) = 0.0_fp
    endif

    !!move previous datum
    waveform_buf(1 : waveform_buf_index_max - sampling_int_use, 1 : nwinch) &
    &  = waveform_buf(sampling_int_use + 1 : waveform_buf_index_max, 1 : nwinch)
    yr(1 : nsec_buf - 1) = yr(2 : nsec_buf)
    mo(1 : nsec_buf - 1) = mo(2 : nsec_buf)
    dy(1 : nsec_buf - 1) = dy(2 : nsec_buf)
    hh(1 : nsec_buf - 1) = hh(2 : nsec_buf)
    mm(1 : nsec_buf - 1) = mm(2 : nsec_buf)
    ss(1 : nsec_buf - 1) = ss(2 : nsec_buf)

    !!read waveforms from stdin, filtering, decimation
    do j = 1, nch
      read(*, *) winch_char, nsample, (waveform_tmp(i), i = 1, nsample)
      read(winch_char, '(z4)') winch_tmp
      winch_index(winch_tmp) = j
      do i = 1, nsampling_int
        if(nsample .eq. sampling_int(i)) exit
      enddo
      waveform_real(1 : nsample) = real(waveform_tmp(1 : nsample), kind = fp) * station_sensitivity(winch_tmp)
      call tandem3(waveform_real(1 : nsample), h(:, i), gn(i), filter_mode, past_uv = uv(:, winch_tmp))
      ndecimate = sampling_int(i) / sampling_int_use
      do i = 1, nsample / ndecimate
        waveform_buf(waveform_buf_index_max - sampling_int_use + i, winch_tmp) = waveform_real(ndecimate * (i - 1) + 1)
      enddo
    enddo

    call pc_clear(iwin)
    call pc_line(iwin, (window_width - width), (window_height - height), (window_width - width), height)
    call pc_line(iwin, (window_width - width), height, width, height)
    call pc_line(iwin, width, height, width, (window_height - height))
    call pc_line(iwin, (window_width - width), height, (window_width - width), (window_height - height))
    call pc_flush(iwin)
    do j = 1, nstation
      maxamp = maxval(waveform_buf(:, station_winch(j)))
      if(maxamp .eq. 0.0_fp) maxamp = 1.0_fp
      plot_x0 = window_width - width
      plot_yref = height - dheight * real(j, kind = sp)
      plot_y0 = plot_yref + real(waveform_buf(1, station_winch(j)) / maxamp, kind = sp) * plotscale
      do i = 1, waveform_buf_index_max - 1
        plot_y1 = plot_yref + real(waveform_buf(i + 1, station_winch(j)) / maxamp, kind = sp) * plotscale
        plot_x1 = plot_x0 + dwidth
        call pc_line(iwin, plot_x0, plot_y0, plot_x1, plot_y1)
        plot_x0 = plot_x1
        plot_y0 = plot_y1
      enddo
    enddo
    call pc_flush(iwin)

    !!calculate correlation between two stations within the array
    narray_success = 0
    do j = 1, ntriangle
      xcorr_flag(j) = .true.
      if(nsta_count(j) .eq. 0) then
        xcorr_flag(j) = .false.
        cycle
      endif
      minval_xcorr(j) = 1.0e+38_fp
      npair_tmp = 0
      do i = 1, nsta_count(j) - 1
        npair_tmp = npair_tmp + i
      enddo
      allocate(lagtime(1 : npair_tmp))
      lagtime(1 : npair_tmp) = 0.0_fp

      i = 1
      do jj = 1, nsta_count(j) - 1
        waveform_fft(1 : ntime_fft, 1) = 0.0_fp
        waveform_fft(1 : nsec_buf * sampling_int_use, 1) &
        &  = waveform_buf(1 : nsec_buf * sampling_int_use, triangle_stationwinch(jj, j)) &
        &  * taper_window(1 : nsec_buf * sampling_int_use) * order
        do ii = jj + 1, nsta_count(j)
          waveform_fft(1 : ntime_fft, 2) = 0.0_fp
          waveform_fft(1 : nsec_buf * sampling_int_use, 2) &
          &  = waveform_buf(1 : nsec_buf * sampling_int_use, triangle_stationwinch(ii, j)) &
          &  * taper_window(1 : nsec_buf * sampling_int_use) * order
          xcorr(-ntime_fft2 + 1 : ntime_fft2) = 0.0_fp
          call correlation_fft(waveform_fft, ntime_fft, xcorr)
          max_xcorr = maxloc(xcorr) - ntime_fft2

          if(xcorr(max_xcorr(1)) .le. minval_xcorr(j)) minval_xcorr(j) = xcorr(max_xcorr(1))
          lagtime(i) = real(max_xcorr(1), kind = fp) * 1.0_fp / real(sampling_int_use, kind = fp)
          i = i + 1
        enddo
      enddo


      !if(.not. (maxval(lagtime) .le. lagtime_max .and. minval(lagtime) .ge. lagtime_min)) then
      !  xcorr_flag(j) = .false.
      !  deallocate(lagtime)
      !  cycle
      !endif
      if(minval_xcorr(j) .le. xcorr_min) then
        xcorr_flag(j) = .false.
        deallocate(lagtime)
        cycle
      endif
      narray_success = narray_success + 1

      slowness(1 : 2, j) = matmul(slowness_matrix(1 : 2, 1 : npair_tmp, j), lagtime(1 : npair_tmp))
      write(0, '(5(f9.4, 1x))') triangle_center(j)%lon, triangle_center(j)%lat, &
      &                         slowness(1, j), slowness(2, j), minval_xcorr(j)

      !do i = 1, 3
      !  call bl2xy(location_sta(triangle_stationwinch(i, j))%lon,     location_sta(triangle_stationwinch(i, j))%lat, &
      !  &          triangle_center(j)%lon,                            triangle_center(j)%lat, &
      !  &          location_sta(triangle_stationwinch(i, j))%y_north, location_sta(triangle_stationwinch(i, jj))%x_east)
      !  location_sta(triangle_stationwinch(i, j))%x_east = location_sta(triangle_stationwinch(i, j))%x_east * 1.0e-3_fp
      !  location_sta(triangle_stationwinch(i, j))%y_north = location_sta(triangle_stationwinch(i, j))%y_north * 1.0e-3_fp
      !enddo
      !write(0, *) slowness(1, j) &
      !&         * (location_sta(triangle_stationwinch(2, j))%x_east - location_sta(triangle_stationwinch(1, j))%x_east) &
      !&         + slowness(2, j) &
      !&         * (location_sta(triangle_stationwinch(2, j))%y_north - location_sta(triangle_stationwinch(1, j))%y_north), &
      !&           lagtime(1)
      !write(0, *) slowness(1, j) &
      !&         * (location_sta(triangle_stationwinch(3, j))%x_east - location_sta(triangle_stationwinch(1, j))%x_east) &
      !&         + slowness(2, j) &
      !&         * (location_sta(triangle_stationwinch(3, j))%y_north - location_sta(triangle_stationwinch(1, j))%y_north), &
      !&           lagtime(2)
      !write(0, *) slowness(1, j) &
      !&         * (location_sta(triangle_stationwinch(3, j))%x_east - location_sta(triangle_stationwinch(2, j))%x_east) &
      !&         + slowness(2, j) &
      !&         * (location_sta(triangle_stationwinch(3, j))%y_north - location_sta(triangle_stationwinch(2, j))%y_north), &
      !&           lagtime(3)

      deallocate(lagtime)
    enddo
    print '(6(i0.2, 1x), i4)', yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf), narray_success
    do i = 1, ntriangle
      if(xcorr_flag(i) .eqv. .true.) then
        print '(i0, 5(1x, f8.4))', &
        &      i, triangle_center(i)%lon, triangle_center(i)%lat, slowness(1, i), slowness(2, i), minval_xcorr(i)
      endif
    enddo


  enddo time_loop

  call pc_plotend(iwin, 1)

  stop
end program AELUMA_shmdump

