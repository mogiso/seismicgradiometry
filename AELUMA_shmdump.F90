program AELUMA_shmdump
  use nrtype, only : fp, sp, dp
  use tandem
  use typedef
  use calc_kernelmatrix, only : calc_slowness_est_matrix_delaunay_shmdump
  use aeluma_parameters
  use lonlat_xy_conv
  use taper
  use correlation

  implicit none
  !!plot using ixp library
  real(kind = sp), parameter   :: wavewindow_width = 350.0_sp, wavewindow_height = 300.0_sp, wavescale = 1.0_sp, &
  &                               wavewidth = 330.0_sp, waveheight = 290.0_sp, waveplotscale = 3.0_sp
  real(kind = sp)              :: plot_x0, plot_y0, plot_x1, plot_y1, plot_yref, dheight, dwidth
  real(kind = fp)              :: maxamp
  character(len = 8)           :: date_c, time_c

  type(location)               :: location_sta(1 : nwinch)
  integer                      :: i, j, k, ii, jj, ios, nsample, nch, nstation, winch_tmp, ndecimate, ntriangle, npair_tmp, &
  &                               recflag, transdelay, mon_order, ad_bit, sensor_amp, narray_success, stack_index, ntimefft
  integer                      :: waveform_tmp (1 : maxval(sampling_int)), max_xcorr(1), maxloc_stack(1), minloc_stack(1)
  integer, allocatable         :: station_winch(:), xcorr_index(:, :)
  character(len = 2)           :: yr(1 : nsec_buf), mo(1 : nsec_buf), dy(1 : nsec_buf), &
  &                               hh(1 : nsec_buf), mm(1 : nsec_buf), ss(1 : nsec_buf)
  real(kind = fp)              :: ad_v_min, sensor_sens, naturalfreq, damp, stlat_tmp, stlon_tmp, stelev_tmp, &
  &                               ptime_cor, stime_cor, sum_abslagtime, sum_lagtime, correlation_checkval
  real(kind = fp)              :: waveform_real(1 : maxval(sampling_int)), station_sensitivity(1 : nwinch)
  real(kind = dp)              :: waveform_fft(1 : ntime_fft, 1 : 2), xcorr(-ntime_fft2 + 1 : ntime_fft2)
  real(kind = fp), allocatable :: slowness_matrix(:, :, :), slowness(:, :), lagtime(:), minval_xcorr(:), waveform_buf(:, :), &
  &                               arrivaltime(:), taper_window(:), waveform_stacked(:)
  character(len = 4)           :: winch_char, comp_tmp
  character(len = 6)           :: stname(1 : nwinch), stname_tmp
  character(len = 255)         :: chtbl, chtbl_line, sensor_unit
  logical                      :: is_usewinch(1 : nwinch)
  logical,         allocatable :: xcorr_flag(:)

  !!delaunay triangulation
  integer,         allocatable :: triangle_stationwinch(:, :), nsta_count(:), tnbr(:, :)
  type(location), allocatable  :: triangle_center(:)

  !!band-pass filter
  integer,         parameter   :: filter_mode = 1
  integer                      :: m(1 : nsampling_int), n(1 : nsampling_int)
  real(kind = fp), allocatable :: h(:, :), uv(:, :)
  real(kind = fp)              :: gn(1 : nsampling_int), c(1 : nsampling_int)

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
  ntimefft = min(nsec_for_fft * sampling_int_use, ntime_fft)
  allocate(taper_window(1 : ntimefft), waveform_stacked(1 : ntimefft))
  call cosine_taper(cos_taper_ratio, ntimefft, taper_window)

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
    !write(0, '(i0, a, 1x, z4, 1x, a)') i, " used station winch and name = ", station_winch(i), trim(stname(station_winch(i)))
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
  allocate(arrivaltime(1 : ntriangle), xcorr_index(1 : maxval(nsta_count), 1 : maxval(nsta_count)))

  !!open waveform canvas
  call pc_plotinit(iwin_wave, "Waveform monitor", 0.0, 0.0, wavewindow_width, wavewindow_height, wavescale)
  call pc_setbkcolor(iwin_wave, 255, 255, 255)
  call pc_setcolor(iwin_wave, 0, 0, 0)
  dheight = (2.0_sp * waveheight - wavewindow_height) / real(nstation + 1, kind = sp)
  dwidth  = (2.0_sp * wavewidth  - wavewindow_width)  / real(waveform_buf_index_max, kind = sp)

  !!read waveforms from standard input, then conduct analysis
  yr(1 : nsec_buf) = "00"
  mo(1 : nsec_buf) = "00"
  dy(1 : nsec_buf) = "00"
  hh(1 : nsec_buf) = "00"
  mm(1 : nsec_buf) = "00"
  ss(1 : nsec_buf) = "00"
  time_loop: do
    !!read date from stdin
    yr(1 : nsec_buf - 1) = yr(2 : nsec_buf)
    mo(1 : nsec_buf - 1) = mo(2 : nsec_buf)
    dy(1 : nsec_buf - 1) = dy(2 : nsec_buf)
    hh(1 : nsec_buf - 1) = hh(2 : nsec_buf)
    mm(1 : nsec_buf - 1) = mm(2 : nsec_buf)
    ss(1 : nsec_buf - 1) = ss(2 : nsec_buf)
    read(*, *, iostat = ios) yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf), nch
    if(ios .ne. 0) error stop
    write(0, '(a, 6(a2, 1x))') "Reading ", &
    &                             yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf)
    write(0, '(a, i0)') "nch = ", nch

    if(.not. allocated(waveform_buf)) then
      allocate(waveform_buf(1 : waveform_buf_index_max, 1 : nwinch))
      waveform_buf(1 : waveform_buf_index_max, 1 : nwinch) = 0.0_fp
    endif
    if(.not. allocated(uv)) then
      allocate(uv(1 : 4 * maxval(m), 1 : nwinch))
      uv(1 : 4 * maxval(m), 1 : nwinch) = 0.0_fp
    endif

    !!move previous datum
    waveform_buf(1 : waveform_buf_index_max - sampling_int_use, 1 : nwinch) &
    &  = waveform_buf(sampling_int_use + 1 : waveform_buf_index_max, 1 : nwinch)

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
        call tandem3(waveform_real(1 : nsample), h(:, i), gn(i), filter_mode, past_uv = uv(:, winch_tmp))
      endif
      ndecimate = sampling_int(i) / sampling_int_use
      do i = 1, nsample / ndecimate
        waveform_buf(waveform_buf_index_max - sampling_int_use + i, winch_tmp) = waveform_real(ndecimate * (i - 1) + 1)
      enddo
    enddo

    call pc_clear(iwin_wave)
    call pc_line(iwin_wave, wavewindow_width - wavewidth, waveheight, wavewidth, waveheight)
    call pc_line(iwin_wave, wavewidth, waveheight, wavewidth, wavewindow_height - waveheight)
    call pc_line(iwin_wave, wavewidth, wavewindow_height - waveheight, &
    &                                  wavewindow_width  - wavewidth,  &
    &                                  wavewindow_height - waveheight)
    call pc_setdash(iwin_wave, 3)
    !!draw time index line
    do i = 1, int(nsec_buf / 60) + 1
      plot_x0 = (wavewindow_width - wavewidth) + 60.0_sp * real(sampling_int_use, kind = sp) * dwidth * real(i - 1, kind = sp)
      plot_x1 = plot_x0
      plot_y0 = wavewindow_height - waveheight
      plot_y1 = waveheight
      call pc_line(iwin_wave, plot_x0, plot_y0, plot_x1, plot_y1)
      !!draw date and time
      plot_y0 = plot_y0 - 2.0_sp
      plot_y1 = plot_y0 - 5.3_sp
      if(i .eq. int(nsec_buf / 60) + 1) then
        date_c = yr(nsec_buf) // "/" // mo(nsec_buf) // "/" // dy(nsec_buf)
        time_c = hh(nsec_buf) // ":" // mm(nsec_buf) // ":" // ss(nsec_buf)
        call pc_text(iwin_wave, plot_x0, plot_y0, 4.0_sp, date_c, 0.0, len(date_c), 5)
        call pc_text(iwin_wave, plot_x0, plot_y1, 4.0_sp, time_c, 0.0, len(time_c), 5)
      else
        date_c = yr(60 * (i - 1) + 1) // "/" // mo(60 * (i - 1) + 1) // "/" // dy(60 * (i - 1) + 1)
        time_c = hh(60 * (i - 1) + 1) // ":" // mm(60 * (i - 1) + 1) // ":" // ss(60 * (i - 1) + 1)
        call pc_text(iwin_wave, plot_x0, plot_y0, 4.0_sp, date_c, 0.0, len(date_c), 4)
        call pc_text(iwin_wave, plot_x0, plot_y1, 4.0_sp, time_c, 0.0, len(time_c), 4)
      endif
    enddo
    do j = 1, nstation
      maxamp = maxval(waveform_buf(:, station_winch(j)))
      if(maxamp .eq. 0.0_fp) maxamp = 1.0_fp
      plot_x0 = wavewindow_width - wavewidth
      plot_yref = waveheight - dheight * real(j, kind = sp)
      plot_y0 = plot_yref + real(waveform_buf(1, station_winch(j)) / maxamp, kind = sp) * waveplotscale
      do i = 1, waveform_buf_index_max - 1
        plot_y1 = plot_yref + real(waveform_buf(i + 1, station_winch(j)) / maxamp, kind = sp) * waveplotscale
        plot_x1 = plot_x0 + dwidth
        call pc_line(iwin_wave, plot_x0, plot_y0, plot_x1, plot_y1)
        plot_x0 = plot_x1
        plot_y0 = plot_y1
      enddo
      !!plot stationname
      plot_x0 = wavewidth + 0.8_sp
      stname_tmp = trim(stname(station_winch(j)))
      call pc_text(iwin_wave, plot_x0, plot_yref, 3.75_sp, stname_tmp, 0.0_sp, len(stname_tmp), 4)
    enddo
    call pc_flush(iwin_wave)

    !!calculate correlation between two stations within the array
    narray_success = 0
    do j = 1, ntriangle
      xcorr_flag(j) = .true.
      if(nsta_count(j) .eq. 0) then
        xcorr_flag(j) = .false.
        cycle
      endif
      npair_tmp = 0
      do i = 1, nsta_count(j) - 1
        npair_tmp = npair_tmp + i
      enddo
      allocate(lagtime(1 : npair_tmp))
      lagtime(1 : npair_tmp) = 0.0_fp

      i = 0
      minval_xcorr(j) = 0.0_fp
      do jj = 1, nsta_count(j) - 1
        waveform_fft(1 : ntime_fft, 1) = 0.0_dp
        do k = 1, ntimefft
          waveform_fft(k, 1) = real(waveform_buf(nsec_buf * sampling_int_use - ntimefft + k, &
          &                                      triangle_stationwinch(jj, j)) * taper_window(k) * order, kind = dp)
        enddo
        do ii = jj + 1, nsta_count(j)
          i = i + 1
          waveform_fft(1 : ntime_fft, 2) = 0.0_fp
          do k = 1, ntimefft
            waveform_fft(k, 2) = real(waveform_buf(nsec_buf * sampling_int_use - ntimefft + k, &
            &                                      triangle_stationwinch(ii, j)) * taper_window(k) * order, kind = dp)
          enddo
          xcorr(-ntime_fft2 + 1 : ntime_fft2) = 0.0_dp
          call correlation_fft(waveform_fft, ntime_fft, xcorr)
          max_xcorr = maxloc(xcorr) - ntime_fft2

          !if(xcorr(max_xcorr(1)) .le. minval_xcorr(j)) minval_xcorr(j) = xcorr(max_xcorr(1))
          minval_xcorr(j) = minval_xcorr(j) + real(xcorr(max_xcorr(1)), kind = fp)
          lagtime(i) = real(max_xcorr(1), kind = fp) * 1.0_fp / real(sampling_int_use, kind = fp)
          xcorr_index(ii, jj) = i
        enddo
      enddo
      minval_xcorr(j) = minval_xcorr(j) / real(i, kind = fp)

      !!check cross-correlation value
      correlation_consistency: do jj = 1, nsta_count(j) - 2
        do ii = jj + 1, nsta_count(j) - 1
          sum_abslagtime = abs(lagtime(xcorr_index(ii, jj))) + abs(lagtime(xcorr_index(ii + 1, jj))) &
          &              + abs(lagtime(xcorr_index(ii + 1, ii)))
          sum_lagtime    = lagtime(xcorr_index(ii, jj)) - lagtime(xcorr_index(ii + 1, jj)) + lagtime(xcorr_index(ii + 1, ii))
          correlation_checkval = 1.0_fp - abs(sum_lagtime) / sum_abslagtime
          if(correlation_checkval .le. lagtime_ratio_threshold) then
            write(0, '(a, i0, 2(a, e15.7))') "cross-correlation consistency error, array num = ", j, &
            &                                " checkvalue = ", correlation_checkval, " minval_xcorr = ", minval_xcorr(j)
            minval_xcorr(j) = xcorr_min
            exit correlation_consistency
          endif
        enddo
      enddo correlation_consistency


      !if(.not. (maxval(lagtime) .le. lagtime_max .and. minval(lagtime) .ge. lagtime_min)) then
      !  xcorr_flag(j) = .false.
      !  deallocate(lagtime)
      !  cycle
      !endif
      !if(minval_xcorr(j) .le. xcorr_min) then
      !  xcorr_flag(j) = .false.
      !  deallocate(lagtime)
      !  cycle
      !endif

      slowness(1 : 2, j) = matmul(slowness_matrix(1 : 2, 1 : npair_tmp, j), lagtime(1 : npair_tmp))
      if(minval_xcorr(j) .le. xcorr_min) then
        xcorr_flag(j) = .false.
        deallocate(lagtime)
        cycle
      endif
      if(slowness(1, j) .eq. 0.0_fp .and. slowness(2, j) .eq. 0.0_fp) then
        xcorr_flag(j) = .false.
        deallocate(lagtime)
        cycle
      endif
      if(slowness(1, j) * slowness(1, j) + slowness(2, j) * slowness(2, j) .gt. max_slowness * max_slowness) then
        xcorr_flag(j) = .false.
        deallocate(lagtime)
        cycle
      endif
      narray_success = narray_success + 1
      waveform_stacked(1 : ntimefft) = 0.0_fp
      do jj = 1, nsta_count(j)
        do ii = 1, ntimefft
          stack_index = &
          &  int((slowness(1, j) * location_sta(triangle_stationwinch(jj, j))%x_east &
          &     + slowness(2, j) * location_sta(triangle_stationwinch(jj, j))%y_north) * real(sampling_int_use, kind = fp)) &
          &     + nsec_buf * sampling_int_use - ntimefft + ii
          if(stack_index .lt. 1 .or. stack_index .gt. nsec_buf * sampling_int_use) cycle
          waveform_stacked(ii) = waveform_stacked(ii) + waveform_buf(stack_index, triangle_stationwinch(jj, j))
        enddo
      enddo
      maxloc_stack = maxloc(waveform_stacked)
      minloc_stack = minloc(waveform_stacked)
      if(abs(waveform_stacked(minloc_stack(1))) .gt. waveform_stacked(maxloc_stack(1))) then
        arrivaltime(j) = real(minloc_stack(1), kind = fp) / real(sampling_int_use, kind = fp)
      else
        arrivaltime(j) = real(maxloc_stack(1), kind = fp) / real(sampling_int_use, kind = fp)
      endif
      !arrivaltime(j) = (real(maxloc_stack(1), kind = fp) + real(minloc_stack(1), kind = fp)) * 0.5_fp &
      !&              / real(sampling_int_use, kind = fp)
      if(abs(arrivaltime(j)) .lt. 1.0_fp / real(sampling_int_use, kind = fp)) xcorr_flag(j) = .false.
      !write(0, '(5(f9.4, 1x))') triangle_center(j)%lon, triangle_center(j)%lat, &
      !&                         slowness(1, j), slowness(2, j), minval_xcorr(j)

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
    print '(6(a2, 1x), 2(i0, 1x))', yr(nsec_buf), mo(nsec_buf), dy(nsec_buf), hh(nsec_buf), mm(nsec_buf), ss(nsec_buf), &
    &                               narray_success, ntriangle
    do i = 1, ntriangle
      if(xcorr_flag(i) .eqv. .true.) then
        print '(i0, 6(1x, f9.4))', &
        &      i, triangle_center(i)%lon, triangle_center(i)%lat, slowness(1, i), slowness(2, i), &
        &         minval_xcorr(i), arrivaltime(i)
        !write(0, '(i0, 6(1x, f9.4))') &
        !&      i, triangle_center(i)%lon, triangle_center(i)%lat, slowness(1, i), slowness(2, i), &
        !&         minval_xcorr(i), arrivaltime(i)
      endif
    enddo


  enddo time_loop

  call pc_plotend(iwin_wave, 1)
  deallocate(taper_window)

  stop
end program AELUMA_shmdump

