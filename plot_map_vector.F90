program plot_map_vector
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  use particle_filter
  use origintime
  use legend
  implicit none


  integer         :: i, ios, ncoastline, mapcount, narray, ntriangle, narray_use, maxloc_likelihood(1)
  integer         :: year, month, day, hr, mi, sc, julianday, sec_from_day
  real(kind = fp) :: slowness_x, slowness_y, origintime_median
  real(kind = fp), allocatable :: appvel_obs(:), az_obs(:), lon_array(:), lat_array(:), min_correlation(:), arrivaltime(:), &
  &                               origintime_candidate(:)
  integer,         allocatable :: arrayindex(:)
  logical,         allocatable :: result_exist(:)

  integer                           :: color(1 : 3)
  real(kind = sp)                   :: plot_x, plot_y, plot_x1, plot_y1, plot_theta
  real(kind = fp)                   :: width_min, width_max, height_min, height_max, dwidth, dheight, maplon, maplat, &
  &                                    map_x, map_y, map_x1, map_y1
  character(len = 255)              :: coastline_txt, mapbuf_tmp
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 2)                :: yr, mo, dy, hh, mm, ss
  character(len = 19)               :: date_txt
  
  integer, allocatable :: seed(:)
  integer              :: seedsize
  real(kind = fp)      :: lon_particle(1 : nparticle), lat_particle(1 : nparticle), &
  &                       likelihood_particle(1 : nparticle), az_weight(1 : int(2.0_fp * pi / daz_weight))
  real(kind = fp)      :: likelihood_tmp

  call random_seed(size = seedsize)
  allocate(seed(1 : seedsize))
  do i = 1, seedsize
    call system_clock(count = seed(i))
  enddo
  call random_seed(put = seed(:))

  !!Read coastline
  call getarg(1, coastline_txt)
  ncoastline = 0
  open(unit = 10, file = trim(coastline_txt))
  do
    read(10, '(a255)', iostat = ios) mapbuf_tmp
    if(ios .ne. 0) exit
    if(mapbuf_tmp(1 : 1) .eq. "#") cycle
    ncoastline = ncoastline + 1
  enddo
  rewind(10)
  allocate(mapbuf(1 : ncoastline))
  i = 1
  do
    read(10, '(a255)', iostat = ios) mapbuf_tmp
    if(ios .ne. 0) exit
    if(mapbuf_tmp(1 : 1) .eq. "#") cycle
    mapbuf(i) = mapbuf_tmp
    i = i + 1
  enddo
  close(10)    
 
  !!Plot legend
  call pc_plotinit(iwin_legend, "Legend", 0.0_sp, -300.0_sp, width / 2, 27.0_sp, scale)
  call plot_legend

  !!Read and plot AELUMA results
  call pc_plotinit(iwin_map, "AELUMA results", 0.0_sp, 0.0_sp, width, height, scale)
  call pc_setbkcolor(iwin_map, 255, 255, 255)
  call mercator(center_lon, lon_w, lat_s, width_min, height_min)
  call mercator(center_lon, lon_e, lat_n, width_max, height_max)
  dwidth = 1.0_fp / (width_max - width_min)
  dheight = 1.0_fp / (height_max - height_min)
  !!read AELUMA results from stdin
  do 
    call pc_clear(iwin_map)
    call pc_setcolor(iwin_map, 0, 0, 0)
    call pc_setline(iwin_map, 1)
    az_weight(1 : int(2.0_fp * pi / daz_weight)) = 0.0_fp
    read(*, *, iostat = ios) yr, mo, dy, hh, mm, ss, narray, ntriangle
    read(yr, *) year; year = year + 2000
    read(mo, *) month
    read(dy, *) day
    read(hh, *) hr
    read(mm, *) mi
    read(ss, *) sc
    call ymd2jday(year, month, day, julianday)
    sec_from_day = hr * 60 * 60 + mi * 60 + sc
    !print *, '(8(i0, 1x))', yr, mo, dy, hh, mm, ss, narray, ntriangle
    print '(5(i0, 1x))', year, julianday, sec_from_day, narray, ntriangle
    if(ios .ne. 0) stop
    if(.not. allocated(arrayindex)) then
      allocate(az_obs(1 : ntriangle), appvel_obs(1 : ntriangle), result_exist(1 : ntriangle), arrayindex(1 : ntriangle), &
      &        lon_array(1 : ntriangle), lat_array(1 : ntriangle), min_correlation(1 : ntriangle), arrivaltime(1 : ntriangle))
    endif

    if(narray .ge. 1) then
      call pc_setline(iwin_map, 4)
      result_exist(1 : ntriangle) = .false.
      !!read and plot slowness vector
      do i = 1, narray
        read(*, *) arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
        &          slowness_x, slowness_y, min_correlation(arrayindex(i)), &
        &          arrivaltime(arrayindex(i))
        !!arrival time: relative time in s from current time
        arrivaltime(arrayindex(i)) = -(real(nsec_buf, kind = fp) - arrivaltime(arrayindex(i)))

        az_obs(arrayindex(i)) = atan2(slowness_x, slowness_y)
        if(az_obs(arrayindex(i)) .lt. 0.0_fp) az_obs(arrayindex(i)) = az_obs(arrayindex(i)) + 2.0_fp * pi
        appvel_obs(arrayindex(i)) = 1.0_fp / sqrt(slowness_x ** 2 + slowness_y ** 2)
        if(min_correlation(arrayindex(i)) .ge. correlation_threshold) then
          az_weight(int(az_obs(arrayindex(i)) / daz_weight) + 1) &
          &  = az_weight(int(az_obs(arrayindex(i)) / daz_weight) + 1) + 1.0_fp
        endif

        result_exist(arrayindex(i)) = .true.
        print '(i0, 6(1x, e15.7))', arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
        &                           az_obs(arrayindex(i)), appvel_obs(arrayindex(i)), min_correlation(arrayindex(i)), &
        &                           arrivaltime(arrayindex(i))
      enddo

      !!estimate location
      if(narray .ge. 5) then
        call particle_filter_search(narray, arrayindex, result_exist, lon_array, lat_array, min_correlation, az_obs, &
        &                           az_weight, lon_particle, lat_particle, likelihood_particle)
         
        !!write particles
        call pc_setline(iwin_map, 1)
        do i = 1, nparticle
          call mercator(center_lon, lon_particle(i), lat_particle(i), map_x, map_y)
          plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
          plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
          likelihood_tmp = likelihood_particle(i) * likelihood_legend_normalize
          if(likelihood_tmp .le. 0.1_fp) then
            color(1 : 3) = [252, 238, 158]
          elseif(likelihood_tmp .gt. 0.1_fp .and. likelihood_tmp .le. 0.3_fp) then
            color(1 : 3) = [238, 179, 87]
          elseif(likelihood_tmp .gt. 0.3_fp .and. likelihood_tmp .le. 1_fp) then
            color(1 : 3) = [222, 117, 79]
          elseif(likelihood_tmp .gt. 1_fp .and. likelihood_tmp .le. 3.0_fp) then
            color(1 : 3) = [149, 66, 62]
          elseif(likelihood_tmp .gt. 3.0_fp .and. likelihood_tmp .le. 10.0_fp) then
            color(1 : 3) = [63, 39, 23]
          elseif(likelihood_tmp .gt. 10.0_fp) then
            color(1 : 3) = [26, 26, 1]
          endif
          call pc_setcolor(iwin_map, color(1), color(2), color(3))
          call pc_symbol(iwin_map, plot_x, plot_y, 3.0_sp, 1, 0)
          call pc_setcolor(iwin_map, 0, 0, 0)
          call pc_symbol(iwin_map, plot_x, plot_y, 3.0_sp, 1, 1)
        enddo
        maxloc_likelihood = maxloc(likelihood_particle)
        call mercator(center_lon, lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), map_x, map_y)
        plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
        plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
        call pc_symbol(iwin_map, plot_x, plot_y, 9.0_sp, 1, 0)
        call pc_setcolor(iwin_map, 255, 255, 255)
        call pc_symbol(iwin_map, plot_x, plot_y, 4.0_sp, 1, 0)

        !!estimate origintime
        if(.not. allocated(origintime_candidate)) allocate(origintime_candidate(1 : ntriangle))
        call ot_search_median(narray, arrayindex, result_exist, lon_array, lat_array, min_correlation, az_obs, &
        &                     appvel_obs, arrivaltime, origintime_candidate, &
        &                     lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), origintime_median, &
        &                     narray_use)
        if(origintime_median .lt. 0.0_fp) then
          origintime_median = origintime_median + 86400.0_fp
          julianday = julianday - 1
        endif
        if(julianday .eq. 0) then
          year = year - 1
          month = 12
          day = 31
        else
          call jday2ymd(julianday, year, month, day)
        endif
        hr = int(origintime_median / 3600.0_fp)
        mi = int((origintime_median - 3600.0_fp * real(hr, kind = fp)) / 60.0_fp)
        sc = int(origintime_median - 3600.0_fp * real(hr, kind = fp) - 60.0_fp * real(mi, kind = fp))
        !write(0, '(4(i0, 1x), f8.4)') year, julianday, sec_from_day, narray_use, origintime_median
        !write(0, '(6(i0, 1x))') year, month, day, hr, mi, sc
 
      endif

    endif

    !!plot slowness vector
    call pc_setline(iwin_map, 4)
    do i = 1, narray
      if(.not. result_exist(arrayindex(i))) cycle
      !theta = atan2(slowness_x, slowness_y) * rad2deg
      plot_theta = 90.0_sp - real(az_obs(arrayindex(i)) * rad2deg, kind = sp)
      call mercator(center_lon, lon_array(arrayindex(i)), lat_array(arrayindex(i)), map_x, map_y)
      plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
      if(min_correlation(arrayindex(i)) .lt. 0.01_fp) then
        color(1 : 3) = [220, 204, 222]
      elseif(min_correlation(arrayindex(i)) .ge. 0.2_fp .and. min_correlation(arrayindex(i)) .lt. 0.4_fp) then
        color(1 : 3) = [212, 156, 189]
      elseif(min_correlation(arrayindex(i)) .ge. 0.4_fp .and. min_correlation(arrayindex(i)) .lt. 0.6_fp) then
        color(1 : 3) = [196, 110, 155]
      elseif(min_correlation(arrayindex(i)) .ge. 0.6_fp .and. min_correlation(arrayindex(i)) .lt. 0.8_fp) then
        color(1 : 3) = [136, 97, 141]
      elseif(min_correlation(arrayindex(i)) .ge. 0.8_fp) then
        color(1 : 3) = [73, 57, 100]
      endif
      call pc_setcolor(iwin_map, color(1), color(2), color(3))
      call pc_vector(iwin_map, plot_x, plot_y, plot_theta, vector_len, vector_width, vector_head1, vector_head2, 1)
    enddo
    call pc_setline(iwin_map, 1)

 
    !!plot map
    call pc_setcolor(iwin_map, 0, 0, 0)
    date_txt = "20" // yr // "/" // mo // "/" // dy // " " // hh // ":" // mm // ":" // ss
    call mercator(center_lon, lon_w, lat_n, map_x, map_y)
    plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
    call pc_text(iwin_map, plot_x, plot_y, 7.0, date_txt, 0.0, len(date_txt), 7)

    mapcount = 0
    do i = 1, ncoastline
      if(mapbuf(i)(1 : 1) .eq. ">") then
        mapcount = 0
        cycle
      endif
      read(mapbuf(i), *) maplon, maplat
      call mercator(center_lon, maplon, maplat, map_x, map_y)
    
      if(mapcount .eq. 0) then
        mapcount = 1
        map_x1 = map_x
        map_y1 = map_y
        cycle
      endif
      plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
      plot_x1 = real((map_x1 - width_min)  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
      plot_y1 = real((map_y1 - height_min) * dheight, kind = sp) * height
      map_x1 = map_x
      map_y1 = map_y
      call pc_line(iwin_map, plot_x, plot_y, plot_x1, plot_y1)
    enddo
    call pc_flush(iwin_map)
  enddo


  call pc_plotend(iwin_map, 1)
  call pc_plotend(iwin_legend, 1)

  stop
end program plot_map_vector

subroutine mercator(center_lon, lon, lat, x_east, y_north)
  use nrtype, only : fp
  use constants, only : r_earth, deg2rad, pi
  implicit none
  real(kind = fp), intent(in)  :: center_lon, lon, lat
  real(kind = fp), intent(out) :: x_east, y_north

  x_east  = r_earth * deg2rad * (lon - center_lon)
  y_north = r_earth * log(tan(pi * 0.25_fp + lat * deg2rad * 0.5_fp))

  return
end subroutine mercator
   
