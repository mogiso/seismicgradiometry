module plotmodule
  private
  public :: plot_particle, plot_particle_maxlikelihood, plot_slowness_vector, read_coastline, plot_coastline, &
  &         plot_currentdate, plot_legend, epicenter2char, plot_eplist

  contains

  subroutine plot_particle(iwin_plot, lon_particle, lat_particle, likelihood_particle, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use aeluma_parameters
    use mapprojection
    implicit none
    integer, intent(in)         :: iwin_plot
    real(kind = fp), intent(in) :: lon_particle(:), lat_particle(:), likelihood_particle(:), &
    &                              width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    real(kind = fp) :: map_x, map_y, likelihood_tmp
    real(kind = sp) :: plot_x, plot_y
    integer :: i, colorindex, maxloc_likelihood(1)

    call pc_setline(iwin_plot, 1)
    do i = 1, nparticle
      if(mod(i, 8) .ne. 0) cycle
      call mercator(center_lon, lon_particle(i), lat_particle(i), map_x, map_y)
      plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
      if(.not. (plot_x .gt. 0.0_sp .and. plot_x .lt. width)) cycle
      if(.not. (plot_y .gt. 0.0_sp .and. plot_y .lt. height)) cycle
      likelihood_tmp = likelihood_particle(i) * likelihood_legend_normalize
      if(likelihood_tmp .le. 0.1_fp) then
        colorindex = 1
      elseif(likelihood_tmp .gt. 0.1_fp .and. likelihood_tmp .le. 0.2_fp) then
        colorindex = 2
      elseif(likelihood_tmp .gt. 0.2_fp .and. likelihood_tmp .le. 0.5_fp) then
        colorindex = 3
      elseif(likelihood_tmp .gt. 0.5_fp .and. likelihood_tmp .le. 1.0_fp) then
        colorindex = 4
      elseif(likelihood_tmp .gt. 1.0_fp .and. likelihood_tmp .le. 2.0_fp) then
        colorindex = 5
      elseif(likelihood_tmp .gt. 2.0_fp .and. likelihood_tmp .le. 5.0_fp) then
        colorindex = 6
      elseif(likelihood_tmp .gt. 5.0_fp .and. likelihood_tmp .le. 10.0_fp) then
        colorindex = 7
      elseif(likelihood_tmp .gt. 10.0_fp .and. likelihood_tmp .le. 20.0_fp) then
        colorindex = 8
      elseif(likelihood_tmp .gt. 20.0_fp .and. likelihood_tmp .le. 50.0_fp) then
        colorindex = 9
      elseif(likelihood_tmp .gt. 50.0_fp) then
        colorindex = 10
      endif
      call pc_setcolor(iwin_plot, &
      &                color_likelihood(1, colorindex), color_likelihood(2, colorindex), color_likelihood(3, colorindex))
      call pc_symbol(iwin_plot, plot_x, plot_y, 3.0_sp, 1, 0)
      call pc_setcolor(iwin_plot, 0, 0, 0)
      call pc_symbol(iwin_plot, plot_x, plot_y, 3.0_sp, 1, 1)
    enddo
    maxloc_likelihood = maxloc(likelihood_particle)
    call mercator(center_lon, lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), map_x, map_y)
    plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
    if((plot_x .gt. 0.0_sp .and. plot_x .lt. width) .and. (plot_y .gt. 0.0_sp .and. plot_y .lt. height)) then 
      call pc_symbol(iwin_map, plot_x, plot_y, 8.0_sp, 1, 0)
      call pc_setcolor(iwin_map, 255, 255, 255)
      call pc_symbol(iwin_map, plot_x, plot_y, 3.5_sp, 1, 0)
    endif

    write(0, '(a, f0.4, a, f0.4, a, e15.7)') " Lon = ", lon_particle(maxloc_likelihood(1)), &
    &                                        " Lat = ", lat_particle(maxloc_likelihood(1)), &
    &                                        " likelihood = ", likelihood_particle(maxloc_likelihood(1))

    return
  end subroutine plot_particle

  subroutine plot_particle_maxlikelihood(iwin_plot, lon_particle, lat_particle, likelihood_particle, &
  &                                      width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use aeluma_parameters
    use mapprojection
    implicit none
    integer,         intent(in) :: iwin_plot
    real(kind = fp), intent(in) :: lon_particle(:), lat_particle(:), likelihood_particle(:), &
    &                              width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    real(kind = fp) :: map_x, map_y
    real(kind = sp) :: plot_x, plot_y
    integer :: maxloc_likelihood(1)

    maxloc_likelihood = maxloc(likelihood_particle)
    call mercator(center_lon, lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), map_x, map_y)
    plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
    call pc_setcolor(iwin_map, 0, 0, 0)
    !call pc_setline(iwin_map, 9)
    call pc_symbol(iwin_plot, plot_x, plot_y, 8.0_sp, 1, 0)
    call pc_setcolor(iwin_plot, 255, 255, 255)
    call pc_symbol(iwin_plot, plot_x, plot_y, 3.5_sp, 1, 0)

    write(0, '(a, f0.4, a, f0.4, a, e15.7)') " Lon = ", lon_particle(maxloc_likelihood(1)), &
    &                                        " Lat = ", lat_particle(maxloc_likelihood(1)), &
    &                                        " likelihood = ", likelihood_particle(maxloc_likelihood(1))

  end subroutine plot_particle_maxlikelihood

  subroutine plot_slowness_vector(iwin_plot, narray, arrayindex, result_exist, lon_array, lat_array, az_obs, appvel_obs, &
  &                               min_correlation, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use constants, only : rad2deg
    use aeluma_parameters
    use mapprojection

    implicit none
    integer,         intent(in) :: iwin_plot, narray, arrayindex(1 : narray)
    logical,         intent(in) :: result_exist(:)
    real(kind = fp), intent(in) :: lon_array(:), lat_array(:), az_obs(:), appvel_obs(:), min_correlation(:), &
    &                              width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight

    integer :: i, colorindex
    real(kind = sp) :: plot_theta, plot_x, plot_y, vector_len
    real(kind = fp) :: map_x, map_y

    call pc_setline(iwin_plot, 4)
    do i = 1, narray
      if(.not. result_exist(arrayindex(i))) cycle
      plot_theta = 90.0_sp - real(az_obs(arrayindex(i)) * rad2deg, kind = sp)
      vector_len = vector_len_ref * real(appvel_obs(arrayindex(i)) / ref_appvelocity, kind = sp)
      if(vector_len .gt. 20.0_sp) vector_len = 20.0_sp
      call mercator(center_lon, lon_array(arrayindex(i)), lat_array(arrayindex(i)), map_x, map_y)
      plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
      if(min_correlation(arrayindex(i)) .lt. 0.1_fp) then
        colorindex = 1
      elseif(min_correlation(arrayindex(i)) .ge. 0.1_fp .and. min_correlation(arrayindex(i)) .lt. 0.2_fp) then
        colorindex = 2
      elseif(min_correlation(arrayindex(i)) .ge. 0.2_fp .and. min_correlation(arrayindex(i)) .lt. 0.3_fp) then
        colorindex = 3 
      elseif(min_correlation(arrayindex(i)) .ge. 0.3_fp .and. min_correlation(arrayindex(i)) .lt. 0.4_fp) then
        colorindex = 4 
      elseif(min_correlation(arrayindex(i)) .ge. 0.4_fp .and. min_correlation(arrayindex(i)) .lt. 0.5_fp) then
        colorindex = 5
      elseif(min_correlation(arrayindex(i)) .ge. 0.5_fp .and. min_correlation(arrayindex(i)) .lt. 0.6_fp) then
        colorindex = 6
      elseif(min_correlation(arrayindex(i)) .ge. 0.6_fp .and. min_correlation(arrayindex(i)) .lt. 0.7_fp) then
        colorindex = 7
      elseif(min_correlation(arrayindex(i)) .ge. 0.7_fp .and. min_correlation(arrayindex(i)) .lt. 0.8_fp) then
        colorindex = 8
      elseif(min_correlation(arrayindex(i)) .ge. 0.8_fp .and. min_correlation(arrayindex(i)) .lt. 0.9_fp) then
        colorindex = 9
      elseif(min_correlation(arrayindex(i)) .ge. 0.9_fp) then
        colorindex = 10
      endif
      call pc_setcolor(iwin_plot, &
      &                color_correlation(1, colorindex), color_correlation(2, colorindex), color_correlation(3, colorindex))
      call pc_vector(iwin_plot, plot_x, plot_y, plot_theta, vector_len, vector_width, vector_head1, vector_head2, 1)
    enddo
  end subroutine plot_slowness_vector

  subroutine read_coastline(coastline_txt, mapbuf)
    implicit none
    character(len = *),              intent(in)  :: coastline_txt
    character(len = *), allocatable, intent(out) :: mapbuf(:)
    character(len = 255)                         :: mapbuf_tmp
    integer :: i, ncoastline, ios
    
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
      read(10, '(a255)', iostat = ios) mapbuf
      if(ios .ne. 0) exit
      if(mapbuf_tmp(1 : 1) .eq. "#") cycle
      mapbuf(i) = mapbuf_tmp
      i = i + 1
    enddo
    close(10)
    return
  end subroutine read_coastline

  subroutine plot_coastline(iwin_plot, ncoastline, mapbuf, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use mapprojection
    use aeluma_parameters
    implicit none

    real(kind = fp),    intent(in) :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    integer,            intent(in) :: iwin_plot, ncoastline
    character(len = *), intent(in) :: mapbuf(1 : ncoastline)
    integer                        :: i, mapcount
    real(kind = fp)                :: map_x, map_y, map_x1, map_y1, maplon, maplat
    real(kind = sp)                :: plot_x, plot_x1, plot_y, plot_y1

    mapcount = 0
    call pc_setcolor(iwin_map, 0, 0, 0)
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
      plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
      plot_x1 = real((map_x1 - width_tmp(1))  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
      plot_y1 = real((map_y1 - height_tmp(1)) * dheight, kind = sp) * height
      map_x1 = map_x
      map_y1 = map_y
      call pc_line(iwin_plot, plot_x, plot_y, plot_x1, plot_y1)
    enddo
  end subroutine plot_coastline

  subroutine plot_currentdate(iwin_plot, yr, mo, dy, hh, mm, ss, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : sp, fp
    use aeluma_parameters
    use mapprojection
    implicit none
    integer,            intent(in) :: iwin_plot
    character(len = 2), intent(in) :: yr, mo, dy, hh, mm, ss
    real(kind = fp),    intent(in) :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    real(kind = sp)                :: plot_x, plot_y
    real(kind = fp)                :: map_x, map_y
    character(len = 19)            :: date_txt

    date_txt = "20" // yr // "/" // mo // "/" // dy // " " // hh // ":" // mm // ":" // ss
    call mercator(center_lon, lon_w, lat_n, map_x, map_y)
    plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
    call pc_rect(iwin_plot, plot_x, plot_y - 10.0_sp, plot_x + 80.0_sp, plot_y - 10.0_sp)
    call pc_setcolor(iwin_plot, 255, 255, 255)
    call pc_text(iwin_plot, plot_x, plot_y, 7.0, date_txt, 0.0, len(date_txt), 7)

    return
  end subroutine plot_currentdate
 

  subroutine plot_legend(iwin_plot)
    use nrtype, only : sp
    use aeluma_parameters

    implicit none
    integer, intent(in) :: iwin_plot
    integer             :: i
    real(kind = sp)     :: plot_x, plot_y, vector_len
    character(len = 6)  :: plottext

    call pc_setbkcolor(iwin_legend, 255, 255, 255)
    vector_len = vector_len_ref * 10.0_sp
    do i = 1, 10
      if(i .eq. 1) then
        plottext = "0.1  "
      elseif(i .eq. 2) then
        plottext = "0.2   "
      elseif(i .eq. 3) then
        plottext = "0.3   "
      elseif(i .eq. 4) then
        plottext = "0.4   "
      elseif(i .eq. 5) then
        plottext = "0.5   "
      elseif(i .eq. 6) then
        plottext = "0.6   "
      elseif(i .eq. 7) then
        plottext = "0.7   "
      elseif(i .eq. 8) then
        plottext = "0.8   "
      elseif(i .eq. 9) then
        plottext = "0.9   "
      elseif(i .eq. 10) then
        plottext = "1.0   "
      endif
      plot_x = real(i, kind = sp) * 14.0_sp - 10.0_sp
      plot_y = 20.0_sp
      call pc_setcolor(iwin_plot, color_correlation(1, i), color_correlation(2, i), color_correlation(3, i))
      call pc_setline(iwin_plot, 4)
      call pc_vector(iwin_plot, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
      call pc_setcolor(iwin_plot, 0, 0, 0)
      plot_x = plot_x + 6.0_sp
      if(i .ne. 10) call pc_text(iwin_plot, plot_x, plot_y, 4.0_sp, trim(plottext), 0.0_sp, len(trim(plottext)), 4) 
    enddo

    plot_x = 4.0_sp
    plot_y = 13.0_sp
    call pc_setcolor(iwin_plot, 0, 0, 0)
    vector_len = vector_len_ref * 10.0_sp
    plottext = "3 km/s"
    call pc_vector(iwin_plot, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
    plot_x = plot_x + vector_len + 2.0_sp
    call pc_text(iwin_plot, plot_x, plot_y, 4.0_sp, trim(plottext), 0.0_sp, len(trim(plottext)), 4) 

    plot_x = plot_x + 20.0_sp
    vector_len = vector_len_ref * 20.0_sp
    plottext = "6 km/s"
    call pc_setline(iwin_plot, 4)
    call pc_vector(iwin_plot, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
    plot_x = plot_x + vector_len + 2.0_sp
    call pc_text(iwin_plot, plot_x, plot_y, 4.0_sp, trim(plottext), 0.0_sp, len(trim(plottext)), 4) 

    plot_x = plot_x + 20.0_sp
    vector_len = vector_len_ref * 30.0_sp
    plottext = "9 km/s"
    call pc_setline(iwin_plot, 4)
    call pc_vector(iwin_plot, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
    plot_x = plot_x + vector_len + 2.0_sp
    call pc_text(iwin_plot, plot_x, plot_y, 4.0_sp, trim(plottext), 0.0_sp, len(trim(plottext)), 4) 

    plot_x = plot_x + 20.0_sp
    vector_len = vector_len_ref * 40.0_sp
    plottext = "12km/s"
    call pc_setline(iwin_plot, 4)
    call pc_vector(iwin_plot, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
    plot_x = plot_x + vector_len + 2.0_sp
    call pc_text(iwin_plot, plot_x, plot_y, 4.0_sp, trim(plottext), 0.0_sp, len(trim(plottext)), 4) 

    do i = 1, 10 
      if(i .eq. 1) then
        plottext = "0.1"
      elseif(i .eq. 2) then
        plottext = "0.2"
      elseif(i .eq. 3) then
        plottext = "0.5"
      elseif(i .eq. 4) then
        plottext = "1.0"
      elseif(i .eq. 5) then
        plottext = "2.0"
      elseif(i .eq. 6) then
        plottext = "5.0"
      elseif(i .eq. 7) then
        plottext = "10.0"
      elseif(i .eq. 8) then
        plottext = "20.0"
      elseif(i .eq. 9) then
        plottext = "50.0"
      elseif(i .eq. 10) then
        plottext = "5.0"
      endif
      plot_x = real(i, kind = sp) * 14.0_sp - 8.0_sp
      plot_y = 5.0_sp
      call pc_setline(iwin_plot, 1)
      call pc_setcolor(iwin_plot, color_likelihood(1, i), color_likelihood(2, i), color_likelihood(3, i))
      call pc_symbol(iwin_plot, plot_x, plot_y, 3.0_sp, 1, 0)
      call pc_setcolor(iwin_plot, 0, 0, 0)
      call pc_symbol(iwin_plot, plot_x, plot_y, 3.0_sp, 1, 1)
      plot_x = plot_x + 7.0_sp
      if(i .ne. 10) call pc_text(iwin_plot, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 5) 
    enddo
    plot_x = plot_x + 3.0_sp
    call pc_text(iwin_plot, plot_x, plot_y, 4.0, trim(likelihood_legend_normalize_c), 0.0, &
    &            len(trim(likelihood_legend_normalize_c)), 5) 
    call pc_flush(iwin_plot)

    return
  end subroutine plot_legend

  subroutine epicenter2char(year, julianday, sec_from_day, lon_particle, lat_particle, origintime, likelihood_particle, &
  &                         epicenter_info, sigma_lon, sigma_lat, sigma_ot, maxval_likelihood, &
  &                         lon_array, lat_array, az_obs, appvel_obs, array_used, min_correlation, array_maxamp, array_lta)
    use nrtype, only : fp, sp
    use aeluma_parameters
    use jday
    implicit none

    integer,            intent(in)            :: year, julianday, sec_from_day
    real(kind = fp),    intent(in)            :: lon_particle(:), lat_particle(:), origintime(:), likelihood_particle(:)
    character(len = *), intent(out)           :: epicenter_info
    real(kind = fp),    intent(out), optional :: sigma_lon(1 : 2), sigma_lat(1 : 2), sigma_ot(1 : 2), maxval_likelihood
    real(kind = fp),    intent(in),  optional :: lon_array(:), lat_array(:)
    real(kind = sp),    intent(in),  optional :: az_obs(:), appvel_obs(:), min_correlation(:), array_maxamp(:), &
    &                                            array_lta(:)
    logical,            intent(in),  optional :: array_used(:)

    integer         :: i, ios, ot_year, ot_julianday, ot_mo, ot_dy, ot_hour, ot_min, ot_sec, maxloc_likelihood_particle(1), &
    &                  sigma_index
    real(kind = fp) :: epicenter_lon, epicenter_lat, ot_list, ot_from_day, sigma_diff, &
    &                  sigma_normalize_lon(1 : 2), sigma_normalize_lat(1 : 2), sigma_normalize_ot(1 : 2)
    logical         :: leap
    character(len = 255) :: outfile

    maxloc_likelihood_particle = maxloc(likelihood_particle(:))
    if(present(maxval_likelihood)) maxval_likelihood = likelihood_particle(maxloc_likelihood_particle(1))
    epicenter_lon = lon_particle(maxloc_likelihood_particle(1))
    epicenter_lat = lat_particle(maxloc_likelihood_particle(1))
    ot_list = origintime(maxloc_likelihood_particle(1))
    ot_from_day = real(sec_from_day, kind = fp) + ot_list
    ot_julianday = julianday
    ot_year = year
    if(ot_from_day .le. 0.0_fp) then
      ot_from_day = ot_from_day + 86400.0_fp
      ot_julianday = ot_julianday - 1
    endif
    if(ot_julianday .le. 0) then
      ot_year = ot_year - 1
      call leapyear(ot_year, leap)
      if(leap) then
        ot_julianday = ot_julianday + 366
      else
        ot_julianday = ot_julianday + 365
      endif
    endif
    ot_hour = int(ot_from_day / 3600.0_fp)
    ot_min  = int((ot_from_day - real(ot_hour, kind = fp) * 3600.0_fp) / 60.0_fp)
    ot_sec  = int(ot_from_day  - real(ot_hour, kind = fp) * 3600.0_fp - real(ot_min, kind = fp) * 60.0_fp)
    call jday2ymd(ot_julianday, ot_year, ot_mo, ot_dy)
    write(epicenter_info, '(i4, a, 5(i2.2, a), 2(a, f0.4))') &
    &  ot_year, "-", ot_mo, "-", ot_dy, "T", ot_hour, ":", ot_min, ":", ot_sec, " ", &
    &  "Lon = ", epicenter_lon, " Lat = ", epicenter_lat

    if(present(sigma_lon) .and. present(sigma_lat) .and. present(sigma_ot)) then
      sigma_lon(1 : 2) = 0.0_fp
      sigma_lat(1 : 2) = 0.0_fp
      sigma_ot(1 : 2) = 0.0_fp
      sigma_normalize_lon(1 : 2) = 0.0_fp
      sigma_normalize_lat(1 : 2) = 0.0_fp
      sigma_normalize_ot(1 : 2) = 0.0_fp
      do i = 1, nparticle
        sigma_diff = epicenter_lon - lon_particle(i)
        if(sigma_diff .ge. 0.0_fp) then
          sigma_index = 1
        else
          sigma_index = 2
        endif
        if(sigma_diff .ne. 0.0_fp) then
          sigma_lon(sigma_index) = sigma_lon(sigma_index) + likelihood_particle(i) * abs(sigma_diff)
          sigma_normalize_lon(sigma_index) = sigma_normalize_lon(sigma_index) + likelihood_particle(i)
        endif

        sigma_diff = epicenter_lat - lat_particle(i)
        if(sigma_diff .ge. 0.0_fp) then
          sigma_index = 1
        else
          sigma_index = 2
        endif
        if(sigma_diff .ne. 0.0_fp) then
          sigma_lat(sigma_index) = sigma_lat(sigma_index) + likelihood_particle(i) * abs(sigma_diff)
          sigma_normalize_lat(sigma_index) = sigma_normalize_lat(sigma_index) + likelihood_particle(i)
        endif
        sigma_ot  = sigma_ot  + likelihood_particle(i) * abs(ot_list - origintime(i))

        sigma_diff = ot_list - origintime(i)
        if(sigma_diff .ge. 0.0_fp) then
          sigma_index = 1
        else
          sigma_index = 2
        endif
        if(sigma_diff .ne. 0.0_fp) then
          sigma_ot(sigma_index) = sigma_ot(sigma_index) + likelihood_particle(i) * abs(sigma_diff)
          sigma_normalize_ot(sigma_index) = sigma_normalize_ot(sigma_index) + likelihood_particle(i)
        endif
      enddo
      do i = 1, 2
        if(sigma_normalize_lon(i) .ne. 0.0_fp) sigma_lon(i) = sigma_lon(i) / sigma_normalize_lon(i)
        if(sigma_normalize_lat(i) .ne. 0.0_fp) sigma_lat(i) = sigma_lat(i) / sigma_normalize_lat(i)
        if(sigma_normalize_ot(i)  .ne. 0.0_fp) sigma_ot(i)  = sigma_ot(i)  / sigma_normalize_ot(i)
      enddo

      write(outfile, '(i4, 5(i2.2), a)') ot_year, ot_mo, ot_dy, ot_hour, ot_min, ot_sec, "_likelihood_particle.dat"
      open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 3)
      ios = 0
      do i = 1, nparticle
        ios = ios + 1
        write(10, rec = ios) real(lon_particle(i), kind = sp), real(lat_particle(i), kind = sp), &
        &                    real(likelihood_particle(i), kind = sp)
      enddo
      close(10)
      if(present(lon_array) .and. present(lat_array) .and. present(az_obs) .and. present(appvel_obs) .and. &
      &  present(array_used) .and. present(min_correlation) .and. present(array_maxamp) .and. present(array_lta)) then
        write(outfile, '(i4, 5(i2.2), a)') ot_year, ot_mo, ot_dy, ot_hour, ot_min, ot_sec, "_array_obs.dat"
        open(unit = 10, file = trim(outfile), form = "unformatted", access = "direct", recl = 4 * 7)
        ios = 0
        do i = 1, ubound(array_used, 1)
          if(array_used(i)) then
            ios = ios + 1
            write(10, rec = ios) real(lon_array(i), kind = sp), real(lat_array(i), kind = sp), az_obs(i), appvel_obs(i), &
            &                    min_correlation(i), array_maxamp(i), array_lta(i)
          endif
        enddo
        close(10)
      endif
    endif
 
    return
  end subroutine epicenter2char

  subroutine plot_eplist(iwin_plot, epicenter_info, plot_x, plot_y)
    use nrtype, only : sp
    use aeluma_parameters
    implicit none

    integer,            intent(in) :: iwin_plot
    character(len = *), intent(in) :: epicenter_info
    real(kind = sp), intent(inout) :: plot_x, plot_y

    call pc_text(iwin_plot, plot_x, plot_y, 6.0, trim(epicenter_info), 0.0, len(trim(epicenter_info)), 7)
    plot_y = plot_y - plot_dy_eplist
 
    return
 end subroutine plot_eplist    

end module plotmodule
