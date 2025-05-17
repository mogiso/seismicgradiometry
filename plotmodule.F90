module plotmodule
  private
  public :: plot_particle, plot_particle_maxlikelihood, plot_slowness_vector, read_coastline, plot_coastline, &
  &         plot_currentdate, plot_legend, epicenter2char, plot_eplist

  contains

  subroutine plot_particle(lon_particle, lat_particle, likelihood_particle, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use aeluma_parameters
    use mapprojection
    implicit none
    real(kind = fp), intent(in) :: lon_particle(:), lat_particle(:), likelihood_particle(:), &
    &                              width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    real(kind = fp) :: map_x, map_y, likelihood_tmp
    real(kind = sp) :: plot_x, plot_y
    integer :: i, color(1 : 3)

    call pc_setline(iwin_map, 1)
    do i = 1, nparticle
      call mercator(center_lon, lon_particle(i), lat_particle(i), map_x, map_y)
      plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
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
    !maxloc_likelihood = maxloc(likelihood_particle)
    !call mercator(center_lon, lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), map_x, map_y)
    !plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
    !plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
    !call pc_symbol(iwin_map, plot_x, plot_y, 9.0_sp, 1, 0)
    !call pc_setcolor(iwin_map, 255, 255, 255)
    !call pc_symbol(iwin_map, plot_x, plot_y, 4.0_sp, 1, 0)

    !write(0, '(a, f0.4, a, f0.4, a, e15.7)') " Lon = ", lon_particle(maxloc_likelihood(1)), &
    !&                                        " Lat = ", lat_particle(maxloc_likelihood(1)), &
    !&                                        " likelihood = ", likelihood_particle(maxloc_likelihood(1))

    return
  end subroutine plot_particle

  subroutine plot_particle_maxlikelihood(lon_particle, lat_particle, likelihood_particle, &
  &                                      width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use aeluma_parameters
    use mapprojection
    implicit none
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
    call pc_symbol(iwin_map, plot_x, plot_y, 8.0_sp, 1, 0)
    call pc_setcolor(iwin_map, 255, 255, 255)
    call pc_symbol(iwin_map, plot_x, plot_y, 3.5_sp, 1, 0)

    write(0, '(a, f0.4, a, f0.4, a, e15.7)') " Lon = ", lon_particle(maxloc_likelihood(1)), &
    &                                        " Lat = ", lat_particle(maxloc_likelihood(1)), &
    &                                        " likelihood = ", likelihood_particle(maxloc_likelihood(1))

  end subroutine plot_particle_maxlikelihood

  subroutine plot_slowness_vector(narray, arrayindex, result_exist, lon_array, lat_array, az_obs, min_correlation, &
  &                               width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use constants, only : rad2deg
    use aeluma_parameters
    use mapprojection

    implicit none
    integer,         intent(in) :: narray, arrayindex(1 : narray)
    logical,         intent(in) :: result_exist(:)
    real(kind = fp), intent(in) :: lon_array(:), lat_array(:), az_obs(:), min_correlation(:), &
    &                              width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight

    integer :: i, color(1 : 3)
    real(kind = sp) :: plot_theta, plot_x, plot_y
    real(kind = fp) :: map_x, map_y

    call pc_setline(iwin_map, 4)
    do i = 1, narray
      if(.not. result_exist(arrayindex(i))) cycle
      !theta = atan2(slowness_x, slowness_y) * rad2deg
      plot_theta = 90.0_sp - real(az_obs(arrayindex(i)) * rad2deg, kind = sp)
      call mercator(center_lon, lon_array(arrayindex(i)), lat_array(arrayindex(i)), map_x, map_y)
      plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
      plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
      if(min_correlation(arrayindex(i)) .lt. 0.1_fp) then
        cycle
      elseif(min_correlation(arrayindex(i)) .ge. 0.1_fp .and. min_correlation(arrayindex(i)) .lt. 0.3_fp) then
        color(1 : 3) = [220, 204, 222]
      elseif(min_correlation(arrayindex(i)) .ge. 0.3_fp .and. min_correlation(arrayindex(i)) .lt. 0.5_fp) then
        color(1 : 3) = [212, 156, 189]
      elseif(min_correlation(arrayindex(i)) .ge. 0.5_fp .and. min_correlation(arrayindex(i)) .lt. 0.7_fp) then
        color(1 : 3) = [196, 110, 155]
      elseif(min_correlation(arrayindex(i)) .ge. 0.7_fp .and. min_correlation(arrayindex(i)) .lt. 0.9_fp) then
        color(1 : 3) = [136, 97, 141]
      elseif(min_correlation(arrayindex(i)) .ge. 0.9_fp) then
        color(1 : 3) = [73, 57, 100]
      endif
      call pc_setcolor(iwin_map, color(1), color(2), color(3))
      call pc_vector(iwin_map, plot_x, plot_y, plot_theta, vector_len, vector_width, vector_head1, vector_head2, 1)
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

  subroutine plot_coastline(ncoastline, mapbuf, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : fp, sp
    use mapprojection
    use aeluma_parameters
    implicit none

    real(kind = fp),    intent(in) :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    integer,            intent(in) :: ncoastline
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
      call pc_line(iwin_map, plot_x, plot_y, plot_x1, plot_y1)
    enddo
  end subroutine plot_coastline

  subroutine plot_currentdate(yr, mo, dy, hh, mm, ss, width_tmp, height_tmp, dwidth, dheight)
    use nrtype, only : sp, fp
    use aeluma_parameters
    use mapprojection
    implicit none
    character(len = 2), intent(in) :: yr, mo, dy, hh, mm, ss
    real(kind = fp),    intent(in) :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
    real(kind = sp)                :: plot_x, plot_y
    real(kind = fp)                :: map_x, map_y
    character(len = 19)            :: date_txt

    date_txt = "20" // yr // "/" // mo // "/" // dy // " " // hh // ":" // mm // ":" // ss
    call mercator(center_lon, lon_w, lat_n, map_x, map_y)
    plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
    call pc_rect(iwin_map, plot_x, plot_y - 10.0_sp, plot_x + 80.0_sp, plot_y - 10.0_sp)
    call pc_setcolor(iwin_map, 255, 255, 255)
    call pc_text(iwin_map, plot_x, plot_y, 7.0, date_txt, 0.0, len(date_txt), 7)

    return
  end subroutine plot_currentdate
 

  subroutine plot_legend
    use nrtype, only : sp
    use aeluma_parameters

    implicit none
    integer            :: i, color(1 : 3)
    real(kind = sp)    :: plot_x, plot_y
    character(len = 6) :: plottext

    call pc_setbkcolor(iwin_legend, 255, 255, 255)
    do i = 1, 5
      if(i .eq. 1) then
        color(1 : 3) = [220, 204, 222]
        plottext = "0.1   "
      elseif(i .eq. 2) then
        color(1 : 3) = [212, 156, 189]
        plottext = "0.3   "
      elseif(i .eq. 3) then
        color(1 : 3) = [196, 110, 155]
        plottext = "0.5   "
      elseif(i .eq. 4) then
        color(1 : 3) = [136, 97, 141]
        plottext = "0.7   "
      elseif(i .eq. 5) then
        color(1 : 3) = [73, 57, 100]
        plottext = "0.9   "
      endif
      plot_x = real(i, kind = sp) * 16.0_sp
      plot_y = 20.0_sp
      call pc_setcolor(iwin_legend, color(1), color(2), color(3))
      call pc_setline(iwin_legend, 4)
      call pc_vector(iwin_legend, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
      call pc_setcolor(iwin_legend, 0, 0, 0)
      plot_x = real(i, kind = sp) * 16.0_sp - 10.0_sp
      call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 4) 
    enddo

    do i = 1, 6
      if(i .eq. 1) then
        color(1 : 3) = [252, 238, 158]
        plottext = "<0.1"
      elseif(i .eq. 2) then
        color(1 : 3) = [238, 179, 87]
        plottext = "0.1"
      elseif(i .eq. 3) then
        color(1 : 3) = [222, 117, 79]
        plottext = "0.3"
      elseif(i .eq. 4) then
        color(1 : 3) = [149, 66, 62]
        plottext = "1.0"
      elseif(i .eq. 5) then
        color(1 : 3) = [63, 39, 23]
        plottext = "3.0"
      elseif(i .eq. 6) then
        color(1 : 3) = [26, 26, 1]
        plottext = "10.0"
      endif
      plot_x = real(i, kind = sp) * 16.0_sp
      plot_y = 10.0_sp
      call pc_setcolor(iwin_legend, color(1), color(2), color(3))
      call pc_symbol(iwin_legend, plot_x, plot_y, 3.0_sp, 1, 0)
      call pc_setcolor(iwin_legend, 0, 0, 0)
      call pc_symbol(iwin_legend, plot_x, plot_y, 3.0_sp, 1, 1)
      plot_x = real(i, kind = sp) * 16.0_sp - 8.0_sp
      call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 5) 
    enddo
    plot_x = real(i, kind = sp) * 16.0_sp - 8.0_sp
    plottext = "10.0<"
    call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 5) 
    plot_x = plot_x + 20.0_sp
    call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(likelihood_legend_normalize_c), 0.0, &
    &            len(trim(likelihood_legend_normalize_c)), 5) 
    call pc_flush(iwin_legend)

    return
  end subroutine plot_legend

  subroutine epicenter2char(year, julianday, sec_from_day, lon_particle, lat_particle, origintime, likelihood_particle, &
  &                         epicenter_info, sigma_lon, sigma_lat, sigma_ot, maxval_likelihood)
    use nrtype, only : fp, sp
    use aeluma_parameters
    use jday
    implicit none

    integer,            intent(in)            :: year, julianday, sec_from_day
    real(kind = fp),    intent(in)            :: lon_particle(:), lat_particle(:), origintime(:), likelihood_particle(:)
    character(len = *), intent(out)           :: epicenter_info
    real(kind = fp),    intent(out), optional :: sigma_lon, sigma_lat, sigma_ot, maxval_likelihood

    integer         :: i, ot_year, ot_julianday, ot_mo, ot_dy, ot_hour, ot_min, ot_sec, maxloc_likelihood_particle(1)
    real(kind = fp) :: epicenter_lon, epicenter_lat, ot_list, ot_from_day
    logical         :: leap

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
    write(epicenter_info, '(i4, a, 5(i2.2, a), 2(a, f0.3))') &
    &  ot_year, "/", ot_mo, "/", ot_dy, " ", ot_hour, ":", ot_min, ":", ot_sec, " ", &
    &  "Lon = ", epicenter_lon, " Lat = ", epicenter_lat

    if(present(sigma_lon) .and. present(sigma_lat) .and. present(sigma_ot)) then
      sigma_lon = 0.0_fp
      sigma_lat = 0.0_fp
      sigma_ot = 0.0_fp
      do i = 1, nparticle
        sigma_lon = sigma_lon + likelihood_particle(i) * (epicenter_lon - lon_particle(i)) ** 2
        sigma_lat = sigma_lat + likelihood_particle(i) * (epicenter_lat - lat_particle(i)) ** 2
        sigma_ot  = sigma_ot  + likelihood_particle(i) * (ot_list - origintime(i)) ** 2
      enddo
      sigma_lon = sqrt(sigma_lon)
      sigma_lat = sqrt(sigma_lat)
      sigma_ot  = sqrt(sigma_ot)
    endif
 
    return
  end subroutine epicenter2char

  subroutine plot_eplist(epicenter_info, plot_x, plot_y)
    use nrtype, only : sp
    use aeluma_parameters
    implicit none

    character(len = *), intent(in) :: epicenter_info
    real(kind = sp), intent(inout) :: plot_x, plot_y

    call pc_text(iwin_eplist, plot_x, plot_y, 5.0, trim(epicenter_info), 0.0, len(trim(epicenter_info)), 7)
    plot_y = plot_y - plot_dy_eplist
 
    return
 end subroutine plot_eplist    

end module plotmodule
