module plotmodule
  private
  public :: plot_particle, plot_slowness_vector, plot_coastline, plot_currentdate

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
    integer :: i, maxloc_likelihood(1), color(1 : 3)

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
    maxloc_likelihood = maxloc(likelihood_particle)
    call mercator(center_lon, lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), map_x, map_y)
    plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
    call pc_symbol(iwin_map, plot_x, plot_y, 9.0_sp, 1, 0)
    call pc_setcolor(iwin_map, 255, 255, 255)
    call pc_symbol(iwin_map, plot_x, plot_y, 4.0_sp, 1, 0)

    write(0, '(a, f0.4, a, f0.4, a, e15.7)') " Lon = ", lon_particle(maxloc_likelihood(1)), &
    &                                        " Lat = ", lat_particle(maxloc_likelihood(1)), &
    &                                        " likelihood = ", likelihood_particle(maxloc_likelihood(1))

    return
  end subroutine plot_particle

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
      if(min_correlation(arrayindex(i)) .lt. 0.5_fp) then
        cycle
      elseif(min_correlation(arrayindex(i)) .ge. 0.5_fp .and. min_correlation(arrayindex(i)) .lt. 0.6_fp) then
        color(1 : 3) = [220, 204, 222]
      elseif(min_correlation(arrayindex(i)) .ge. 0.6_fp .and. min_correlation(arrayindex(i)) .lt. 0.7_fp) then
        color(1 : 3) = [212, 156, 189]
      elseif(min_correlation(arrayindex(i)) .ge. 0.7_fp .and. min_correlation(arrayindex(i)) .lt. 0.8_fp) then
        color(1 : 3) = [196, 110, 155]
      elseif(min_correlation(arrayindex(i)) .ge. 0.8_fp .and. min_correlation(arrayindex(i)) .lt. 0.9_fp) then
        color(1 : 3) = [136, 97, 141]
      elseif(min_correlation(arrayindex(i)) .ge. 0.9_fp) then
        color(1 : 3) = [73, 57, 100]
      endif
      call pc_setcolor(iwin_map, color(1), color(2), color(3))
      call pc_vector(iwin_map, plot_x, plot_y, plot_theta, vector_len, vector_width, vector_head1, vector_head2, 1)
    enddo
  end subroutine plot_slowness_vector

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
    call pc_text(iwin_map, plot_x, plot_y, 7.0, date_txt, 0.0, len(date_txt), 7)

    return
  end subroutine plot_currentdate
 

end module plotmodule
