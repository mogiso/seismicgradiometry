program plot_map_vector
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  use particle_filter
  use legend
  use random_number
  use mapprojection
  use plotmodule
  implicit none


  integer         :: i, j, k, ios, ncoastline, narray, ntriangle
  integer         :: year, month, day, hr, mi, sc, julianday, sec_from_day
  real(kind = fp) :: slowness_x, slowness_y, origintime_median
  real(kind = fp), allocatable :: appvel_obs(:), az_obs(:), lon_array(:), lat_array(:), min_correlation(:), arrivaltime(:), &
  &                               origintime_candidate(:)
  integer,         allocatable :: arrayindex(:)
  logical,         allocatable :: result_exist(:, :)

  real(kind = fp)                   :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
  character(len = 255)              :: coastline_txt
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 2)                :: yr, mo, dy, hh, mm, ss
  
  real(kind = fp)      :: lon_particle(1 : nparticle), lat_particle(1 : nparticle), &
  &                       lon_particle_list(1 : nparticle, 1 : nepicenter), &
  &                       lat_particle_list(1 : nparticle, 1 : nepicenter), &
  &                       likelihood_particle_list(1 : nparticle, 1 : nepicenter), &
  &                       origintime_list(1 : nepicenter), &
  &                       max_likelihood(1 : nepicenter), &
  &                       likelihood_particle(1 : nparticle), az_weight(1 : int(2.0_fp * pi / daz_weight), 1 : nepicenter)
  logical              :: epicenter_exist(1 : nepicenter)

  !!initiate random number generator
  call make_seed(seed)
  call random_generator_init(random_status, seed)

  !!Read coastline
  call getarg(1, coastline_txt)
  call read_coastline(coastline_txt, mapbuf)
  ncoastline = ubound(mapbuf, 1)

  !!initiaalize event list
  epicenter_exist(1 : nepicenter) = .false.
 
  !!Plot legend
  call pc_plotinit(iwin_legend, "Legend", 0.0_sp, -300.0_sp, width / 2, 27.0_sp, scale)
  call plot_legend

  !!Read and plot AELUMA results
  call pc_plotinit(iwin_map, "AELUMA results", 0.0_sp, 0.0_sp, width, height, scale)
  call pc_setbkcolor(iwin_map, 255, 255, 255)
  call mercator(center_lon, lon_w, lat_s, width_tmp(1), height_tmp(1))
  call mercator(center_lon, lon_e, lat_n, width_tmp(2), height_tmp(2))
  dwidth = 1.0_fp / (width_tmp(2) - width_tmp(1))
  dheight = 1.0_fp / (height_tmp(2) - height_tmp(1))
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
      allocate(az_obs(1 : ntriangle), appvel_obs(1 : ntriangle), result_exist(1 : ntriangle, 1 : nepicenter), &
      &        arrayindex(1 : ntriangle), lon_array(1 : ntriangle), lat_array(1 : ntriangle), &
      &        min_correlation(1 : ntriangle), arrivaltime(1 : ntriangle))
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
        if(appvel_obs(arrayindex(i)) .lt. phasevelocity) appvel_obs(arrayindex(i)) = phasevelocity
        az_weight(int(az_obs(arrayindex(i)) / daz_weight) + 1) = az_weight(int(az_obs(arrayindex(i)) / daz_weight) + 1) + 1.0_fp
        result_exist(arrayindex(i), 1 : nepicenter) = .true.
      enddo

      !!estimate location
      narray_use(1 : nepicenter) = narray
      do k = 1, nepicenter
        if(epicenter_exist(j)) then
          do j = 1, narray
            likelihood_tmp = 0.0_fp
            do i = 1, nparticle
              call greatcircle_dist(lat_array(arrayindex(j)), lon_array(arrayindex(j)), &
              &                     lat_particle_list(i, k),  lon_particle_list(i, k),  &
              &                     distance = dist_tmp,      azimuth = az_tmp)
              ttime_diff = (origintime_list(i, k) - dist_tmp * appvel_obs(j))
              az_tmp = az_tmp * pi
              if(az_tmp .ge. 2.0_fp * pi) az_tmp = az_tmp - 2.0_fp * pi
              daz = az_obs(j) - az_tmp
              likelihood_tmp = likelihood_tmp &
              &              + likelihood_particle(i, k) * 0.5_fp / (pi * sigma_traveltimediff * daz_weight) &
              &              * exp(-0.5_fp * ((ttime_diff * ttime_diff / sigma_traveltimediff / sigma_traveltimediff) &
              &                            +  (daz        * daz        / daz_weight2)))
            enddo
            if(likelihood_tmp .lt. min_likelihood_eqobs) then
              result_exist(arrayindex(j), k) = .false.
              narray_use(k) = narray_use(k) - 1
            else
              result_exist(arrayindex(j), k + 1 : nepicenter) = .false.
              narray_use(k + 1 : nepicenter) = narray_use(k * 1 : nepicenter) - 1
            endif
          enddo
          if(narray_use .ge. narray_use_min) then
            
        else
          call particle_filter_init(random_status, lon_particle_list(:, k), lat_particle_list(:, k))
        endif

        call particle_filter_search(narray, arrayindex, result_exist(:, k), lon_array, lat_array, min_correlation, az_obs, &
        &                           az_weight, random_status, &
        &                           lon_particle_list(:, k), lat_particle_list(:, k), likelihood_particle_list(:, k), &
        &                           appvel = appvel_obs, arrivaltime = arrivaltime, origintime = origintime_list(k))

        !!write particles
        if(epicenter_exist(k) &
        &  call plot_particle(lon_particle_list(:, k), lat_particle_list(:, k), likelihood_particle_list(:, k), &
        &                     width_tmp, height_tmp, dwidth, dheight)

    endif

    !!plot slowness vector
    call plot_slowness_vector(narray, arrayindex, result_exist, lon_array, lat_array, az_obs, min_correlation, &
    &                         width_tmp, height_tmp, dwidth, dheight)
    !!plot map
    call pc_setline(iwin_map, 1)
    call pc_setcolor(iwin_map, 0, 0, 0)
    call plot_currentdate(yr, mo, dy, hh, mm, ss, width_tmp, height_tmp, dwidth, dheight)
    call plot_coastline(ncoastline, mapbuf, width_tmp, height_tmp, dwidth, dheight)

    call pc_flush(iwin_map)
  enddo


  call pc_plotend(iwin_map, 1)
  call pc_plotend(iwin_legend, 1)

  stop
end program plot_map_vector

   
