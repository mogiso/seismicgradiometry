program plot_map_vector
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  use particlefilter
  use random_number
  use mapprojection
  use plotmodule
  use particlefilter_functions
  implicit none

  real(kind = sp) :: plot_x_tmp, plot_y_tmp
  integer         :: i, j, k, ios, ncoastline, narray, ntriangle, integer_tmp
  integer         :: year, month, day, hr, mi, sc, julianday, sec_from_day
  real(kind = fp) :: slowness_x, slowness_y, daz, az_tmp, dist_tmp, likelihood_tmp, ot_diff, kahan_val1, kahan_val2, &
  &                  error_lon, error_lat, error_ot, maxval_likelihood, appvel_median
  logical         :: no_associated_arrayuse
  real(kind = fp), allocatable :: appvel_obs(:), az_obs(:), lon_array(:), lat_array(:), min_correlation(:), arrivaltime(:)
  real(kind = sp), allocatable :: az_obs_used(:, :), appvel_obs_used(:, :)
  integer,         allocatable :: arrayindex(:)
  logical,         allocatable :: result_exist(:, :), result_exist_org(:), result_exist_tmp(:), array_used_list(:, :)

  real(kind = fp)                   :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
  character(len = 255)              :: coastline_txt, epicenter_info
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 2)                :: yr, mo, dy, hh, mm, ss
  character(len = 10)               :: text_tmp
  
  integer              :: narray_use(1 : nepicenter), narray_use_list(1 : nepicenter), &
  &                       epicenter_acceptcount(1 : nepicenter)
  real(kind = fp)      :: lon_particle_list(1 : nparticle, 1 : nepicenter),        &
  &                       lat_particle_list(1 : nparticle, 1 : nepicenter),        &
  &                       origintime_list(1 : nparticle, 1 : nepicenter),          &
  &                       maxval_likelihood_particle_list(1 : nepicenter),         &
  &                       appvel_median_list(1 : nepicenter),                      &
  &                       likelihood_particle_list(1 : nparticle, 1 : nepicenter), &
  &                       lon_particle(1 : nparticle), lat_particle(1 : nparticle), likelihood_particle(1 : nparticle), &
  &                       origintime(1 : nparticle), az_weight(1 : int(2.0_fp * pi / daz_weight) + 1)
  logical              :: epicenter_exist(1 : nepicenter)

  !!initiate random number generator
  call make_seed(seed)
  call random_generator_init(random_status, seed)
  print '(a, i0)', "seed = ", seed

  !!Read coastline
  call getarg(1, coastline_txt)
  call read_coastline(coastline_txt, mapbuf)
  ncoastline = ubound(mapbuf, 1)

  !!initiaalize event list
  epicenter_exist(1 : nepicenter) = .false.
 
  !!Plot legend
  call pc_plotinit(iwin_legend, "Legend", 0.0_sp, -300.0_sp, width / 2, 27.0_sp, scale)
  call plot_legend(iwin_legend)

  !!Open epicenter window
  call pc_plotinit(iwin_eplist, "Epicenter list", 0.0_sp, -300.0_sp, width / 2 + 60.0_sp, 27.0_sp, scale)

  !!Read and plot AELUMA results
  call pc_plotinit(iwin_map, "AELUMA results", 0.0_sp, 0.0_sp, width, height, scale)
  call pc_setbkcolor(iwin_map, 255, 255, 255)
  call mercator(center_lon, lon_w, lat_s, width_tmp(1), height_tmp(1))
  call mercator(center_lon, lon_e, lat_n, width_tmp(2), height_tmp(2))
  dwidth = 1.0_fp / (width_tmp(2) - width_tmp(1))
  dheight = 1.0_fp / (height_tmp(2) - height_tmp(1))

  epicenter_acceptcount(1 : nepicenter) = 0
  !!read AELUMA results from stdin
  do 
    plot_x_tmp = plot_x_eplist
    plot_y_tmp = plot_y_eplist
    call pc_clear(iwin_map)
    call pc_setcolor(iwin_map, 0, 0, 0)
    call pc_clear(iwin_eplist)
    call pc_setcolor(iwin_eplist, 0, 0, 0)
    call pc_setline(iwin_eplist, 2)
    !az_weight(1 : int(2.0_fp * pi / daz_weight)) = 0.0_fp
    read(*, *, iostat = ios) yr, mo, dy, hh, mm, ss, narray, ntriangle
    if(ios .ne. 0) error stop
    read(yr, *) year; year = year + 2000
    read(mo, *) month
    read(dy, *) day
    read(hh, *) hr
    read(mm, *) mi
    read(ss, *) sc
    call ymd2jday(year, month, day, julianday)
    sec_from_day = hr * 60 * 60 + mi * 60 + sc
    write(0, '(8(i0, 1x))') year, month, day, hr, mi, sc, narray, ntriangle
    !print '(5(i0, 1x))', year, julianday, sec_from_day, narray, ntriangle

    if(.not. allocated(arrayindex)) then
      allocate(az_obs(1 : ntriangle), appvel_obs(1 : ntriangle), result_exist(1 : ntriangle, 1 : nepicenter), &
      &        result_exist_org(1 : ntriangle), result_exist_tmp(1 : ntriangle), &
      &        arrayindex(1 : ntriangle), lon_array(1 : ntriangle), lat_array(1 : ntriangle), &
      &        min_correlation(1 : ntriangle), arrivaltime(1 : ntriangle))
      allocate(az_obs_used(1 : ntriangle, 1 : nepicenter), appvel_obs_used(1 : ntriangle, 1 : nepicenter), &
      &        array_used_list(1 : ntriangle, 1 : nepicenter))
    endif

    result_exist_org(1 : ntriangle)             = .false.
    result_exist(1 : ntriangle, 1 : nepicenter) = .false.
    call pc_setline(iwin_map, 4)

    !!read and plot slowness vector
    narray_use(1 : nepicenter) = narray
    do i = 1, narray
      read(*, *) arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
      &          slowness_x, slowness_y, min_correlation(arrayindex(i)), arrivaltime(arrayindex(i))
      result_exist_org(arrayindex(i))  = .true.
      do j = 1, nepicenter
        result_exist(arrayindex(i), j) = .true.
        if(min_correlation(arrayindex(i)) .le. correlation_threshold) then
          result_exist(arrayindex(i), j) = .false.
          narray_use(j) = narray_use(j) - 1
        endif
      enddo

      !!arrival time: relative time in s from current time
      arrivaltime(arrayindex(i)) = -(real(nsec_buf, kind = fp) - arrivaltime(arrayindex(i)))
      az_obs(arrayindex(i)) = atan2(slowness_x, slowness_y)
      if(az_obs(arrayindex(i)) .lt. 0.0_fp) az_obs(arrayindex(i)) = az_obs(arrayindex(i)) + 2.0_fp * pi
      appvel_obs(arrayindex(i)) = 1.0_fp / sqrt(slowness_x ** 2 + slowness_y ** 2)
    enddo

    !!associate observation and events
    do k = 1, nepicenter
      origintime_list(1 : nparticle, k) = origintime_list(1 : nparticle, k) - dtimestep

      if(epicenter_exist(k)) then
        do j = 1, narray
          if(.not. result_exist(arrayindex(j), k)) cycle 
          likelihood_tmp = 0.0_fp
          kahan_val1 = 0.0_fp
          do i = 1, nparticle
            call greatcircle_dist(lat_array(arrayindex(j)), lon_array(arrayindex(j)), &
            &                     lat_particle_list(i, k),  lon_particle_list(i, k),  &
            &                     distance = dist_tmp,      azimuth = az_tmp)
            az_tmp = az_tmp + pi
            if(az_tmp .ge. 2.0_fp * pi) az_tmp = az_tmp - 2.0_fp * pi
            daz = delta_az(az_obs(arrayindex(j)), az_tmp)
            ot_diff = origintime_list(i, k) - origintime_cal(arrivaltime(arrayindex(j)), dist_tmp, appvel_obs(arrayindex(j)))
            kahan_val1 = kahan_val1 &
            &          + likelihood_particle_list(i, k) * 0.5_fp / (pi * sigma_otdiff * daz_weight) &
            &          * exp(-0.5_fp * ((ot_diff * ot_diff / sigma_otdiff2) &
            &                        +  (daz     * daz     / daz_weight2)))
            kahan_val2 = likelihood_tmp
            likelihood_tmp = likelihood_tmp + kahan_val1
            kahan_val2 = likelihood_tmp - kahan_val2
            kahan_val1 = kahan_val1 - kahan_val2
          enddo
          if(likelihood_tmp .lt. min_likelihood_eqobs) then
            result_exist(arrayindex(j), k) = .false.
            narray_use(k) = narray_use(k) - 1
          else
            do i = 1, nepicenter
              if(i .eq. k) cycle
              if(result_exist(arrayindex(j), i)) then
                result_exist(arrayindex(j), i) = .false.
                narray_use(i) = narray_use(i) - 1
              endif
            enddo
          endif
        enddo
      endif
    enddo

    no_associated_arrayuse = .false.
    do i = 1, nepicenter
      if(.not. epicenter_exist(i)) then
        if(.not. no_associated_arrayuse) then
          no_associated_arrayuse = .true.
        else
          narray_use(i) = 0
        endif
      endif
    enddo

    do i = 1, nepicenter
      if(epicenter_exist(i)) then
        if(narray_use(i) .lt. 1) then
          if(epicenter_acceptcount(i) .ge. epicenter_acceptcount_threshold) then
            call epicenter2char(year, julianday, sec_from_day, lon_particle_list(:, i), lat_particle_list(:, i), &
            &                   origintime_list(:, i), likelihood_particle_list(:, i), epicenter_info, &
            &                   sigma_lon = error_lon, sigma_lat = error_lat, sigma_ot = error_ot, &
            &                   maxval_likelihood = maxval_likelihood, &
            &                   lon_array = lon_array, lat_array = lat_array, &
            &                   az_obs = az_obs_used(:, i), appvel_obs = appvel_obs_used(:, i), &
            &                   array_used = array_used_list(:, i))
            print '(a, 4(1x, e15.7), 2(1x, i0), 1x, f0.3)', trim(epicenter_info), error_lon, error_lat, error_ot, &
            &                                               maxval_likelihood, narray_use_list(i), epicenter_acceptcount(i), &
            &                                               appvel_median_list(i)

          endif
          epicenter_exist(i) = .false.
          epicenter_acceptcount(i) = 0
          cycle
        endif
        !!plot particles
        call plot_particle(iwin_map, lon_particle_list(:, i), lat_particle_list(:, i), likelihood_particle_list(:, i), &
                           width_tmp, height_tmp, dwidth, dheight)
        call plot_particle_maxlikelihood(iwin_map, lon_particle_list(:, i), lat_particle_list(:, i), &
        &                                likelihood_particle_list(:, i), width_tmp, height_tmp, dwidth, dheight)
        call epicenter2char(year, julianday, sec_from_day, lon_particle_list(:, i), lat_particle_list(:, i), &
        &                   origintime_list(:, i), likelihood_particle_list(:, i), epicenter_info)
        write(text_tmp, '(f4.1)') appvel_median_list(i)
        epicenter_info = trim(epicenter_info) // " " // trim(text_tmp) // "km/s"
        write(text_tmp, '(i0)') narray_use(i)
        epicenter_info = trim(epicenter_info) // " " // trim(text_tmp)
        call plot_eplist(iwin_eplist, epicenter_info, plot_x_tmp, plot_y_tmp)
        if(narray_use(i) .ge. narray_use_min) then
          epicenter_acceptcount(i) = epicenter_acceptcount(i) + 1
          write(0, '(i4.4, 5(a, i2.2), 2a)') year, "-", month, "-", day, "T", hr, ":", mi, ":", sc, " ", trim(epicenter_info)
        endif
      endif
    enddo
    call pc_flush(iwin_eplist)

    !!swap the order of epicenter list
    do j = 1, nepicenter - 1
      if(epicenter_exist(j)) cycle
      do i = j + 1, nepicenter
        if(epicenter_exist(i)) then
          result_exist_tmp(1 : ntriangle) = result_exist(1 : ntriangle, i)
          result_exist(1 : ntriangle, i) = result_exist(1 : ntriangle, j)
          result_exist(1 : ntriangle, j) = result_exist_tmp(1 : ntriangle)
          lon_particle_list(1 : nparticle, j) = lon_particle_list(1 : nparticle, i)
          lat_particle_list(1 : nparticle, j) = lat_particle_list(1 : nparticle, i)
          likelihood_particle_list(1 : nparticle, j) = likelihood_particle_list(1 : nparticle, i)
          origintime_list(1 : nparticle, j) = origintime_list(1 : nparticle, i)
          maxval_likelihood_particle_list(j) = maxval_likelihood_particle_list(i)
          appvel_median_list(j) = appvel_median_list(i)
          narray_use_list(j) = narray_use_list(i)
          integer_tmp = narray_use(j)
          narray_use(j) = narray_use(i)
          narray_use(i) = integer_tmp
          epicenter_acceptcount(j) = epicenter_acceptcount(i)
          az_obs_used(1 : ntriangle, j) = az_obs_used(1 : ntriangle, i)
          appvel_obs_used(1 : ntriangle, j) = appvel_obs_used(1 : ntriangle, i)
          array_used_list(1 : ntriangle, j) = array_used_list(1 : ntriangle, i)
          epicenter_exist(j) = .true.
          epicenter_exist(i) = .false.
          exit
        endif
      enddo
    enddo
 
    !!estimate event epicenters
    do i = 1, nepicenter
      if(narray_use(i) .lt. narray_use_min) cycle
      if(epicenter_exist(i)) then
        if(narray_use(i) .le. narray_use_list(i)) cycle
      endif

      !!calculate azimuthal weighting array
      az_weight(1 : int(2.0_fp * pi / daz_weight) + 1) = 0.0_fp
      do j = 1, narray
        if(result_exist(arrayindex(j), i)) then
          az_weight(int(az_obs(arrayindex(j)) / daz_weight) + 1) &
          &  = az_weight(int(az_obs(arrayindex(j)) / daz_weight) + 1) + 1.0_fp
        endif
      enddo

      !!do particle filter
      !!initialize location of each particle
      if(.not. epicenter_exist(i)) then
        call particlefilter_init(random_status, lon_particle, lat_particle)
      else
        lon_particle(1 : nparticle) = lon_particle_list(1 : nparticle, i)
        lat_particle(1 : nparticle) = lat_particle_list(1 : nparticle, i)
      endif

      call particlefilter_search(narray, arrayindex, result_exist(:, i), lon_array, lat_array, az_obs, &
      &                           az_weight, random_status, lon_particle, lat_particle, likelihood_particle, &
      &                           appvel = appvel_obs, arrivaltime = arrivaltime, origintime = origintime, &
      &                           appvel_median = appvel_median)
      maxval_likelihood = maxval(likelihood_particle)

      !!renew epicenter parameters
      if(.not. epicenter_exist(i)) then
        epicenter_exist(i) = .true.
        lon_particle_list       (1 : nparticle, i) = lon_particle       (1 : nparticle)
        lat_particle_list       (1 : nparticle, i) = lat_particle       (1 : nparticle)
        origintime_list         (1 : nparticle, i) = origintime         (1 : nparticle)
        likelihood_particle_list(1 : nparticle, i) = likelihood_particle(1 : nparticle)
        maxval_likelihood_particle_list(i) = maxval_likelihood
        appvel_median_list(i) = appvel_median
        narray_use_list(i) = narray_use(i)
        array_used_list(1 : ntriangle, i) = result_exist(1 : ntriangle, i)
        az_obs_used(1 : ntriangle, i) = real(az_obs(1 : ntriangle) * rad2deg, kind = sp)
        appvel_obs_used(1 : ntriangle, i) = real(appvel_obs(1 : ntriangle), kind = sp)
      else
        if(maxval_likelihood .ge. maxval_likelihood_particle_list(i)) then
          lon_particle_list       (1 : nparticle, i) = lon_particle       (1 : nparticle)
          lat_particle_list       (1 : nparticle, i) = lat_particle       (1 : nparticle)
          origintime_list         (1 : nparticle, i) = origintime         (1 : nparticle)
          likelihood_particle_list(1 : nparticle, i) = likelihood_particle(1 : nparticle)
          maxval_likelihood_particle_list(i) = maxval_likelihood
          appvel_median_list(i) = appvel_median
          narray_use_list(i) = narray_use(i)
          array_used_list(1 : ntriangle, i) = result_exist(1 : ntriangle, i)
          az_obs_used(1 : ntriangle, i) = real(az_obs(1 : ntriangle) * rad2deg, kind = sp)
          appvel_obs_used(1 : ntriangle, i) = real(appvel_obs(1 : ntriangle), kind = sp)
        endif
      endif
    enddo

    !!plot slowness vector
    call plot_slowness_vector(iwin_map, narray, arrayindex, result_exist_org, lon_array, lat_array, az_obs, appvel_obs, &
    &                         min_correlation, width_tmp, height_tmp, dwidth, dheight)
    !!plot map
    call pc_setline(iwin_map, 1)
    call pc_setcolor(iwin_map, 0, 0, 0)
    call plot_currentdate(iwin_map, yr, mo, dy, hh, mm, ss, width_tmp, height_tmp, dwidth, dheight)
    call plot_coastline(iwin_map, ncoastline, mapbuf, width_tmp, height_tmp, dwidth, dheight)

    call pc_flush(iwin_map)
  enddo

  call pc_plotend(iwin_map, 1)
  call pc_plotend(iwin_legend, 1)

  stop
end program plot_map_vector

