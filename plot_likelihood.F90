program plot_likelihood
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  use random_number
  use mapprojection
  use plotmodule
  use particlefilter_functions
  use sort
  implicit none

  integer, parameter :: iwin_likelihood(1 : 3) = [0, 1, 2]
  integer, parameter :: nparticle_plot = 3000000
  integer, parameter :: color(1 : 3, 1 : 10) = reshape([255, 255, 153, &
  &                                                     255, 255,  51, &
  &                                                     255, 221,   0, &
  &                                                     255, 153,   0, &
  &                                                     255,  85,   0, &
  &                                                     255,  17,   0, &
  &                                                     204,   0,   0, &
  &                                                     136,   0,   0, &
  &                                                      68,   0,   0, &
  &                                                       0,   0,   0], [3, 10])
  real(kind = fp), parameter :: coefficient = 1.0e+2;

  real(kind = sp) :: plot_x, plot_y
  integer         :: i, j, ios, ncoastline, narray, ntriangle, narray_use_tmp, az_weight_index, color_index, narray_tmp
  integer         :: year, month, day, hr, mi, sc, julianday, sec_from_day
  real(kind = fp) :: slowness_x, slowness_y, daz, az_tmp, likelihood_tmp, ot_diff, origintime, &
  &                  likelihood_azweight, likelihood_distweight, map_x, map_y, rnd
  real(kind = fp), allocatable :: appvel_obs(:), az_obs(:), lon_array(:), lat_array(:), min_correlation(:), arrivaltime(:), &
  &                               ot_est(:), dist_tmp(:)
  integer,         allocatable :: arrayindex(:)
  logical,         allocatable :: result_exist_org(:)

  real(kind = fp)                   :: width_tmp(1 : 2), height_tmp(1 : 2), dwidth, dheight
  character(len = 255)              :: coastline_txt
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 2)                :: yr, mo, dy, hh, mm, ss
  
  
  real(kind = fp)      :: lon_particle(1 : nparticle_plot), lat_particle(1 : nparticle_plot), &
  &                       likelihood(1 : nparticle_plot, 1 : 3), &
  &                       az_weight(1 : int(2.0_fp * pi / daz_weight) + 1), &
  &                       kahan_val1(1 : 3), kahan_val2(1 : 3), sum_likelihood(1 : 3)

  !!initiate random number generator
  call make_seed(seed)
  call random_generator_init(random_status, seed)
  print '(a, i0)', "seed = ", seed

  !!Read coastline
  call getarg(1, coastline_txt)
  call read_coastline(coastline_txt, mapbuf)
  ncoastline = ubound(mapbuf, 1)

  !!Plot legend
  !call pc_plotinit(iwin_legend, "Legend", 0.0_sp, -300.0_sp, width / 2, 27.0_sp, scale)
  !call plot_legend

  !!Open epicenter window

  !!Read and plot AELUMA results
  call pc_plotinit(iwin_likelihood(1), "Total likelihoood distribution",      0.0_sp, 0.0_sp, width, height, scale)
  call pc_plotinit(iwin_likelihood(2), "Azimuth likelihoood distribution",    0.0_sp, 0.0_sp, width, height, scale)
  call pc_plotinit(iwin_likelihood(3), "Origintime likelihoood distribution", 0.0_sp, 0.0_sp, width, height, scale)
  do i = 1, 3
    call pc_setbkcolor(iwin_likelihood(i), 255, 255, 255)
  enddo
  call mercator(center_lon, lon_w, lat_s, width_tmp(1), height_tmp(1))
  call mercator(center_lon, lon_e, lat_n, width_tmp(2), height_tmp(2))
  dwidth = 1.0_fp / (width_tmp(2) - width_tmp(1))
  dheight = 1.0_fp / (height_tmp(2) - height_tmp(1))

  !!read AELUMA results from stdin
  do 
    do i = 1, 3
      call pc_clear(iwin_likelihood(i))
      call pc_setcolor(iwin_likelihood(i), 0, 0, 0)
    enddo
 
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
    narray_tmp = narray

    if(.not. allocated(arrayindex)) then
      allocate(az_obs(1 : ntriangle), appvel_obs(1 : ntriangle), result_exist_org(1 : ntriangle), &
      &        arrayindex(1 : ntriangle), lon_array(1 : ntriangle), lat_array(1 : ntriangle), &
      &        min_correlation(1 : ntriangle), arrivaltime(1 : ntriangle), dist_tmp(1 : ntriangle))
    endif


    do i = 1, narray
      read(*, *) arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
      &          slowness_x, slowness_y, min_correlation(arrayindex(i)), arrivaltime(arrayindex(i))
      result_exist_org(arrayindex(i))  = .true.
      if(min_correlation(arrayindex(i)) .le. correlation_threshold) narray_tmp = narray_tmp - 1

      !!arrival time: relative time in s from current time
      arrivaltime(arrayindex(i)) = -(real(nsec_buf, kind = fp) - arrivaltime(arrayindex(i)))
      az_obs(arrayindex(i)) = atan2(slowness_x, slowness_y)
      if(az_obs(arrayindex(i)) .lt. 0.0_fp) az_obs(arrayindex(i)) = az_obs(arrayindex(i)) + 2.0_fp * pi
      appvel_obs(arrayindex(i)) = 1.0_fp / sqrt(slowness_x ** 2 + slowness_y ** 2)
    enddo
    if(narray_tmp .eq. 0) cycle

    !!calculate azimuthal weighting array
    az_weight(1 : int(2.0_fp * pi / daz_weight) + 1) = 0.0_fp
    do i = 1, narray
      if(min_correlation(arrayindex(i)) .gt. correlation_threshold) then
        az_weight(int(az_obs(arrayindex(i)) / daz_weight) + 1) = az_weight(int(az_obs(arrayindex(i)) / daz_weight) + 1) + 1.0_fp
      endif
    enddo


    kahan_val1(1 : 3) = 0.0_fp
    sum_likelihood(1 : 3) = 0.0_fp
    allocate(ot_est(1 : narray))
    do j = 1, nparticle
      call gen_random_number(random_status, rnd)
      lon_particle(j) = lon_w + (lon_e - lon_w) * rnd
      call gen_random_number(random_status, rnd)
      lat_particle(j) = lat_s + (lat_n - lat_s) * rnd

      likelihood(j, 2) = 0.0_fp !!az
      likelihood(j, 3) = 0.0_fp !!ot
      narray_use_tmp = 0
      do i = 1, narray
        if(min_correlation(arrayindex(i)) .le. correlation_threshold) cycle

        call greatcircle_dist(lat_array(arrayindex(i)),           lon_array(arrayindex(i)), &
        &                     lat_particle(j),                    lon_particle(j), &
        &                     distance = dist_tmp(arrayindex(i)), azimuth = az_tmp)
        az_tmp = az_tmp + pi
        if(az_tmp .ge. 2.0_fp * pi) az_tmp = az_tmp - 2.0_fp * pi
        daz = delta_az(az_obs(arrayindex(i)), az_tmp)
        az_weight_index = int(az_obs(arrayindex(i)) / daz_weight) + 1
        likelihood_azweight = likelihood_weight(azweight_coef, sameaz_num2, az_weight(az_weight_index))
        likelihood_tmp      = likelihood_modified(daz, daz_weight2, likelihood_azweight)
        likelihood(j, 2)    = likelihood_renew(likelihood(j, 2), likelihood_tmp)
 
        narray_use_tmp = narray_use_tmp + 1
        ot_est(narray_use_tmp) = origintime_cal(arrivaltime(arrayindex(i)), dist_tmp(i), appvel_obs(arrayindex(i)))
      enddo

      !!calculate origin time
      call bubblesort(ot_est(1 : narray_use_tmp))
      call pickup_medianval(ot_est(1 : narray_use_tmp), origintime)
      narray_use_tmp = 0
      do i = 1, narray
        if(min_correlation(arrayindex(i)) .le. correlation_threshold) cycle
        narray_use_tmp = narray_use_tmp + 1
        ot_diff = origintime - ot_est(narray_use_tmp)
        likelihood_distweight = likelihood_weight(ot_coef, sigma_dist2, dist_tmp(i))
        likelihood_tmp        = likelihood_modified(ot_diff, sigma_otdiff2, likelihood_distweight)
        likelihood(j, 3)      = likelihood_renew(likelihood(j, 3), likelihood_tmp)
      enddo

      likelihood(j, 1) = likelihood(j, 2) * likelihood(j, 3)
        
      do i = 1, 3
        kahan_val1(i) = kahan_val1(i) + likelihood(j, i)
        kahan_val2(i) = sum_likelihood(i)
        sum_likelihood(i) = sum_likelihood(i) + kahan_val1(i)
        kahan_val2(i) = sum_likelihood(i) - kahan_val2(i)
        kahan_val1(i) = kahan_val1(i) - kahan_val2(i)
      enddo

    enddo
    deallocate(ot_est)

    do j = 1, 3
      do i = 1, nparticle
        likelihood(i, j) = likelihood(i, j) / sum_likelihood(j) * coefficient
        call mercator(center_lon, lon_particle(i), lat_particle(i), map_x, map_y)
        plot_x  = real((map_x  - width_tmp(1))  * dwidth,  kind = sp) * width
        plot_y  = real((map_y  - height_tmp(1)) * dheight, kind = sp) * height
        if(likelihood(i, j) .le. 0.02_fp) then
          color_index = 1
        elseif(likelihood(i, j) .gt. 0.02_fp .and. likelihood(i, j) .lt. 0.05_fp) then
          color_index = 2
        elseif(likelihood(i, j) .gt. 0.05_fp .and. likelihood(i, j) .lt. 0.1_fp) then
          color_index = 3 
        elseif(likelihood(i, j) .gt. 0.1_fp .and. likelihood(i, j) .lt. 0.2_fp) then
          color_index = 4
        elseif(likelihood(i, j) .gt. 0.2_fp .and. likelihood(i, j) .lt. 0.5_fp) then
          color_index = 5 
        elseif(likelihood(i, j) .gt. 0.5_fp .and. likelihood(i, j) .lt. 1.0_fp) then
          color_index = 6 
        elseif(likelihood(i, j) .gt. 1.0_fp .and. likelihood(i, j) .lt. 2.0_fp) then
          color_index = 7
        elseif(likelihood(i, j) .gt. 2.0_fp .and. likelihood(i, j) .lt. 5.0_fp) then
          color_index = 8
        elseif(likelihood(i, j) .gt. 5.0_fp .and. likelihood(i, j) .lt. 10.0_fp) then
          color_index = 9
        elseif(likelihood(i, j) .gt. 10.0_fp) then
          color_index = 10
        endif
        call pc_setcolor(iwin_likelihood(j), color(1, color_index), color(2, color_index), color(3, color_index))
        call pc_symbol(iwin_likelihood(j), plot_x, plot_y, 3.0_sp, 1, 0)
        call pc_setcolor(iwin_likelihood(j), 0, 0, 0)
        call pc_symbol(iwin_likelihood(j), plot_x, plot_y, 3.0_sp, 1, 1)
      enddo
      !!plot slowness vector
      call plot_slowness_vector(iwin_likelihood(j), narray, arrayindex, result_exist_org, lon_array, lat_array, az_obs, &
      &                         min_correlation, width_tmp, height_tmp, dwidth, dheight)
      !!plot map
      call pc_setline(iwin_likelihood(j), 1)
      call pc_setcolor(iwin_likelihood(j), 0, 0, 0)
      call plot_currentdate(iwin_likelihood(j), yr, mo, dy, hh, mm, ss, width_tmp, height_tmp, dwidth, dheight)
      call pc_setcolor(iwin_likelihood(j), 0, 0, 0)
      call plot_coastline(iwin_likelihood(j), ncoastline, mapbuf, width_tmp, height_tmp, dwidth, dheight)
      call pc_flush(iwin_likelihood(j))
    enddo
      
  enddo

  do i = 1, 3
    call pc_plotend(iwin_likelihood(i), 1)
  enddo

  stop
end program plot_likelihood

