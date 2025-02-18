program plot_map_vector
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  implicit none

  real(kind = sp), parameter :: likelihood_legend_normalize = 1.0e+3

  integer         :: i, j, k, ios, ncoastline, mapcount, narray, ntriangle, color(1 : 3), maxloc_likelihood(1), az_weight_index
  integer         :: year, month, day, hr, mi, sc, julianday, sec_from_day
  real(kind = fp) :: width_min, width_max, height_min, height_max, dwidth, dheight, maplon, maplat, map_x, map_y, &
  &                  map_x1, map_y1
  real(kind = fp), allocatable :: slowness_x(:), slowness_y(:), lon_array(:), lat_array(:), min_correlation(:), arrivaltime(:)
  integer,         allocatable :: arrayindex(:)
  logical, allocatable :: result_exist(:)
  real(kind = sp) :: plot_x, plot_y, plot_x1, plot_y1, plot_theta
  character(len = 255) :: coastline_txt, mapbuf_tmp
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 19) :: date_txt
  character(len = 6)  :: plottext
  
  integer, allocatable :: seed(:)
  integer              :: seedsize
  real(kind = fp)      :: lon_particle(1 : nparticle), lat_particle(1 : nparticle), &
  &                       likelihood_particle(1 : nparticle), &
  &                       lon_particle_new(1 : nparticle), lat_particle_new(1 : nparticle), particle_probability(1 : nparticle), &
  &                       az_weight(1 : int(2.0_fp * pi / daz_weight))
  real(kind = fp)      :: rnd, rnd1, rnd2, normalize_likelihood, maxval_likelihood_particle, daz, likelihood_tmp, &
  &                       kahan_val1, kahan_val2, sum_likelihood, likelihood_azweight
  real(kind = fp), allocatable :: az(:), az_obs(:)

  call random_seed(size = seedsize)
  allocate(seed(1 : seedsize))
  do i = 1, seedsize
    call system_clock(count = seed(i))
  enddo
  call random_seed(put = seed(:))

  !!read AELUMA results from stdin
  do 
    az_weight(1 : int(2.0_fp * pi / daz_weight)) = 0.0_fp
    read(*, *, iostat = ios) yr_tmp, jday_tmp, sec_from_day_tmp, narray, ntriangle
    if(ios .ne. 0) stop
    if(.not. allocated(slowness_x)) then
      allocate(propagation_azimuth(1 : ntriangle, 1 : ntrg), app_velocity(1 : ntriangle, 1 : ntrg), &
      &        result_exist(1 : ntriangle, 1 : ntrg), arrayindex(1 : ntriangle), &
      &        trg_sec_from_current(1 : ntriangle, 1 : ntrg))
    endif

    if(i .lt. 1) cycle
    do i = 1, narray
      read(*, *) arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
      &          slowness_x_tmp, slowness_y_tmp, min_correlation_tmp, 


    if(narray .ge. 1) then
      call pc_setline(iwin, 4)
      result_exist(1 : ntriangle) = .false.
      !!read and plot slowness vector
      do i = 1, narray
        read(*, *) arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
        &          slowness_x(arrayindex(i)), slowness_y(arrayindex(i)), min_correlation(arrayindex(i)), &
        &          arrivaltime(arrayindex(i))
        if(.not. (slowness_x(arrayindex(i)) .ne. 0.0_fp .and. slowness_y(arrayindex(i)) .ne. 0.0_fp)) cycle
        if(min_correlation(arrayindex(i)) .ge. correlation_threshold) then
          az_obs(i) = atan2(slowness_x(arrayindex(i)), slowness_y(arrayindex(i)))
          if(az_obs(i) .lt. 0.0_fp) az_obs(i) = az_obs(i) + 2.0_fp * pi
          az_weight(int(az_obs(i) / daz_weight) + 1) = az_weight(int(az_obs(i) / daz_weight) + 1) + 1.0_fp
        endif
        result_exist(arrayindex(i)) = .true.
        print '(i0, 6(1x, e15.7))', arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
        &                           slowness_x(arrayindex(i)), slowness_y(arrayindex(i)), min_correlation(arrayindex(i)), &
        &                           arrivaltime(arrayindex(i))
      enddo

      !!estimate location
      if(narray .ge. 5) then
        !!initial particle
        do i = 1, nparticle
          call random_number(rnd)
          lon_particle(i) = lon_w + (lon_e - lon_w) * rnd   
          call random_number(rnd)
          lat_particle(i) = lat_s + (lat_n - lat_s) * rnd
        enddo
        maxval_likelihood_particle = 0.0_fp
        particlefilter: do k = 1, niter
          !write(0, '(a, i0)') "iter num = ", k
          sum_likelihood = 0.0_fp
          kahan_val1 = 0.0_fp
          particleloop: do j = 1, nparticle
            likelihood_particle(j) = 0.0_fp
            do i = 1, narray
              if(.not. result_exist(arrayindex(i))) cycle
              if(min_correlation(arrayindex(i)) .lt. correlation_threshold) cycle
              call greatcircle_dist(lat_array(arrayindex(i)), lon_array(arrayindex(i)), &
              &                     lat_particle(j), lon_particle(j), azimuth = az(i))
              az(i) = az(i) + pi
              if(az(i) .ge. 2.0_fp * pi) az(i) = az(i) - 2.0_fp * pi
              daz = az_obs(i) - az(i)
              if(daz .gt.  pi) daz = 2.0_fp * pi - daz
              if(daz .lt. -pi) daz = 2.0_fp * pi + daz
              az_weight_index = int(az_obs(i) / daz_weight) + 1
              likelihood_azweight = 1.0_fp - 0.9_fp * exp(-(real(az_weight(az_weight_index), kind = fp) ** 2) / sameaz_num2)
              likelihood_tmp = exp(-(daz ** 2) * 0.5_fp / daz_weight2)
              likelihood_tmp = (1.0_fp - likelihood_azweight) * likelihood_tmp + likelihood_azweight
              if(likelihood_particle(j) .eq. 0.0_fp) then
                likelihood_particle(j) = likelihood_tmp
              else
                likelihood_particle(j) = likelihood_particle(j) * likelihood_tmp
              endif
            enddo
            kahan_val1 = kahan_val1 + likelihood_particle(j)
            kahan_val2 = sum_likelihood
            sum_likelihood = sum_likelihood + kahan_val1
            kahan_val2 = sum_likelihood - kahan_val2
            kahan_val1 = kahan_val1 - kahan_val2
          enddo particleloop
          if(sum_likelihood .le. 1.0e-100_fp) exit particlefilter
          normalize_likelihood = 1.0_fp / sum_likelihood
          !print *, maxval(likelihood_particle)
          if(k .eq. niter) exit particlefilter
          !print *, "iter num = ", k, maxval_likelihood_particle, maxval(likelihood_particle)
          if(maxval(likelihood_particle) .le. maxval_likelihood_particle) exit particlefilter
          maxval_likelihood_particle = maxval(likelihood_particle)

          particle_probability(1) = likelihood_particle(1) * normalize_likelihood
          do i = 2, nparticle
            particle_probability(i) = particle_probability(i - 1) + (likelihood_particle(i) * normalize_likelihood)
          enddo


          !!redistribute particle
          do j = 1, nparticle
            call random_number(rnd)
            do i = 1, nparticle
              if(rnd .le. particle_probability(i)) exit
            enddo
            if(i .gt. nparticle) then
              do i = 1, nparticle
                write(0, *) i, particle_probability(i), likelihood_particle(i), normalize_likelihood, sum_likelihood
              enddo
            endif
            call random_number(rnd1)
            call random_number(rnd2)
            rnd = sqrt(-2.0_fp * log(rnd1)) * cos(2.0_fp * pi * rnd2)
            lon_particle_new(j) = lon_particle(i) + rnd * sigma_particle
            call random_number(rnd1)
            call random_number(rnd2)
            rnd = sqrt(-2.0_fp * log(rnd1)) * sin(2.0_fp * pi * rnd2)
            lat_particle_new(j) = lat_particle(i) + rnd * sigma_particle
          enddo
          lon_particle(1 : nparticle) = lon_particle_new(1 : nparticle)
          lat_particle(1 : nparticle) = lat_particle_new(1 : nparticle)
        enddo particlefilter
         
        !!write circles
        call pc_setline(iwin, 1)
        do i = 1, nparticle
          call mercator(center_lon, lon_particle(i), lat_particle(i), map_x, map_y)
          plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
          plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
          likelihood_tmp = likelihood_particle(i) * normalize_likelihood * likelihood_legend_normalize
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
          call pc_setcolor(iwin, color(1), color(2), color(3))
          call pc_symbol(iwin, plot_x, plot_y, 3.0_sp, 1, 0)
          call pc_setcolor(iwin, 0, 0, 0)
          call pc_symbol(iwin, plot_x, plot_y, 3.0_sp, 1, 1)
        enddo
        maxloc_likelihood = maxloc(likelihood_particle)
        call mercator(center_lon, lon_particle(maxloc_likelihood(1)), lat_particle(maxloc_likelihood(1)), map_x, map_y)
        plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
        plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
        call pc_symbol(iwin, plot_x, plot_y, 9.0_sp, 1, 0)
        call pc_setcolor(iwin, 255, 255, 255)
        call pc_symbol(iwin, plot_x, plot_y, 4.0_sp, 1, 0)
      endif
    endif

    !!plot slowness vector
    call pc_setline(iwin, 4)
    do i = 1, narray
      if(.not. result_exist(arrayindex(i))) cycle
      !theta = atan2(slowness_x, slowness_y) * rad2deg
      plot_theta = real(atan2(slowness_y(arrayindex(i)), slowness_x(arrayindex(i))) * rad2deg, kind = sp)
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
      call pc_setcolor(iwin, color(1), color(2), color(3))
      call pc_vector(iwin, plot_x, plot_y, plot_theta, vector_len, vector_width, vector_head1, vector_head2, 1)
    enddo
    call pc_setline(iwin, 1)

 
    !!plot map
    call pc_setcolor(iwin, 0, 0, 0)
    date_txt = "20" // yr // "/" // mo // "/" // dy // " " // hh // ":" // mm // ":" // ss
    call mercator(center_lon, lon_w, lat_n, map_x, map_y)
    plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
    call pc_text(iwin, plot_x, plot_y, 7.0, date_txt, 0.0, len(date_txt), 7)

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
      call pc_line(iwin, plot_x, plot_y, plot_x1, plot_y1)
    enddo
    call pc_flush(iwin)
  enddo


  call pc_plotend(iwin, 1)
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
   
