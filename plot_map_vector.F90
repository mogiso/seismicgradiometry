program plot_map_vector
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  implicit none

  real(kind = sp), parameter :: width = 300.0_sp, height = 300.0_sp, scale = 1.0_sp
  real(kind = sp), parameter :: vector_len = 5.0, vector_width = 1.5, vector_head1 = 2.5, vector_head2 = 4.0
  real(kind = fp), parameter :: lon_w = 120.0_fp, lon_e = 149.0_fp, lat_s = 22.5_fp, lat_n = 48.0_fp, center_lon = 135.0_fp
  integer, parameter :: iwin = 0, iwin_legend = 1

  integer         :: i, j, k, ios, ncoastline, mapcount, narray, ntriangle, color(1 : 3), max_similarity(1)
  real(kind = fp) :: width_min, width_max, height_min, height_max, dwidth, dheight, maplon, maplat, map_x, map_y, &
  &                  map_x1, map_y1
  real(kind = fp), allocatable :: slowness_x(:), slowness_y(:), lon_array(:), lat_array(:), min_correlation(:)
  integer,         allocatable :: arrayindex(:)
  logical, allocatable :: result_exist(:)
  real(kind = sp) :: plot_x, plot_y, plot_x1, plot_y1, plot_theta
  character(len = 255) :: coastline_txt, mapbuf_tmp
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 2) :: yr, mo, dy, hh, mm, ss
  character(len = 19) :: date_txt
  character(len = 6)  :: plottext
  
  !!location estimation
  integer, parameter   :: nparticle = 1000, niter = 6
  real(kind = fp), parameter :: daz_weight = 10.0_fp * deg2rad
  real(kind = fp), parameter :: sigma_particle = 0.1_fp
  real(kind = fp), parameter :: cos_similarity_accept_degree = 5.0_fp * deg2rad
  real(kind = fp), parameter :: correlation_threshold = 0.5_fp
  integer, allocatable :: seed(:)
  integer              :: seedsize
  real(kind = fp)      :: cos_similarity(1 : nparticle), lon_particle(1 : nparticle), lat_particle(1 : nparticle), &
  &                       lon_particle_new(1 : nparticle), lat_particle_new(1 : nparticle), particle_probability(1 : nparticle), &
  &                       az_weight(1 : int(2.0_fp * pi / daz_weight))
  real(kind = fp)      :: rnd, rnd1, rnd2, normalize_cos_similarity, cos_similarity_tmp
  real(kind = fp), allocatable :: az(:)

  call random_seed(size = seedsize)
  allocate(seed(1 : seedsize))
  do i = 1, seedsize
    call system_clock(count = seed(i))
  enddo
  call random_seed(put = seed(:))

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
 

  call pc_plotinit(iwin, "AELUMA results", 0.0_sp, 0.0_sp, width, height, scale)
  call pc_setbkcolor(iwin, 255, 255, 255)
  call mercator(center_lon, lon_w, lat_s, width_min, height_min)
  call mercator(center_lon, lon_e, lat_n, width_max, height_max)
  dwidth = 1.0_fp / (width_max - width_min)
  dheight = 1.0_fp / (height_max - height_min)

  call pc_plotinit(iwin_legend, "Legend", 0.0_sp, -300.0_sp, width / 2, 27.0_sp, scale)
  call pc_setbkcolor(iwin_legend, 255, 255, 255)
  do i = 1, 5
    if(i .eq. 1) then
      color(1 : 3) = [220, 204, 222]
      plottext = "0.0  "
    elseif(i .eq. 2) then
      color(1 : 3) = [212, 156, 189]
      plottext = "0.2  "
    elseif(i .eq. 3) then
      color(1 : 3) = [196, 110, 155]
      plottext = "0.4  "
    elseif(i .eq. 4) then
      color(1 : 3) = [136, 97, 141]
      plottext = "0.6  "
    elseif(i .eq. 5) then
      color(1 : 3) = [73, 57, 100]
      plottext = "0.8  "
    endif
    plot_x = real(i, kind = sp) * 16.0_sp
    plot_y = 20.0_sp
    call pc_setcolor(iwin_legend, color(1), color(2), color(3))
    call pc_setline(iwin_legend, 4)
    call pc_vector(iwin_legend, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
    call pc_setcolor(iwin_legend, 0, 0, 0)
    plot_x = real(i, kind = sp) * 16.0_sp - 10.0_sp
    call pc_text(iwin_legend, plot_x, plot_y, 4.5, trim(plottext), 0.0, len(trim(plottext)), 4) 
  enddo
  plot_x = real(i, kind = sp) * 16.0_sp - 10.0_sp
  plottext = "1.0  "
  call pc_text(iwin_legend, plot_x, plot_y, 4.5, trim(plottext), 0.0, len(trim(plottext)), 4) 

  do i = 1, 6
    if(i .eq. 1) then
      color(1 : 3) = [252, 238, 158]
      plottext = "0.0  "
    elseif(i .eq. 2) then
      color(1 : 3) = [238, 179, 87]
      plottext = "0.2  "
    elseif(i .eq. 3) then
      color(1 : 3) = [222, 117, 79]
      plottext = "0.4  "
    elseif(i .eq. 4) then
      color(1 : 3) = [149, 66, 62]
      plottext = "0.6  "
    elseif(i .eq. 5) then
      color(1 : 3) = [63, 39, 23]
      plottext = "0.8  "
    elseif(i .eq. 6) then
      color(1 : 3) = [26, 26, 1]
      plottext = "1.0  "
    endif
    plot_x = real(i, kind = sp) * 16.0_sp
    plot_y = 10.0_sp
    call pc_setcolor(iwin_legend, color(1), color(2), color(3))
    call pc_symbol(iwin_legend, plot_x, plot_y, 3.0_sp, 1, 0)
    call pc_setcolor(iwin_legend, 0, 0, 0)
    call pc_symbol(iwin_legend, plot_x, plot_y, 3.0_sp, 1, 1)
    plot_x = real(i, kind = sp) * 16.0_sp - 8.0_sp
    call pc_text(iwin_legend, plot_x, plot_y, 4.5, trim(plottext), 0.0, len(trim(plottext)), 5) 
  enddo
  plot_x = real(i, kind = sp) * 16.0_sp - 8.0_sp
  plottext = "1.0< "
  call pc_text(iwin_legend, plot_x, plot_y, 4.5, trim(plottext), 0.0, len(trim(plottext)), 5) 
  plot_x = plot_x + 20.0_sp
  plottext = "x3e+2"
  call pc_text(iwin_legend, plot_x, plot_y, 4.5, trim(plottext), 0.0, len(trim(plottext)), 5) 
  call pc_flush(iwin_legend)
    
  

  !!read AELUMA results from stdin
  do 
    call pc_clear(iwin)
    call pc_setcolor(iwin, 0, 0, 0)
    call pc_setline(iwin, 1)
    read(*, *, iostat = ios) yr, mo, dy, hh, mm, ss, narray, ntriangle
    if(ios .ne. 0) stop
    if(.not. allocated(slowness_x)) then
      allocate(slowness_x(1 : ntriangle), slowness_y(1 : ntriangle), result_exist(1 : ntriangle), arrayindex(1 : ntriangle), &
      &        lon_array(1 : ntriangle), lat_array(1 : ntriangle), min_correlation(1 : ntriangle), az(1 : ntriangle))
    endif
    date_txt = "20" // yr // "/" // mo // "/" // dy // " " // hh // ":" // mm // ":" // ss
    call mercator(center_lon, lon_w, lat_n, map_x, map_y)
    plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
    plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
    call pc_text(iwin, plot_x, plot_y, 7.0, date_txt, 0.0, len(date_txt), 7)


    !!writing map
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

    if(narray .ge. 1) then
      call pc_setline(iwin, 4)
      result_exist(1 : ntriangle) = .false.
      !!read and plot slowness vector
      do i = 1, narray
        read(*, *) arrayindex(i), lon_array(arrayindex(i)), lat_array(arrayindex(i)), &
        &          slowness_x(arrayindex(i)), slowness_y(arrayindex(i)), min_correlation(arrayindex(i))
        if(.not. (slowness_x(arrayindex(i)) .ne. 0.0_fp .and. slowness_y(arrayindex(i)) .ne. 0.0_fp)) cycle
        result_exist(arrayindex(i)) = .true.
        !theta = atan2(slowness_x, slowness_y) * rad2deg
        plot_theta = real(atan2(slowness_y(arrayindex(i)), slowness_x(arrayindex(i))) * rad2deg, kind = sp)
        call mercator(center_lon, lon_array(arrayindex(i)), lat_array(arrayindex(i)), map_x, map_y)
        plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
        plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
        if(min_correlation(arrayindex(i)) .lt. 0.2_fp) then
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
      call pc_setcolor(iwin, 0, 0, 0)

      !!estimate location
      if(narray .ge. 5) then
        !!initial particle
        do i = 1, nparticle
          call random_number(rnd)
          lon_particle(i) = lon_w + (lon_e - lon_w) * rnd   
          call random_number(rnd)
          lat_particle(i) = lat_s + (lat_n - lat_s) * rnd
        enddo
        particlefilter: do k = 1, niter
          !!calculate cosine-similarity
          do j = 1, nparticle
            cos_similarity(j) = 0.0_fp
            particle_probability(j) = 0.0_fp
            az_weight(1 : int(2.0_fp * pi / daz_weight)) = 0.0_fp
            do i = 1, narray
              if(min_correlation(i) .lt. correlation_threshold) cycle
              call greatcircle_dist(lat_array(arrayindex(i)), lon_array(arrayindex(i)), &
              &                     lat_particle(j), lon_particle(j), azimuth = az(i))
              az(i) = az(i) + pi
              if(az(i) .ge. 2.0_fp * pi) az(i) = az(i) - 2.0_fp * pi
              az_weight(int(az(i) / daz_weight) + 1) = az_weight(int(az(i) / daz_weight) + 1) + 1.0_fp
            enddo
            do i = 1, narray
              if(min_correlation(i) .lt. correlation_threshold) cycle
              cos_similarity_tmp = (slowness_x(arrayindex(i)) * sin(az(i)) + slowness_y(arrayindex(i)) * cos(az(i))) &
              &              / sqrt(slowness_x(arrayindex(i)) ** 2 + slowness_y(arrayindex(i)) ** 2)
              if(cos_similarity_tmp .le. cos(cos_similarity_accept_degree)) cos_similarity_tmp = 0.0_fp
              cos_similarity(j) = cos_similarity(j) &
              &                 + cos_similarity_tmp / az_weight(int(az(i) / daz_weight) + 1) * min_correlation(arrayindex(i))
            enddo
          enddo
          if(sum(cos_similarity) .eq. 0.0_fp) exit
          normalize_cos_similarity = 1.0_fp / sum(cos_similarity(1 : nparticle))
          if(k .eq. niter) exit particlefilter
          particle_probability(1) = cos_similarity(1) * normalize_cos_similarity
          do i = 2, nparticle
            particle_probability(i) = particle_probability(i - 1) + cos_similarity(i) * normalize_cos_similarity
          enddo
          !!check
          !write(0, *) "max_particle_prob = ", particle_probability(nparticle)

          !!redistribute particle
          do j = 1, nparticle
            call random_number(rnd)
            do i = 1, nparticle
              if(rnd .le. particle_probability(i)) exit
            enddo
            if(i .gt. nparticle) then
              write(0, *) i, rnd, particle_probability(i - 1), normalize_cos_similarity
              do i = 1, narray
                write(0, *) az(i), min_correlation(arrayindex(i))
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
          cos_similarity_tmp = cos_similarity(i) * normalize_cos_similarity * 3.0e+2_fp
          if(cos_similarity_tmp .le. 0.2_fp) then
            color(1 : 3) = [252, 238, 158]
          elseif(cos_similarity_tmp .gt. 0.2_fp .and. cos_similarity_tmp .le. 0.4_fp) then
            color(1 : 3) = [238, 179, 87]
          elseif(cos_similarity_tmp .gt. 0.4_fp .and. cos_similarity_tmp .le. 0.6_fp) then
            color(1 : 3) = [222, 117, 79]
          elseif(cos_similarity_tmp .gt. 0.6_fp .and. cos_similarity_tmp .le. 0.8_fp) then
            color(1 : 3) = [149, 66, 62]
          elseif(cos_similarity_tmp .gt. 0.8_fp .and. cos_similarity_tmp .le. 1.0_fp) then
            color(1 : 3) = [63, 39, 23]
          elseif(cos_similarity_tmp .gt. 1.0_fp) then
            color(1 : 3) = [26, 26, 1]
          endif
          call pc_setcolor(iwin, color(1), color(2), color(3))
          call pc_symbol(iwin, plot_x, plot_y, 3.0_sp, 1, 0)
          call pc_setcolor(iwin, 0, 0, 0)
          call pc_symbol(iwin, plot_x, plot_y, 3.0_sp, 1, 1)
        enddo
        max_similarity = maxloc(cos_similarity)
        call mercator(center_lon, lon_particle(max_similarity(1)), lat_particle(max_similarity(1)), map_x, map_y)
        plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
        plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
        call pc_symbol(iwin, plot_x, plot_y, 9.0_sp, 1, 0)
        call pc_setcolor(iwin, 255, 255, 255)
        call pc_symbol(iwin, plot_x, plot_y, 4.0_sp, 1, 0)
        
     
      endif
    endif
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
   
