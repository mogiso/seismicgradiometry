program plot_map_vector
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  implicit none

  real(kind = sp), parameter :: width = 300.0_sp, height = 300.0_sp, scale = 1.0_sp
  real(kind = sp), parameter :: vector_len = 5.0, vector_width = 1.5, vector_head1 = 2.5, vector_head2 = 4.0
  real(kind = fp), parameter :: lon_w = 120.0_fp, lon_e = 149.0_fp, lat_s = 22.5_fp, lat_n = 48.0_fp, center_lon = 135.0_fp
  integer, parameter :: iwin = 0

  integer         :: i, j, k, ios, ncoastline, mapcount, narray, arrayindex, ntriangle, color(1 : 3)
  real(kind = fp) :: width_min, width_max, height_min, height_max, dwidth, dheight, maplon, maplat, map_x, map_y, &
  &                  map_x1, map_y1
  real(kind = fp), allocatable :: slowness_x(:), slowness_y(:), lon_array(:), lat_array(:), min_correlation(:)
  logical, allocatable :: result_exist(:)
  real(kind = sp) :: plot_x, plot_y, plot_x1, plot_y1, plot_theta
  character(len = 255) :: coastline_txt, mapbuf_tmp
  character(len = 255), allocatable :: mapbuf(:)
  character(len = 2) :: yr, mo, dy, hh, mm, ss
  character(len = 19) :: date_txt
  
  !!location estimation
  integer, parameter   :: nparticle = 1000, niter = 5
  real(kind = fp), parameter :: sigma_particle = 0.1_fp
  integer, allocatable :: seed(:)
  integer              :: seedsize
  real(kind = fp)      :: cos_similarity(1 : nparticle), lon_particle(1 : nparticle), lat_particle(1 : nparticle), &
  &                       lon_particle_new(1 : nparticle), lat_particle_new(1 : nparticle), particle_probability(1 : nparticle)
  real(kind = fp)      :: rnd, rnd1, rnd2, az, normalize_cos_similarity, cos_similarity_tmp

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
 

  call pc_plotinit(iwin, "AELUMA results", 0.0, 0.0, width, height, scale)
  call pc_setbkcolor(iwin, 255, 255, 255)
  call mercator(center_lon, lon_w, lat_s, width_min, height_min)
  call mercator(center_lon, lon_e, lat_n, width_max, height_max)
  dwidth = 1.0_fp / (width_max - width_min)
  dheight = 1.0_fp / (height_max - height_min)


  !!read AELUMA results from stdin
  do 
    call pc_clear(iwin)
    call pc_setcolor(iwin, 0, 0, 0)
    call pc_setline(iwin, 1)
    read(*, *, iostat = ios) yr, mo, dy, hh, mm, ss, narray, ntriangle
    if(ios .ne. 0) stop
    if(.not. allocated(slowness_x)) then
      allocate(slowness_x(1 : ntriangle), slowness_y(1 : ntriangle), result_exist(1 : ntriangle), &
      &        lon_array(1 : ntriangle), lat_array(1 : ntriangle), min_correlation(1 : ntriangle))
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
      do i = 1, narray
        read(*, *) arrayindex, lon_array(arrayindex), lat_array(arrayindex), &
        &          slowness_x(arrayindex), slowness_y(arrayindex), min_correlation(arrayindex)
        if(.not. (slowness_x(arrayindex) .ne. 0.0_fp .and. slowness_y(arrayindex) .ne. 0.0_fp)) cycle
        result_exist(arrayindex) = .true.
        !theta = atan2(slowness_x, slowness_y) * rad2deg
        plot_theta = real(atan2(slowness_y(arrayindex), slowness_x(arrayindex)) * rad2deg, kind = sp)
        call mercator(center_lon, lon_array(arrayindex), lat_array(arrayindex), map_x, map_y)
        plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
        plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
        if(min_correlation(arrayindex) .lt. 0.2_fp) then
          color(1 : 3) = [220, 204, 222]
        elseif(min_correlation(arrayindex) .ge. 0.2_fp .and. min_correlation(arrayindex) .lt. 0.4_fp) then
          color(1 : 3) = [212, 156, 189]
        elseif(min_correlation(arrayindex) .ge. 0.4_fp .and. min_correlation(arrayindex) .lt. 0.6_fp) then
          color(1 : 3) = [196, 110, 155]
        elseif(min_correlation(arrayindex) .ge. 0.6_fp .and. min_correlation(arrayindex) .lt. 0.8_fp) then
          color(1 : 3) = [136, 97, 141]
        elseif(min_correlation(arrayindex) .ge. 0.8_fp) then
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
        do k = 1, niter
          !!calculate cosine-similarity
          do j = 1, nparticle
            cos_similarity(j) = 0.0_fp
            particle_probability(j) = 0.0_fp
            do i = 1, ntriangle
              if(result_exist(i) .eqv. .false.) cycle
              call greatcircle_dist(lat_array(i), lon_array(i), lat_particle(j), lon_particle(j), azimuth = az)
              az = az + pi
              cos_similarity_tmp =(slowness_x(i) * sin(az) + slowness_y(i) * cos(az)) &
              &                  / sqrt(slowness_x(i) ** 2 + slowness_y(i) ** 2)
              if(cos_similarity_tmp .le. cos(pi * 0.25_fp)) cos_similarity_tmp = 0.0_fp
              cos_similarity(j) = cos_similarity(j) + cos_similarity_tmp * min_correlation(i) * min_correlation(i)
            enddo
          enddo
          if(sum(cos_similarity) .eq. 0.0_fp) exit
          normalize_cos_similarity = 1.0_fp / sum(cos_similarity)
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
            call random_number(rnd1)
            call random_number(rnd2)
            rnd = sqrt(-2.0_fp * log(rnd1)) * cos(2.0_fp * pi * rnd2)
            lon_particle_new(j) = lon_particle(i) + rnd * sigma_particle
            rnd = sqrt(-2.0_fp * log(rnd1)) * sin(2.0_fp * pi * rnd2)
            lat_particle_new(j) = lat_particle(i) + rnd * sigma_particle
          enddo
          lon_particle(1 : nparticle) = lon_particle_new(1 : nparticle)
          lat_particle(1 : nparticle) = lat_particle_new(1 : nparticle)
        enddo
          
         
        !!write circles
        call pc_setline(iwin, 1)
        do i = 1, nparticle
          call mercator(center_lon, lon_particle(i), lat_particle(i), map_x, map_y)
          plot_x  = real((map_x  - width_min)  * dwidth,  kind = sp) * width
          plot_y  = real((map_y  - height_min) * dheight, kind = sp) * height
          call pc_symbol(iwin, plot_x, plot_y, 3.0_sp, 1, 1)
        enddo
        
     
      endif
    endif
    call pc_flush(iwin)
  enddo


  call pc_plotend(iwin, 1)

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
   
