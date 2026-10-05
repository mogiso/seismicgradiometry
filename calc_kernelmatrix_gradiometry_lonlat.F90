module calc_kernelmatrix_lonlat
  private
  public :: cartesian_dist, calc_kernelmatrix_delaunay2, calculate_kernelmatrix  
  
  contains
  
  subroutine cartesian_dist(x_east2, x_east1, y_north2, y_north1, dist_x, dist_y, distance)
    use nrtype, only : fp
    implicit none
    real(kind = fp), intent(in)            :: x_east2, x_east1, y_north2, y_north1
    real(kind = fp), intent(out), optional :: dist_x, dist_y, distance
    real(kind = fp) :: dist_x_tmp, dist_y_tmp
  
    dist_x_tmp = x_east2  - x_east1
    dist_y_tmp = y_north2 - y_north1
  
    if(present(dist_x)) dist_x = dist_x_tmp 
    if(present(dist_y)) dist_y = dist_y_tmp 
    if(present(distance)) distance = sqrt(dist_x_tmp * dist_x_tmp + dist_y_tmp * dist_y_tmp)
    
    return
  end subroutine cartesian_dist
  
  
  subroutine calc_kernelmatrix_delaunay2(paramfile, location_grid, location_sta, &
  &                                      grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix, error_matrix)
    use nrtype, only : fp
    use constants, only : pi
    use typedef
    use greatcircle
    use m_readini
    implicit none
  
    character(len = 129), intent(in)            :: paramfile
    type(location),       intent(in)            :: location_grid(:, :), location_sta(:)
    logical,              intent(out)           :: grid_enough_sta(:, :)
    integer,              intent(out)           :: nsta_count(:, :), grid_stationindex(:, :, :)
    real(kind = fp),      intent(out), optional :: kernel_matrix(:, :, :, :)
    real(kind = fp),      intent(out), optional :: error_matrix(:, :, :)
  
    real(kind = fp) :: interstationdistance_min, cutoff_dist
    integer         :: nsta_grid_min, nsta_grid_max, naddstation_array
    integer         :: unitnum, i, j, ii, jj, kk, info, ntriangle, nsta_use
    integer         :: ngrid_lon, ngrid_lat, nsta
    logical         :: is_inside
    real(kind = fp) :: point_tmp(1 : 2), triangle_vertix_tmp(1 : 2, 1 : 3), dist_tmp
    type(location), allocatable  :: location_sta_tmp(:)
    real(kind = fp), allocatable :: add_station_index(:), add_station_distance(:), vertices(:, :)
    integer,         allocatable :: vertix_index(:), triangle_indices(:, :), tnbr(:, :), index_org(:)
    logical,         allocatable :: is_usestation(:), used_station(:)
   
    open(newunit = unitnum, file = trim(paramfile))
    call readini__strict_mode(.true.)
    call readini(unitnum, "interstationdistance_min", interstationdistance_min)
    call readini(unitnum, "nsta_grid_min", nsta_grid_min)
    call readini(unitnum, "nsta_grid_max", nsta_grid_max)
    call readini(unitnum, "naddstation_array", naddstation_array)
    call readini(unitnum, "cutoff_dist", cutoff_dist)
    close(unitnum)
  
    ngrid_lon = ubound(location_grid, 1)
    ngrid_lat = ubound(location_grid, 2)
    nsta = ubound(location_sta, 1)
  
    nsta_use = nsta
    allocate(is_usestation(1 : nsta), index_org(1 : nsta), used_station(1 : nsta))
    is_usestation(1 : nsta) = .true.
    if(present(error_matrix)) error_matrix(1 : 3, 1 : ngrid_lon, 1 : ngrid_lat) = 0.0_fp
  
    !!check interstation distance
    do j = 1, nsta - 1
      if(is_usestation(j) .eqv. .false.) cycle
      do i = j + 1, nsta
        if(is_usestation(i) .eqv. .false.) cycle
         call greatcircle_dist(location_sta(i)%lat, location_sta(i)%lon, &
         &                     location_sta(j)%lat, location_sta(j)%lon, distance = dist_tmp)
        if(dist_tmp .le. interstationdistance_min) then
          is_usestation(i) = .false.
          nsta_use = nsta_use - 1
        endif
      enddo
    enddo
  
    !!Do delaunay triangulation
    allocate(vertix_index(1 : nsta_use), vertices(1 : 2, 1 : nsta_use), triangle_indices(1 : 3, 1 : 2 * nsta_use), &
    &        tnbr(1 : 3, 1 : 2 * nsta_use))
    j = 1
    do i = 1, nsta
      if(is_usestation(i) .eqv. .false.) cycle
      vertices(1, j) = location_sta(i)%x_east
      vertices(2, j) = location_sta(i)%y_north
      vertix_index(j) = j 
      index_org(j) = i
      j = j + 1
    enddo
    call dtris2(nsta_use, vertices, vertix_index, ntriangle, triangle_indices, tnbr, info)
  
    !!Select stations at each grid
    !open(newunit = unitnum, file = "stationlist_grid.txt")
    do kk = 1, ngrid_lat
      do jj = 1, ngrid_lon
        nsta_count(jj, kk) = 0
        grid_enough_sta(jj, kk) = .false.
        grid_stationindex(1 : nsta_grid_max, jj, kk) = 0
  
        used_station(1 : nsta) = .false.
        !!find triangle that contains the grid
        point_tmp(1 : 2) = [location_grid(jj, kk)%x_east, location_grid(jj, kk)%y_north]
        do j = 1, ntriangle
          do i = 1, 3
            triangle_vertix_tmp(1 : 2, i) = [location_sta(index_org(triangle_indices(i, j)))%x_east, &
            &                                location_sta(index_org(triangle_indices(i, j)))%y_north]
          enddo
          call triangle_contains_point_2d_3(triangle_vertix_tmp, point_tmp, is_inside)
          if(is_inside .eqv. .true.) then
            grid_enough_sta(jj, kk) = .true.
            nsta_count(jj, kk) = 3
            do i = 1, 3
              grid_stationindex(i, jj, kk) = index_org(triangle_indices(i, j))
              used_station(grid_stationindex(i, jj, kk)) = .true.
            enddo
            exit
          endif
        enddo
        if(grid_enough_sta(jj, kk) .eqv. .false.) cycle 
  
        !!find nadd_station additional stations based on the distance between grid and station
        if(naddstation_array .ge. 1) then
          allocate(add_station_distance(1 : naddstation_array), add_station_index(1 : naddstation_array))
          add_station_distance(1 : naddstation_array) = huge(1.0_fp)
          add_station_index(1 : naddstation_array) = 0
          do ii = 1, nsta
            if(is_usestation(ii) .eqv. .false.) cycle
            if(used_station(ii)  .eqv. .true.) cycle
            call greatcircle_dist(location_sta(ii)%lat, location_sta(ii)%lon, &
            &                     location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
            &                     distance = dist_tmp)
            if(dist_tmp .gt. cutoff_dist) cycle
            do j = 1, naddstation_array
              if(dist_tmp .le. add_station_distance(j)) then
                do i = naddstation_array, j + 1, -1
                  add_station_distance(i) = add_station_distance(i - 1)
                  add_station_index(i) = add_station_index(i - 1)
                enddo
                add_station_distance(j) = dist_tmp
                add_station_index(j) = ii
                exit
              endif
            enddo
          enddo
          do i = 1, naddstation_array
            nsta_count(jj, kk) = nsta_count(jj, kk) + 1
            grid_stationindex(nsta_count(jj, kk), jj, kk) = add_station_index(i)
          enddo
          deallocate(add_station_distance, add_station_index)
        endif
        if(nsta_count(jj, kk) .lt. nsta_grid_min) grid_enough_sta(jj, kk) = .false.
        if(nsta_count(jj, kk) .gt. nsta_grid_max) nsta_count(jj, kk) = nsta_grid_max
  
        !!check distance between grid and stations (vertices)
        do i = 1, nsta_count(jj, kk)
          call greatcircle_dist(location_sta(grid_stationindex(i, jj, kk))%lat, location_sta(grid_stationindex(i, jj, kk))%lon, &
          &                     location_grid(jj, kk)%lat,                      location_grid(jj, kk)%lon, &
          &                     distance = dist_tmp)
          if(dist_tmp .gt. cutoff_dist) grid_enough_sta(jj, kk) = .false.
        enddo
  
        if(grid_enough_sta(jj, kk) .eqv. .false.) cycle
  
        if(present(kernel_matrix) .and. present(error_matrix)) then
          location_sta_tmp(1 : nsta_count(jj, kk)) = location_sta(grid_stationindex(1 : nsta_count(jj, kk), jj, kk))
          call calculate_kernelmatrix(nsta_count(jj, kk), location_sta_tmp(1 : nsta_count(jj, kk)), location_grid(ii, jj), &
          &                           cutoff_dist, info, kernel_matrix(:, :, jj, kk), error_matrix(:, jj, kk))
          if(info .ne. 0) grid_enough_sta(jj, kk) = .false.
          deallocate(location_sta_tmp)
        endif
         
      enddo
    enddo
  
    deallocate(is_usestation, used_station, index_org, vertix_index, vertices, triangle_indices, tnbr)
  
    return
  end subroutine calc_kernelmatrix_delaunay2
  
  subroutine calculate_kernelmatrix(nsta_count, location_sta, location_grid, cutoff_dist, info, kernel_matrix, error_matrix)
    use nrtype, only : fp
    use typedef
#ifdef MKL
    use lapack95
#else
    use f95_lapack
#endif
    implicit none
  
    integer,         intent(in)  :: nsta_count
    type(location),  intent(in)  :: location_sta(1 : nsta_count), location_grid
    real(kind = fp), intent(in)  :: cutoff_dist
    integer,         intent(out) :: info
    real(kind = fp), intent(out) :: kernel_matrix(1 : 3, 1 : nsta_count)
    real(kind = fp), intent(out) :: error_matrix(1 : 3)
  
    integer :: i, ipiv(1 : 3)
    real(kind = fp) :: g(1 : nsta_count, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_count, 1 : nsta_count)
  
    g(1 : nsta_count, 1 : 3) = 0.0_fp
    weight(1 : nsta_count, 1 : nsta_count) = 0.0_fp
    do i = 1, nsta_count
      g(i, 1) = 1.0_fp
      call cartesian_dist(location_sta(i)%x_east,  location_grid%x_east, &
      &                   location_sta(i)%y_north, location_grid%y_north, &
      &                   dist_x = g(i, 2), dist_y = g(i, 3))
      weight(i, i) = exp(-(g(i, 2) ** 2 + g(i, 3) ** 2) / (cutoff_dist ** 2))
    enddo
    g_tmp = matmul(transpose(g), matmul(weight, g))
  
#ifdef MKL
    call getrf(g_tmp, ipiv = ipiv, info = info); call getri(g_tmp, ipiv, info = info)
#else
    call LA_GETRF(g_tmp, ipiv, info = info);     call LA_GETRI(g_tmp, ipiv, info = info)
#endif
    if(info .ne. 0) return
  
    kernel_matrix(1 : 3, 1 : nsta_count) = matmul(matmul(g_tmp, transpose(g)), weight)
    do i = 1, 3
       error_matrix(i) = g_tmp(i, i)
    enddo
    return
  end subroutine calculate_kernelmatrix

end module calc_kernelmatrix_lonlat 
  
