! Copyright 2024 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

module calc_kernelmatrix

  private
  public :: calc_slowness_est_matrix_delaunay, calc_slowness_est_matrix_delaunay_shmdump

contains

subroutine calc_slowness_est_matrix_delaunay(location_sta, nadd_station, ntriangle, &
&                                            triangle_center, slowness_matrix,      &
&                                            triangle_stationindex, nsta_count, tnbr)
  use nrtype, only : fp
  use constants, only : pi
  use typedef
  use aeluma_parameters
  use lonlat_xy_conv, only : xy2bl
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif
  implicit none

  type(location),  intent(inout) :: location_sta(:)
  integer,         intent(in)    :: nadd_station
  integer,         intent(out)   :: ntriangle
  type(location),  intent(out), allocatable :: triangle_center(:)
  real(kind = fp), intent(out), allocatable :: slowness_matrix(:, :, :)
  integer,         intent(out), allocatable :: triangle_stationindex(:, :), nsta_count(:), tnbr(:, :)

  integer         :: i, j, ii, jj, info, nsta, nsta_use, npair, npair_tmp
  integer         :: ipiv(1 : 3)
  real(kind = fp) :: dist_tmp
  real(kind = fp), allocatable :: vertices(:, :), add_station_distance(:), g2(:, :), g_tmp2(:, :)
  integer,         allocatable :: vertix_index(:), triangle_indices(:, :), index_org(:), add_station_index(:)
  logical,         allocatable :: is_usestation(:), used_station(:)
 
  nsta = size(location_sta)
  nsta_use = nsta
  allocate(is_usestation(1 : nsta))
  is_usestation(1 : nsta) = .true.

  !!check interstation distance
  do j = 1, nsta - 1
    if(is_usestation(j) .eqv. .false.) cycle
    do i = j + 1, nsta
      if(is_usestation(i) .eqv. .false.) cycle
      call cartesian_dist(location_sta(i)%x_east,  location_sta(j)%x_east, &
      &                   location_sta(i)%y_north, location_sta(j)%y_north, &
      &                   distance = dist_tmp)
      if(dist_tmp .le. interstationdistance_min) then
        is_usestation(i) = .false.
        nsta_use = nsta_use - 1
      endif
    enddo
  enddo
  

  !!Do delaunay triangulation
  allocate(vertix_index(1 : nsta_use),                &
  &        vertices(1 : 2, 1 : nsta_use),             &
  &        triangle_indices(1 : 3, 1 : 2 * nsta_use), &
  &        tnbr(1 : 3, 1 : 2 * nsta_use),             &
  &        index_org(1 : nsta),                       &
  &        used_station(1 : nsta))

  j = 1
  do i = 1, nsta
    if(is_usestation(i) .eqv. .false.) cycle
    vertices(1, j) = location_sta(i)%x_east
    vertices(2, j) = location_sta(i)%y_north
    vertix_index(j) =j 
    index_org(j) = i
    j = j + 1
  enddo
  call dtris2(nsta_use, vertices, vertix_index, ntriangle, triangle_indices, tnbr, info)

  open(unit = 10, file = "station_triangle.txt")
  do j = 1, ntriangle
    do i = 1, 3
      write(10, '(2(e15.7, 1x))') location_sta(index_org(triangle_indices(i, j)))%lon, &
      &                           location_sta(index_org(triangle_indices(i, j)))%lat
    enddo
    write(10, '(a)') ">"
  enddo
  close(10)

  allocate(triangle_center(1 : ntriangle),                        &
  &        nsta_count(1 : ntriangle),                             &
  &        triangle_stationindex(1 : 3 + nadd_station, 1 : ntriangle))
  do j = 1, ntriangle
    used_station(1 : nsta) = .false.
    triangle_center(j)%x_east  = 0.0_fp
    triangle_center(j)%y_north = 0.0_fp
    nsta_count(j) = 0
    do i = 1, 3
      used_station(index_org(triangle_indices(i, j))) = .true.
      triangle_center(j)%x_east = triangle_center(j)%x_east + location_sta(index_org(triangle_indices(i, j)))%x_east
      triangle_center(j)%y_north = triangle_center(j)%y_north + location_sta(index_org(triangle_indices(i, j)))%y_north
      triangle_stationindex(i, j) = index_org(triangle_indices(i, j))
      nsta_count(j) = nsta_count(j) + 1
    enddo
    triangle_center(j)%x_east  = triangle_center(j)%x_east  / 3.0_fp
    triangle_center(j)%y_north = triangle_center(j)%y_north / 3.0_fp
    call xy2bl(triangle_center(j)%y_north * 1000.0_fp, triangle_center(j)%x_east * 1000.0_fp, &
    &          center_lon_aeluma, center_lat_aeluma, &
    &          triangle_center(j)%lon, triangle_center(j)%lat)

    !!find nadd_station additional stations based on the distance between grid and station
    if(nadd_station .ge. 1) then
      allocate(add_station_distance(1 : nadd_station), add_station_index(1 : nadd_station))
      add_station_distance(1 : nadd_station) = 1.0e+38
      add_station_index(1 : nadd_station) = 0
      do jj = 1, nsta
        if(is_usestation(jj) .eqv. .false.) cycle
        if(used_station(jj)  .eqv. .true. ) cycle
        call cartesian_dist(location_sta(jj)%x_east,  triangle_center(j)%x_east,  &
        &                   location_sta(jj)%y_north, triangle_center(j)%y_north, &
        &                   distance = dist_tmp)
        do ii = 1, nadd_station
          if(dist_tmp .le. add_station_distance(ii)) then
            do i = nadd_station, ii + 1, -1
              add_station_distance(i) = add_station_distance(i - 1)
              add_station_index(i) = add_station_index(i - 1)
            enddo
            add_station_distance(ii) = dist_tmp
            add_station_index(ii) = jj
            exit
          endif
        enddo
      enddo
      do i = 1, nadd_station
        nsta_count(j) = nsta_count(j) + 1
        triangle_stationindex(nsta_count(j), j) = add_station_index(i)
      enddo
      deallocate(add_station_distance, add_station_index)
    endif

    !!check distance between the center of triangle and stations (vertices)
    do i = 1, nsta_count(j)
      call cartesian_dist(location_sta(triangle_stationindex(i, j))%x_east,  triangle_center(j)%x_east,  &
      &                   location_sta(triangle_stationindex(i, j))%y_north, triangle_center(j)%y_north, &
      &                   distance = dist_tmp)
      if(dist_tmp .gt. cutoff_dist) then
        nsta_count(j) = 0
        exit
      endif
    enddo
  enddo
  npair = 0
  do i = 1, maxval(nsta_count) - 1
    npair = npair + i
  enddo

  allocate(slowness_matrix(1 : 2, 1 : npair, 1 : ntriangle))
  do jj = 1, ntriangle
    npair_tmp = 0
    do ii = 1, nsta_count(jj) - 1
      npair_tmp = npair_tmp + ii
    enddo
    slowness_matrix(1 : 2, 1 : npair_tmp, jj) = 0.0_fp
    if(nsta_count(jj) .eq. 0) cycle
    allocate(g2(1 : npair_tmp, 1 : 2), g_tmp2(1 : 2, 1 : 2))
    ii = 1
    do j = 1, nsta_count(jj) - 1
      do i = j + 1, nsta_count(jj)
        call cartesian_dist(location_sta(triangle_stationindex(i, jj))%x_east,  &
        &                   location_sta(triangle_stationindex(j, jj))%x_east,  &
        &                   location_sta(triangle_stationindex(i, jj))%y_north, &
        &                   location_sta(triangle_stationindex(j, jj))%y_north, &
        &                   dist_x = g2(ii, 1), dist_y = g2(ii, 2))
        ii = ii + 1
      enddo
    enddo

    g_tmp2 = matmul(transpose(g2), g2)
#ifdef MKL
    call getrf(g_tmp2, ipiv = ipiv(1 : 2), info = info)
    call getri(g_tmp2, ipiv(1 : 2), info = info)
#else
    call LA_GETRF(g_tmp2, ipiv(1 : 2), info = info)
    call LA_GETRI(g_tmp2, ipiv(1 : 2), info = info)
#endif
    slowness_matrix(1 : 2, 1 : npair_tmp, jj) = matmul(g_tmp2, transpose(g2))
    deallocate(g_tmp2, g2)
  enddo

  deallocate(is_usestation, index_org, used_station, vertix_index, vertices, triangle_indices)
  return
end subroutine calc_slowness_est_matrix_delaunay

subroutine calc_slowness_est_matrix_delaunay_shmdump(location_sta, station_winch, nadd_station, &
&                                                    ntriangle, triangle_center, slowness_matrix, &
&                                                    triangle_stationwinch, nsta_count, tnbr)
  use nrtype, only : fp
  use constants, only : pi
  use typedef
  use lonlat_xy_conv, only : bl2xy
  use aeluma_parameters
  use greatcircle
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif
  implicit none

  type(location),  intent(inout)            :: location_sta(:)
  integer,         intent(in)               :: station_winch(:)
  integer,         intent(in)               :: nadd_station
  integer,         intent(out)              :: ntriangle
  type(location),  intent(out), allocatable :: triangle_center(:)
  real(kind = fp), intent(out), allocatable :: slowness_matrix(:, :, :)
  integer,         intent(out), allocatable :: triangle_stationwinch(:, :), nsta_count(:), tnbr(:, :)

  integer                      :: i, j, ii, jj, info, nsta, nsta_use, npair, npair_tmp
  integer                      :: ipiv(1 : 3)
  real(kind = fp)              :: dist_tmp
  real(kind = fp), allocatable :: vertices(:, :), add_station_distance(:), g2(:, :), g_tmp2(:, :)
  integer,         allocatable :: vertix_index(:), triangle_indices(:, :), index_org(:), add_station_index(:)
  logical,         allocatable :: is_usestation(:), used_station(:)
 
  nsta = ubound(station_winch, 1)
  nsta_use = nsta
  allocate(is_usestation(1 : nsta))
  is_usestation(1 : nsta) = .true.

  !!check interstation distance
  do j = 1, nsta - 1
    do i = j + 1, nsta
      if(is_usestation(i) .eqv. .false.) cycle
      call greatcircle_dist(location_sta(station_winch(j))%lat, location_sta(station_winch(j))%lon, &
      &                     location_sta(station_winch(i))%lat, location_sta(station_winch(i))%lon, &
      &                     distance = dist_tmp)
      if(dist_tmp .le. interstationdistance_min) then
        is_usestation(i) = .false.
        nsta_use = nsta_use - 1
      endif
    enddo
  enddo
  

  !!Do delaunay triangulation
  allocate(vertix_index(1 : nsta_use),                &
  &        vertices(1 : 2, 1 : nsta_use),             &
  &        triangle_indices(1 : 3, 1 : 2 * nsta_use), &
  &        tnbr(1 : 3, 1 : 2 * nsta_use),             &
  &        index_org(1 : nsta),                       &
  &        used_station(1 : nsta))
  i = 1
  do j = 1, nsta
    if(is_usestation(j) .eqv. .false.) cycle
    vertices(1 : 2, i) = [location_sta(station_winch(j))%x_east, location_sta(station_winch(j))%y_north]
    vertix_index(i) = i
    index_org(i) = j
    i = i + 1
  enddo
  call dtris2(nsta_use, vertices, vertix_index, ntriangle, triangle_indices, tnbr, info)

  allocate(triangle_center(1 : ntriangle), nsta_count(1 : ntriangle), &
  &        triangle_stationwinch(1 : 3 + nadd_station, 1 : ntriangle))
  do j = 1, ntriangle
    used_station(1 : nsta) = .false.
    triangle_center(j)%lon = 0.0_fp
    triangle_center(j)%lat = 0.0_fp
    nsta_count(j) = 0
    do i = 1, 3
      triangle_stationwinch(i, j) = station_winch(index_org(triangle_indices(i, j)))
      used_station(index_org(triangle_indices(i, j))) = .true.
      triangle_center(j)%lon = triangle_center(j)%lon + location_sta(triangle_stationwinch(i, j))%lon
      triangle_center(j)%lat = triangle_center(j)%lat + location_sta(triangle_stationwinch(i, j))%lat
      nsta_count(j) = nsta_count(j) + 1
    enddo
    triangle_center(j)%lon = triangle_center(j)%lon / 3.0_fp
    triangle_center(j)%lat = triangle_center(j)%lat / 3.0_fp

    !!find nadd_station additional stations based on the distance between grid and station
    if(nadd_station .ge. 1) then
      allocate(add_station_distance(1 : nadd_station), add_station_index(1 : nadd_station))
      add_station_distance(1 : nadd_station) = 1.0e+38_fp
      add_station_index(1 : nadd_station) = 0
      do jj = 1, nsta
        if(is_usestation(jj) .eqv. .false.) cycle
        if(used_station(jj)  .eqv. .true. ) cycle
        call greatcircle_dist(triangle_center(j)%lat,              triangle_center(j)%lon, &
        &                     location_sta(station_winch(jj))%lat, location_sta(station_winch(jj))%lon, &
        &                     distance = dist_tmp)
        do ii = 1, nadd_station
          if(dist_tmp .le. add_station_distance(ii)) then
            do i = nadd_station, ii + 1, -1
              add_station_distance(i) = add_station_distance(i - 1)
              add_station_index(i)    = add_station_index(i - 1)
            enddo
            add_station_distance(ii) = dist_tmp
            add_station_index(ii)    = jj
            exit
          endif
        enddo
      enddo
      do i = 1, nadd_station
        nsta_count(j) = nsta_count(j) + 1
        triangle_stationwinch(nsta_count(j), j) = station_winch(add_station_index(i))
      enddo
      deallocate(add_station_distance, add_station_index)

      triangle_center(j)%lon = 0.0_fp
      triangle_center(j)%lat = 0.0_fp
      do i = 1, nsta_count(j)
        triangle_center(j)%lon = triangle_center(j)%lon + location_sta(triangle_stationwinch(i, j))%lon
        triangle_center(j)%lat = triangle_center(j)%lat + location_sta(triangle_stationwinch(i, j))%lat
      enddo
      triangle_center(j)%lon = triangle_center(j)%lon / real(nsta_count(j), kind = fp)
      triangle_center(j)%lat = triangle_center(j)%lat / real(nsta_count(j), kind = fp)
    endif
    
    !!check distance between the center of triangle and stations (vertices)
    do i = 1, nsta_count(j)
      call greatcircle_dist(triangle_center(j)%lat, triangle_center(j)%lon, &
      &                     location_sta(triangle_stationwinch(i, j))%lat, location_sta(triangle_stationwinch(i, j))%lon, &
      &                     distance = dist_tmp)
      if(dist_tmp .gt. cutoff_dist) then
        nsta_count(j) = 0
        exit
      endif
    enddo
  enddo
  npair = 0
  do i = 1, maxval(nsta_count) - 1
    npair = npair + i
  enddo

  allocate(slowness_matrix(1 : 2, 1 : npair, 1 : ntriangle))
  do jj = 1, ntriangle
    if(nsta_count(jj) .eq. 0) cycle
    npair_tmp = 0
    do ii = 1, nsta_count(jj) - 1
      npair_tmp = npair_tmp + ii
    enddo
    slowness_matrix(1 : 2, 1 : npair_tmp, jj) = 0.0_fp
    allocate(g2(1 : npair_tmp, 1 : 2), g_tmp2(1 : 2, 1 : 2))
    ii = 1
    do j = 1, nsta_count(jj) - 1
      call bl2xy(location_sta(triangle_stationwinch(j, jj))%lon,     location_sta(triangle_stationwinch(j, jj))%lat, &
      &          triangle_center(jj)%lon,                            triangle_center(jj)%lat, &
      &          location_sta(triangle_stationwinch(j, jj))%y_north, location_sta(triangle_stationwinch(j, jj))%x_east)
      location_sta(triangle_stationwinch(j, jj))%x_east  = location_sta(triangle_stationwinch(j, jj))%x_east  * 1.0e-3_fp
      location_sta(triangle_stationwinch(j, jj))%y_north = location_sta(triangle_stationwinch(j, jj))%y_north * 1.0e-3_fp

      do i = j + 1, nsta_count(jj)
        call bl2xy(location_sta(triangle_stationwinch(i, jj))%lon,     location_sta(triangle_stationwinch(i, jj))%lat, &
        &          triangle_center(jj)%lon,                            triangle_center(jj)%lat, &
        &          location_sta(triangle_stationwinch(i, jj))%y_north, location_sta(triangle_stationwinch(i, jj))%x_east)
        location_sta(triangle_stationwinch(i, jj))%x_east  = location_sta(triangle_stationwinch(i, jj))%x_east  * 1.0e-3_fp
        location_sta(triangle_stationwinch(i, jj))%y_north = location_sta(triangle_stationwinch(i, jj))%y_north * 1.0e-3_fp

        call cartesian_dist(location_sta(triangle_stationwinch(i, jj))%x_east,  &
        &                   location_sta(triangle_stationwinch(j, jj))%x_east,  &
        &                   location_sta(triangle_stationwinch(i, jj))%y_north, &
        &                   location_sta(triangle_stationwinch(j, jj))%y_north, &
        &                   dist_x = g2(ii, 1), dist_y = g2(ii, 2))
        ii = ii + 1
      enddo
    enddo
    g_tmp2 = matmul(transpose(g2), g2)
#ifdef MKL
    call getrf(g_tmp2, ipiv = ipiv(1 : 2), info = info)
    call getri(g_tmp2, ipiv(1 : 2),        info = info)
#else
    call LA_GETRF(g_tmp2, ipiv(1 : 2), info = info)
    call LA_GETRI(g_tmp2, ipiv(1 : 2), info = info)
#endif
    slowness_matrix(1 : 2, 1 : npair_tmp, jj) = matmul(g_tmp2, transpose(g2))
    deallocate(g_tmp2, g2)
  enddo

  open(unit = 10, file = "station_aelumaarray.txt")
  do j = 1, ntriangle
    if(nsta_count(j) .lt. 3 + nadd_station) cycle
    do i = 1, nsta_count(j)
      write(10, '(4(e15.7, 1x), z4)') location_sta(triangle_stationwinch(i, j))%lon, &
      &                               location_sta(triangle_stationwinch(i, j))%lat, &
      &                               triangle_center(j)%lon, triangle_center(j)%lat, &
      &                               triangle_stationwinch(i, j)
    enddo
    write(10, '(2(e15.7, 1x))') location_sta(triangle_stationwinch(1, j))%lon, &
    &                           location_sta(triangle_stationwinch(1, j))%lat
    write(10, '(a)') ">"
  enddo
  close(10)

  deallocate(is_usestation, index_org, used_station, vertix_index, vertices, triangle_indices)
  return
end subroutine calc_slowness_est_matrix_delaunay_shmdump

    
subroutine cartesian_dist(x_east2, x_east1, y_north2, y_north1, theta, dist_x, dist_y, distance)
  use nrtype, only : fp
  implicit none
  real(kind = fp), intent(in) :: x_east2, x_east1, y_north2, y_north1
  real(kind = fp), intent(in), optional :: theta          !!in radian, measured in clockwise direction from north 
  real(kind = fp), intent(out), optional :: dist_x, dist_y, distance
  real(kind = fp) :: dist_x_tmp, dist_y_tmp

  dist_x_tmp = x_east2 - x_east1
  dist_y_tmp = y_north2 - y_north1

  if(present(dist_x)) then
    if(present(theta)) then
      dist_x = dist_x_tmp * cos(theta) + dist_y_tmp * sin(theta)
    else
      dist_x = dist_x_tmp 
    endif
  endif
  if(present(dist_y)) then
    if(present(theta)) then
      dist_y = -dist_x_tmp * sin(theta) + dist_y_tmp * cos(theta)
    else
      dist_y = dist_y_tmp 
    endif
  endif

  if(present(distance)) distance = sqrt(dist_x_tmp * dist_x_tmp + dist_y_tmp * dist_y_tmp)
  
  return
end subroutine cartesian_dist

end module calc_kernelmatrix
