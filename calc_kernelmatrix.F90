! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

module calc_kernelmatrix

  private
  public :: calc_kernelmatrix_circle, calc_kernelmatrix_delaunay, calc_kernelmatrix_delaunay2

contains

subroutine calc_kernelmatrix_circle(location_grid, location_sta, &
&                                   grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)
  use nrtype, only : fp
  use constants, only : pi
  use typedef
  use gradiometry_parameters
  use greatcircle, only : greatcircle_dist
  use sort, only : bubblesort

#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif
  implicit none

  type(location),  intent(in)  :: location_grid(:, :), location_sta(:)
  logical,         intent(out) :: grid_enough_sta(:, :)
  integer,         intent(out) :: nsta_count(:, :), grid_stationindex(:, :, :)
  real(kind = fp), intent(out) :: kernel_matrix(:, :, :, :)

  integer         :: i, j, ii, jj, kk, info, nsta
  integer         :: ipiv(3)
  real(kind = fp) :: propagation_direction, dist_x_tmp, dist_y_tmp, dist_tmp, az_diff_tmp
  real(kind = fp) :: azimuth_grid_sta(1 : nsta_grid_max), dist_grid_sta(1 : nsta_grid_max), &
                     g(1 : nsta_grid_max, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_grid_max, 1 : nsta_grid_max)

  nsta = size(location_sta)

  do kk = 1, ngrid_y
    do jj = 1, ngrid_x
      nsta_count(jj, kk) = 0
      grid_enough_sta(jj, kk) = .false.
      grid_stationindex(1 : nsta_grid_max, jj, kk) = 0
      azimuth_grid_sta(1 : nsta_grid_max) = 1.0e+15

      !!count usable stations for each grid
#ifdef ELLIPSE
      dist_grid_sta(1 : nsta_grid_max) = 0.0_fp
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      do ii = 1, nsta
        call cartesian_dist(location_sta(ii)%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(ii)%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)
        dist_tmp = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
        if(dist_tmp .lt. exp(-2.0_fp)) cycle

        nsta_count(jj, kk) = nsta_count(jj, kk) + 1
        do j = 1, nsta_grid_max
          if(dist_tmp .ge. dist_grid_sta(j)) then
            do i = nsta_grid_max, j + 1, -1
              dist_grid_sta(i) = dist_grid_sta(i - 1)
              grid_stationindex(i, jj, kk) = grid_stationindex(i - 1, jj, kk)
              azimuth_grid_sta(i) = azimuth_grid_sta(i - 1)
            enddo
            dist_grid_sta(j) = dist_tmp
            grid_stationindex(j, jj, kk) = ii
            call greatcircle_dist(location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
            &                     location_sta(ii)%lat,      location_sta(ii)%lon, &
            &                     azimuth = azimuth_grid_sta(j))
            exit
          endif
        enddo
      enddo
      if(nsta_count(jj, kk) .ge. nsta_grid_min) grid_enough_sta(jj, kk) = .true.
      if(nsta_count(jj, kk) .gt. nsta_grid_max) nsta_count(jj, kk) = nsta_grid_max
#else
      dist_grid_sta(1 : nsta_grid_max) = 1.0e+15
      do ii = 1, nsta
        call cartesian_dist(location_sta(ii)%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(ii)%y_north, location_grid(jj, kk)%y_north, &
        &                   distance = dist_tmp)
        if(dist_tmp .gt. cutoff_dist) cycle
        nsta_count(jj, kk) = nsta_count(jj, kk) + 1
        do j = 1, nsta_grid_max
          if(dist_tmp .le. dist_grid_sta(j)) then
            do i = nsta_grid_max, j + 1, -1
              dist_grid_sta(i) = dist_grid_sta(i - 1)
              grid_stationindex(i, jj, kk) = grid_stationindex(i - 1, jj, kk)
              azimuth_grid_sta(i) = azimuth_grid_sta(i - 1)
            enddo
            dist_grid_sta(j) = dist_tmp
            grid_stationindex(j, jj, kk) = ii
            call greatcircle_dist(location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
            &                     location_sta(ii)%lat, location_sta(ii)%lon, &
            &                     azimuth = azimuth_grid_sta(j))
            exit
          endif
        enddo
      enddo
      if(nsta_count(jj, kk) .ge. nsta_grid_min) grid_enough_sta(jj, kk) = .true.
      if(nsta_count(jj, kk) .gt. nsta_grid_max) nsta_count(jj, kk) = nsta_grid_max
#endif

      call bubblesort(azimuth_grid_sta)
      az_diff_tmp = 2.0_fp * pi - (azimuth_grid_sta(1) - azimuth_grid_sta(nsta_count(jj, kk)))
      do i = 2, nsta_count(jj, kk)
        if(azimuth_grid_sta(i - 1) - azimuth_grid_sta(i) .gt. az_diff_tmp) then
          az_diff_tmp = azimuth_grid_sta(i - 1) - azimuth_grid_sta(i)
        endif
      enddo
      if(az_diff_tmp .gt. az_diff_max) grid_enough_sta(jj, kk) = .false.

      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle

      g(1 : nsta_grid_max, 1 : 3) = 0.0_fp
      weight(1 : nsta_grid_max, 1 : nsta_grid_max) = 0.0_fp
      do i = 1, nsta_count(jj, kk)
        g(i, 1) = 1.0_fp
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   dist_x = g(i, 2), dist_y = g(i, 3))

#ifdef ELLIPSE
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)

        weight(i, i) = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
#else
        weight(i, i) = exp(-(g(i, 2) ** 2 + g(i, 3) ** 2) / (cutoff_dist ** 2))
#endif

      enddo
      g_tmp = matmul(transpose(g), matmul(weight, g))

#ifdef MKL
      call getrf(g_tmp, ipiv = ipiv, info = info)
      call getri(g_tmp, ipiv, info = info)
#else
      call LA_GETRF(g_tmp, ipiv, info = info)
      call LA_GETRI(g_tmp, ipiv, info = info)
#endif

      if(info .ne. 0) then
        grid_enough_sta(jj, kk) = .false.
        cycle
      endif

      kernel_matrix(1 : 3, 1 : nsta_grid_max, jj, kk) = matmul(g_tmp, matmul(transpose(g), weight))
      !if(grid_enough_sta(jj, kk) .eqv. .true.) print *, jj, kk, nsta_count(jj, kk)

    enddo
  enddo

  return
end subroutine calc_kernelmatrix_circle
    

subroutine calc_kernelmatrix_delaunay(location_grid, location_sta, &
&                                   grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)
  use nrtype, only : fp
  use constants, only : pi
  use typedef
  use gradiometry_parameters
  use greatcircle, only : greatcircle_dist
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif
  implicit none

  type(location),  intent(in)  :: location_grid(:, :), location_sta(:)
  logical,         intent(out) :: grid_enough_sta(:, :)
  integer,         intent(out) :: nsta_count(:, :), grid_stationindex(:, :, :)
  real(kind = fp), intent(out) :: kernel_matrix(:, :, :, :)

  integer         :: i, j, ii, jj, kk, info, nsta, ntriangle, nsta_use
  integer         :: ipiv(1 : 3)
  logical         :: is_inside
  real(kind = fp) :: propagation_direction, dist_x_tmp, dist_y_tmp, dist_tmp, az_diff_tmp
  real(kind = fp) :: g(1 : nsta_grid_max, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_grid_max, 1 : nsta_grid_max), &
  &                  point_tmp(1 : 2), triangle_vertix_tmp(1 : 2, 1 : 3)
  integer, allocatable :: vertix_index(:), triangle_indices(:, :), tnbr(:, :), index_org(:)
  real(kind = fp), allocatable :: vertices(:, :)
  logical, allocatable :: is_usestation(:)
 
  nsta = size(location_sta)
  nsta_use = nsta
  allocate(is_usestation(1 : nsta), index_org(1 : nsta))
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
  allocate(vertix_index(1 : nsta_use), vertices(1 : 2, 1 : nsta_use), triangle_indices(1 : 3, 1 : 2 * nsta_use), &
  &        tnbr(1 : 3, 1 : 2 * nsta_use))
  ii = 1
  do i = 1, nsta
    if(is_usestation(i) .eqv. .false.) cycle
    vertices(1, ii) = location_sta(i)%x_east
    vertices(2, ii) = location_sta(i)%y_north
    vertix_index(ii) = ii
    index_org(ii) = i
    ii = ii + 1
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

  do kk = 1, ngrid_y
    do jj = 1, ngrid_x
      nsta_count(jj, kk) = 0
      grid_enough_sta(jj, kk) = .false.
      grid_stationindex(1 : nsta_grid_max, jj, kk) = 0

      !!find triangle that contains the grid
      point_tmp(1) = location_grid(jj, kk)%x_east
      point_tmp(2) = location_grid(jj, kk)%y_north
      do j = 1, ntriangle
        do i = 1, 3
          triangle_vertix_tmp(1, i) = location_sta(index_org(triangle_indices(i, j)))%x_east
          triangle_vertix_tmp(2, i) = location_sta(index_org(triangle_indices(i, j)))%y_north
        enddo
        call triangle_contains_point_2d_3(triangle_vertix_tmp, point_tmp, is_inside)
        if(is_inside .eqv. .true.) then
          grid_enough_sta(jj, kk) = .true.
          nsta_count(jj, kk) = 3
          do i = 1, 3
            grid_stationindex(i, jj, kk) = index_org(triangle_indices(i, j))
          enddo
          exit
        endif
      enddo
      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle 

#ifdef ELLIPSE
      !!Check whether all stations (vertices) are within the ellipse
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      do ii = 1, nsta_count(jj, kk)
        call cartesian_dist(location_sta(grid_stationindex(ii, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(ii, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)
        dist_tmp = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
        if(dist_tmp .lt. exp(-2.0_fp)) then
          grid_enough_sta(jj, kk) = .false.
          exit
        endif
      enddo
#else
      !!check distance between grid and stations (vertices)
      do ii = 1, nsta_count(jj, kk)
        call cartesian_dist(location_sta(grid_stationindex(ii, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(ii, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   distance = dist_tmp)
        if(dist_tmp .gt. cutoff_dist) then
          grid_enough_sta(jj, kk) = .false.
          exit
        endif
      enddo
#endif

      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle

      !!calculate kernel matrix for interpolation
      g(1 : nsta_grid_max, 1 : 3) = 0.0_fp
      weight(1 : nsta_grid_max, 1 : nsta_grid_max) = 0.0_fp
      do i = 1, nsta_count(jj, kk)
        g(i, 1) = 1.0_fp
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   dist_x = g(i, 2), dist_y = g(i, 3))

#ifdef ELLIPSE
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)
        weight(i, i) = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
#else
        weight(i, i) = exp(-(g(i, 2) ** 2 + g(i, 3) ** 2) / (cutoff_dist ** 2))
#endif

      enddo
      g_tmp = matmul(transpose(g), matmul(weight, g))

#ifdef MKL
      call getrf(g_tmp, ipiv = ipiv, info = info)
      call getri(g_tmp, ipiv, info = info)
#else
      call LA_GETRF(g_tmp, ipiv, info = info)
      call LA_GETRI(g_tmp, ipiv, info = info)
#endif

      if(info .ne. 0) then
        grid_enough_sta(jj, kk) = .false.
        cycle
      endif

      kernel_matrix(1 : 3, 1 : nsta_grid_max, jj, kk) = matmul(g_tmp, matmul(transpose(g), weight))

    enddo
  enddo

  deallocate(vertix_index, vertices, triangle_indices, tnbr, is_usestation)

  return
end subroutine calc_kernelmatrix_delaunay
    

subroutine calc_kernelmatrix_delaunay2(location_grid, location_sta, nadd_station, &
&                                   grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix, slowness_est_matrix)
  use nrtype, only : fp
  use constants, only : pi
  use typedef
  use gradiometry_parameters
  use greatcircle, only : greatcircle_dist
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif

  type(location),  intent(in)  :: location_grid(:, :), location_sta(:)
  integer,         intent(in)  :: nadd_station
  logical,         intent(out) :: grid_enough_sta(:, :)
  integer,         intent(out) :: nsta_count(:, :), grid_stationindex(:, :, :)
  real(kind = fp), intent(out) :: kernel_matrix(:, :, :, :)
  real(kind = fp), intent(out), optional :: slowness_est_matrix(:, :, :, :)

  integer         :: i, j, ii, jj, kk, info, nsta, ntriangle, nsta_use, ndata_slowness_max
  integer         :: ipiv(1 : 3)
  logical         :: is_inside
  real(kind = fp) :: propagation_direction, dist_x_tmp, dist_y_tmp, dist_tmp, az_diff_tmp
  real(kind = fp) :: g(1 : nsta_grid_max, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_grid_max, 1 : nsta_grid_max), &
  &                  point_tmp(1 : 2), triangle_vertix_tmp(1 : 2, 1 : 3)
  integer, allocatable :: vertix_index(:), triangle_indices(:, :), tnbr(:, :), index_org(:), add_station_index(:)
  real(kind = fp), allocatable :: vertices(:, :), add_station_distance(:), g2(:, :), g_tmp2(:, :)
  logical, allocatable :: is_usestation(:), used_station(:)
 
  nsta = size(location_sta)
  nsta_use = nsta
  allocate(is_usestation(1 : nsta), index_org(1 : nsta))
  is_usestation(1 : nsta) = .true.
  if(present(slowness_est_matrix)) ndata_slowness_max = size(slowness_est_matrix) / (2 * ngrid_x * ngrid_y)

  allocate(add_station_index(1 : nadd_station), add_station_distance(1 : nadd_station), used_station(1 : nsta))

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
  allocate(vertix_index(1 : nsta_use), vertices(1 : 2, 1 : nsta_use), triangle_indices(1 : 3, 1 : 2 * nsta_use), &
  &        tnbr(1 : 3, 1 : 2 * nsta_use))
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

  do kk = 1, ngrid_y
    do jj = 1, ngrid_x
      nsta_count(jj, kk) = 0
      grid_enough_sta(jj, kk) = .false.
      grid_stationindex(1 : nsta_grid_max, jj, kk) = 0

      used_station(1 : nsta) = .false.
      !!find triangle that contains the grid
      point_tmp(1) = location_grid(jj, kk)%x_east
      point_tmp(2) = location_grid(jj, kk)%y_north
      do j = 1, ntriangle
        do i = 1, 3
          triangle_vertix_tmp(1, i) = location_sta(index_org(triangle_indices(i, j)))%x_east
          triangle_vertix_tmp(2, i) = location_sta(index_org(triangle_indices(i, j)))%y_north
        enddo
        call triangle_contains_point_2d_3(triangle_vertix_tmp, point_tmp, is_inside)
        if(is_inside .eqv. .true.) then
          grid_enough_sta(jj, kk) = .true.
          nsta_count(jj, kk) = 3
          do i = 1, 3
            grid_stationindex(i, jj, kk) = index_org(triangle_indices(i, j))
            used_station(index_org(triangle_indices(i, j))) = .true.
          enddo
          exit
        endif
      enddo
      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle 

      !!find nadd_station additional stations based on the distance between grid and station
      add_station_distance(1 : nadd_station) = 1.0e+38
      add_station_index(1 : nadd_station) = 0
      do ii = 1, nsta
        if(is_usestation(ii) .eqv. .false.) cycle
        if(used_station(ii) .eqv. .false.) cycle
        call cartesian_dist(location_sta(ii)%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(ii)%y_north, location_grid(jj, kk)%y_north, &
        &                   distance = dist_tmp)
        do j = 1, nadd_station
          if(dist_tmp .le. add_station_distance(j)) then
            do i = nadd_station, j + 1, -1
              add_station_distance(i) = add_station_distance(i - 1)
              add_station_index(i) = add_station_index(i - 1)
            enddo
            add_station_distance(j) = dist_tmp
            add_station_index(j) = ii
            exit
          endif
        enddo
      enddo
      do i = 1, nadd_station
        nsta_count(jj, kk) = nsta_count(jj, kk) + 1
        grid_stationindex(nsta_count(jj, kk), jj, kk) = add_station_index(i)
      enddo

#ifdef ELLIPSE
      !!Check whether all stations (vertices) are within the ellipse
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      do i = 1, nsta_count(jj, kk)
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)
        dist_tmp = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
        if(dist_tmp .lt. exp(-2.0_fp)) then
          grid_enough_sta(jj, kk) = .false.
          exit
        endif
      enddo
#else
      !!check distance between grid and stations (vertices)
      do i = 1, nsta_count(jj, kk)
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   distance = dist_tmp)
        if(dist_tmp .gt. cutoff_dist) then
          grid_enough_sta(jj, kk) = .false.
          exit
        endif
      enddo
#endif

      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle

      !!calculate kernel matrix for interpolation
      g(1 : nsta_grid_max, 1 : 3) = 0.0_fp
      weight(1 : nsta_grid_max, 1 : nsta_grid_max) = 0.0_fp
      do i = 1, nsta_count(jj, kk)
        g(i, 1) = 1.0_fp
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   dist_x = g(i, 2), dist_y = g(i, 3))

#ifdef ELLIPSE
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)
        weight(i, i) = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
#else
        weight(i, i) = exp(-(g(i, 2) ** 2 + g(i, 3) ** 2) / (cutoff_dist ** 2))
#endif

      enddo
      g_tmp = matmul(transpose(g), matmul(weight, g))

#ifdef MKL
      call getrf(g_tmp, ipiv = ipiv, info = info)
      call getri(g_tmp, ipiv, info = info)
#else
      call LA_GETRF(g_tmp, ipiv, info = info)
      call LA_GETRI(g_tmp, ipiv, info = info)
#endif

      if(info .ne. 0) then
        grid_enough_sta(jj, kk) = .false.
        cycle
      endif

      kernel_matrix(1 : 3, 1 : nsta_grid_max, jj, kk) = matmul(g_tmp, matmul(transpose(g), weight))

      if(present(slowness_est_matrix)) then
        allocate(g2(1 : ndata_slowness_max, 1 : 2), g_tmp2(1 : 2, 1 : 2))
        slowness_est_matrix(1 : 2, 1 : ndata_slowness_max, jj, kk) = 0.0_fp
        ii = 1
        do i = 1, nsta_count(jj, kk) - 1
          do j = i + 1, nsta_count(jj, kk)
            call cartesian_dist(location_sta(grid_stationindex(j, jj, kk))%x_east, &
            &                   location_sta(grid_stationindex(i, jj, kk))%x_east, &
            &                   location_sta(grid_stationindex(j, jj, kk))%y_north, &
            &                   location_sta(grid_stationindex(i, jj, kk))%y_north, &
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
        slowness_est_matrix(1 : 2, 1 : ndata_slowness_max, jj, kk) = matmul(g_tmp2, transpose(g2))
        deallocate(g_tmp2, g2)
      endif


    enddo
  enddo

  deallocate(vertix_index, vertices, triangle_indices, tnbr, is_usestation)
  deallocate(add_station_distance, add_station_index, used_station)

  return
end subroutine calc_kernelmatrix_delaunay2
    
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
