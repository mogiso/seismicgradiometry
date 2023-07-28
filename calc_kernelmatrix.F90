! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

module calc_kernelmatrix

  private
  public :: calc_kernelmatrix_circle, calc_kernelmatrix_delaunay, &
  &         calc_kernelmatrix_delaunay2, calc_slowness_est_matrix_delaunay

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
#ifdef ELLIPSE
  real(kind = fp) :: propagation_direction, dist_x_tmp, dist_y_tmp, az_diff_tmp
#endif
  real(kind = fp) :: g(1 : nsta_grid_max, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_grid_max, 1 : nsta_grid_max), &
  &                  point_tmp(1 : 2), triangle_vertix_tmp(1 : 2, 1 : 3), dist_tmp
  integer, allocatable :: vertix_index(:), triangle_indices(:, :), tnbr(:, :), index_org(:)
  real(kind = fp), allocatable :: vertices(:, :)
  logical, allocatable :: is_usestation(:), used_triangle(:)
 
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
  allocate(used_triangle(1 : ntriangle))
  used_triangle(1 : ntriangle) = .false.

  open(unit = 10, file = "station_triangle.txt")
  !do j = 1, ntriangle
  !  do i = 1, 3
  !    write(10, '(2(e15.7, 1x))') location_sta(index_org(triangle_indices(i, j)))%x_east, &
  !    &                           location_sta(index_org(triangle_indices(i, j)))%y_north
  !  enddo
  !  write(10, '(a)') ">"
  !enddo
  !close(10)

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

      if(used_triangle(j) .eqv. .false.) then
        do i = 1, 3
          write(10, '(2(e15.7, 1x))') location_sta(index_org(triangle_indices(i, j)))%x_east, &
          &                           location_sta(index_org(triangle_indices(i, j)))%y_north
        enddo
        write(10, '(a)') ">"
        used_triangle(j) = .true.
      endif

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
  close(10)

  deallocate(vertix_index, vertices, triangle_indices, tnbr, is_usestation, used_triangle)

  return
end subroutine calc_kernelmatrix_delaunay
    

subroutine calc_kernelmatrix_delaunay2(location_grid, location_sta, nadd_station, &
&                                      grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix, &
&                                      nsta_correlation, slowness_est_matrix)
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
  integer,         intent(in)  :: nadd_station
  logical,         intent(out) :: grid_enough_sta(:, :)
  integer,         intent(out) :: nsta_count(:, :), grid_stationindex(:, :, :)
  real(kind = fp), intent(out) :: kernel_matrix(:, :, :, :)
  integer,         intent(out), optional :: nsta_correlation(:, :)
  real(kind = fp), intent(out), allocatable, optional :: slowness_est_matrix(:, :, :, :)

  integer         :: i, j, ii, jj, kk, info, nsta, ntriangle, nsta_use, nsta_count_tmp
  integer         :: ipiv(1 : 3)
  logical         :: is_inside
#ifdef ELLIPSE
  real(kind = fp) :: propagation_direction, dist_x_tmp, dist_y_tmp, az_diff_tmp
#endif
  real(kind = fp) :: g(1 : nsta_grid_max, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_grid_max, 1 : nsta_grid_max), &
  &                  point_tmp(1 : 2), triangle_vertix_tmp(1 : 2, 1 : 3), dist_tmp
  real(kind = fp), allocatable :: vertices(:, :), add_station_distance(:), g2(:, :), g_tmp2(:, :)
  integer,         allocatable :: vertix_index(:), triangle_indices(:, :), tnbr(:, :), index_org(:), add_station_index(:)
  logical,         allocatable :: is_usestation(:), used_station(:)
 
  nsta = size(location_sta)
  nsta_use = nsta
  allocate(is_usestation(1 : nsta), index_org(1 : nsta), used_station(1 : nsta))
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


  !!Select stations at each grid
  open(unit = 10, file = "stationlist_grid.txt")
  do kk = 1, ngrid_y
    do jj = 1, ngrid_x
      nsta_count(jj, kk) = 0
      grid_enough_sta(jj, kk) = .false.
      grid_stationindex(1 : nsta_grid_max, jj, kk) = 0

      used_station(1 : nsta) = .false.
      !!find triangle that contains the grid
      point_tmp(1 : 2) = (/location_grid(jj, kk)%x_east, location_grid(jj, kk)%y_north/)
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
      if(nadd_station .ge. 1) then
        allocate(add_station_distance(1 : nadd_station), add_station_index(1 : nadd_station))
        add_station_distance(1 : nadd_station) = 1.0e+38
        add_station_index(1 : nadd_station) = 0
        do ii = 1, nsta
          if(is_usestation(ii) .eqv. .false.) cycle
          if(used_station(ii) .eqv. .true.) cycle
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
        deallocate(add_station_distance, add_station_index)
      endif

      nsta_count_tmp = nsta_count(jj, kk)
#ifdef ELLIPSE
      !!Check whether all stations (vertices) are within the ellipse
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      do i = 1, nsta_count_tmp
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   theta = propagation_direction, dist_x = dist_x_tmp, dist_y = dist_y_tmp)
        dist_tmp = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
        if(dist_tmp .lt. exp(-2.0_fp)) then
          grid_enough_sta(jj, kk) = .false.
        endif
      enddo
#else
      !!check distance between grid and stations (vertices)
      do i = 1, nsta_count_tmp
        call cartesian_dist(location_sta(grid_stationindex(i, jj, kk))%x_east,  location_grid(jj, kk)%x_east, &
        &                   location_sta(grid_stationindex(i, jj, kk))%y_north, location_grid(jj, kk)%y_north, &
        &                   distance = dist_tmp)
        if(dist_tmp .gt. cutoff_dist) then
          grid_enough_sta(jj, kk) = .false.
        endif
      enddo
#endif

      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle

      do i = 1, nsta_count(jj, kk)
        write(10, '(2(e15.7, 1x))') location_sta(grid_stationindex(i, jj, kk))%x_east, &
        &                           location_sta(grid_stationindex(i, jj, kk))%y_north
      enddo
      write(10, '(a)') ">"

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

      if(present(slowness_est_matrix) .and. present(nsta_correlation)) then
        nsta_correlation(jj, kk) = 0
        do i = 1, nsta_count(jj, kk) - 1
          do j = i + 1, nsta_count(jj, kk)
            nsta_correlation(jj, kk) = nsta_correlation(jj, kk) + 1
          enddo
        enddo
      endif

    enddo
  enddo
  close(10)

  deallocate(is_usestation, used_station, index_org, vertix_index, vertices, triangle_indices, tnbr)

  if(present(slowness_est_matrix) .and. present(nsta_correlation)) then
    allocate(slowness_est_matrix(1 : 2, 1 : maxval(nsta_correlation), 1 : ngrid_x, 1 : ngrid_y))
    do kk = 1, ngrid_y
      do jj = 1, ngrid_x
        if(grid_enough_sta(jj, kk) .eqv. .false.) cycle
        allocate(g2(1 : nsta_correlation(jj, kk), 1 : 2), g_tmp2(1 : 2, 1 : 2))
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

        slowness_est_matrix(1 : 2, 1 : nsta_correlation(jj, kk), jj, kk) = 0.0_fp
        g_tmp2 = matmul(transpose(g2), g2)
#ifdef MKL
        call getrf(g_tmp2, ipiv = ipiv(1 : 2), info = info)
        call getri(g_tmp2, ipiv(1 : 2), info = info)
#else
        call LA_GETRF(g_tmp2, ipiv(1 : 2), info = info)
        call LA_GETRI(g_tmp2, ipiv(1 : 2), info = info)
#endif
        slowness_est_matrix(1 : 2, 1 : nsta_correlation(jj, kk), jj, kk) = matmul(g_tmp2, transpose(g2))
        deallocate(g_tmp2, g2)
      enddo
    enddo
  endif


  return
end subroutine calc_kernelmatrix_delaunay2


subroutine calc_slowness_est_matrix_delaunay(location_sta, nadd_station, ntriangle, &
&                                            triangle_center, slowness_matrix,      &
&                                            triangle_stationindex, nsta_count, tnbr)
  use nrtype, only : fp
  use constants, only : pi
  use typedef
  use gradiometry_parameters
  use lonlat_xy_conv, only : xy2bl
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif
  implicit none

  type(location),  intent(in)  :: location_sta(:)
  integer,         intent(in)  :: nadd_station
  integer,         intent(out) :: ntriangle
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
    call xy2bl(triangle_center(j)%y_north * 1000.0_fp, triangle_center(j)%x_east * 1000.0_fp, center_lon, center_lat, &
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
  npair = 1
  do i = 1, maxval(nsta_count) - 1
    npair = npair * i
  enddo

  allocate(slowness_matrix(1 : 2, 1 : npair, 1 : ntriangle))
  do jj = 1, ntriangle
    npair_tmp = 1
    do ii = 1, nsta_count(jj) - 1
      npair_tmp = npair_tmp * ii
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
