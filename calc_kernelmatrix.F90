! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

module calc_kernelmatrix

  private
  public :: calc_kernelmatrix_circle, calc_kernelmatrix_delaunay

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

#ifdef ELLIPSE
      dist_grid_sta(1 : nsta_grid_max) = 0.0_fp
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      propagation_direction = propagation_direction - pi
      !print *, location_grid(jj, kk)%lon, location_grid(jj, kk)%lat, propagation_direction * rad2deg
      do ii = 1, nsta
        dist_x_tmp =   (location_sta(ii)%y_north - location_grid(jj, kk)%y_north) * cos(propagation_direction) &
        &            + (location_sta(ii)%x_east  - location_grid(jj, kk)%x_east)  * sin(propagation_direction)
        dist_y_tmp = - (location_sta(ii)%y_north - location_grid(jj, kk)%y_north) * sin(propagation_direction) &
        &            + (location_sta(ii)%x_east  - location_grid(jj, kk)%x_east)  * cos(propagation_direction)
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
            &                     location_sta(ii)%lat, location_sta(ii)%lon, &
            &                     azimuth = azimuth_grid_sta(j))
            exit
          endif
        enddo
      enddo
      if(nsta_count(jj, kk) .ge. nsta_grid_min) grid_enough_sta(jj, kk) = .true.
      if(nsta_count(jj, kk) .gt. nsta_grid_max) nsta_count(jj, kk) = nsta_grid_max
#else
      !!count usable stations for each grid
      dist_grid_sta(1 : nsta_grid_max) = 1.0e+15
      do ii = 1, nsta
        dist_tmp = sqrt((location_sta(ii)%x_east  - location_grid(jj, kk)%x_east)  ** 2 &
        &             + (location_sta(ii)%y_north - location_grid(jj, kk)%y_north) ** 2)
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
        g(i, 2) = location_sta(grid_stationindex(i, jj, kk))%x_east - location_grid(jj, kk)%x_east
        g(i, 3) = location_sta(grid_stationindex(i, jj, kk))%y_north - location_grid(jj, kk)%y_north

#ifdef ELLIPSE
        dist_x_tmp =   (location_sta(grid_stationindex(i, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * cos(propagation_direction) &
        &            + (location_sta(grid_stationindex(i, jj, kk))%x_east  &
        &            -  location_grid(jj, kk)%x_east)  * sin(propagation_direction)
        dist_y_tmp = - (location_sta(grid_stationindex(i, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * sin(propagation_direction) &
        &            + (location_sta(grid_stationindex(i, jj, kk))%x_east  &
        &            -  location_grid(jj, kk)%x_east)  * cos(propagation_direction)

        weight(i, i) = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
#else
        weight(i, i) = exp(-(g(i, 2) ** 2 + g(i, 3) ** 2) / (cutoff_dist ** 2))
#endif

      enddo
      g_tmp = matmul(transpose(g), matmul(weight, g))

#ifdef MKL
      call getrf(g_tmp, ipiv = ipiv, info = info)
      !write(0, '(a, i0)') "GETRF info = ", info
      call getri(g_tmp, ipiv, info = info)
      !write(0, '(a, i0)') "GETRI info = ", info
#else
      call LA_GETRF(g_tmp, ipiv, info = info)
      !write(0, '(a, i0)') "LA_GETRF info = ", info
      call LA_GETRI(g_tmp, ipiv, info = info)
      !write(0, '(a, i0)') "LA_GETRI info = ", info
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
  use sort, only : bubblesort
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif

  type(location),  intent(in)  :: location_grid(:, :), location_sta(:)
  logical,         intent(out) :: grid_enough_sta(:, :)
  integer,         intent(out) :: nsta_count(:, :), grid_stationindex(:, :, :)
  real(kind = fp), intent(out) :: kernel_matrix(:, :, :, :)

  integer         :: i, j, ii, jj, kk, info, nsta, ntriangle
  integer         :: ipiv(1 : 3)
  logical         :: is_inside
  real(kind = fp) :: propagation_direction, dist_x_tmp, dist_y_tmp, dist_tmp, az_diff_tmp
  real(kind = fp) :: azimuth_grid_sta(1 : nsta_grid_max), &
  &                  g(1 : nsta_grid_max, 1 : 3), g_tmp(1 : 3, 1 : 3), weight(1 : nsta_grid_max, 1 : nsta_grid_max), &
  &                  point_tmp(1 : 2), triangle_vertix_tmp(1 : 2, 1 : 3)
  integer, allocatable :: vertix_index(:), triangle_indices(:, :), tnbr(:, :)
  real(kind = fp), allocatable :: vertices(:, :)
 
  nsta = size(location_sta)

  !!Do delaunay triangulation
  allocate(vertix_index(1 : nsta), vertices(1 : 2, 1 : nsta), triangle_indices(1 : 3, 1 : 2 * nsta), &
  &        tnbr(1 : 3, 1 : 2 * nsta))
  do i = 1, nsta
    vertices(1, i) = location_sta(i)%x_east
    vertices(2, i) = location_sta(i)%y_north
    vertix_index(i) = i
  enddo
  call dtris2(nsta, vertices, vertix_index, ntriangle, triangle_indices, tnbr, info)

  do kk = 1, ngrid_y
    do jj = 1, ngrid_x
      nsta_count(jj, kk) = 0
      grid_enough_sta(jj, kk) = .false.
      grid_stationindex(1 : nsta_grid_max, jj, kk) = 0
      azimuth_grid_sta(1 : nsta_grid_max) = 1.0e+15

      !!find triangle that contains the grid
      point_tmp(1) = location_grid(jj, kk)%x_east
      point_tmp(2) = location_grid(jj, kk)%y_north
      do j = 1, ntriangle
        do i = 1, 3
          triangle_vertix_tmp(1, i) = location_sta(triangle_indices(i, j))%x_east
          triangle_vertix_tmp(2, i) = location_sta(triangle_indices(i, j))%y_north
        enddo
        call triangle_contains_point_2d_3(triangle_vertix_tmp, point_tmp, is_inside)
        if(is_inside .eqv. .true.) then
          grid_enough_sta(jj, kk) = .true.
          nsta_count(jj, kk) = 3
          do i = 1, nsta_count(jj, kk)
            grid_stationindex(i, jj, kk) = triangle_indices(i, j)
            call greatcircle_dist(location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
            &                     location_sta(grid_stationindex(ii, jj, kk))%lat, &
            &                     location_sta(grid_stationindex(ii, jj, kk))%lon, &
            &                     azimuth = azimuth_grid_sta(ii))
          enddo
          exit
        endif
      enddo
      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle 

#ifdef ELLIPSE
      !!Check whether all stations (vertices) are within the ellipse
      call greatcircle_dist(evlat, evlon, location_grid(jj, kk)%lat, location_grid(jj, kk)%lon, &
      &                     backazimuth=propagation_direction)
      propagation_direction = propagation_direction - pi
      do ii = 1, nsta_count(jj, kk)
        dist_x_tmp =   (location_sta(grid_stationindex(ii, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * cos(propagation_direction) &
        &            + (location_sta(grid_stationindex(ii, jj, kk))%x_east &
        &            -  location_grid(jj, kk)%x_east)  * sin(propagation_direction)
        dist_y_tmp = - (location_sta(grid_stationindex(ii, jj, kk))%y_north &
        &            - location_grid(jj, kk)%y_north) * sin(propagation_direction) &
        &            + (location_sta(grid_stationindex(ii, jj, kk))%x_east &
        &            -  location_grid(jj, kk)%x_east)  * cos(propagation_direction)
        dist_tmp = exp(-0.5_fp * ((dist_x_tmp / sigma_x) ** 2 + (dist_y_tmp / sigma_y) ** 2))
        if(dist_tmp .lt. exp(-2.0_fp)) then
          grid_enough_sta(jj, kk) = .false.
          exit
        endif
      enddo
#else
      !!check distance between grid and stations (vertices)
      do ii = 1, nsta_count(jj, kk)
        dist_tmp = sqrt((location_sta(grid_stationindex(ii, jj, kk))%x_east  - location_grid(jj, kk)%x_east)  ** 2 &
        &             + (location_sta(grid_stationindex(ii, jj, kk))%y_north - location_grid(jj, kk)%y_north) ** 2)
        if(dist_tmp .gt. cutoff_dist) then
          grid_enough_sta(jj, kk) = .false.
          exit
        endif
      enddo
#endif

      if(grid_enough_sta(jj, kk) .eqv. .false.) cycle

      g(1 : nsta_grid_max, 1 : 3) = 0.0_fp
      weight(1 : nsta_grid_max, 1 : nsta_grid_max) = 0.0_fp
      do i = 1, nsta_count(jj, kk)
        g(i, 1) = 1.0_fp
        g(i, 2) = location_sta(grid_stationindex(i, jj, kk))%x_east - location_grid(jj, kk)%x_east
        g(i, 3) = location_sta(grid_stationindex(i, jj, kk))%y_north - location_grid(jj, kk)%y_north

#ifdef ELLIPSE
        dist_x_tmp =   (location_sta(grid_stationindex(i, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * cos(propagation_direction) &
        &            + (location_sta(grid_stationindex(i, jj, kk))%x_east  &
        &            -  location_grid(jj, kk)%x_east)  * sin(propagation_direction)
        dist_y_tmp = - (location_sta(grid_stationindex(i, jj, kk))%y_north &
        &            -  location_grid(jj, kk)%y_north) * sin(propagation_direction) &
        &            + (location_sta(grid_stationindex(i, jj, kk))%x_east  &
        &            -  location_grid(jj, kk)%x_east)  * cos(propagation_direction)
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

  deallocate(vertix_index, vertices, triangle_indices, tnbr)

  return
end subroutine calc_kernelmatrix_delaunay
    

end module calc_kernelmatrix
