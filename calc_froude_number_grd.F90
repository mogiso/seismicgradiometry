program calc_froude_number_grd
  use nrtype, only : fp, sp
  use constants
  use gradiometry_parameters
  use grdfile_io
  use lonlat_xy_conv
  use calc_kernelmatrix
  use typedef

  implicit none

  character(len = 129) :: outfile, velocity_atm_t, direction_atm_deg_t, slownessfile, stationlocation_file
  real(kind = fp) :: x_east, y_north, direction_tsunami, amp_geospread, amp_radterm, delta_direction, weight_sum, depth_mean
  real(kind = fp) :: velocity_atm, velocity_atm_proj, direction_atm_deg, direction_atm_rad, velocity_tsunami
  real(kind = fp) :: dist_x, dist_y
  real(kind = fp) :: froude_number(1 : ngrid_x, 1 : ngrid_y), slowness_x(1 : ngrid_x, 1 : ngrid_y), &
  &                  slowness_y(1 : ngrid_x, 1 : ngrid_y), ampterm_x(1 : ngrid_x, 1 : ngrid_y), &
  &                  ampterm_y(1 : ngrid_x, 1 : ngrid_y), kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                  weight(1 : nsta_grid_max)
  real(kind = sp) :: x_tmp, y_tmp, sx_tmp, sy_tmp, sigma_sx_tmp, sigma_sy_tmp, ampterm_x_tmp, ampterm_y_tmp, &
  &                  sigma_ampterm_x_tmp, sigma_ampterm_y_tmp
  integer :: i, j, k, x_index, y_index, icount, ios, nsta
  integer :: grid_stationindex(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), nsta_count(1 : ngrid_x, 1 : ngrid_y)
  logical :: grid_enough_sta(1 : ngrid_x, 1 : ngrid_y)

  type(location), allocatable :: location_sta(:)
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)


  call getarg(1, slownessfile)
  call getarg(2, outfile)
  call getarg(3, velocity_atm_t); read(velocity_atm_t, *) velocity_atm
  call getarg(4, direction_atm_deg_t); read(direction_atm_deg_t, *) direction_atm_deg
  call getarg(5, stationlocation_file)

  if(direction_atm_deg .gt. 180.0_fp) direction_atm_deg = direction_atm_deg - 360.0_fp
  direction_atm_rad = direction_atm_deg * deg2rad


  open(unit = 10, file = stationlocation_file)
  nsta = 0
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nsta = nsta + 1
  enddo
  rewind(10)
  allocate(location_sta(1 : nsta))
  do i = 1, nsta
    read(10, *) location_sta(i)%x_east, location_sta(i)%y_north, location_sta(i)%lon, location_sta(i)%lat, &
    &           location_sta(i)%depth
  enddo
  close(10)
  !!set grid location
  do j = 1, ngrid_y
    do i = 1, ngrid_x
      location_grid(i, j)%x_east  = x_start + dgrid_x * real(i - 1, kind = fp)
      location_grid(i, j)%y_north = y_start + dgrid_y * real(j - 1, kind = fp)
      call xy2bl(location_grid(i, j)%y_north * 1000.0_fp, location_grid(i, j)%x_east * 1000.0_fp, &
      &          center_lon, center_lat, location_grid(i, j)%lon, location_grid(i, j)%lat)
    enddo
  enddo
  !!make kernel matrix for each grid
  !call calc_kernelmatrix_circle(location_grid, location_sta, grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)
  call calc_kernelmatrix_delaunay2(location_grid, location_sta, naddstation_array, grid_enough_sta, &
  &                                nsta_count, grid_stationindex, kernel_matrix)


  open(unit = 10, file = trim(slownessfile), form = "unformatted", access = "direct", recl = 4 * 10)
  icount = 1
  slowness_x(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  slowness_y(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  ampterm_x(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  ampterm_y(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  do
    read(10, rec = icount, iostat = ios) x_tmp, y_tmp, sx_tmp, sy_tmp, sigma_sx_tmp, sigma_sy_tmp, &
    &                                    ampterm_x_tmp, ampterm_y_tmp, sigma_ampterm_x_tmp, sigma_ampterm_y_tmp
    if(ios .ne. 0) exit
    icount = icount + 1
    x_index = nint((real(x_tmp, kind = fp) - x_start) / dgrid_x) + 1
    y_index = nint((real(y_tmp, kind = fp) - y_start) / dgrid_y) + 1
    slowness_x(x_index, y_index) = real(sx_tmp, kind = fp)
    slowness_y(x_index, y_index) = real(sy_tmp, kind = fp)
    ampterm_x(x_index, y_index) = real(ampterm_x_tmp, kind = fp)
    ampterm_y(x_index, y_index) = real(ampterm_y_tmp, kind = fp)
  enddo
  close(10)
  

  froude_number(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  do j = 1, ngrid_y
    y_north = y_start + dgrid_y * real(j - 1, kind = fp)
    do i = 1, ngrid_x
      x_east = x_start + dgrid_x * real(i - 1, kind = fp)
      if(y_north .eq. -340.0_fp) then
        !print '(2(i0, 1x), 2(e15.7, 1x), l)', i, j, slowness_x(i, j), slowness_y(i, j), grid_enough_sta(i, j)
      endif
      if(.not. (slowness_x(i, j) .ne. 0.0_fp .and. slowness_y(i, j) .ne. 0.0_fp)) cycle
      if(grid_enough_sta(i, j) .eqv. .false.) cycle

      weight_sum = 0.0_fp
      do k = 1, nsta_count(i, j)
        dist_x = location_sta(grid_stationindex(k, i, j))%x_east - location_grid(i, j)%x_east
        dist_y = location_sta(grid_stationindex(k, i, j))%y_north - location_grid(i, j)%y_north
        weight(k) = exp(-(dist_x ** 2 + dist_y ** 2) / (cutoff_dist ** 2))
        weight_sum = weight_sum + weight(k)
      enddo

      depth_mean = 0.0_fp
      do k = 1, nsta_count(i, j)
        depth_mean = depth_mean + location_sta(grid_stationindex(k, i, j))%depth * weight(k) / weight_sum
      enddo
      depth_mean = depth_mean * 1000.0_fp

      direction_tsunami = atan2(slowness_x(i, j), slowness_y(i, j))
      delta_direction = direction_atm_rad - direction_tsunami
      !delta_direction = 0.0_fp
      velocity_atm_proj = velocity_atm * cos(delta_direction)
      !velocity_tsunami = sqrt(1.0_fp / (slowness_x(i, j) ** 2 + slowness_y(i, j) ** 2)) * 1000.0_fp 
      velocity_tsunami = sqrt(9.80655_fp * depth_mean)
      froude_number(i, j) = velocity_atm_proj / velocity_tsunami
      amp_geospread = ampterm_x(i, j) * sin(direction_tsunami) + ampterm_y(i, j) * cos(direction_tsunami)
      amp_geospread = amp_geospread * 100.0_fp
      amp_radterm = ampterm_x(i, j) * cos(direction_tsunami) - ampterm_y(i, j) * sin(direction_tsunami)
      amp_radterm = amp_radterm * deg2rad * 10000.0_fp
      print '(2(i0, 1x), 3(f6.1, 1x), 4(e15.7, 1x))', i, j, x_east, y_north, depth_mean, froude_number(i, j), amp_geospread, &
      &                                               delta_direction * rad2deg, amp_radterm
      print *, velocity_tsunami, velocity_atm, velocity_atm_proj
      
    enddo
  enddo
  call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, froude_number, outfile, nanval = 0.0_fp)

  stop
end program calc_froude_number_grd
      
