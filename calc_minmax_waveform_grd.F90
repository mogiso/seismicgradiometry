!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT

program calc_minmax_waveform_grd
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use typedef
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy, xy2bl
  use gradiometry_parameters
  use calc_kernelmatrix
  use tandem

  implicit none

  integer, parameter :: minmax_begin = 250, minmax_end = 320

  integer :: nsta, npts_tmp, i, j, ii, jj, m, n
  type(location) :: location_grid(1 : ngrid_x, 1 : ngrid_y)
  type(location),  allocatable :: location_sta(:)
  real(kind = fp), allocatable :: begin(:)
  real(kind = fp), allocatable :: waveform_obs(:, :)
  real(kind = fp), allocatable :: minmaxval_sta(:, :)
  real(kind = fp), allocatable :: h(:), uv(:, :)             !!For time-domain recursive filter 
  real(kind = fp)              :: dt, c, gn, &
  &                               obsvector(1 : nsta_grid_max), &
  &                               minmaxval_grd(1 : ngrid_x, 1 : ngrid_y, 1 : 3, 1 : 2), &
  &                               kernel_matrix(1 : 3, 1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y)
  integer                      :: grid_stationindex(1 : nsta_grid_max, 1 : ngrid_x, 1 : ngrid_y), &
  &                               nsta_count(1 : ngrid_x, 1 : ngrid_y)
  logical                      :: grid_enough_sta(1 : ngrid_x, 1 : ngrid_y)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 8),   allocatable :: stname(:)
  character(len = 129) :: outfile

  nsta = command_argument_count()

  allocate(waveform_obs(1 : ntime, 1 : nsta))
  allocate(location_sta(1 : nsta), sacfile(1 : nsta), begin(1 : nsta), stname(1 : nsta))

  do i = 1, nsta
    call get_command_argument(i, value = sacfile(i))
  enddo

  !!read sac-formatted waveforms
  do i = 1, nsta
    call read_sachdr(sacfile(i), begin = begin(i), delta = dt, npts = npts_tmp, stname = stname(i), &
    &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
    call read_sacdata(sacfile(i), ntime, waveform_obs(:, i))
#ifdef ENVELOPE
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order * order
#else
    waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * order
#endif
#ifdef GREEN_CORRECTION
    if(location_sta(i)%depth .ne. -12.345_fp) &
    &  waveform_obs(1 : ntime, i) = waveform_obs(1 : ntime, i) * (location_sta(i)%depth / depth_ref) ** 0.25_fp
#endif
  enddo

#ifdef ENVELOPE
#else
  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(1 : 4 * m), uv(1 : 4 * m, 1 : nsta))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  do i = 1, nsta
    call tandem3(waveform_obs(:, i), h, gn, 1, past_uv = uv(:, i))
  enddo
  !!inverse filtering
  uv(1 : 4 * m, 1 : nsta) = 0.0_fp
  do i = 1, nsta
    call tandem3(waveform_obs(:, i), h, gn, -1, past_uv = uv(:, i))
  enddo
  deallocate(h, uv)
#endif

  !!find min/max value
  allocate(minmaxval_sta(1 : nsta, 1 : 2))
  do i = 1, nsta
    minmaxval_sta(i, 1) = minval(waveform_obs((minmax_begin - 1) * ntimestep + 1: (minmax_end - 1) * ntimestep + 1, i))
    minmaxval_sta(i, 2) = maxval(waveform_obs((minmax_begin - 1) * ntimestep + 1: (minmax_end - 1) * ntimestep + 1, i))
  enddo

  !!set grid location
  do j = 1, ngrid_y
    do i = 1, ngrid_x
      location_grid(i, j)%x_east  = x_start + dgrid_x * real(i - 1, kind = fp)
      location_grid(i, j)%y_north = y_start + dgrid_y * real(j - 1, kind = fp)
      call xy2bl(location_grid(i, j)%y_north * 1000.0_fp, location_grid(i, j)%x_east * 1000.0_fp, &
      &          center_lon, center_lat, location_grid(i, j)%lon, location_grid(i, j)%lat)
    enddo
  enddo
  !!convert station longitude/latitude to x_east/y_north
  open(unit = 12, file = "station_location_minmaxamp.txt")
  do i = 1, nsta
    call bl2xy(location_sta(i)%lon, location_sta(i)%lat, center_lon, center_lat, &
    &          location_sta(i)%y_north, location_sta(i)%x_east)
    location_sta(i)%y_north = location_sta(i)%y_north / 1000.0_fp
    location_sta(i)%x_east  = location_sta(i)%x_east / 1000.0_fp
    write(12, '(7(e15.7, 1x), a)') location_sta(i)%x_east, location_sta(i)%y_north, location_sta(i)%lon, location_sta(i)%lat, &
    &                           location_sta(i)%depth, minmaxval_sta(i, 1), minmaxval_sta(i, 2), trim(stname(i))
  enddo
  close(12)

  !!make kernel matrix for each grid
  !call calc_kernelmatrix_circle(location_grid, location_sta, grid_enough_sta, nsta_count, grid_stationindex, kernel_matrix)
  call calc_kernelmatrix_delaunay2(location_grid, location_sta, naddstation_array, grid_enough_sta, &
  &                                nsta_count, grid_stationindex, kernel_matrix)

  !!calculate amplitude and its spatial derivatives at each grid
  do jj = 1, ngrid_y
    do ii = 1, ngrid_x
      if(grid_enough_sta(ii, jj) .eqv. .false.) cycle
      do i = 1, 2
        obsvector(1 : nsta_count(ii, jj)) = minmaxval_sta(grid_stationindex(1 : nsta_count(ii, jj), ii, jj), i)
        minmaxval_grd(ii, jj, 1 : 3, i) &
        &  = matmul(kernel_matrix(1 : 3, 1 : nsta_count(ii, jj), ii, jj), obsvector(1 : nsta_count(ii, jj)))
      enddo
    enddo
  enddo

  outfile = "minval_waveform.grd"
  outfile = trim(outfile)
  call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, &
  &                        minmaxval_grd(:, :, 1, 1), outfile, nanval = 0.0_fp)
  outfile = "maxval_waveform.grd"
  outfile = trim(outfile)
  call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, &
  &                        minmaxval_grd(:, :, 1, 2), outfile, nanval = 0.0_fp)
  

  stop
end program calc_minmax_waveform_grd

