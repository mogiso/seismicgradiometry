program calc_semblance
  use nrtype, only : fp, sp
  use constants, only : pi, deg2rad, rad2deg
  use read_sacfile, only : read_sachdr, read_sacdata
  use grdfile_io, only : write_grdfile_fp_2d
  use lonlat_xy_conv, only : bl2xy

  implicit none


  type location
    real(kind = fp) :: lon, lat, x_east, y_north, depth
  end type location

  real(kind = fp), parameter :: sx_min = -10.0_fp, sx_max = 10.0_fp, sy_min = -10.0_fp, sy_max = 10.0_fp
  real(kind = fp), parameter :: dgrid_x = 0.2_fp, dgrid_y = 0.2_fp
  real(kind = fp), parameter :: center_lon = 135.75_fp, center_lat = 33.0_fp   !!origin of array
  real(kind = fp), parameter :: order = 1.0_fp

  !real(kind = fp), parameter :: fl = 1.0_fp / 100.0_fp, fh = 1.0_fp / 50.0_fp, fs = 1.0_fp / 20.0_fp, &
  !&                             ap = 0.5_fp, as = 5.0_fp
  !real(kind = fp), parameter :: fl = 1.0_fp / (60.0_fp * 60.0_fp), fh = 1.0_fp / (20.0_fp * 60.0_fp), &
  !&                             fs = 1.0_fp / (10.0_fp * 60.0_fp), ap = 0.5_fp, as = 5.0_fp  !!DONET long-period
  !integer, parameter :: ntime_slowness = 241                                                 !!DONET long-period
  real(kind = fp), parameter :: fl = 1.0_fp / (20.0_fp * 60.0_fp), fh = 1.0_fp / (6.0_fp * 60.0_fp), &
  &                             fs = 1.0_fp / (3.0_fp * 60.0_fp), ap = 0.5_fp, as = 10.0_fp  !!DONET short-period test
  integer, parameter :: ntime_slowness = 121                                                 !!DONET short-period

  integer, parameter :: ntime_slowness2 = (ntime_slowness - 1) / 2
  integer, parameter :: ngrid_x = int((sx_max - sx_min) / real(dgrid_x, kind = fp)) + 1
  integer, parameter :: ngrid_y = int((sy_max - sy_min) / real(dgrid_y, kind = fp)) + 1
  integer, parameter :: ntime = 630
  integer, parameter :: ntime_decimate = 4


  integer :: nsta, npts_tmp, i, j, k, ii, jj, m, n, delta_t_int, nsta_count
  type(location), allocatable :: location_sta(:)
  real(kind = fp), allocatable :: waveform_obs(:, :)         !!(1 : ntime, 1 : nsta)
  real(kind = fp), allocatable :: h(:)
  real(kind = fp) :: semblance(1 : ngrid_x, 1 : ngrid_y), semblance_maxval, slowness_x, slowness_y, begin, dt, c, gn, &
  &                  numerator, denominator, numerator_tmp, denominator_tmp, delta_t, velocity, azimuth
  integer :: semblance_maxloc(2)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  character(len = 4) :: time_index


  nsta = iargc()

  allocate(waveform_obs(1 : ntime * ntime_decimate, 1 : nsta), location_sta(1 : nsta), sacfile(1 : nsta))

  do i = 1, nsta
    call getarg(i, sacfile(i))
  enddo

  !!read sac-formatted waveforms
  do i = 1, nsta
    call read_sachdr(sacfile(i), begin = begin, delta = dt, npts = npts_tmp, &
    &                stlon = location_sta(i)%lon, stlat = location_sta(i)%lat, stdp = location_sta(i)%depth)
    call read_sacdata(sacfile(i), ntime * ntime_decimate, waveform_obs(:, i))
    waveform_obs(1 : ntime * ntime_decimate, i) = waveform_obs(1 : ntime * ntime_decimate, i) * order
  enddo


  !!calculate filter parameter
  call calc_bpf_order(fl, fh, fs, ap, as, dt, m, n, c)
  allocate(h(4 * m))
  call calc_bpf_coef(fl, fh, dt, m, n, h, c, gn)
  !apply filter
  do i = 1, nsta
    call tandem1(waveform_obs(:, i), waveform_obs(:, i), ntime * ntime_decimate, h, m, 1)
    waveform_obs(1 : ntime * ntime_decimate, i) = waveform_obs(1 : ntime * ntime_decimate, i) * gn
    call tandem1(waveform_obs(:, i), waveform_obs(:, i), ntime * ntime_decimate, h, m, -1)
    waveform_obs(1 : ntime * ntime_decimate, i) = waveform_obs(1 : ntime * ntime_decimate, i) * gn
  enddo


  !!convert station longitude/latitude to x_east/y_north
  open(unit = 12, file = "station_location_array.txt")
  do i = 1, nsta
    call bl2xy(location_sta(i)%lon, location_sta(i)%lat, center_lon, center_lat, &
    &          location_sta(i)%y_north, location_sta(i)%x_east)
    location_sta(i)%y_north = location_sta(i)%y_north / 1000.0_fp
    location_sta(i)%x_east = location_sta(i)%x_east / 1000.0_fp
    write(12, '(5(e15.7, 1x))') location_sta(i)%x_east, location_sta(i)%y_north, location_sta(i)%lon, location_sta(i)%lat, &
    &                           location_sta(i)%depth
  enddo
  close(12)

  !!calculate semblance at each grid
  open(unit = 10, file = "semblance_maximum.txt")
  write(10, '(a)') "# time slowness_x slowness_y semblance"
  do k = 1, ntime
    write(time_index, '(i4)') k
    do i = 1, 4
      if(time_index(i : i) .eq. " ") time_index(i : i) = "0"
    enddo
    write(0, '(a)') "calculate semblance distribution for time " // time_index
    semblance(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
    do jj = 1, ngrid_y
      slowness_y = sy_min + dgrid_y * real(jj - 1, kind = fp)
      do ii = 1, ngrid_x
        slowness_x = sx_min + dgrid_x * real(ii - 1, kind = fp)
        numerator = 0.0_fp
        denominator = 0.0_fp
        do j = -ntime_slowness2, ntime_slowness2
          numerator_tmp = 0.0_fp
          denominator_tmp = 0.0_fp
          do i = 1, nsta
            delta_t = location_sta(i)%x_east * slowness_x + location_sta(i)%y_north * slowness_y
            delta_t_int = -int(delta_t / dt)
            if(ntime_decimate * (k - 1) + 1 + j + delta_t_int .ge. 1 &
            &  .and. ntime_decimate * (k - 1) + 1 + j + delta_t_int .le. ntime * ntime_decimate) then
              numerator_tmp = numerator_tmp + waveform_obs(ntime_decimate * (k - 1) + 1 + j + delta_t_int, i)
              denominator_tmp = denominator_tmp + waveform_obs(ntime_decimate * (k - 1) + 1 + j + delta_t_int, i) ** 2
            endif
          enddo
          numerator = numerator + numerator_tmp ** 2
          denominator = denominator + denominator_tmp
        enddo
        if(denominator .ne. 0.0_fp) then
          semblance(ii, jj) = numerator / (real(nsta, kind = fp) * denominator)
        endif
      enddo
    enddo
    outfile = "semblance_distribution_" // trim(time_index) // ".grd"
    call write_grdfile_fp_2d(sx_min, sy_min, dgrid_x, dgrid_y, ngrid_x, ngrid_y, semblance, outfile)

    semblance_maxval = maxval(semblance)
    semblance_maxloc = maxloc(semblance)
    slowness_x = sx_min + dgrid_x * real(semblance_maxloc(1) - 1, kind = fp)
    slowness_y = sy_min + dgrid_y * real(semblance_maxloc(2) - 1, kind = fp)
    velocity = 1.0_fp / sqrt(slowness_x ** 2 + slowness_y ** 2)
    azimuth = (pi / 2.0_fp - atan2(slowness_y, slowness_x)) * rad2deg
    if(azimuth .gt. 180.0_fp) azimuth = azimuth - 360.0_fp
    if(azimuth .lt. -180.0_fp) azimuth = azimuth + 360.0_fp
    write(0, '(4(e15.7, 1x))') begin + real((k - 1) * ntime_decimate + 1, kind = fp) * dt, &
    &                          velocity, azimuth, semblance_maxval
    write(10, '(4(e15.7, 1x))') begin + real((k - 1) * ntime_decimate + 1, kind = fp) * dt, &
    &                          velocity, azimuth, semblance_maxval
  enddo 
  close(10)


  stop
end program calc_semblance

