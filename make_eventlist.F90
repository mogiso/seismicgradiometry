program make_event_list
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  implicit none

  type(trigger_list), allocatable :: triglist

  call random_seed(size = seedsize)
  allocate(seed(1 : seedsize))
  do i = 1, seedsize
    call system_clock(count = seed(i))
  enddo
  call random_seed(put = seed(:))

  !!read AELUMA results from stdin
  do 
    read(*, *, iostat = ios) yr_tmp, jday_tmp, sec_from_day_tmp, narray, ntriangle
    if(ios .ne. 0) error stop
    if(.not. allocated(triglist)) then
      allocate(trglist(1 : ntrig_max, 1 : ntriangle))
    endif

    do j = 1, narray
      read(*, *) arrayindex(i), lon_array_tmp, lat_array_tmp, az_obs_tmp, appvel_obs_tmp, min_correlation_tmp, arrivaltime_tmp
      arrivaltime_tmp = sec_from_day_tmp - (real(nsec_buf, kind = fp) + arrivaltime_tmp)
      do i = 1, ntrig_max
      

        



end program make_eventlist

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
   
