program make_eventlist
  use nrtype, only : sp, fp
  use constants, only : rad2deg, deg2rad, pi
  use greatcircle, only : greatcircle_dist
  use aeluma_parameters
  use jday
  use particle_filter
  use origintime
  use typedef
  implicit none

  type(trigger_list), allocatable :: triglist(:, :)
  integer                         :: i, j, ios, yr_tmp, jday_tmp, sec_from_day_tmp, narray, ntriangle, &
  &                                  arrayindex
  real(kind = fp)                 :: lon_array_tmp, lat_array_tmp, az_obs_tmp, appvel_obs_tmp, min_correlation_tmp, &
  &                                  arrivaltime_tmp
  real(kind = fp),    allocatable :: lon_array(:), lat_array(:)

  integer, allocatable :: seed(:)
  integer              :: seedsize

  call random_seed(size = seedsize)
  allocate(seed(1 : seedsize))
  do i = 1, seedsize
    call system_clock(count = seed(i))
  enddo
  call random_seed(put = seed(:))

  !!read AELUMA results from stdin
  infinite_loop: do 
    read(*, *, iostat = ios) yr_tmp, jday_tmp, sec_from_day_tmp, narray, ntriangle
    if(ios .ne. 0) error stop
    if(.not. allocated(triglist)) then
      allocate(triglist(1 : ntrig_max, 1 : ntriangle), lon_array(1 : ntriangle), lat_array(1 : ntriangle))
      triglist(1 : ntrig_max, 1 : ntriangle)%trig_exist = .false.
    endif

    !!read array results and renew trigger list
    array_loop: do j = 1, narray
      read(*, *) arrayindex, lon_array(arrayindex), lat_array(arrayindex), &
      &          az_obs_tmp, appvel_obs_tmp, min_correlation_tmp, arrivaltime_tmp
      if(min_correlation_tmp .lt. correlation_threshold) cycle
      triglist_loop: do i = 1, ntrig_max
        if(triglist(i, arrayindex)%trig_exist) then
          if(abs(arrivaltime_tmp - triglist(i, arrayindex)%arrivaltime) .lt. arrivaltime_diff_threshold) then
            if(min_correlation_tmp .ge. triglist(i, arrayindex)%min_correlation) then
              triglist(i, arrayindex)%az_obs          = az_obs_tmp
              triglist(i, arrayindex)%appvel_obs      = appvel_obs_tmp
              triglist(i, arrayindex)%min_correlation = min_correlation_tmp
              triglist(i, arrayindex)%arrivaltime     = arrivaltime_tmp
            endif
            cycle array_loop
          endif
        else
          triglist(i, arrayindex)%az_obs          = az_obs_tmp
          triglist(i, arrayindex)%appvel_obs      = appvel_obs_tmp
          triglist(i, arrayindex)%min_correlation = min_correlation_tmp
          triglist(i, arrayindex)%arrivaltime     = arrivaltime_tmp
          triglist(i, arrayindex)%trig_exist      = .true.
          cycle array_loop
        endif
      enddo triglist_loop
    enddo array_loop

    do j = 1, ntriangle
      do i = 1, ntrig_max
        if(triglist(i, j)%trig_exist) then
          print *, i, j, lon_array(j), lat_array(j), triglist(i, j)%min_correlation, &
          &              triglist(i, j)%arrivaltime
        endif
      enddo
    enddo



  enddo infinite_loop
           
              
                
        

        



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

