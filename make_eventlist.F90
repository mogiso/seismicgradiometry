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
      allocate(triglist(1 : ntrig_max, 1 : ntriangle))
      triglist(1 : ntrig_max, 1 : ntriangle)%trig_exist = .false.
    endif

    !!read array results and renew trigger list
    do j = 1, narray
      print *, yr_tmp, jday_tmp, j, narray
      read(*, *) arrayindex, lon_array_tmp, lat_array_tmp, az_obs_tmp, appvel_obs_tmp, min_correlation_tmp, arrivaltime_tmp
      if(min_correlation_tmp .ge. correlation_threshold) then
        triglist_loop: do i = 1, ntrig_max
          if(triglist(i, arrayindex)%trig_exist) then
            if(abs(arrivaltime_tmp - triglist(i, arrayindex)%arrivaltime) .lt. arrivaltime_diff_threshold) then
              !!case that arrival times are almost the same
              if(min_correlation_tmp .ge. triglist(i, arrayindex)%min_correlation) then
                triglist(i, arrayindex)%lon_array       = lon_array_tmp
                triglist(i, arrayindex)%lat_array       = lat_array_tmp
                triglist(i, arrayindex)%az_obs          = az_obs_tmp
                triglist(i, arrayindex)%appvel_obs      = appvel_obs_tmp
                triglist(i, arrayindex)%min_correlation = min_correlation_tmp
                triglist(i, arrayindex)%arrivaltime     = arrivaltime_tmp
                exit triglist_loop
              endif
            else
              cycle triglist_loop
            endif
          else
            triglist(i, arrayindex)%lon_array       = lon_array_tmp
            triglist(i, arrayindex)%lat_array       = lat_array_tmp
            triglist(i, arrayindex)%az_obs          = az_obs_tmp
            triglist(i, arrayindex)%appvel_obs      = appvel_obs_tmp
            triglist(i, arrayindex)%min_correlation = min_correlation_tmp
            triglist(i, arrayindex)%arrivaltime     = arrivaltime_tmp
            triglist(i, arrayindex)%trig_exist      = .true.
            exit triglist_loop
          endif
        enddo triglist_loop
      endif
    enddo

    do j = 1, ntriangle
      do i = 1, ntrig_max
        if(triglist(i, j)%trig_exist) then
          print *, i, j, triglist(i, j)%lon_array, triglist(i, j)%lat_array, triglist(i, j)%min_correlation, &
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

