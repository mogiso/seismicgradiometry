module origintime
  public :: ot_search_median

  contains

  subroutine ot_search_median(narray, arrayindex, result_exist, lon_array, lat_array, min_correlation, az_obs, appvel_obs, &
  &                           arrivaltime, origintime_candidate, lon_epicenter, lat_epicenter, origintime_median, narray_use)
    use nrtype, only : fp
    use constants
    use AELUMA_parameters
    use greatcircle
    use sort

    implicit none
    integer,         intent(in)  :: narray, arrayindex(:)
    logical,         intent(in)  :: result_exist(:)
    real(kind = fp), intent(in)  :: lon_array(:), lat_array(:), min_correlation(:), az_obs(:), appvel_obs(:), arrivaltime(:), &
    &                               lon_epicenter, lat_epicenter
    real(kind = fp), intent(out) :: origintime_candidate(:), origintime_median
    integer,         intent(out) :: narray_use
    integer                      :: i
    real(kind = fp)              :: az, dist_tmp, daz
    real(kind = fp), parameter   :: groupvelocity = 3.0_fp

    origintime_candidate(1 : size(origintime_candidate)) = 0.0_fp
    narray_use = 0
    do i = 1, narray
      if(.not. result_exist(arrayindex(i))) cycle
      if(min_correlation(arrayindex(i)) .lt. correlation_threshold) cycle
      call greatcircle_dist(lat_array(arrayindex(i)), lon_array(arrayindex(i)), &
      &                     lat_epicenter,            lon_epicenter,            &
      &                     azimuth = az, distance = dist_tmp)
      az = az + pi
      if(az .ge. 2.0_fp * pi) az = az - 2.0_fp * pi
      
      daz = az_obs(arrayindex(i)) - az
      if(daz .gt.  pi) daz = 2.0_fp * pi - daz
      if(daz .lt. -pi) daz = 2.0_fp * pi + daz
      if(abs(daz) .gt. daz_weight) cycle
      narray_use = narray_use + 1
      if(appvel_obs(arrayindex(i)) .lt. 3.5_fp) then
        origintime_candidate(narray_use) = arrivaltime(arrayindex(i)) - dist_tmp / groupvelocity
      else
        origintime_candidate(narray_use) = arrivaltime(arrayindex(i)) - dist_tmp / appvel_obs(arrayindex(i))
      endif
    enddo
    if(narray_use .eq. 0) then
      origintime_median = 0.0_fp
    elseif(narray_use .eq. 1) then
      origintime_median = origintime_candidate(narray_use)
    else
      call bubblesort(origintime_candidate(1 : narray_use))
      if(mod(narray_use, 2) .eq. 1) then
        origintime_median = origintime_candidate(int(narray_use / 2) + 1)
      else
        origintime_median = 0.5_fp * (origintime_candidate(int(narray_use / 2)) &
        &                           + origintime_candidate(int(narray_use / 2) + 1))
      endif
    endif
    return
  end subroutine ot_search_median
end module origintime

