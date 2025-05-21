module particlefilter_functions
  private
  public :: delta_az, likelihood_weight, likelihood_modified, origintime_cal, pickup_medianval

  contains

  function delta_az(az_obs, az_cal)
    use nrtype, only    : fp
    use constants, only : pi
    implicit none
    real(kind = fp)     :: delta_az, az_obs, az_cal

    delta_az = az_obs - az_cal
    if(delta_az .gt.  pi) delta_az = 2.0_fp * pi - delta_az
    if(delta_az .lt. -pi) delta_az = 2.0_fp * pi + delta_az

    return
  end function delta_az

  function likelihood_weight(c0, c1, k0)
    !!calculate 1.0 - c0 * exp(-(k0 ** 2) / (2.0 * c1))
    !!c1 should be already squared
    use nrtype, only : fp
    implicit none
    real(kind = fp) :: c0, c1, k0, likelihood_weight

    likelihood_weight = 1.0_fp - c0 * exp(-0.5_fp * k0 * k0 / c1)

  end function likelihood_weight

  function likelihood_modified(delta_obs, sigma2, likelihood_weight)
    !!calculate (1.0 - likelihood_weight) * exp(-(delta_obs ** 2 / (2 * sigma2))) + likelihood_weight
    !!sigma2 should be squared variance (sigma2 = sigma ** 2)
    use nrtype, only : fp
    implicit none
    real(kind = fp) :: delta_obs, sigma2, likelihood_weight, likelihood_modified

    likelihood_modified = (1.0_fp - likelihood_weight) * exp(-0.5_fp * delta_obs * delta_obs / sigma2) + likelihood_weight

    return
  end function likelihood_modified

  function origintime_cal(arrivaltime, distance, apparentvel)
    use nrtype, only : fp
    implicit none
    real(kind = fp) :: arrivaltime, distance, apparentvel, origintime_cal

    origintime_cal = arrivaltime - distance / apparentvel
    return
  end function origintime_cal

  subroutine pickup_medianval(array, val)
    !!pickup median value of array(1 : n), array should be already sorted
    use nrtype, only : fp
    implicit none
    real(kind = fp), intent(in)  :: array(:)
    real(kind = fp), intent(out) :: val
    integer :: maxindex

    maxindex = ubound(array, 1)
    if(mod(maxindex, 2) .eq. 0) then
      val = 0.5_fp * (array(maxindex / 2) + array(maxindex / 2 + 1))
    else
      val = array((maxindex + 1) / 2)
    endif
    return
  end subroutine pickup_medianval


end module particlefilter_functions

