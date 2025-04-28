module typedef

  use nrtype, only : fp
  implicit none

  type location
    real(kind = fp) :: lon, lat, x_east, y_north, depth
  end type location

  type trigger_list
    real(kind = fp) :: az_obs, appvel_obs, min_correlation, arrivaltime
    integer         :: ref_yr, ref_jday, ref_secfromday
    logical         :: trig_exist, trig_used
  end type trigger_list

end module typedef
