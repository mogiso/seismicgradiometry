module typedef

  use nrtype, only : fp
  implicit none

  type location
    real(kind = fp) :: lon, lat, x_east, y_north, depth
  end type location

end module typedef
