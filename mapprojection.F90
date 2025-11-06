module mapprojection
  private
  public :: mercator

  contains

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
end module mapprojection

