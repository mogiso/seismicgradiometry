program conv_xy2lonlat
  use nrtype, only : fp
  use lonlat_xy_conv, only : xy2bl

  implicit none
  real(kind = fp), parameter :: order = 1.0e-3_fp
  real(kind = fp) :: lon, lat, x_east, y_north, center_lon, center_lat
  character(len = 129) :: x_east_t, y_north_t, center_lon_t, center_lat_t

  call getarg(1, x_east_t); read(x_east_t, *) x_east
  call getarg(1, y_north_t); read(y_north_t, *) y_north
  call getarg(3, center_lon_t); read(center_lon_t, *) center_lon
  call getarg(4, center_lat_t); read(center_lat_t, *) center_lat
  call xy2bl(y_north / order, x_east / order, center_lon, center_lat, lon, lat)

  print '(a)', "lon lat x_east y_north"
  print '(2(f9.4, 1x), 2(e15.7, 1x))', lon, lat, x_east, y_north

  stop
end program conv_xy2lonlat

