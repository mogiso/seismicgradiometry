! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

!!convert geographical-coordinate coastline data to cartesian-coordinate using Gauss-Kruger transformation 
program conv_coastline
  use nrtype, only : fp
  use lonlat_xy_conv, only : bl2xy

  implicit none
  real(kind = fp), parameter :: order = 1.0e-3_fp
  real(kind = fp) :: lon, lat, x_east, y_north, center_lon, center_lat
  character(len = 129) :: infile, outfile, center_lon_t, center_lat_t, buf
  integer :: ios

  call getarg(1, infile)
  call getarg(2, outfile)
  call getarg(3, center_lon_t); read(center_lon_t, *) center_lon
  call getarg(4, center_lat_t); read(center_lat_t, *) center_lat

  open(unit = 10, file = trim(infile))
  open(unit = 20, file = trim(outfile))
  do
    read(10, '(a128)', iostat = ios) buf
    if(ios .ne. 0) exit
    if(buf(1 : 1) .eq. ">") then
      write(20, '(a128)') buf
    else
      read(buf, *) lon, lat
      call bl2xy(lon, lat, center_lon, center_lat, y_north, x_east)
      y_north = y_north * order
      x_east = x_east * order
      write(20, '(e15.7, 1x, e15.7)') x_east, y_north
    endif
  enddo
  close(10)
  close(20)

  stop
end program conv_coastline

