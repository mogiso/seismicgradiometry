program calc_ampradterm_diff
  use nrtype, only : fp, sp
  use constants
  use grdfile_io
  use gradiometry_parameters

  integer :: i, j, icount, ios, x_index, y_index
  real(kind = fp) :: slowness_x(1 : 2, 1 : ngrid_x, 1 : ngrid_y), slowness_y(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                  ampterm_x(1 : 2, 1 : ngrid_x, 1 : ngrid_y), ampterm_y(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                  ampterm_x_diff(1 : ngrid_x, 1 : ngrid_y), ampterm_y_diff(1 : ngrid_x, 1 : ngrid_y), &
  &                  amp_geospread_diff(1 : ngrid_x, 1 : ngrid_y), amp_radterm_diff(1 : ngrid_x, 1 : ngrid_y), &
  &                  direction(1 : 2), amp_geospread(1 : 2), amp_radterm(1 : 2)
  real(kind = sp) :: x_tmp, y_tmp, sx_tmp, sy_tmp, sigma_sx_tmp, sigma_sy_tmp, ampterm_x_tmp, ampterm_y_tmp, &
  &                  sigma_ampterm_x_tmp, sigma_ampterm_y_tmp
  character(len = 129) :: infile(1 : 2), geospread_diff_grd, radterm_diff_grd



  call getarg(1, infile(1))
  call getarg(2, infile(2))
  call getarg(3, geospread_diff_grd)
  call getarg(4, radterm_diff_grd)


  slowness_x(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  slowness_y(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  ampterm_x(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  ampterm_y(1 : 2, 1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  do i = 1, 2
    open(unit = 10, file = trim(infile(i)), form = "unformatted", access = "direct", recl = 4 * 10)
    icount = 1
    do
      read(10, rec = icount, iostat = ios) x_tmp, y_tmp, sx_tmp, sy_tmp, sigma_sx_tmp, sigma_sy_tmp, &
      &                                    ampterm_x_tmp, ampterm_y_tmp, sigma_ampterm_x_tmp, sigma_ampterm_y_tmp
      if(ios .ne. 0) exit
      icount = icount + 1
      x_index = int((real(x_tmp, kind = fp) - x_start) / dgrid_x) + 1
      y_index = int((real(y_tmp, kind = fp) - y_start) / dgrid_y) + 1
      slowness_x(i, x_index, y_index) = real(sx_tmp, kind = fp)
      slowness_y(i, x_index, y_index) = real(sy_tmp, kind = fp)
      ampterm_x(i, x_index, y_index) = real(ampterm_x_tmp, kind = fp)
      ampterm_y(i, x_index, y_index) = real(ampterm_y_tmp, kind = fp)
    enddo
    close(10)
  enddo

  amp_geospread_diff(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  amp_radterm_diff(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  do j = 1, ngrid_y
    do i = 1, ngrid_x
      if(slowness_x(1, i, j) .ne. 0.0_fp .and. slowness_y(1, i, j) .ne. 0.0_fp .and. &
      &  slowness_x(2, i, j) .ne. 0.0_fp .and. slowness_y(2, i, j) .ne. 0.0_fp) then
        direction(1 : 2) = atan2(slowness_x(1 : 2, i, j), slowness_y(1 : 2, i, j))
        ampterm_x_diff(i, j) = ampterm_x(2, i, j) - ampterm_x(1, i, j)
        ampterm_y_diff(i, j) = ampterm_y(2, i, j) - ampterm_y(1, i, j)
        amp_geospread(1 : 2) = ampterm_x(1 : 2, i, j) * sin(direction(1 : 2)) + ampterm_y(1 : 2, i, j) * cos(direction(1 : 2))
        amp_geospread(1 : 2) = amp_geospread(1 : 2) * 100.0_fp
        amp_radterm(1 : 2) = ampterm_x(1 : 2, i, j) * cos(direction(1 : 2)) - ampterm_y(1 : 2, i, j) * sin(direction(1 : 2))
        amp_radterm(1 : 2) = amp_radterm(1 : 2) * 100.0_fp
        amp_geospread_diff(i, j) = amp_geospread(2) - amp_geospread(1)
        amp_radterm_diff(i, j) = amp_radterm(2) - amp_radterm(1)
        !amp_geospread_diff(i, j) = ampterm_x_diff(i, j) * sin(direction(2)) + ampterm_y_diff(i, j) * cos(direction(2))
        !amp_geospread_diff(i, j) = amp_geospread_diff(i, j) * 100.0_fp
        !amp_radterm_diff(i, j) = ampterm_x_diff(i, j) * cos(direction(2)) - ampterm_y_diff(i, j) * sin(direction(2))
        !amp_radterm_diff(i, j) = amp_radterm_diff(i, j) * 100.0_fp
      endif
    enddo
  enddo

  call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, amp_geospread_diff, geospread_diff_grd, &
  &                        nanval = 0.0_fp)
  call write_grdfile_fp_2d(x_start, y_start, dgrid_x, dgrid_y, ngrid_x, ngrid_y, amp_radterm_diff, radterm_diff_grd, &
  &                        nanval = 0.0_fp)


  stop
end program calc_ampradterm_diff
