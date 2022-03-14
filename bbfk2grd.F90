program bbfk2grd
  use nrtype, only : fp, sp
  use grdfile_io, only : write_grdfile_fp_2d

  implicit none

  character(len = 129) :: infile, outfile
  integer :: nxsize, nysize, i, j, icount
  real(kind = sp) :: buf
  real(kind = fp) :: xmin, xmax, ymin, ymax, dx, dy, zmax
  real(kind = fp), allocatable :: data(:, :)

  call getarg(1, infile)
  call getarg(2, outfile)

  open(unit = 10, file = trim(infile), form = "unformatted", access = "direct", recl = 4)
  read(10, rec = 60) buf; xmin = real(buf, kind = fp)
  read(10, rec = 61) buf; xmax = real(buf, kind = fp)
  read(10, rec = 62) buf; ymin = real(buf, kind = fp)
  read(10, rec = 63) buf; ymax = real(buf, kind = fp)
  read(10, rec = 83) nxsize
  read(10, rec = 84) nysize
  allocate(data(nxsize, nysize))
  zmax = 0.0_fp
  do j = 1, nysize
    do i = 1, nxsize
      icount = nxsize * (j - 1) + i
      read(10, rec = 158 + icount) buf; data(i, j) = real(buf, kind = fp)
      if(abs(data(i, j)) .ge. zmax) zmax = abs(data(i, j))
    enddo
  enddo
  close(10)

  dx = (xmax - xmin) / real(nxsize - 1, kind = fp)
  dy = (ymax - ymin) / real(nysize - 1, kind = fp)
  
  data(1 : nxsize, 1 : nysize) = data(1 : nxsize, 1 : nysize) / zmax
  call write_grdfile_fp_2d(xmin, ymin, dx, dy, nxsize, nysize, data(:, :), outfile)

  stop
end program bbfk2grd

