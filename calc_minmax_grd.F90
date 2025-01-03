program calc_minmax_grd
  use nrtype, only : fp
  use gradiometry_parameters
  use grdfile_io

  integer :: i, j, k, nfile 
  character(len = 129), allocatable :: grdfile(:)

  real(kind = fp) :: xval(1 : ngrid_x), yval(1 : ngrid_y), &
  &                  zval_min(1 : ngrid_x, 1 : ngrid_y), zval_max(1 : ngrid_x, 1 : ngrid_y)
  real(kind = fp), allocatable :: xval_tmp(:), yval_tmp(:), zval_tmp(:, :)

  nfile = command_argument_count()
  allocate(grdfile(1 : nfile))
  do i = 1, nfile
    call get_command_argument(i, value = grdfile(i))
  enddo

  zval_min(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp
  zval_max(1 : ngrid_x, 1 : ngrid_y) = 0.0_fp

  do k = 1, nfile
    
    call read_grdfile_2d(grdfile(k), xval_tmp, yval_tmp, zval_tmp)
    xval(1 : ngrid_x) = xval_tmp(1 : ngrid_x)
    yval(1 : ngrid_y) = yval_tmp(1 : ngrid_y)

    do j = 1, ngrid_y
      do i = 1, ngrid_x
        if(zval_tmp(i, j) .ne. zval_tmp(i, j)) zval_tmp(i, j) = 0.0_fp

        if(zval_tmp(i, j) .le. zval_min(i, j)) zval_min(i, j) = zval_tmp(i, j)
        if(zval_tmp(i, j) .ge. zval_max(i, j)) zval_max(i, j) = zval_tmp(i, j)

      enddo
    enddo
    deallocate(xval_tmp, yval_tmp, zval_tmp)
  enddo

  do j = 1, ngrid_y
    do i = 1, ngrid_x
      print '(2(f6.1, 1x), 2(i0, 1x), 2(e15.8, 1x))', xval(i), yval(j), i, j, zval_min(i, j), zval_max(i, j)
    enddo
  enddo
  deallocate(grdfile)

  stop
end program calc_minmax_grd
