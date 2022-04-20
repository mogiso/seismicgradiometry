program calc_lpf_test
  implicit none
  real(8) :: fl, fh, fs, ap, as, sample, c, gn, data(360000), data_f(360000)
  real(8), allocatable :: h(:)
  real(4) :: tmp_h(200), data_r(360000), gn1, data_rf(360000)
  integer(4) :: m, n, i, buf, m1

  fl = 0.1d0
  fh = 0.5d0
  fs = 1.0d0
  ap = 0.5d0
  as = 5.0d0
  sample = 0.01d0

  call calc_bpf_order(fl, fh, fs, ap, as, sample, m, n, c)
  allocate(h(4 * m))
  call calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)

  print *, m, n, c, gn

  do i = 1, m
    print '(I3, 4F13.8)', i, h(4 * i - 3), h(4 * i - 2), h(4 * i - 1), h(4 * i)
  enddo
  print *, ' '

  call butpas(tmp_h, m1, gn1, n, real(fl * sample), real(fh * sample), real(fs * sample), real(ap), real(as))
  print *, m1, n, gn1
  do i = 1, m
    print '(I3, 4F13.8)', i, tmp_h(4 * i - 3), tmp_h(4 * i - 2), tmp_h(4 * i - 1), tmp_h(4 * i)
  enddo

  !open(unit = 3, file = "pns_cmg_ew.dat", form = 'unformatted', access = 'direct', recl = 4)
  !do i = 1, 360000
  !  read(3, rec = i) buf
  !  data(i) = dble(buf)
  !  data_r(i) = real(buf)
  !enddo
  !close(3)

  !call tandem1(data, data_f, 360000, h, m, 1)
  !call tandem(data_r, data_rf, 360000, tmp_h, m1, 1)

  !do i = 1, 360000
  !  print '(I6, 3F15.8)', i, data(i), data_f(i) * gn, data_rf(i) * gn1
  !enddo


  stop
end program calc_lpf_test
