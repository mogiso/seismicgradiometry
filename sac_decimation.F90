program sac_decimation
  use nrtype, only : fp, dp, sp
  use constants, only : pi
  use read_sacfile, only : read_sachdr, read_sacdata

  implicit none

  real(kind = fp), parameter :: ap = 0.5_fp, as = 15.0_fp
  !integer, parameter :: nmean = 100 * 15
  
  integer :: i, j, nfile, npts, decimate
  character(len = 4) :: sachdr_char(1 : 158)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  real(kind = fp), allocatable :: waveform_org(:), waveform_decimate(:)
  real(kind = fp) :: mean, sampling, sampling_new

  real(kind = fp), allocatable :: h(:)
  real(kind = fp) :: fpass, fstop, c, gn
  integer :: m, n
  character(len = 65) :: decimate_c

  nfile = iargc()
  nfile = nfile - 1
  allocate(sacfile(1 : nfile))
  call getarg(1, decimate_c)
  read(decimate_c, *) decimate
  do i = 1, nfile
    call getarg(i + 1, sacfile(i))
  enddo

  do j = 1, nfile
    open(unit = 10, file = sacfile(j), form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      read(10, rec = i) sachdr_char(i)
    enddo
    close(10)

    call read_sachdr(sacfile(j), delta = sampling, npts = npts)
    allocate(waveform_org(1 : npts))
    call read_sacdata(sacfile(j), npts, waveform_org)

    !mean = 0.0_dp
    !do i = 1, nmean
    !  mean = mean + waveform_org(i)
    !enddo
    !mean = mean / real(nmean, kind = dp)
    !waveform_org(1 : npts) = waveform_org(1 : npts) - mean

    write(0, '(a)') "Appling low-pass filter"
    !!new sampling period
    sampling_new = sampling * real(decimate, kind = fp)
    !!fp and fs are determined from new sampling period
    fpass = 1.0_fp / (sampling_new * 3.0_fp)
    fstop = 1.0_fp / (sampling_new * 2.5_fp)
    write(0, '(a, 2(f6.3, 1x))') "old and new sampling period = ", sampling, sampling_new
    write(0, '(a, e15.7)') "new nyquist frequency (Hz) = ", 1.0_fp / (sampling_new * 2.0_fp)
    write(0, '(a, 2(e15.7, 1x))') "parameter fpass and fstop (Hz) = ", fpass, fstop
    call calc_lpf_order(fpass, fstop, ap, as, sampling, m, n, c)
    allocate(h(4 * m))
    call calc_lpf_coef(m, n, h, c, gn)
    call tandem1(waveform_org, waveform_org, npts, h, m, 1)
    waveform_org(1 : npts) = waveform_org(1 : npts) * gn
    call tandem1(waveform_org, waveform_org, npts, h, m, -1)
    waveform_org(1 : npts) = waveform_org(1 : npts) * gn

    allocate(waveform_decimate(npts / decimate))
    do i = 1, npts / decimate
      waveform_decimate(i) = waveform_org(decimate * (i - 1) + 1)
    enddo

    outfile = "decimate_" // trim(sacfile(j))
    write(0, '(2a)') "output ", trim(outfile)
    open(unit = 20, file = outfile, form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      write(20, rec = i) sachdr_char(i)
    enddo
    write(20, rec = 1) real(sampling_new, kind = sp)
    write(20, rec = 80) int(npts / decimate)
    do i = 1, npts / decimate
      write(20, rec = 158 + i) real(waveform_decimate(i), kind = sp)
    enddo
    close(20)

    deallocate(h, waveform_org, waveform_decimate)
  enddo

  stop
end program sac_decimation
    
