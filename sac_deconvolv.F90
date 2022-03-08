program sac_deconvolv
  use nrtype, only : fp, dp, sp
  use constants, only : pi
  use read_sacfile, only : read_sachdr, read_sacdata

  implicit none

  real(kind = fp), parameter :: damping_constant = 0.7_fp, natural_freq = 1.0_fp
  real(kind = fp), parameter :: fl = 0.02_fp, fh = 0.05_fp, fs = 0.1_fp, ap = 0.5_fp, as = 5.0_fp
  integer, parameter :: decimate = 100
  integer, parameter :: nmean = 100 * 15
  
  integer :: i, j, nfile, npts
  character(len = 4) :: sachdr_char(1 : 158)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile
  real(kind = fp), allocatable :: waveform_org(:), waveform_acc(:), waveform_vel(:), waveform_decimate(:)
  real(kind = fp) :: mean, sampling, sampling_new

  real(kind = fp), allocatable :: h(:)
  real(kind = fp) :: c, gn
  integer :: m, n

  nfile = iargc()
  allocate(sacfile(1 : nfile))
  do i = 1, nfile
    call getarg(i, sacfile(i))
  enddo

  do j = 1, nfile
    open(unit = 10, file = sacfile(j), form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      read(10, rec = i) sachdr_char(i)
    enddo
    close(10)

    call read_sachdr(sacfile(j), delta = sampling, npts = npts)
    allocate(waveform_org(1 : npts), waveform_acc(1 : npts), waveform_vel(1 : npts))
    call read_sacdata(sacfile(j), npts, waveform_org)

    mean = 0.0_dp
    do i = 1, nmean
      mean = mean + waveform_org(i)
    enddo
    mean = mean / real(nmean, kind = dp)
    waveform_org(1 : npts) = waveform_org(1 : npts) - mean

    write(0, '(a)') "Doing deconvolution"
    call deconvolution(waveform_org, npts, sampling, 1.0_fp, damping_constant, natural_freq, waveform_acc)

    mean = 0.0_dp
    do i = 1, nmean
      mean = mean + waveform_acc(i)
    enddo
    mean = mean / real(nmean, kind = dp)
    waveform_acc(1 : npts) = waveform_acc(1 : npts) - mean

    write(0, '(a)') "Appling band-pass filter"
    call calc_bpf_order(fl, fh, fs, ap, as, sampling, m, n, c)
    allocate(h(4 * m))
    call calc_bpf_coef(fl, fh, sampling, m, n, h, c, gn)
    call tandem1(waveform_acc, waveform_acc, npts, h, m, 1)
    waveform_acc(1 : npts) = waveform_acc(1 : npts) * gn
    call tandem1(waveform_acc, waveform_acc, npts, h, m, -1)
    waveform_acc(1 : npts) = waveform_acc(1 : npts) * gn

    waveform_vel(1) = 0.0_dp
    do i = 2, npts
      waveform_vel(i) = waveform_vel(i - 1) + (waveform_acc(i) + waveform_acc(i - 1)) * sampling * 0.5_dp
    enddo
    allocate(waveform_decimate(npts / decimate))
    do i = 1, npts / decimate
      waveform_decimate(i) = waveform_vel(decimate * (i - 1) + 1)
    enddo
    sampling_new = sampling * real(decimate, kind = fp)

    outfile = "decimate_decon_" // trim(sacfile(j))
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

    deallocate(waveform_org, waveform_acc, waveform_vel, waveform_decimate)
    deallocate(h)
  enddo

  stop
end program sac_deconvolv
    
