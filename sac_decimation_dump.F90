program sac_decimation_dump
  use nrtype, only : fp
  use constants, only : pi
  use read_sacfile, only : read_sachdr, read_sacdata
  use tandem, only : tandem3
  use jday

  implicit none

  real(kind = fp), parameter :: ap = 0.1_fp, as = 15.0_fp, order = 1.0e+6_fp
  integer, parameter :: nmean = 100 * 15
  
  integer :: i, j, nfile, npts, decimate, julianday, jday_tmp, yr, mo, dy, hh, mm
  character(len = 4) :: sachdr_char(1 : 158), beginswitch
  character(len = 129) :: outfile, sacfile
  real(kind = fp), allocatable :: waveform_org(:)
  real(kind = fp) :: mean, sampling, sampling_new, waveform_decimate, ss

  real(kind = fp), allocatable :: h(:)
  real(kind = fp) :: c, gn, fl, fh, fs, begintime, sec_from_begin
  integer :: m, n 
  character(len = 65) :: decimate_c, fl_c, fh_c, fs_c

  call getarg(1, decimate_c); read(decimate_c, *) decimate
  call getarg(2, fl_c); read(fl_c, *) fl
  call getarg(3, fh_c); read(fh_c, *) fh
  call getarg(4, fs_c); read(fs_c, *) fs
  call getarg(5, sacfile)
  call getarg(6, beginswitch)

  open(unit = 10, file = trim(sacfile), form = "unformatted", access = "direct", recl = 4)

  call read_sachdr(trim(sacfile), delta = sampling, npts = npts, begin = begintime, year = yr, julianday = julianday)
  allocate(waveform_org(1 : npts))
  call read_sacdata(trim(sacfile), npts, waveform_org)
  !waveform_org(1 : npts) = waveform_org(1 : npts) * order

  mean = 0.0_fp
  do i = 1, nmean
    mean = mean + waveform_org(i)
  enddo
  mean = mean / real(nmean, kind = fp)
  waveform_org(1 : npts) = waveform_org(1 : npts) - mean

  write(0, '(a)') "Appling band-pass filter"
  !!new sampling period
  sampling_new = sampling * real(decimate, kind = fp)
  !!fp and fs are determined from new sampling period
  write(0, '(a, 2(f6.3, 1x))') "old and new sampling period = ", sampling, sampling_new
  write(0, '(a, e15.7)') "new nyquist frequency (Hz) = ", 1.0_fp / (sampling_new * 2.0_fp)
  write(0, '(a, 3(e15.7, 1x))') "parameter fl, fh, fs (Hz) = ", fl, fh, fs
  call calc_bpf_order(fl, fh, fs, ap, as, sampling, m, n, c)
  allocate(h(4 * m))
  call calc_bpf_coef(fl, fh, sampling, m, n, h, c, gn)
  call tandem3(waveform_org, h, gn, 1)
  !call tandem3(waveform_org, h, gn, -1)

  do i = 1, npts / decimate
    waveform_decimate = waveform_org(decimate * (i - 1) + 1)
    sec_from_begin = begintime + sampling * real(decimate * (i - 1), kind = fp)
    if(trim(beginswitch) .eq. "1") then
      print '(2(e15.7, 1x))', sec_from_begin, waveform_decimate
    else
      jday_tmp = julianday
      if(sec_from_begin .lt. 0.0_fp) then
        do
          sec_from_begin = sec_from_begin + 86400.0_fp
          jday_tmp = jday_tmp - 1
          if(sec_from_begin .ge. 0.0_fp) exit
        enddo
      endif
      if(sec_from_begin .ge. 86400.0_fp) then
        do
          sec_from_begin = sec_from_begin - 86400.0_fp
          jday_tmp = jday_tmp + 1
          if(sec_from_begin .lt. 86400.0_fp) exit
        enddo
      endif
      call jday2ymd(jday_tmp, yr, mo, dy)
      hh = int(sec_from_begin / (60.0_fp * 60.0_fp))
      mm = int((sec_from_begin - 60.0_fp * 60.0_fp * real(hh, kind = fp)) / 60.0_fp)
      ss = sec_from_begin - 60.0_fp * 60.0_fp * real(hh, kind = fp) - 60.0_fp * real(mm, kind = fp)
      print '(i4, 4(a, i2.2), a, f0.2, 1x, e15.7)', yr, "-", mo, "-", dy, "T", hh, ":", mm, ":", ss, waveform_decimate
    endif
  enddo


  deallocate(h, waveform_org)

  stop
end program sac_decimation_dump
    
