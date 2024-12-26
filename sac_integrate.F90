! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

!!Integrate sac-formatted waveform
program sac_integrate
  use nrtype, only : fp, dp, sp
  use constants, only : pi
  use read_sacfile, only : read_sachdr, read_sacdata
  use tandem, only : tandem3

  implicit none

  real(kind = fp), parameter :: ap = 0.5_fp, as = 5.0_fp
  integer, parameter :: ivel = 7, idisp = 6
  
  integer :: i, j, nfile, npts, decimate
  character(len = 4) :: sachdr_char(1 : 158)
  character(len = 129), allocatable :: sacfile(:)
  character(len = 129) :: outfile, carg
  real(kind = fp), allocatable :: waveform_org(:), waveform_vel(:), waveform_disp(:), waveform_decimate(:), time_mean(:)
  real(kind = fp) :: time_tmp, sampling, sampling_new, a0, a1

  real(kind = fp), allocatable :: h(:)
  real(kind = fp) :: c, gn, fl, fh, fs
  integer :: m, n

  call getarg(1, carg); read(carg, *) fl
  call getarg(2, carg); read(carg, *) fh
  call getarg(3, carg); read(carg, *) fs
  call getarg(4, carg); read(carg, *) sampling_new
  nfile = iargc() - 4
  allocate(sacfile(1 : nfile))
  do i = 1, nfile
    call getarg(i + 4, sacfile(i))
  enddo

  do j = 1, nfile
    open(unit = 10, file = sacfile(j), form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      read(10, rec = i) sachdr_char(i)
    enddo
    close(10)

    call read_sachdr(sacfile(j), delta = sampling, npts = npts)
    allocate(waveform_org(1 : npts), waveform_vel(1 : npts), waveform_disp(1 : npts), time_mean(1 : npts))
    call read_sacdata(sacfile(j), npts, waveform_org)

    do i = 1, npts
      time_mean(i) = real(i - 1, kind = fp) * sampling
    enddo
    call line_fit(time_mean, waveform_org, npts, a0, a1)
    do i = 1, npts
      time_tmp = real(i - 1, kind = fp) * sampling
      waveform_org(i) = waveform_org(i) - (a0 + a1 * time_tmp)
    enddo

    !!integrate
    write(0, '(a)') "Integrate once/twice"
    waveform_vel(1) = 0.0_fp
    waveform_disp(1 : 2) = 0.0_fp
    do i = 2, npts
      waveform_vel(i) = waveform_vel(i - 1) + (waveform_org(i) + waveform_org(i - 1)) * sampling * 0.5_fp
      if(i .gt. 2) waveform_disp(i) = sampling * sampling * 0.25_fp &
      &                             * (waveform_org(i) + 2.0_fp * waveform_org(i - 1) + waveform_org(i - 2)) &
      &                             + 2.0_fp * waveform_disp(i - 1) - waveform_disp(i - 2)
    enddo

    write(0, '(a)') "Appling band-pass filter"
    call calc_bpf_order(fl, fh, fs, ap, as, sampling, m, n, c)
    allocate(h(4 * m))
    call calc_bpf_coef(fl, fh, sampling, m, n, h, c, gn)
    call tandem3(waveform_vel, h, gn, 1)
    call tandem3(waveform_disp, h, gn, 1)

    decimate = int(sampling_new / sampling)
    allocate(waveform_decimate(npts / decimate))
    do i = 1, npts / decimate
      waveform_decimate(i) = waveform_vel(decimate * (i - 1) + 1)
    enddo

    outfile = "decimate_vel_" // trim(sacfile(j))
    write(0, '(2a)') "output ", trim(outfile)
    open(unit = 20, file = outfile, form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      write(20, rec = i) sachdr_char(i)
    enddo
    write(20, rec = 1) real(sampling_new, kind = sp)
    write(20, rec = 80) int(npts / decimate)
    write(20, rec = 87) ivel
    do i = 1, npts / decimate
      write(20, rec = 158 + i) real(waveform_decimate(i), kind = sp)
    enddo
    close(20)


    do i = 1, npts / decimate
      waveform_decimate(i) = waveform_disp(decimate * (i - 1) + 1)
    enddo

    outfile = "decimate_disp_" // trim(sacfile(j))
    write(0, '(2a)') "output ", trim(outfile)
    open(unit = 20, file = outfile, form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      write(20, rec = i) sachdr_char(i)
    enddo
    write(20, rec = 1) real(sampling_new, kind = sp)
    write(20, rec = 80) int(npts / decimate)
    write(20, rec = 87) idisp
    do i = 1, npts / decimate
      write(20, rec = 158 + i) real(waveform_decimate(i), kind = sp)
    enddo
    close(20)

    deallocate(waveform_org, waveform_vel, waveform_disp, time_mean, waveform_decimate)
    deallocate(h)
  enddo

  stop
end program sac_integrate
    
