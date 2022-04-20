subroutine calc_lpf_order(fpass, fstop, ap, as, sample, m, n, c)
  use nrtype, only : fp
  use constants, only : pi
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Butterworth low-pass filter by Saito(1978)
  !!! fp: pass band frequency (Hz) not normalized
  !!! fs: stop band frequency (Hz) not normalized
  !!! Ap: maximum attenuation in fp (normally 0.5)
  !!! As: minimum attenuation in fs (normally 5.0)
  !!! sample: sampling interval (sec)
  !!! m: order of filter
  !!! n: order of butterworth function
  !!! c: filter coefficient
  !!! 
  !!! Usage:
  !!! 1. call calc_lpf_order(fp, fs, ap, as, sample, m, n, c)
  !!! 2. allocate filter coefficient allocate(h(4 * m))
  !!! 3. call calc_lpf_coef(m, n, h, c, gn)
  !!! 4. call tandem(data, data, ndata, h, m, nml)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind = fp), intent(IN) :: fpass, fstop, ap, as, sample
  integer, intent(OUT) :: m, n
  real(kind = fp), intent(OUT) :: c

  real(kind = fp) :: sigma_fp, sigma_fs

  sigma_fp = tan(pi * fpass * sample)
  sigma_fs = tan(pi * fstop * sample)

  n = max(2, int(log(as / ap) / log(sigma_fs / sigma_fp) + 0.5_fp))
  m = ceiling(real(n, kind = fp) / 2.0_fp)
  c = sqrt(exp(log(ap * as) / real(n, kind = fp)) / (sigma_fp * sigma_fs))

  return
end subroutine calc_lpf_order
  
