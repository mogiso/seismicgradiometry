subroutine calc_bpf_order(fl, fh, fs, ap, as, sample, m, n, c)
  use nrtype, only : fp
  use constants, only : pi
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Butterworth band-pass filter by Saito(1978)                
  !!! fl: low frequency cut-off (Hz) not normalized              
  !!! fh: high frequency cut-off (Hz) not normalized             
  !!! fs: stopped frequency (Hz) not normalized                  
  !!! Ap: maximum attenuation in fp (usually 0.5)               
  !!! As: minimum attenuation in fs (usually 5.0)               
  !!! sample: sampling interval (sec)                            
  !!! m: order of filter                                         
  !!! n: order of butterworth function                           
  !!! c: filter coefficient                                      
  !!!                                                            
  !!! Usage:                                                     
  !!! 1. call calc_bpf_order(fl, fh, fs, ap, as, sample, m, n, c)
  !!! 2. allocate filter coefficient allocate(h(4 * m))          
  !!! 3. call calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)      
  !!! 4. call tandem(data, data, ndata, h, m, nml)               
  !!!
  !!! Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!! Reference: Saito, M. (1978) An automatic design algorithm for band selective recursive
  !!!                             digital filters, Geophysical Exploration, 31(4), 240-263 (In Japanese)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind = fp), intent(IN) :: fl,fh, fs, ap, as, sample
  integer, intent(OUT) :: m, n
  real(kind = fp), intent(OUT) :: c

  real(kind = fp) :: sigma_fl, sigma_fh, sigma_fs, op, os

  sigma_fl = pi * fl * sample
  sigma_fh = pi * fh * sample
  sigma_fs = abs(fs * sample) * pi

  op = sin(sigma_fh - sigma_fl) / (cos(sigma_fl) * cos(sigma_fh))
  os = abs(tan(sigma_fs) - tan(sigma_fh) * tan(sigma_fl) / tan(sigma_fs))

  n = max(2, int(abs(real(log(ap / as) / log(op / os), kind = fp)) + 0.5_fp))
  m = n
  c = sqrt(exp(log(ap * as) / real(n, kind = fp)) / (op * os))

  return
end subroutine calc_bpf_order

