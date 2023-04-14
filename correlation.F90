!! Module that contains subroutine to calculate correlation
!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT
!! use fftsg.f
!! calculation: cross_spectrum(\omega) = conjg(Wave1(\omega)) * Wave2(\omega)
!! lag index: maxloc(corr) - ndata / 2
!! correlation becomes maximum between wave1(i) and wave2(i + (lag index))

module correlation

  private
  public :: correlation_fft

contains

subroutine correlation_fft(waveform, ndata, corr)
  use nrtype, only : dp
  implicit none
  integer, intent(in) :: ndata
  real(kind = dp), intent(in) :: waveform(1 : ndata, 1 : 2)
  real(kind = dp), intent(out) :: corr(-ndata / 2 + 1 : ndata / 2)

  real(kind = dp) :: avg(1 : 2), w(0 : ndata / 2 - 1), &
  &                  cross_spectrum(0 : 2 * ndata - 1), waveform_cmplx(0 : 2 * ndata - 1, 1 : 2)
  integer         :: ip(0 : 2 + int(sqrt(real(ndata, kind = dp)) + 0.5_dp))
  integer :: i

  do i = 1, ndata
    waveform_cmplx(2 * (i - 1), 1)     = waveform(i, 1)
    waveform_cmplx(2 * (i - 1) + 1, 1) = 0.0_dp
    waveform_cmplx(2 * (i - 1), 2)     = waveform(i, 2)
    waveform_cmplx(2 * (i - 1) + 1, 2) = 0.0_dp
  enddo

  ip(0) = 0
  do i = 1, 2
    call cdft(2 * ndata, -1, waveform_cmplx(:, i), ip(:), w)
  enddo

  avg(1 : 2) = 0.0_dp
  do i = 0, ndata - 1
    cross_spectrum(2 * i) = waveform_cmplx(2 * i, 1)     * waveform_cmplx(2 * i, 2) &
    &                     + waveform_cmplx(2 * i + 1, 1) * waveform_cmplx(2 * i + 1, 2)
    cross_spectrum(2 * i + 1) = -waveform_cmplx(2 * i + 1, 1) * waveform_cmplx(2 * i, 2) &
    &                         +  waveform_cmplx(2 * i, 1)     * waveform_cmplx(2 * i + 1, 2)
    avg(1) = avg(1) + (waveform_cmplx(2 * i, 1) ** 2 + waveform_cmplx(2 * i + 1, 1) ** 2)
    avg(2) = avg(2) + (waveform_cmplx(2 * i, 2) ** 2 + waveform_cmplx(2 * i + 1, 2) ** 2)
  enddo
  if(avg(1) .eq. 0.0_dp .or. avg(2) .eq. 0.0_dp) then
    corr(-ndata / 2 + 1 : ndata / 2) = 0.0_dp
    return
  endif

  !!omit multiplying real(ndata, kind = dp) after ifft of cross spectrum, for it is implicitly included in normalization
  call cdft(2 * ndata, 1, cross_spectrum, ip(:), w)
  cross_spectrum(0 : 2 * ndata - 1) = cross_spectrum(0 : 2 * ndata - 1) / sqrt(avg(1) * avg(2))

  corr(0) = cross_spectrum(0)
  corr(ndata / 2) = cross_spectrum(ndata)
  do i = 1, ndata / 2 - 1
    corr(i) = cross_spectrum(2 * i)
    corr(-ndata / 2 + i) = cross_spectrum(2 * (ndata / 2 + i))
  enddo

  return
end subroutine correlation_fft

end module correlation
