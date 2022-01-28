subroutine deconvolution(velocity, npts, sampling, gain, h, freq, acc)
  use nrtype, only : fp
  use constants, only : pi
  implicit none
  !!Deconvolv seismometer response in time domain (Kanamori and Livera(2008), eq. 4)
  !!This subroutine can be used for both broadband and short-period seismometer
  !!Input time series: velocity, output time series: acceleration
  !!Input variables: velocity time series, number of wavedata, sampling period (second), 
  !!                 gain factor between digits and physical quantities (e.g., m/s),
  !!                 damping factor (h), and natural frequency

  integer, intent(in) :: npts
  real(kind = fp), intent(in) :: velocity(npts) 
  real(kind = fp), intent(in) :: sampling, gain, h, freq
  real(kind = fp), intent(out) :: acc(npts)

  integer :: i
  real(kind = fp) :: omega, c0, c1, c2

  !!mapping frequency range (-inf, inf) to (-Nyquist, Nyquist) ref.: Maeda et al.(2011, JGR) eq.(4)
  omega = 2.0_fp * pi * (2.0_fp / sampling * atan2(sampling * freq, 2.0_fp))

  c0 = 1.0_fp / (gain * sampling)
  c1 = -2.0_fp * (1.0_fp + h * omega * sampling) / (gain * sampling)
  c2 = (1.0_fp + 2.0_fp * h * omega * sampling + sampling ** 2 * omega ** 2)/ (gain * sampling)

  acc(1) = 0.0_fp
  acc(2) = 0.0_fp

  do i = 3, npts
    acc(i) = acc(i - 1) + c2 * velocity(i) * gain + c1 * velocity(i - 1) * gain + c0 * velocity(i - 2) * gain
  enddo


  return
end subroutine deconvolution
