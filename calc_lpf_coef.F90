subroutine calc_lpf_coef(m, n, h, c, gn)
  use nrtype, only : fp
  use constants, only : pi
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Butterworth low-pass filter by Saito(1978)
  !!! Input
  !!! m, n: filter order (calculated by calc_lpf_order.f90)
  !!! h(4 * m) : filter coefficient
  !!! gn : gain factor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer, intent(IN) :: m, n
  real(kind = fp), intent(IN) :: c
  real(kind = fp), intent(OUT) :: h(4 * m)
  real(kind = fp), intent(OUT) :: gn

  real(kind = fp) :: g, fj, c2, sj, tj, a, dp
  integer :: i

  dp = pi / 2.0_fp / real(n, kind = fp)

  g = 1.0_fp
  fj = 1.0_fp
  c2 = 2.0_fp * (1.0_fp - c) * (1.0_fp + c)

  do i = 1, int(n / 2)
    sj = cos(dp * fj) ** 2
    tj = sin(dp * fj)
    fj = fj + 2.0_fp
    a = 1.0_fp / ((c + tj) ** 2 + sj)
    g = g * a
    h(4 * i - 3) = 2.0_fp
    h(4 * i - 2) = 1.0_fp
    h(4 * i - 1) = c2 * a
    h(4 * i) = ((c - tj) ** 2 + sj) * a
  enddo
  gn = g
  if(mod(n, 2) .ne. 0) then
    gn = g / (1.0_fp + c)
    h(4 * m - 3) = 1.0_fp
    h(4 * m - 2) = 0.0_fp
    h(4 * m - 1) = (1.0_fp - c) / (1.0_fp + c)
    h(4 * m) = 0.0_fp
  endif

  return
end subroutine calc_lpf_coef

