subroutine line_fit(x, y, data, a0, a1)
  implicit none

  !! Fit y = a0 + a1 * x to data (x(i), y(i)) (i = 1, m)
  !! in least-square sense
  !! Input variables ::
  !! data: number of data
  !! x(data), y(data): data used for fitting
  !! Output variables ::
  !! a0, a1: line parameter (y = a0 + a1 * x)

  integer, intent(IN) :: data
  real*8, intent(IN) :: x(data), y(data)
  real*8, intent(OUT) :: a0, a1

  integer :: i
  real*8 :: sigma_x, sigma_y, sigma_xx, sigma_xy


  sigma_xy = 0.0d0
  sigma_x = 0.0d0
  sigma_y = 0.0d0
  sigma_xx = 0.0d0

  do i = 1, data
    sigma_xy = sigma_xy + x(i) * y(i)
    sigma_x = sigma_x + x(i)
    sigma_y = sigma_y + y(i)
    sigma_xx = sigma_xx + x(i) ** 2
  enddo

  a1 = (dble(data) * sigma_xy - sigma_x * sigma_y) / (dble(data) * sigma_xx - sigma_x ** 2)
  a0 = (sigma_xx * sigma_y - sigma_xy * sigma_x) / (dble(data) * sigma_xx - sigma_x ** 2)

  return
end subroutine line_fit
