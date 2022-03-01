module lonlat_xy_conv
  use nrtype, only : fp
  implicit none
  private
  public :: bl2xy, xy2bl

  real(kind = fp), parameter :: inv_f        = 298.2572222_fp   !!GRS1980
  real(kind = fp), parameter :: major_radius = 6378137.0_fp     !!GRS1980
  real(kind = fp), parameter :: n            = 1.0_fp / (2.0_fp * inv_f - 1.0_fp)
  real(kind = fp), parameter :: m0           = 0.9999_fp

contains

  function atanh(x)
    real(kind = fp), intent(in) :: x
    real(kind = fp) :: atanh

    atanh = 0.5_fp * log((1.0_fp + x) / (1.0_fp - x))
  end function atanh

  subroutine bl2xy(lon, lat, lon0, lat0, x_north, y_east)
    use constants, only : deg2rad
    implicit none

    real(kind = fp), parameter :: alpha(1 : 5) = [ &
    &                              n * (0.5_fp + n * (-2.0_fp / 3.0_fp &
    &                            + n * (5.0_fp / 16.0_fp + n * (41.0_fp / 180.0_fp - n * 127.0_fp / 288.0_fp)))), &
    &                              n * n * (13.0_fp / 48.0_fp + n * (-0.6_fp &
    &                            + n * (557.0_fp / 1440.0_fp + n * 281.0_fp / 630.0_fp))), &
    &                              n * n * n * (61.0_fp / 240.0_fp + n * (-103.0_fp / 140.0_fp &
    &                            + n * 15061.0_fp / 26880.0_fp)), &
    &                              n * n * n * n * (49561.0_fp / 161280.0_fp + n * (-179.0_fp / 168.0_fp)), &
    &                              n * n * n * n * n * 34729.0_fp / 80640.0_fp ]
    real(kind = fp), parameter :: a(0 : 5) = [ &
    &                              1.0_fp + n * n * (0.25_fp + n * n / 64.0_fp), &
    &                             -1.5_fp * n * (1.0_fp + n * n * (-1.0_fp / 8.0_fp - n * n / 64.0_fp)), &
    &                              15.0_fp / 16.0_fp * n * n * (1.0_fp - 0.25_fp * n * n), &
    &                             -35.0_fp / 48.0_fp * n * n * n * (1.0_fp - n * n * 5.0_fp / 16.0_fp), &
    &                              315.0_fp / 512.0_fp * n * n * n * n, &
    &                             -693.0_fp / 1280.0_fp * n * n * n * n * n]

    real(kind = fp), intent(in) :: lon, lat         !!longitude and ratitude (degree)
    real(kind = fp), intent(in) :: lon0, lat0       !!origin longitude and latitude (degree)
    real(kind = fp), intent(out) :: x_north, y_east !!output location in xy-coordinate (meters)
   
    integer :: i
    real(kind = fp) :: phi, lambda_c, lambda_s, t, tbar, zeta, eta, abar, sigma, tau, s_phi0, ncount

    phi = lat * deg2rad

    lambda_c = cos((lon - lon0) * deg2rad)
    lambda_s = sin((lon - lon0) * deg2rad)
    t = sinh(atanh(sin(phi)) - 2.0_fp * sqrt(n) / (1.0_fp + n) * atanh(2.0_fp * sqrt(n) / (1.0_fp + n) * sin(phi)))
    tbar = sqrt(1.0_fp + t * t)
    zeta = atan2(t, lambda_c)
    eta = atanh(lambda_s / tbar)

    abar = m0 * major_radius / (1.0_fp + n) * a(0)
    sigma = 1.0_fp
    tau = 0.0_fp
    s_phi0 = 0.0_fp
    x_north = zeta
    y_east = eta
    do i = 1, 5
      ncount  = 2.0_fp * real(i, kind = fp)
      sigma   = sigma + ncount * alpha(i) * cos(ncount * zeta) * cosh(ncount * eta)
      tau     = tau + ncount * alpha(i) * sin(ncount * zeta) * sinh(ncount * eta)
      s_phi0  = s_phi0 + a(i) * sin(ncount * (lat0 * deg2rad))
      x_north = x_north + alpha(i) * sin(ncount * zeta) * cosh(ncount * eta)
      y_east  = y_east + alpha(i) * cos(ncount * zeta) * sinh(ncount * eta)
    enddo
    s_phi0 = (s_phi0 + a(0) * lat0 * deg2rad) * m0 * major_radius / (1.0_fp + n)
    x_north = abar * x_north - s_phi0
    y_east = abar * y_east

    return
  end subroutine bl2xy

  subroutine xy2bl(x_north, y_east, lon0, lat0, lon, lat)
    use constants, only : deg2rad, rad2deg
    implicit none

    real(kind = fp), parameter :: a(0 : 5) = [ &
    &                              1.0_fp + n * n * (0.25_fp + n * n / 64.0_fp), &
    &                             -1.5_fp * n * (1.0_fp + n * n * (-1.0_fp / 8.0_fp - n * n / 64.0_fp)), &
    &                              15.0_fp / 16.0_fp * n * n * (1.0_fp - 0.25_fp * n * n), &
    &                             -35.0_fp / 48.0_fp * n * n * n * (1.0_fp - n * n * 5.0_fp / 16.0_fp), &
    &                              315.0_fp / 512.0_fp * n * n * n * n, &
    &                             -693.0_fp / 1280.0_fp * n * n * n * n * n ]
    real(kind = fp), parameter :: beta(1 : 5) = [ &
    &                              n * (0.5_fp + n * (-2.0_fp / 3.0_fp + n * (37.0_fp / 96.0_fp + n * (1.0_fp / 360.0_fp &
    &                            + n * (-81.0_fp / 512.0_fp))))), &
    &                              n * n * (1.0_fp / 48.0_fp + n * (1.0_fp / 15.0_fp + n * (-437.0_fp / 1440.0_fp &
    &                            + n * 46.0_fp / 105.0_fp))), &
    &                              n * n * n * (17.0_fp / 480.0_fp + n * (-37.0_fp / 840.0_fp + n * (-209.0_fp / 4480.0_fp))), &
    &                              n * n * n * n * (4397.0_fp / 161280.0_fp + n * (-11.0_fp / 504.0_fp)), &
    &                              n * n * n * n * n * 4583.0_fp / 161280.0_fp ]
    real(kind = fp), parameter :: delta(1 : 6) = [ &
    &                              n * (2.0_fp + n * (-2.0_fp / 3.0_fp + n * (-2.0_fp + n * (116.0_fp / 45.0_fp &
    &                            + n * (26.0_fp / 45.0_fp + n * (-2854.0_fp / 675.0_fp)))))), &
    &                              n * n * (7.0_fp / 3.0_fp + n * (-8.0_fp / 5.0_fp + n * (-227.0_fp / 45.0_fp &
    &                            + n * (2704.0_fp / 315.0_fp + n * 2323.0_fp / 945.0_fp)))), &
    &                              n * n * n * (56.0_fp / 15.0_fp + n * (-136.0_fp / 35.0_fp + n * (1262.0_fp / 105.0_fp &
    &                            + n * 73814.0_fp / 2835.0_fp))), &
    &                              n * n * n * n * (4279.0_fp / 630.0_fp + n * (-332.0_fp / 35.0_fp &
    &                            + n * (-399572.0_fp / 14175.0_fp))), &
    &                              n * n * n * n * n * (4174.0_fp / 315.0_fp + n * (-144838.0_fp / 6237.0_fp)), &
    &                              n * n * n * n * n * n * 601676.0_fp / 22275.0_fp ]

    real(kind = fp), intent(in) :: x_north, y_east !!(x, y)-value in meters
    real(kind = fp), intent(in) :: lon0, lat0      !!longitude and latitude of origin (degree)
    real(kind = fp), intent(out) :: lon, lat       !!(x, y) -> (lon, lat) (degree)

    integer :: i
    real(kind = fp) :: abar, s_phi0, ncount, zeta, eta, zeta_d, eta_d, sigma, tau, kai

    abar = m0 * major_radius / (1.0_fp + n) * a(0)
    s_phi0 = 0.0_fp
    do i = 1, 5
      ncount = 2.0_fp * real(i, kind = fp)
      s_phi0 = s_phi0 + a(i) * sin(ncount * lat0 * deg2rad)
    enddo
    s_phi0 = (s_phi0 + a(0) * lat0 * deg2rad) * m0 * major_radius / (1.0_fp + n)
    
    zeta = (x_north + s_phi0) / abar
    eta = y_east / abar
    zeta_d = 0.0_fp
    eta_d = 0.0_fp 
    sigma = 0.0_fp
    tau = 0.0_fp
    do i = 1, 5
      ncount = 2.0_fp * real(i, kind = fp)
      zeta_d = zeta_d + beta(i) * sin(ncount * zeta) * cosh(ncount * eta)
      eta_d = eta_d + beta(i) * cos(ncount * zeta) * sinh(ncount * eta)
      sigma = sigma + ncount * beta(i) * cos(ncount * zeta) * cosh(ncount * eta)
      tau = tau + ncount * beta(i) * sin(ncount * zeta) * sinh(ncount * eta)
    enddo
    zeta_d = zeta - zeta_d
    eta_d = eta - eta_d
    sigma = 1.0_fp - sigma
    kai = asin(sin(zeta_d) / cosh(eta_d))
    
    lat = 0.0_fp
    do i = 1, 6
      ncount = 2.0_fp * real(i, kind = fp)
      lat = lat + delta(i) * sin(ncount * kai)
    enddo
    lat = rad2deg * (kai + lat)
    lon = lon0 + rad2deg * atan2(sinh(eta_d), cos(zeta_d))

    return
  end subroutine xy2bl
    
end module lonlat_xy_conv
