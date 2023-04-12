module taper
  private
  public :: cosine_taper

contains

subroutine cosine_taper(alpha, ndata, window)
  use nrtype, only : fp
  use constants, only : pi
  implicit none
  integer, intent(IN) :: ndata
  real(kind = fp), intent(IN) :: alpha
  real(kind = fp), intent(OUT) :: window(ndata)

  integer :: i, m

  m = int((ndata - 2) * alpha)
  do i = 1, ndata
    if(i .ge. 1 .and. i .le. m + 1) then
      window(i) = 0.5_fp * (1.0_fp - cos(pi * real(i - 1, kind = fp) / real(m + 1, kind = fp)))
    elseif(i .gt. m + 1 .and. i .le. ndata - m - 1) then
      window(i) = 1.0_fp
    else
      window(i) = 0.5_fp * (1.0_fp - cos(pi * real(ndata - i, kind = fp) / real(m + 1, kind = fp)))
    endif
  enddo

  return
end subroutine cosine_taper

end module taper
