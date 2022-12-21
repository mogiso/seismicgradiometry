! Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
! Released under the MIT license.
! see https://opensource.org/licenses/MIT

module sort
  use nrtype, only : fp
  private
  public :: bubblesort

contains

  subroutine bubblesort(data)
    implicit none
    real(kind = fp), intent(inout)  :: data(:)
    integer :: i, j, n
    real(kind = fp) :: tmp
 
    n = size(data)
    do j = 1, n - 1
      do i = 2, n - j + 1
        if(data(i) .gt. data(i - 1)) then
          tmp = data(i)
          data(i) = data(i - 1)
          data(i - 1) = tmp
        endif
      enddo
    enddo
  
    return 
  end subroutine bubblesort

end module sort

