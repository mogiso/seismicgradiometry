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

  subroutine combsort(data)
    implicit none
    real(kind = fp), intent(inout)  :: data(:)
    integer :: space, i, swap, n
    real(kind = fp) :: tmp

    n = size(data)
    space = int(real(n, kind = fp) / 1.3_fp)
    do
      if(space .lt. 1) exit
      swap = 0
      i = 1
      do
        if(space .eq. 9 .or. space .eq. 10) space = 11
        if(i + space .gt. n) exit
          if(data(i) .gt. data(i + space)) then
          tmp = data(i)
          data(i) = data(i + space)
          data(i + space) = tmp
          swap = swap + 1
        endif
        i = i + 1
      enddo
      if(space .eq. 1) then
        if(swap .eq. 0) exit
      else
        space = int(real(space, kind = fp) / 1.3_fp)
      endif
    enddo

    return 
  end subroutine combsort

end module sort

