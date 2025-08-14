module jday

public :: leapyear, ymd2jday, jday2ymd

contains
  subroutine leapyear(year, leap)
  !!Return leapyear or not
  !!Input: year
  !!Output: leap (true if year is leapyear)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: leap

    if(mod(year, 400) .eq. 0) then
      leap = .true.
    elseif(mod(year, 100) .eq. 0) then
      leap = .false.
    elseif(mod(year, 4) .eq. 0) then
      leap = .true.
    else
      leap = .false.
    endif
    return
  end subroutine leapyear

  subroutine ymd2jday(yy, mm, dd, jday)
  !!Return Julian day from yy/1/1
  !!Input:: yy(year), mm(month), dd(day)
  !!output :: jday
    implicit none
    integer, intent(in) :: yy, mm, dd
    integer, intent(out) :: jday
    integer :: dom(11), i
    logical :: leap

    dom(1) = 31
    dom(3) = 31
    dom(4) = 30
    dom(5) = 31
    dom(6) = 30
    dom(7) = 31
    dom(8) = 31
    dom(9) = 30
    dom(10) = 31
    dom(11) = 30

    call leapyear(yy, leap)
    if(leap) then
      dom(2) = 29
    else
      dom(2) = 28
    endif

    jday = dd
    if(mm .ne. 1) then
      do i = 1, mm - 1
        jday = jday + dom(i)
      enddo
    endif

    return
  end subroutine ymd2jday

  subroutine jday2ymd(jday, yy, mm, dd)
  !!Return yy/mm/dd from Julian day yy/jday
  !!Input:: jday, yy(year)
  !!output :: mm(month), dd(day) 

    implicit none
    integer, intent(in) :: yy, jday
    integer, intent(out) :: mm, dd
    integer :: dom(12), i, jday_tmp
    logical :: leap

    dom(1) = 31
    dom(3) = 31
    dom(4) = 30
    dom(5) = 31
    dom(6) = 30
    dom(7) = 31
    dom(8) = 31
    dom(9) = 30
    dom(10) = 31
    dom(11) = 30
    dom(12) = 31

    call leapyear(yy, leap)
    if(leap) then
      dom(2) = 29
    else
      dom(2) = 28
    endif

    jday_tmp = dom(1)
    do i = 2, 13
      if(jday .le. jday_tmp) then
        mm = i - 1
        exit
      endif
      jday_tmp = jday_tmp + dom(i)
    enddo
    dd = jday - (jday_tmp - dom(mm))

    return
  end subroutine jday2ymd
end module jday

