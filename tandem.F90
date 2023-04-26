!! Subroutine for applying time domain recursive filter
!! Copyright 2023 Masashi Ogiso (masashi.ogiso@gmail.com)
!! Released under the MIT license.
!! see https://opensource.org/licenses/MIT
!! Reference: Saito, M. (1978) An automatic design algorithm for band selective recursive
!!                             digital filters, Geophysical Exploration, 31(4), 240-263 (In Japanese)

module tandem
  private
  public :: tandem3

contains

  subroutine tandem3(data_inout, h, gn, nml, past_uv)
    use nrtype, only : fp
    implicit none
    real(kind = fp), intent(inout) :: data_inout(:)
    real(kind = fp), intent(in)    :: h(:)
    real(kind = fp), intent(in)    :: gn
    integer,         intent(in)    :: nml !!>0: normal filtering, <0: reverse filtering
    real(kind = fp), intent(inout), optional :: past_uv(:)

    integer :: i, j, ii, jj, kk, ndata, filter_order
    real(kind = fp) :: data_tmp, past(1 : 4)

    ndata = size(data_inout)
    filter_order = size(h) / 4
    if(present(past_uv)) then
      if(size(past_uv) .ne. 4 * filter_order) then
        write(0, '(2(a, i0), a)') "order mismatch between h(", (4 * filter_order), " and past_uv(", size(past_uv), ")"
        error stop
      endif
    endif

    do j = 1, filter_order
      if(present(past_uv)) then
        past(1 : 4) = past_uv(4 * (j - 1) + 1 : 4 * (j - 1) + 4)
      else
        past(1 : 4) = 0.0_fp
      endif
      if(nml .gt. 0) then
        ii = 1
        jj = ndata
        kk = 1
      elseif(nml .lt. 0) then
        ii = ndata
        jj = 1
        kk = -1
      else
         write(0, '(a)') "parameter nml should be >0 or <0"
         error stop
      endif
      do i = ii, jj, kk
        data_tmp = data_inout(i) &
        &        + h(4 * (j - 1) + 1) * past(1) + h(4 * (j - 1) + 2) * past(2) &
        &        - h(4 * (j - 1) + 3) * past(3) - h(4 * (j - 1) + 4) * past(4)
        past(2) = past(1)
        past(1) = data_inout(i)
        past(4) = past(3)
        past(3) = data_tmp
        data_inout(i) = data_tmp
      enddo
      if(present(past_uv)) past_uv(4 * (j - 1) + 1 : 4 * (j - 1) + 4) = past(1 : 4)
    enddo
    data_inout(1 : ndata) = data_inout(1 : ndata) * gn
 
    return
  end subroutine tandem3
end module tandem

