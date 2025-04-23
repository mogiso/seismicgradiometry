module legend
  private
  public :: plot_legend

  contains
  subroutine plot_legend
    use nrtype, only : sp
    use aeluma_parameters

    implicit none

    integer, parameter :: iwin_legend = 1
    integer            :: i, color(1 : 3)
    real(kind = sp)    :: plot_x, plot_y
    character(len = 6) :: plottext
    character(len = 5),   parameter   :: likelihood_legend_normalize_c = "x1e-3"


    call pc_plotinit(iwin_legend, "Legend", 0.0_sp, -300.0_sp, width / 2, 27.0_sp, scale)
    call pc_setbkcolor(iwin_legend, 255, 255, 255)
    do i = 1, 5
      if(i .eq. 1) then
        color(1 : 3) = [220, 204, 222]
        plottext = "0.0   "
      elseif(i .eq. 2) then
        color(1 : 3) = [212, 156, 189]
        plottext = "0.2   "
      elseif(i .eq. 3) then
        color(1 : 3) = [196, 110, 155]
        plottext = "0.4   "
      elseif(i .eq. 4) then
        color(1 : 3) = [136, 97, 141]
        plottext = "0.6   "
      elseif(i .eq. 5) then
        color(1 : 3) = [73, 57, 100]
        plottext = "0.8   "
      endif
      plot_x = real(i, kind = sp) * 16.0_sp
      plot_y = 20.0_sp
      call pc_setcolor(iwin_legend, color(1), color(2), color(3))
      call pc_setline(iwin_legend, 4)
      call pc_vector(iwin_legend, plot_x, plot_y, 0.0_sp, vector_len, vector_width, vector_head1, vector_head2, 1)
      call pc_setcolor(iwin_legend, 0, 0, 0)
      plot_x = real(i, kind = sp) * 16.0_sp - 10.0_sp
      call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 4) 
    enddo
    plot_x = real(i, kind = sp) * 16.0_sp - 10.0_sp
    plottext = "1.0  "
    call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 4) 

    do i = 1, 6
      if(i .eq. 1) then
        color(1 : 3) = [252, 238, 158]
        plottext = "<0.1"
      elseif(i .eq. 2) then
        color(1 : 3) = [238, 179, 87]
        plottext = "0.1"
      elseif(i .eq. 3) then
        color(1 : 3) = [222, 117, 79]
        plottext = "0.3"
      elseif(i .eq. 4) then
        color(1 : 3) = [149, 66, 62]
        plottext = "1.0"
      elseif(i .eq. 5) then
        color(1 : 3) = [63, 39, 23]
        plottext = "3.0"
      elseif(i .eq. 6) then
        color(1 : 3) = [26, 26, 1]
        plottext = "10.0"
      endif
      plot_x = real(i, kind = sp) * 16.0_sp
      plot_y = 10.0_sp
      call pc_setcolor(iwin_legend, color(1), color(2), color(3))
      call pc_symbol(iwin_legend, plot_x, plot_y, 3.0_sp, 1, 0)
      call pc_setcolor(iwin_legend, 0, 0, 0)
      call pc_symbol(iwin_legend, plot_x, plot_y, 3.0_sp, 1, 1)
      plot_x = real(i, kind = sp) * 16.0_sp - 8.0_sp
      call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 5) 
    enddo
    plot_x = real(i, kind = sp) * 16.0_sp - 8.0_sp
    plottext = "10.0<"
    call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(plottext), 0.0, len(trim(plottext)), 5) 
    plot_x = plot_x + 20.0_sp
    call pc_text(iwin_legend, plot_x, plot_y, 4.0, trim(likelihood_legend_normalize_c), 0.0, &
    &            len(trim(likelihood_legend_normalize_c)), 5) 
    call pc_flush(iwin_legend)

    return
  end subroutine plot_legend
end module legend


