program trace_ampterm
  use nrtype, only : fp, sp
  use gradiometry_parameters

  implicit none

  real(kind = fp), parameter :: dtimestep = 30.0_fp
  real(kind = fp) :: slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y), &
  &                  sigma_slowness(1 : 2, 1 : ngrid_x, 1 : ngrid_y), sigma_ampterm(1 : 2, 1 : ngrid_x, 1 : ngrid_y)
  real(kind = fp) :: ray_xini, ray_yini, ray_xtmp, ray_ytmp, ray_az, ray_slowness, ray_length, ray_geom_spread
  real(kind = sp) :: xloc, yloc, slowness_x_tmp, slowness_y_tmp, sigma_slowness_x_tmp, sigma_slowness_y_tmp, &
  &                  ampterm_x_tmp, ampterm_y_tmp, sigma_ampterm_x_tmp, sigma_ampterm_y_tmp
  integer         :: xindex, yindex, irec, ios
  character(len = 129) :: slowness_gradiometry_file, ray_xini_t, ray_yini_t
 
  slowness(1 : 2, 1 : ngrid_x, 1 :ngrid_y) = 0.0_fp
  sigma_slowness(1 : 2, 1 : ngrid_x, 1 :ngrid_y) = 0.0_fp
  ampterm(1 : 2, 1 : ngrid_x, 1 :ngrid_y) = 0.0_fp
  sigma_ampterm(1 : 2, 1 : ngrid_x, 1 :ngrid_y) = 0.0_fp

  call get_command_argument(1, value = slowness_gradiometry_file)
  call get_command_argument(2, value = ray_xini_t); read(ray_xini_t, *) ray_xini
  call get_command_argument(3, value = ray_yini_t); read(ray_yini_t, *) ray_yini

  open(unit = 10, file = slowness_gradiometry_file, form = "unformatted", access = "direct", recl = 4 * 10)
  irec = 1
  do 
    read(10, rec = irec, iostat = ios) xloc, yloc, &
    &                                  slowness_x_tmp, slowness_y_tmp, sigma_slowness_x_tmp, sigma_slowness_y_tmp, &
    &                                  ampterm_x_tmp,  ampterm_y_tmp,  sigma_ampterm_x_tmp,  sigma_ampterm_y_tmp
    if(ios .ne. 0) exit
    xindex = (real(xloc, kind = fp) - x_start) / dgrid_x + 1
    yindex = (real(yloc, kind = fp) - y_start) / dgrid_y + 1
    slowness(1 : 2, xindex, yindex) = [real(slowness_x_tmp, kind = fp), real(slowness_y_tmp, kind = fp)]
    ampterm(1 : 2, xindex, yindex)  = [real(ampterm_x_tmp, kind = fp),  real(ampterm_y_tmp, kind = fp)]
    sigma_slowness(1 : 2, xindex, yindex) = [real(sigma_slowness_x_tmp, kind = fp), real(sigma_slowness_y_tmp, kind = fp)]
    sigma_ampterm(1 : 2, xindex, yindex)  = [real(sigma_ampterm_x_tmp, kind = fp),  real(sigma_slowness_y_tmp, kind = fp)]
    irec = irec + 1
  enddo
  close(10)

  !!do ray trace
  ray_xtmp = ray_xini
  ray_ytmp = ray_yini
  ray_geom_spread = 0.0_fp
  do
    xindex = (ray_xtmp - x_start) / dgrid_x + 1
    yindex = (ray_ytmp - y_start) / dgrid_y + 1


    if(.not. (xindex .ge. 1 .and. xindex .le. ngrid_x .and. yindex .ge. 1 .and. yindex .le. ngrid_y)) exit
    if(slowness(2, xindex, yindex) .eq. 0.0_fp) exit
    print '(2(f6.1, 1x), 2(i0, 1x), e15.8)', ray_xtmp, ray_ytmp, xindex, yindex, ray_geom_spread


    ray_az = atan2(slowness(1, xindex, yindex), slowness(2, xindex, yindex))
    ray_slowness = sqrt(slowness(1, xindex, yindex) ** 2 + slowness(2, xindex, yindex) ** 2)
    ray_length = dtimestep / ray_slowness
    ray_xtmp = ray_xtmp + ray_length * sin(ray_az)
    ray_ytmp = ray_ytmp + ray_length * cos(ray_az)
    ray_geom_spread = ray_geom_spread &
    &               + ray_length * (ampterm(1, xindex, yindex) * sin(ray_az) + ampterm(2, xindex, yindex) * cos(ray_az))
  enddo

  stop
end program trace_ampterm
