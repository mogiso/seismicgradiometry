!!random_number.F90
!! interface for random number generator
!! Last update: 2019-11-29 16:39:22 slightly modified

module random_number
  use xorshift1024star
  use nrtype, only : fp

  implicit none
  private
  public :: make_seed, random_generator_init, random_generator_jump, gen_random_number

contains

!!make seed for random number
  subroutine make_seed(seed)
    integer, intent(out) :: seed
    call system_clock(count = seed)
    write(0, '(a, i15)') "seed = ", seed
    return
  end subroutine make_seed

!!initialize random number generator
  subroutine random_generator_init(status, seed)
    integer, intent(in) :: seed
    type(xorshift1024star_state), intent(out) :: status

    call state_init(status, seed)

    return
  end subroutine random_generator_init

!!jump random_number_status for parallel computing
  subroutine random_generator_jump(status, ntime)
    type(xorshift1024star_state), intent(inout) :: status
    integer, intent(in) :: ntime
    integer :: i

    do i = 1, ntime
      call state_jump(status)
    enddo
    return
  end subroutine random_generator_jump

!!return random number
  subroutine gen_random_number(status, return_num)
    type(xorshift1024star_state), intent(inout) :: status

    real(fp), intent(out) :: return_num    

    do
      return_num = real(draw_uniform(status), kind = fp)
      if(return_num .ne. 1.0_fp) exit
    enddo

    return
  end subroutine gen_random_number
    
end module random_number

