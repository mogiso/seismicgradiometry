module particle_filter
  public :: particle_filter_search, particle_filter_init

contains

  subroutine particle_filter_search(narray, arrayindex, result_exist, lon_array, lat_array, min_correlation, az_obs, az_weight, &
  &                                 randomnumber_status, &
  &                                 lon_particle, lat_particle, likelihood_particle, appvel, arrivaltime, origintime)
    use nrtype, only : fp
    use constants
    use aeluma_parameters
    use greatcircle
    use random_number
    use xorshift1024star

    implicit none
    integer,         intent(in)            :: narray, arrayindex(:)
    logical,         intent(in)            :: result_exist(:)
    real(kind = fp), intent(in)            :: min_correlation(:), lon_array(:), lat_array(:), az_obs(:), az_weight(:)
    type(xorshift1024star_state), intent(inout) :: randomnumber_status
    real(kind = fp), intent(inout)         :: lon_particle(:), lat_particle(:), likelihood_particle(:)
    real(kind = fp), intent(in),  optional :: appvel(:), arrivaltime(:)
    real(kind = fp), intent(out), optional :: origintime
    integer                                :: i, j, ntriangle, arrivaltimeindex_min, az_weight_index, particlefilter_count
    real(kind = fp)                        :: rnd, rnd1, rnd2, maxval_likelihood_particle, daz, likelihood_tmp, &
    &                                         likelihood_azweight, kahan_val1, kahan_val2, sum_likelihood, &
    &                                         normalize_likelihood, arrivaltime_min, origintime_tmp, arrivaltime_diff, &
    &                                         likelihood_distweight
    real(kind = fp)                        :: particle_probability(1 : nparticle), &
    &                                         lon_particle_new(1 : nparticle), lat_particle_new(1 : nparticle)
    real(kind = fp), allocatable           :: az(:), dist(:)


    ntriangle = size(result_exist)
    allocate(az(1 : narray), dist(1 : narray))

    if(present(arrivaltime) .and. present(appvel)) then
      arrivaltime_min = 10.0_fp
      arrivaltimeindex_min = 0
      do i = 1, narray
        if(min_correlation(arrayindex(i)) .lt. correlation_threshold) cycle
        if(arrivaltime(arrayindex(i)) .le. arrivaltime_min) then
          arrivaltime_min = arrivaltime(arrayindex(i))
          arrivaltimeindex_min = i
        endif
      enddo
      if(arrivaltimeindex_min .eq. 0) error stop
    endif

    maxval_likelihood_particle = 0.0_fp
    particlefilter_count = 0
    particlefilter: do
      write(0, *) "particlefilter_count = ", particlefilter_count
      sum_likelihood = 0.0_fp
      kahan_val1 = 0.0_fp
      particleloop: do j = 1, nparticle
        likelihood_particle(j) = 0.0_fp
        do i = 1, narray
          if(.not. result_exist(arrayindex(i))) cycle
          if(min_correlation(arrayindex(i)) .lt. correlation_threshold) cycle
          call greatcircle_dist(lat_array(arrayindex(i)), lon_array(arrayindex(i)), &
          &                     lat_particle(j), lon_particle(j), distance = dist(i), azimuth = az(i))
          az(i) = az(i) + pi
          if(az(i) .ge. 2.0_fp * pi) az(i) = az(i) - 2.0_fp * pi
          daz = az_obs(arrayindex(i)) - az(i)
          if(daz .gt.  pi) daz = 2.0_fp * pi - daz
          if(daz .lt. -pi) daz = 2.0_fp * pi + daz
          az_weight_index = int(az_obs(arrayindex(i)) / daz_weight) + 1
          likelihood_azweight = 1.0_fp - azweight_coef * exp(-(az_weight(az_weight_index) ** 2) / sameaz_num2)
          likelihood_tmp = exp(-(daz ** 2) * 0.5_fp / daz_weight2)
          likelihood_tmp = (1.0_fp - likelihood_azweight) * likelihood_tmp + likelihood_azweight
          if(likelihood_particle(j) .eq. 0.0_fp) then
            likelihood_particle(j) = likelihood_tmp
          else
            likelihood_particle(j) = likelihood_particle(j) * likelihood_tmp
          endif
        enddo

        if(present(arrivaltime) .and. present(appvel)) then
          origintime_tmp = arrivaltime_min - dist(arrivaltimeindex_min) * appvel(arrayindex(arrivaltimeindex_min))
          do i = 1, narray
            arrivaltime_diff = arrivaltime(arrayindex(i)) - (dist(i) * appvel(arrayindex(i)))
            likelihood_distweight = 1.0_fp - ttime_coef * exp(-(dist(i) ** 2) / sigma_dist2)
            likelihood_tmp = exp(-(arrivaltime_diff ** 2) / sigma_arrivaltimediff2)
            likelihood_tmp = (1.0_fp - likelihood_distweight) * likelihood_tmp + likelihood_distweight
            if(likelihood_particle(j) .eq. 0.0_fp) then
              likelihood_particle(j) = likelihood_tmp
            else
              likelihood_particle(j) = likelihood_particle(j) * likelihood_tmp
            endif
          enddo
        endif

        kahan_val1 = kahan_val1 + likelihood_particle(j)
        kahan_val2 = sum_likelihood
        sum_likelihood = sum_likelihood + kahan_val1
        kahan_val2 = sum_likelihood - kahan_val2
        kahan_val1 = kahan_val1 - kahan_val2
      enddo particleloop
      if(sum_likelihood .le. 1.0e-100_fp) exit particlefilter
      normalize_likelihood = 1.0_fp / sum_likelihood

      !!Rule for exit particlefilter loop
      if(maxval(likelihood_particle) .le. maxval_likelihood_particle) then
        particlefilter_count = particlefilter_count + 1
      else
        maxval_likelihood_particle = maxval(likelihood_particle)
      endif
      if(particlefilter_count .gt. niter) exit particlefilter

      particle_probability(1) = likelihood_particle(1) * normalize_likelihood
      do i = 2, nparticle
        particle_probability(i) = particle_probability(i - 1) + likelihood_particle(i) * normalize_likelihood
      enddo

      !!Redistribute particle
      do j = 1, nparticle
        call gen_random_number(randomnumber_status, rnd)
        do i = 1, nparticle
          if(rnd .le. particle_probability(i)) exit
        enddo
        if(i .gt. nparticle) then
          do i = 1, nparticle
            write(0, *) i, particle_probability(i), likelihood_particle(i), normalize_likelihood, sum_likelihood
          enddo
        endif
        call gen_random_number(randomnumber_status, rnd1)
        call gen_random_number(randomnumber_status, rnd2)
        rnd = sqrt(-2.0_fp * log(rnd1)) * cos(2.0_fp * pi * rnd2)
        lon_particle_new(j) = lon_particle(i) + rnd * sigma_particle
        rnd = sqrt(-2.0_fp * log(rnd1)) * sin(2.0_fp * pi * rnd2)
        lat_particle_new(j) = lat_particle(i) + rnd * sigma_particle
      enddo
      lon_particle(1 : nparticle) = lon_particle_new(1 : nparticle)
      lat_particle(1 : nparticle) = lat_particle_new(1 : nparticle)
    enddo particlefilter  

    likelihood_particle(1 : nparticle) = likelihood_particle(1 : nparticle) * normalize_likelihood
    deallocate(az)
    return
  end subroutine particle_filter_search

  subroutine particle_filter_init(randomnumber_status, lon_particle, lat_particle)
    use nrtype, only : fp
    use aeluma_parameters
    use random_number
    use xorshift1024star

    implicit none
    type(xorshift1024star_state), intent(inout) :: randomnumber_status
    real(kind = fp),              intent(out)   :: lon_particle(:), lat_particle(:)
    integer                                     :: i
    real(kind = fp)                             :: rnd

    do i = 1, nparticle
      call gen_random_number(randomnumber_status, rnd)
      lon_particle(i) = lon_w + (lon_e - lon_w) * rnd
      call gen_random_number(randomnumber_status, rnd)
      lat_particle(i) = lat_s + (lat_n - lat_s) * rnd
    enddo

    return
  end subroutine particle_filter_init



end module particle_filter
