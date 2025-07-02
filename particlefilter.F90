module particlefilter
  public :: particlefilter_search, particlefilter_init

contains

  subroutine particlefilter_search(narray, arrayindex, result_exist, lon_array, lat_array, az_obs, az_weight, &
  &                                 randomnumber_status, &
  &                                 lon_particle, lat_particle, likelihood_particle, appvel, arrivaltime, origintime, &
  &                                 appvel_median)
    use nrtype, only : fp
    use constants
    use aeluma_parameters
    use greatcircle
    use random_number
    use xorshift1024star
    use sort
    use particlefilter_functions
    !$ use omp_lib

    implicit none
    integer,         intent(in)            :: narray, arrayindex(:)
    logical,         intent(in)            :: result_exist(:)
    real(kind = fp), intent(in)            :: lon_array(:), lat_array(:), az_obs(:), az_weight(:)
    type(xorshift1024star_state), intent(inout) :: randomnumber_status(:)
    real(kind = fp), intent(inout)         :: lon_particle(1 : nparticle), lat_particle(1 : nparticle)
    real(kind = fp), intent(out)           :: likelihood_particle(1 : nparticle)
    real(kind = fp), intent(in),  optional :: appvel(:), arrivaltime(:)
    real(kind = fp), intent(out), optional :: origintime(1 : nparticle), appvel_median
    integer                                :: i, j, narray_use_tmp, az_weight_index, particlefilter_count, iteration_count, &
    &                                         nthread, threadnum, parallelindex
    real(kind = fp)                        :: rnd, rnd1, rnd2, maxval_likelihood_particle, az_diff, likelihood_tmp, &
    &                                         likelihood_azweight, kahan_val1, kahan_val2, likelihood_sum, likelihood_sum_tmp, &
    &                                         normalize_likelihood, ot_diff, likelihood_distweight, dist_tmp, az
    real(kind = fp)                        :: particle_probability(1 : nparticle), &
    &                                         lon_particle_new(1 : nparticle), lat_particle_new(1 : nparticle)
    real(kind = fp), allocatable           :: dist(:), ot_est(:), appvel_tmp(:)
    integer,         allocatable           :: parallelindex_start(:), parallelindex_end(:)


    nthread = ubound(randomnumber_status, 1)
    allocate(parallelindex_start(1 : nthread), parallelindex_end(1 : nthread))
    do i = 1, nthread
      parallelindex_start(i) = nparticle / nthread * (i - 1) + 1
      if(i .ge. 2) parallelindex_end(i - 1) = parallelindex_start(i) - 1
    enddo
    parallelindex_end(nthread) = nparticle

    allocate(dist(1 : narray), ot_est(1 : narray))

    maxval_likelihood_particle = 0.0_fp
    particlefilter_count = 0
    iteration_count = 0

    !$omp parallel
    particlefilter: do
      !$omp single
      iteration_count = iteration_count + 1
      !$omp end single
      likelihood_sum = 0.0_fp

      !$omp do private(j, kahan_val1, likelihood_sum_tmp, az, az_diff, az_weight_index, likelihood_azweight, likelihood_tmp, &
      !$omp&          ot_est, dist, narray_use_tmp, ot_diff, dist_tmp, likelihood_distweight, kahan_val2), &
      !$omp&   reduction(+:likelihood_sum)
      do parallelindex = 1, nthread
        kahan_val1 = 0.0_fp
        likelihood_sum_tmp = 0.0_fp
        particleloop: do j = parallelindex_start(parallelindex), parallelindex_end(parallelindex)
          likelihood_particle(j) = 0.0_fp
          narray_use_tmp = 0
          do i = 1, narray
            if(.not. result_exist(arrayindex(i))) cycle
            narray_use_tmp = narray_use_tmp + 1

            call greatcircle_dist(lat_array(arrayindex(i)), lon_array(arrayindex(i)), &
            &                     lat_particle(j), lon_particle(j), distance = dist(i), azimuth = az)
            az = az + pi
            if(az .ge. 2.0_fp * pi) az = az - 2.0_fp * pi
            az_diff = delta_az(az_obs(arrayindex(i)), az)
            az_weight_index = int(az_obs(arrayindex(i)) / sameaz_num) + 1
            likelihood_azweight = likelihood_weight(azweight_coef, sameaz_num2, az_weight(az_weight_index))
            likelihood_tmp = likelihood_modified(az_diff, sigma_azdiff2, likelihood_azweight)
            likelihood_particle(j) = likelihood_renew(likelihood_particle(j), likelihood_tmp)
            if(present(arrivaltime) .and. present(appvel) .and. present(origintime)) then
              ot_est(narray_use_tmp) = origintime_cal(arrivaltime(arrayindex(i)), dist(i), appvel(arrayindex(i)))
            endif
          enddo

          if(present(arrivaltime) .and. present(appvel) .and. present(origintime)) then
            call bubblesort(ot_est(1 : narray_use_tmp))
            call pickup_medianval(ot_est(1 : narray_use_tmp), origintime(j))
          
            narray_use_tmp = 0
            do i = 1, narray
              if(.not. result_exist(arrayindex(i))) cycle
              narray_use_tmp = narray_use_tmp + 1
  
              ot_diff = origintime(j) - ot_est(narray_use_tmp)
              dist_tmp = log(dist(i))
              if(dist_tmp .lt. 1.0_fp) dist_tmp = 1.0_fp
              likelihood_distweight = likelihood_weight(ot_coef, sigma_dist2, dist_tmp)
              likelihood_tmp = likelihood_modified(ot_diff, sigma_otdiff2, likelihood_distweight)
              likelihood_particle(j) = likelihood_renew(likelihood_particle(j), likelihood_tmp)
              if(likelihood_particle(j) * 0.0_fp .ne. 0.0_fp) then
                print *, i, j, lon_particle(j), lat_particle(j)
                error stop
              endif
            enddo
          endif

          kahan_val2 = likelihood_sum_tmp + likelihood_particle(j)
          if(abs(likelihood_sum_tmp) .ge. abs(likelihood_particle(j))) then
            kahan_val1 = kahan_val1 + (likelihood_sum_tmp - kahan_val2) + likelihood_particle(j)
          else
            kahan_val1 = kahan_val1 + (likelihood_particle(j) - kahan_val2) + likelihood_sum_tmp
          endif
          likelihood_sum_tmp = kahan_val2
        enddo particleloop
        likelihood_sum = likelihood_sum + (likelihood_sum_tmp + kahan_val1)
      enddo
      !$omp end do
      !$omp barrier

      if(likelihood_sum .le. 1.0e-38_fp) exit particlefilter
      normalize_likelihood = 1.0_fp / likelihood_sum

      !!Rule for exit particlefilter loop
      !$omp single
      if(maxval(likelihood_particle) .le. maxval_likelihood_particle) then
        particlefilter_count = particlefilter_count + 1
      else
        maxval_likelihood_particle = maxval(likelihood_particle)
      endif
      !$omp end single
      if(particlefilter_count .gt. niter) exit particlefilter
      if(iteration_count .gt. iteration_count_max) exit

      particle_probability(1) = likelihood_particle(1) * normalize_likelihood
      do i = 2, nparticle
        particle_probability(i) = particle_probability(i - 1) + likelihood_particle(i) * normalize_likelihood
      enddo

      !!Redistribute particle
      !$omp do private(threadnum, rnd, rnd1, rnd2)
      do j = 1, nparticle
        threadnum = 1
        !$ threadnum = omp_get_thread_num() + 1
        call gen_random_number(randomnumber_status(threadnum), rnd)
        do i = 1, nparticle
          if(rnd .le. particle_probability(i)) exit
        enddo
        if(i .gt. nparticle) i = nparticle
        call gen_random_number(randomnumber_status(threadnum), rnd1)
        call gen_random_number(randomnumber_status(threadnum), rnd2)
        rnd = sqrt(-2.0_fp * log(rnd1)) * cos(2.0_fp * pi * rnd2)
        lon_particle_new(j) = lon_particle(i) + rnd * sigma_particle
        rnd = sqrt(-2.0_fp * log(rnd1)) * sin(2.0_fp * pi * rnd2)
        lat_particle_new(j) = lat_particle(i) + rnd * sigma_particle
      enddo
      !$omp end do
      !$omp barrier
      lon_particle(1 : nparticle) = lon_particle_new(1 : nparticle)
      lat_particle(1 : nparticle) = lat_particle_new(1 : nparticle)
    enddo particlefilter  
    !$omp end parallel

    likelihood_particle(1 : nparticle) = likelihood_particle(1 : nparticle) * normalize_likelihood

    if(present(appvel_median) .and. present(appvel) .and. present(arrivaltime)) then
      allocate(appvel_tmp(1 : narray))
      appvel_tmp(1 : narray) = 0.0_fp
      narray_use_tmp = 0
      do i = 1, narray
        if(.not. result_exist(arrayindex(i))) cycle
        narray_use_tmp = narray_use_tmp + 1
        appvel_tmp(narray_use_tmp) = appvel(arrayindex(i))
      enddo
      call bubblesort(appvel_tmp(1 : narray_use_tmp))
      call pickup_medianval(appvel_tmp(1 : narray_use_tmp), appvel_median)
      deallocate(appvel_tmp)
    endif

    deallocate(dist, ot_est)
    return
  end subroutine particlefilter_search

  subroutine particlefilter_init(randomnumber_status, lon_particle, lat_particle)
    use nrtype, only : fp
    use aeluma_parameters
    use random_number
    use xorshift1024star
    !$ use omp_lib

    implicit none
    type(xorshift1024star_state), intent(inout) :: randomnumber_status(:)
    real(kind = fp),              intent(out)   :: lon_particle(1 : nparticle), lat_particle(1 : nparticle)
    integer                                     :: i, threadnum
    real(kind = fp)                             :: rnd

    threadnum = 1
    !$ threadnum = omp_get_thread_num() + 1
    do i = 1, nparticle
      call gen_random_number(randomnumber_status(threadnum), rnd)
      lon_particle(i) = lon_w + (lon_e - lon_w) * rnd
      call gen_random_number(randomnumber_status(threadnum), rnd)
      lat_particle(i) = lat_s + (lat_n - lat_s) * rnd
    enddo

    return
  end subroutine particlefilter_init



end module particlefilter
