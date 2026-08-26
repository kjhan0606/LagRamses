program test_sidm_estimator
  implicit none

  integer::n_particles,n_pairs
  real(kind(1.0d0))::probability,expected,sigma_over_m
  real(kind(1.0d0))::mass_one,mass_two,v_rel,dt_phys,cell_volume
  external sidm_pair_probability

  sigma_over_m = 3.0d0
  mass_one = 2.0d35
  mass_two = 3.0d35
  v_rel = 2.0d7
  dt_phys = 1.0d12
  cell_volume = 1.0d69
  do n_particles=2,129
     n_pairs = n_particles/2
     call sidm_pair_probability(sigma_over_m,mass_one,mass_two, &
          & v_rel,dt_phys,cell_volume,n_particles,n_pairs,probability)
     expected = sigma_over_m*0.5d0*(mass_one+mass_two)*v_rel*dt_phys &
          /cell_volume*dble(n_particles-1)*dble(n_particles) &
          /dble(2*n_pairs)
     if(abs(probability-expected)>epsilon(1.0d0)*abs(expected)) &
          error stop 1
     write(*,'(2I6,2ES24.16)') n_particles,n_pairs,probability,expected
  end do
end program test_sidm_estimator
