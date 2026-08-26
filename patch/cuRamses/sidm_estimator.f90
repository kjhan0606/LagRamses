!################################################################
! SIDM Monte Carlo estimator helper
!################################################################
subroutine sidm_pair_probability(sigma_over_m,mass_one,mass_two, &
     & v_rel,dt_phys,cell_volume,n_particles,n_pairs,probability)
  implicit none

  real(kind(1.0d0)),intent(in)::sigma_over_m,mass_one,mass_two
  real(kind(1.0d0)),intent(in)::v_rel,dt_phys,cell_volume
  integer,intent(in)::n_particles,n_pairs
  real(kind(1.0d0)),intent(out)::probability
  real(kind(1.0d0))::representative_mass

  representative_mass = 0.5d0*(mass_one+mass_two)
  probability = sigma_over_m*representative_mass*v_rel*dt_phys &
       /cell_volume*dble(n_particles-1)*dble(n_particles) &
       /dble(2*n_pairs)
end subroutine sidm_pair_probability
