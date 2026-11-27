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

! Conservative upper bound on the per-sampled-pair probability rate for a
! cell, using the largest particle mass and the diameter of its velocity box.
subroutine sidm_probability_rate_bound(sigma_over_m,mass_max,v_min,v_max, &
     & cell_volume,n_particles,rate_bound,invalid)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(kind(1.0d0)),intent(in)::sigma_over_m,mass_max,cell_volume
  real(kind(1.0d0)),dimension(1:3),intent(in)::v_min,v_max
  integer,intent(in)::n_particles
  real(kind(1.0d0)),intent(out)::rate_bound
  integer,intent(out)::invalid
  integer::n_pairs,dim
  real(kind(1.0d0))::v_diameter,probability
  real(kind(1.0d0))::v_span(1:3)

  rate_bound=0.0d0
  invalid=0
  if(.not.ieee_is_finite(sigma_over_m) .or. &
       .not.ieee_is_finite(mass_max) .or. &
       .not.ieee_is_finite(cell_volume) .or. &
       any(.not.ieee_is_finite(v_min)) .or. &
       any(.not.ieee_is_finite(v_max))) then
     invalid=1
     return
  end if
  if(n_particles<2 .or. sigma_over_m<0.0d0 .or. mass_max<=0.0d0 .or. &
       cell_volume<=0.0d0) then
     invalid=1
     return
  end if
  do dim=1,3
     if(v_max(dim)<v_min(dim)) then
        invalid=1
        return
     end if
  end do

  v_span=v_max-v_min
  v_diameter=sqrt(sum(v_span**2))
  if(.not.ieee_is_finite(v_diameter)) then
     invalid=1
     return
  end if
  n_pairs=n_particles/2
  call sidm_pair_probability(sigma_over_m,mass_max,mass_max,v_diameter, &
       & 1.0d0,cell_volume,n_particles,n_pairs,probability)
  if(.not.ieee_is_finite(probability)) then
     invalid=1
     return
  end if
  if(probability<0.0d0) then
     invalid=1
     return
  end if
  rate_bound=probability*(1.0d0+64.0d0*epsilon(1.0d0))
  if(.not.ieee_is_finite(rate_bound)) then
     invalid=1
     rate_bound=0.0d0
  end if
end subroutine sidm_probability_rate_bound

! Code-unit timestep that targets half the hard probability ceiling.
subroutine sidm_probability_timestep(rate_bound,scale_t,courant,dt_code,invalid)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(kind(1.0d0)),intent(in)::rate_bound,scale_t,courant
  real(kind(1.0d0)),intent(out)::dt_code
  integer,intent(out)::invalid
  real(kind(1.0d0))::p_target

  dt_code=0.0d0
  invalid=0
  if(any(.not.ieee_is_finite((/rate_bound,scale_t,courant/)))) then
     invalid=1
     return
  end if
  if(rate_bound<=0.0d0 .or. scale_t<=0.0d0 .or. &
       courant<=0.0d0 .or. courant>1.0d0) then
     invalid=1
     return
  end if
  p_target=0.5d0*min(0.1d0,courant)
  dt_code=p_target/(rate_bound*scale_t)
  if(.not.ieee_is_finite(dt_code)) then
     invalid=1
     dt_code=0.0d0
     return
  end if
  if(dt_code<=0.0d0) then
     invalid=1
     dt_code=0.0d0
  end if
end subroutine sidm_probability_timestep

! Choose one collision-operator substep without changing the gravity step.
! Intervals admitted by the strict legacy ceiling retain one-sweep behavior;
! larger intervals are split toward half the hard probability ceiling.
subroutine sidm_choose_collision_substep(dt_remaining_code,rate_bound, &
     & scale_t,courant,dt_sub_code,p_remaining,p_substep,invalid)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(kind(1.0d0)),intent(in)::dt_remaining_code,rate_bound,scale_t,courant
  real(kind(1.0d0)),intent(out)::dt_sub_code,p_remaining,p_substep
  integer,intent(out)::invalid
  real(kind(1.0d0))::p_limit,p_limit_safe,p_target,dt_target

  dt_sub_code=0.0d0
  p_remaining=0.0d0
  p_substep=0.0d0
  invalid=0
  if(any(.not.ieee_is_finite((/dt_remaining_code,rate_bound,scale_t,courant/)))) then
     invalid=1
     return
  end if
  if(dt_remaining_code<0.0d0 .or. rate_bound<0.0d0 .or. &
       scale_t<=0.0d0 .or. courant<=0.0d0 .or. courant>1.0d0) then
     invalid=1
     return
  end if

  p_limit=min(0.1d0,courant)
  p_limit_safe=p_limit*(1.0d0-64.0d0*epsilon(1.0d0))
  dt_sub_code=dt_remaining_code
  if(dt_remaining_code<=0.0d0 .or. rate_bound<=0.0d0) return

  p_remaining=rate_bound*dt_remaining_code*scale_t
  if(.not.ieee_is_finite(p_remaining)) then
     invalid=1
     dt_sub_code=0.0d0
     return
  end if
  if(p_remaining<0.0d0) then
     invalid=1
     dt_sub_code=0.0d0
     return
  end if

  if(p_remaining>p_limit_safe) then
     p_target=0.5d0*p_limit
     dt_target=p_target/(rate_bound*scale_t)
     if(.not.ieee_is_finite(dt_target)) then
        invalid=1
        dt_sub_code=0.0d0
        return
     end if
     if(dt_target<=0.0d0) then
        invalid=1
        dt_sub_code=0.0d0
        return
     end if
     dt_sub_code=min(dt_remaining_code,dt_target)
  end if

  if(.not.ieee_is_finite(dt_sub_code)) then
     invalid=1
     dt_sub_code=0.0d0
     return
  end if
  if(dt_sub_code<=0.0d0) then
     invalid=1
     dt_sub_code=0.0d0
     return
  end if
  p_substep=rate_bound*dt_sub_code*scale_t
  if(.not.ieee_is_finite(p_substep)) then
     invalid=1
     dt_sub_code=0.0d0
     return
  end if
  if(p_substep<0.0d0 .or. p_substep>p_limit_safe) then
     invalid=1
     dt_sub_code=0.0d0
  end if
end subroutine sidm_choose_collision_substep

! Refuse rank-dependent endpoints or scale factors before entering a sweep.
subroutine sidm_collision_schedule_agrees(dt_min,dt_max,scale_t_min, &
     & scale_t_max,courant_min,courant_max,consistent)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(kind(1.0d0)),intent(in)::dt_min,dt_max,scale_t_min,scale_t_max
  real(kind(1.0d0)),intent(in)::courant_min,courant_max
  integer,intent(out)::consistent

  consistent=0
  if(any(.not.ieee_is_finite((/dt_min,dt_max,scale_t_min,scale_t_max, &
       & courant_min,courant_max/)))) return
  if(dt_min<0.0d0 .or. scale_t_min<=0.0d0 .or. &
       courant_min<=0.0d0 .or. courant_max>1.0d0) return
  if(dt_max<dt_min .or. scale_t_max<scale_t_min .or. &
       courant_max<courant_min) return
  if(dt_max>dt_min .or. scale_t_max>scale_t_min .or. &
       courant_max>courant_min) return
  consistent=1
end subroutine sidm_collision_schedule_agrees
