module sidm_kinetic_transport_math
  implicit none
  private
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: i8=selected_int_kind(18)
  public :: sidm_kinetic_crossing_range
  public :: sidm_kinetic_mode_flux
  public :: sidm_kinetic_event_source
  public :: sidm_kinetic_flux_decomposition
  public :: sidm_kinetic_bootstrap_allowed
  public :: sidm_kinetic_interval_matches
contains

  pure subroutine sidm_kinetic_crossing_range(x0,x1,dx,kfirst,klast,direction)
    real(dp), intent(in) :: x0,x1,dx
    integer(i8), intent(out) :: kfirst,klast
    integer, intent(out) :: direction

    kfirst=0_i8
    klast=-1_i8
    direction=0
    if(dx<=0.0d0) return
    if(x1>x0) then
       kfirst=floor(x0/dx,kind=i8)+1_i8
       klast=floor(x1/dx,kind=i8)
       direction=1
    else if(x1<x0) then
       kfirst=floor(x0/dx,kind=i8)
       klast=floor(x1/dx,kind=i8)+1_i8
       direction=-1
    end if
  end subroutine sidm_kinetic_crossing_range

  pure subroutine sidm_kinetic_mode_flux(x0,x1,velocity_x,energy,dt, &
       & box_length,wavenumber,flux,flux_abs)
    real(dp), intent(in) :: x0,x1,velocity_x,energy,dt,box_length,wavenumber
    real(dp), intent(out) :: flux(2),flux_abs
    real(dp) :: argument,sinc_value,amplitude,phase

    argument=0.5d0*wavenumber*(x1-x0)
    if(abs(argument)<1.0d-4) then
       sinc_value=1.0d0-argument**2/6.0d0+argument**4/120.0d0
    else
       sinc_value=sin(argument)/argument
    end if
    amplitude=energy*velocity_x*dt/box_length**3*sinc_value
    phase=0.5d0*wavenumber*(x0+x1)
    flux(1)=amplitude*cos(phase)
    flux(2)=-amplitude*sin(phase)
    flux_abs=abs(energy*velocity_x*dt)/box_length**3
  end subroutine sidm_kinetic_mode_flux

  pure subroutine sidm_kinetic_event_source(x1,x2,delta_k1,delta_k2, &
       & volume,wavenumber,source,source_abs)
    real(dp), intent(in) :: x1,x2,delta_k1,delta_k2,volume,wavenumber
    real(dp), intent(out) :: source(2),source_abs
    real(dp) :: phase1,phase2

    phase1=wavenumber*x1
    phase2=wavenumber*x2
    source(1)=(delta_k1*cos(phase1)+delta_k2*cos(phase2))/volume
    source(2)=-(delta_k1*sin(phase1)+delta_k2*sin(phase2))/volume
    source_abs=(abs(delta_k1)+abs(delta_k2))/volume
  end subroutine sidm_kinetic_event_source

  pure subroutine sidm_kinetic_flux_decomposition(raw,volume,velocity, &
       & centered,heat_flux,advective,pressure_work,total_flux)
    real(dp), intent(in) :: raw(12),volume
    real(dp), intent(out) :: velocity(3),centered(6)
    real(dp), intent(out) :: heat_flux,advective,pressure_work,total_flux
    real(dp) :: mass,q_moment

    mass=raw(1)
    velocity=0.0d0
    centered=0.0d0
    heat_flux=0.0d0
    advective=0.0d0
    pressure_work=0.0d0
    total_flux=0.0d0
    if(mass<=0.0d0 .or. volume<=0.0d0) return

    velocity=raw(2:4)/mass
    centered(1)=raw(5)-raw(2)**2/mass
    centered(2)=raw(6)-raw(3)**2/mass
    centered(3)=raw(7)-raw(4)**2/mass
    centered(4)=raw(8)-raw(2)*raw(3)/mass
    centered(5)=raw(9)-raw(2)*raw(4)/mass
    centered(6)=raw(10)-raw(3)*raw(4)/mass
    q_moment=raw(11)-2.0d0*velocity(1)*raw(12) &
         & -2.0d0*(velocity(1)*raw(5)+velocity(2)*raw(8)+ &
         & velocity(3)*raw(9))+2.0d0*mass*velocity(1)*sum(velocity**2)
    heat_flux=0.5d0*q_moment/volume
    advective=0.5d0*(mass*sum(velocity**2)+sum(centered(1:3)))/ &
         & volume*velocity(1)
    pressure_work=(centered(1)*velocity(1)+centered(4)*velocity(2)+ &
         & centered(5)*velocity(3))/volume
    total_flux=raw(11)/(2.0d0*volume)
  end subroutine sidm_kinetic_flux_decomposition

  pure logical function sidm_kinetic_bootstrap_allowed(dt,time_now, &
       & initial_written,previous_drift_valid,bootstrap_seen,scatter_calls)
    real(dp), intent(in) :: dt,time_now
    logical, intent(in) :: initial_written,previous_drift_valid,bootstrap_seen
    integer(i8), intent(in) :: scatter_calls

    sidm_kinetic_bootstrap_allowed=dt==0.0d0 .and. time_now==0.0d0 .and. &
         & initial_written .and. .not.previous_drift_valid .and. &
         & .not.bootstrap_seen .and. scatter_calls==0_i8
  end function sidm_kinetic_bootstrap_allowed

  pure logical function sidm_kinetic_interval_matches(dt_min,dt_max, &
       & drift_min,drift_max,rtol)
    real(dp), intent(in) :: dt_min,dt_max,drift_min,drift_max,rtol
    real(dp) :: scale

    scale=max(abs(dt_min),abs(dt_max),abs(drift_min),abs(drift_max), &
         & tiny(1.0d0))
    sidm_kinetic_interval_matches= &
         & abs(dt_min-drift_min)<=rtol*scale .and. &
         & abs(dt_max-drift_max)<=rtol*scale
  end function sidm_kinetic_interval_matches

end module sidm_kinetic_transport_math
