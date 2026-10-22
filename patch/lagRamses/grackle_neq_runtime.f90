module grackle_neq_runtime
  use iso_c_binding
  use amr_parameters, only: dp,ndim
  use hydro_parameters, only: igrackle
  implicit none
  private
  public :: grackle_neq_initialize,grackle_neq_seed,grackle_neq_step,grackle_neq_thermal_energy
  interface
     integer(c_int) function c_init(path,uv,shield,h2,cosmo,d,l,t,a,z_reion) bind(C,name='lag_grackle_neq_init')
       import
       character(c_char),intent(in)::path(*)
       integer(c_int),value::uv,shield,h2,cosmo
       real(c_double),value::d,l,t,a,z_reion
     end function
     integer(c_int) function c_step(n,rho,u,s,z,len,cool,d,l,t,a,dt) bind(C,name='lag_grackle_neq_step')
       import
       integer(c_int),value::n
       real(c_double),intent(inout)::rho(*),u(*),s(*),z(*),len(*)
       integer(c_int),intent(in)::cool(*)
       real(c_double),value::d,l,t,a,dt
     end function
  end interface
contains
  pure real(dp) function grackle_neq_thermal_energy(rho,momentum,total,nonthermal)
    integer,parameter::qp=selected_real_kind(33,4931)
    real(dp),intent(in)::rho,momentum(ndim),total,nonthermal
    real(dp)::kinetic_dp,thermal_dp,scale_dp
    real(qp)::kinetic_qp,thermal_qp
    kinetic_dp=0.5d0*sum(momentum**2)/rho
    thermal_dp=total-kinetic_dp-nonthermal
    scale_dp=max(abs(total),abs(kinetic_dp+nonthermal))
    grackle_neq_thermal_energy=thermal_dp
    ! In cold, kinetic-energy-dominated gas, E-K may round to zero even
    ! though the stored binary64 conserved fields imply a positive residual.
    ! Grackle receives this residual as a standalone specific-energy value,
    ! where the positive value is representable; evaluate before subtracting
    ! in binary64 instead of manufacturing a floor or changing total energy.
    if(abs(thermal_dp)<=8d0*spacing(scale_dp))then
       kinetic_qp=0.5_qp*sum(real(momentum,qp)**2)/real(rho,qp)
       thermal_qp=real(total,qp)-kinetic_qp-real(nonthermal,qp)
       grackle_neq_thermal_energy=real(thermal_qp,dp)
    endif
  end function grackle_neq_thermal_energy

  subroutine grackle_neq_initialize(path,uv,shield,h2,cosmo,d,l,t,a,z_reion,ierr)
    character(*),intent(in)::path
    logical,intent(in)::uv,shield,cosmo
    integer,intent(in)::h2
    real(dp),intent(in)::d,l,t,a,z_reion
    integer,intent(out)::ierr
    ierr=c_init(trim(path)//c_null_char,merge(1,0,uv),merge(3,0,shield),h2,merge(1,0,cosmo), &
         d,l,t,a,z_reion)
  end subroutine
  subroutine grackle_neq_seed(h,he,xion,xmol,s)
    real(dp),intent(in)::h,he,xion,xmol
    real(dp),intent(out)::s(9)
    s=0;s(1)=h*(1-xion-xmol);s(2)=h*xion;s(4)=he;s(7)=h*xmol;s(9)=s(2)
  end subroutine
  subroutine grackle_neq_step(n,rho,u,s,z,length,cool,d,l,t,a,dt,ierr)
    integer,intent(in)::n
    real(dp),intent(inout)::rho(n),u(n),s(n,9),z(n),length(n)
    integer(c_int),intent(in)::cool(n)
    real(dp),intent(in)::d,l,t,a,dt
    integer,intent(out)::ierr
    ierr=c_step(n,rho,u,s,z,length,cool,d,l,t,a,dt)
  end subroutine
end module

subroutine initialize_grackle_neq_runtime
  use amr_commons
  use hydro_commons
  use grackle_neq_runtime
  implicit none
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2
  integer::ierr
  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
  ! RAMSES uses supercomoving time in cosmological runs (scale_t ~ a^2).
  ! Grackle's comoving unit contract requires a time unit fixed in time.
  if(cosmo)scale_t=scale_t/dble(aexp)**2
  call grackle_neq_initialize(grackle_neq_data_file,haardt_madau,self_shielding, &
       grackle_neq_h2_shielding,cosmo,scale_d,scale_l,scale_t,dble(aexp),dble(z_reion),ierr)
  if(ierr/=0)then
     write(*,*)'ERROR: Grackle NEQ initialization failed: ',ierr
     call clean_stop
  endif
end subroutine
