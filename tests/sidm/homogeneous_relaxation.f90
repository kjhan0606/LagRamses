! Controlled environment for the complete production sub_sidm_scatter routine.
! Uniform periodic leaves only. No gravity, AMR tree or MPI decomposition.
module amr_parameters
  use iso_fortran_env, only: real64
  implicit none
  integer,parameter::dp=real64
end module

module amr_commons
  use amr_parameters
  implicit none
  integer,parameter::twotondim=8,icoarse_min=0,icoarse_max=0
  ! ncputot is only the RNG bank's leading dimension here, not an MPI run.
  integer::ncpu=256
  integer::myid=1,sidm_npart_min=2,sidm_nstates=2
  integer,allocatable::son(:),next(:)
  real(dp)::boxlen=1d0,aexp=1d0,dtnew(8)=0d0
  real(dp)::sidm_epsilon=0.01d0,sidm_vrel_max=0d0
  real(dp)::sidm_cross_section=1d0,sidm_v0=1000d0,sidm_power=-1d0
  real(dp)::sidm_a_transition=0.5d0,sidm_sigma_ratio=1d0,sidm_a_width=0.05d0
  real(dp)::sidm_mchi=1d0,sidm_energy(0:9)=0d0,sidm_fdiss=0d0
  logical::sidm_inelastic=.false.
  character(16)::sidm_type='power_law',sidm_a_type='none',sidm_angular='isotropic'
contains
  integer function ICELL_OF(grid,subcell)
    integer,intent(in)::grid,subcell
    ICELL_OF=8*(grid-1)+subcell
  end function
end module

module pm_commons
  use amr_parameters
  implicit none
  integer,parameter::PTYPE_DM=0,PTYPE_STAR=3,PTYPE_SINK=4
  integer,allocatable::numbp(:),headp(:),nextp(:),idp(:),ptypep(:),particle_cell(:)
  real(dp),allocatable::mp(:),vp(:,:),xp(:,:)
contains
  integer function ptype_to_isidm_state(ptype)
    integer,intent(in)::ptype
    ptype_to_isidm_state=ptype
  end function
  integer function isidm_state_to_ptype(state)
    integer,intent(in)::state
    isidm_state_to_ptype=state
  end function
end module

module relaxation_diagnostics
  use amr_parameters
  implicit none
  integer::sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down,sidm_n_vrel_reject
  integer::sidm_n_leaf_occupied,sidm_n_leaf_active,sidm_n_dm_leaf
  real(dp)::sidm_dp(3),sidm_dEk,sidm_dp_max,sidm_dEk_max
  real(dp)::sidm_Pmax_level,sidm_dEdiss,sidm_expected_rate
  common /sidm_diag/ sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down,sidm_n_vrel_reject
  common /sidm_cons/ sidm_dp,sidm_dEk,sidm_dp_max,sidm_dEk_max
  common /sidm_pmax_c/ sidm_Pmax_level
  common /sidm_diss/ sidm_dEdiss
  common /sidm_estimator_counts/ sidm_n_leaf_occupied,sidm_n_leaf_active,sidm_n_dm_leaf
  common /sidm_estimator_rate/ sidm_expected_rate
contains
  subroutine reset_diagnostics()
    sidm_n_scatter=0; sidm_n_pairs=0; sidm_n_up=0; sidm_n_down=0
    sidm_n_vrel_reject=0; sidm_n_leaf_occupied=0; sidm_n_leaf_active=0; sidm_n_dm_leaf=0
    sidm_dp=0d0; sidm_dEk=0d0; sidm_dp_max=0d0; sidm_dEk_max=0d0
    sidm_Pmax_level=0d0; sidm_dEdiss=0d0; sidm_expected_rate=0d0
  end subroutine
end module

program homogeneous_relaxation
  use amr_commons
  use pm_commons
  use random
  use relaxation_diagnostics
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  integer::npart,level,nside,ngrid,ncell,nstep,seed_id,i,j,step
  integer,allocatable::seed_bank(:,:)
  integer::thread_seeds(IRandNumSize,0:0),coordinate(3),grid,subcell
  integer,allocatable::occupancy(:)
  integer(kind=8)::total_scatter=0,total_pairs=0
  real(dp)::dt,tend,shift,q0,e0,q,e,momentum(3),energy_ledger=0d0
  real(dp)::p_peak=0d0,pair_density_peak=0d0,low_occupancy_peak=0d0
  real(dp)::uniform,meanv(3),pair_density,low_occupancy,ledger_error
  character(80)::arg
  call get_command_argument(1,arg); read(arg,*)npart
  call get_command_argument(2,arg); read(arg,*)level
  call get_command_argument(3,arg); read(arg,*)dt
  call get_command_argument(4,arg); read(arg,*)tend
  call get_command_argument(5,arg); read(arg,*)seed_id
  call get_command_argument(6,sidm_angular)
  call get_command_argument(7,arg); read(arg,*)sidm_epsilon
  call get_command_argument(8,arg); read(arg,*)sidm_cross_section
  call get_command_argument(9,arg); read(arg,*)sidm_fdiss
  call get_command_argument(10,arg); read(arg,*)shift
  if(npart<2.or.level<1.or.level>4.or.dt<=0d0.or.tend<=0d0)error stop 1
  nstep=nint(tend/dt)
  if(nstep<1.or.abs(nstep*dt-tend)>1d-10.or.mod(nstep,10)/=0)error stop 2
  nside=2**level; ngrid=(nside/2)**3; ncell=nside**3
  allocate(son(ncell),next(ngrid),numbp(ngrid),headp(ngrid),occupancy(ncell))
  allocate(nextp(npart),idp(npart),ptypep(npart),particle_cell(npart))
  allocate(mp(npart),vp(npart,3),xp(npart,3))
  son=0; mp=1d0/npart; ptypep=PTYPE_DM
  do i=1,ngrid
    next(i)=i+1
  enddo
  next(ngrid)=0
  ! Use valid, separated production RNG substreams, not an oversized base-4096
  ! digit as StartVal. Each seed digit must stay in [0,4095].
  if(seed_id<1.or.seed_id>ncpu)error stop 5
  allocate(seed_bank(ncpu,IRandNumSize))
  call rans(ncpu,0,seed_bank)
  if(any(seed_bank<0).or.any(seed_bank>4095))error stop 6
  thread_seeds(:,0)=seed_bank(seed_id,:)
  do i=1,npart
    idp(i)=i
    do j=1,3
      call ranf(thread_seeds(:,0),uniform)
      if(uniform<0d0.or.uniform>=1d0)error stop 7
      xp(i,j)=uniform
      call gaussdev(thread_seeds(:,0),vp(i,j))
    enddo
  enddo
  vp(:,2:3)=0.5d0*vp(:,2:3)
  do j=1,3
    meanv(j)=sum(vp(:,j))/npart
    vp(:,j)=vp(:,j)-meanv(j)
  enddo
  q0=sum(mp*(vp(:,1)**2-0.5d0*(vp(:,2)**2+vp(:,3)**2)))
  e0=0.5d0*sum(mp*(vp(:,1)**2+vp(:,2)**2+vp(:,3)**2))
  dtnew(level)=dt
  write(*,'(A)') '# step time Q/Q0 E/E0 momentum_error ledger_error P_peak pair_density_peak lowN_mass_peak Q0 events pairs'
  call print_row(0)
  do step=1,nstep
    ! Periodic free streaming, with velocity and time units consistent with units().
    ! Fixed translation tests the sensitivity to the uniform cell boundaries.
    xp=modulo(xp+0.5d0*dt*vp,1d0)
    headp=0; numbp=0; occupancy=0
    do i=1,npart
      coordinate=int(modulo(xp(i,:)+shift,1d0)*nside)
      grid=1+coordinate(1)/2+(nside/2)*(coordinate(2)/2)+(nside/2)**2*(coordinate(3)/2)
      subcell=1+mod(coordinate(1),2)+2*mod(coordinate(2),2)+4*mod(coordinate(3),2)
      particle_cell(i)=subcell
      occupancy(ICELL_OF(grid,subcell))=occupancy(ICELL_OF(grid,subcell))+1
      nextp(i)=headp(grid); headp(grid)=i; numbp(grid)=numbp(grid)+1
    enddo
    pair_density=dble(ncell)*sum(dble(occupancy)*dble(occupancy-1))/dble(npart)**2
    low_occupancy=dble(sum(occupancy,mask=occupancy<2))/npart
    pair_density_peak=max(pair_density_peak,pair_density)
    low_occupancy_peak=max(low_occupancy_peak,low_occupancy)
    call reset_diagnostics()
    call sub_sidm_scatter(level,1,1,ngrid,thread_seeds,0,1)
    total_scatter=total_scatter+sidm_n_scatter; total_pairs=total_pairs+sidm_n_pairs
    p_peak=max(p_peak,sidm_Pmax_level)
    ! Production dEdiss is in physical energy units. m=1/N and v_scale=1e8.
    energy_ledger=energy_ledger+sidm_dEdiss/1d16
    if(sidm_n_vrel_reject/=0.or.p_peak>1d0)error stop 3
    xp=modulo(xp+0.5d0*dt*vp,1d0)
    if(mod(step,nstep/10)==0)call print_row(step)
  enddo
contains
  subroutine print_row(index)
    integer,intent(in)::index
    q=sum(mp*(vp(:,1)**2-0.5d0*(vp(:,2)**2+vp(:,3)**2)))
    e=0.5d0*sum(mp*(vp(:,1)**2+vp(:,2)**2+vp(:,3)**2))
    do j=1,3
      momentum(j)=sum(mp*vp(:,j))
    enddo
    ledger_error=abs(e+energy_ledger-e0)/e0
    if(.not.all(ieee_is_finite(vp)).or.ledger_error>2d-11)error stop 4
    write(*,'(I8,9ES24.15,2I14)')index,index*dt,q/q0,e/e0, &
      maxval(abs(momentum))/sqrt(2d0*e0),ledger_error,p_peak,pair_density_peak, &
      low_occupancy_peak,q0,total_scatter,total_pairs
  end subroutine
end program

subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_parameters
  implicit none
  real(dp),intent(out)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  scale_l=1d0; scale_t=1d-8; scale_d=1d0; scale_v=1d8
  scale_nH=1d0; scale_T2=1d0
end subroutine

integer function cell_index_from_part(particle,grid,level)
  use pm_commons,only:particle_cell
  implicit none
  integer,intent(in)::particle,grid,level
  cell_index_from_part=particle_cell(particle)
end function
