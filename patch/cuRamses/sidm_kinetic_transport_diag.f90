module sidm_kinetic_transport_diag
  use amr_commons
  use pm_commons
  use sidm_kinetic_transport_math, only: sidm_kinetic_crossing_range, &
       & sidm_kinetic_mode_flux,sidm_kinetic_event_source, &
       & sidm_kinetic_flux_decomposition,sidm_kinetic_bootstrap_allowed, &
       & sidm_kinetic_interval_matches
  use poisson_commons, only: gravity_type,gravity_params
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int8,int32
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  implicit none
  private

  integer, parameter :: nslab=256,ntarget=21
  real(dp), parameter :: target_fraction(ntarget)=(/ &
       & 0.0d0,0.5d0,1.0d0,1.5d0,2.0d0,2.5d0,3.0d0, &
       & 3.5d0,4.0d0,4.5d0,5.0d0,5.5d0,6.0d0,6.5d0, &
       & 7.0d0,7.5d0,8.0d0,8.5d0,9.0d0,9.5d0,10.0d0/)
  real(dp), parameter :: fourier_k=0.25d0,fourier_rtol=1.0d-10
  real(dp), parameter :: closure_rtol=1.0d-10
  real(dp), save :: tau0_code=0.0d0
  real(dp), save :: source_interval(nslab)=0.0d0
  real(dp), save :: source_total(nslab)=0.0d0
  real(dp), save :: flux_interval(nslab)=0.0d0
  real(dp), save :: flux_total(nslab)=0.0d0
  real(dp), save :: energy_previous(nslab)=0.0d0
  real(dp), save :: raw_cache_mass(nslab)=0.0d0
  real(dp), save :: cache_velocity(3,nslab)=0.0d0
  real(dp), save :: cache_pressure(6,nslab)=0.0d0
  real(dp), save :: cache_heat_flux(nslab)=0.0d0
  real(dp), save :: cache_flux_advective(nslab)=0.0d0
  real(dp), save :: cache_flux_pressure_work(nslab)=0.0d0
  real(dp), save :: cache_flux_total(nslab)=0.0d0
  integer(kind=8), save :: events_interval=0_8,events_total=0_8
  real(dp), save :: mode_previous(2)=0.0d0
  real(dp), save :: source_mode_interval(2)=0.0d0
  real(dp), save :: source_mode_total(2)=0.0d0
  real(dp), save :: flux_mode_interval(2)=0.0d0
  real(dp), save :: flux_mode_total(2)=0.0d0
  real(dp), save :: source_mode_abs_interval=0.0d0
  real(dp), save :: flux_mode_abs_interval=0.0d0
  real(dp), save :: collision_energy_residual_max=0.0d0
  real(dp), save :: collision_momentum_residual_max=0.0d0
  real(dp), save :: collision_energy_scale_interval=0.0d0
  real(dp), save :: momentum_current(3)=0.0d0,momentum_previous(3)=0.0d0
  real(dp), save :: momentum_abs_current(3)=0.0d0,momentum_abs_previous(3)=0.0d0
  real(dp), save :: step_probability_max=0.0d0
  real(dp), save :: step_force_cell_max=0.0d0
  real(dp), save :: step_force_max(2)=0.0d0
  real(dp), save :: step_kick_max(2)=0.0d0
  real(dp), save :: step_velocity_delta_max(2)=0.0d0
  real(dp), save :: step_drift_exposure=0.0d0
  real(dp), save :: step_scatter_dt_min=huge(1.0d0),step_scatter_dt_max=0.0d0
  real(dp), save :: step_drift_dt_min=huge(1.0d0),step_drift_dt_max=0.0d0
  real(dp), save :: previous_drift_dt_min=0.0d0,previous_drift_dt_max=0.0d0
  real(dp), save :: step_preceding_drift_dt_min=0.0d0
  real(dp), save :: step_preceding_drift_dt_max=0.0d0
  real(dp), save :: expected_total_mass=0.0d0,current_total_mass=0.0d0
  real(dp), save :: expected_particle_mass=0.0d0
  real(dp), save :: dt_target_code=0.0d0
  real(dp), save :: current_particle_mass_min=huge(1.0d0)
  real(dp), save :: current_particle_mass_max=0.0d0
  integer(kind=8), save :: step_candidates=0_8,step_events=0_8
  integer(kind=8), save :: step_probability_bad=0_8
  integer(kind=8), save :: step_force_cells=0_8,step_force_bad=0_8
  integer(kind=8), save :: step_errors=0_8
  integer(kind=8), save :: step_scatter_errors=0_8,step_kick_errors=0_8
  integer(kind=8), save :: step_kick_count(2)=0_8,step_drift_count=0_8
  integer(kind=8), save :: step_scatter_calls=0_8
  integer(kind=8), save :: expected_particle_count=0_8
  integer(kind=8), save :: configured_particle_count=0_8
  integer, save :: target_index=1,nerror=0
  integer, save :: configured_level=0
  logical, save :: initialized=.false.,enabled=.false.,initial_written=.false.
  logical, save :: target_control_enabled=.false.,control_only=.false.
  logical, save :: state_export_enabled=.false.
  logical, save :: previous_drift_valid=.false.
  logical, save :: bootstrap_zero_dt_seen=.false.
  logical, save :: step_bootstrap_zero_dt=.false.
  integer, save :: output_unit=-1,fourier_unit=-1,step_unit=-1

  interface
     function sidm_kinetic_diag_rename(old_path,new_path) &
          & bind(C,name='rename') result(status)
       import :: c_char,c_int
       character(kind=c_char), intent(in) :: old_path(*)
       character(kind=c_char), intent(in) :: new_path(*)
       integer(c_int) :: status
     end function sidm_kinetic_diag_rename
  end interface

  public :: sidm_kinetic_diag_initial
  public :: sidm_kinetic_diag_export_input_state
  public :: sidm_kinetic_diag_clip_dt
  public :: sidm_kinetic_diag_scatter
  public :: sidm_kinetic_diag_drift
  public :: sidm_kinetic_diag_target
  public :: sidm_kinetic_diag_is_enabled
  public :: sidm_kinetic_diag_should_stop
  public :: sidm_kinetic_diag_probability
  public :: sidm_kinetic_diag_begin_scatter
  public :: sidm_kinetic_diag_forcecell
  public :: sidm_kinetic_diag_gravitykick
  public :: sidm_kinetic_diag_step_complete
  public :: sidm_kinetic_diag_orderly_stop

contains

  logical function sidm_kinetic_diag_is_enabled()
    call sidm_kinetic_diag_init()
    sidm_kinetic_diag_is_enabled=enabled
  end function sidm_kinetic_diag_is_enabled

  logical function sidm_kinetic_diag_should_stop()
    call sidm_kinetic_diag_init()
    sidm_kinetic_diag_should_stop=target_control_enabled .and. &
         & target_index>ntarget
  end function sidm_kinetic_diag_should_stop

  subroutine sidm_kinetic_diag_export_input_state()
    character(len=16) :: audit_flag
    character(len=16) :: byte_order
    character(len=80) :: filename
    integer :: status,length,ios,info,ipart,state_unit,manifest_unit
    integer(kind=8) :: local_count,global_count,file_bytes
    integer(int32) :: endian_word
    integer(int8) :: endian_bytes(4)
    logical, save :: audit_done=.false.
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    if(audit_done) return
    audit_done=.true.
    audit_flag=' '
    call get_environment_variable('SIDM_PARTICLE_READER_AUDIT',audit_flag, &
         & length=length,status=status)
    if(status==1) return
    if(status/=0 .or. length==0) &
         & call sidm_kinetic_diag_abort('invalid particle-reader audit flag')
    if(trim(audit_flag)=='0') return
    if(trim(audit_flag)/='1') &
         & call sidm_kinetic_diag_abort('invalid particle-reader audit flag')

    write(filename,'("sidm_reader_state_rank",I4.4,".bin")') myid
    open(newunit=state_unit,file=trim(filename),status='new',action='write', &
         & access='stream',form='unformatted',iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('create particle-reader state stream')
    local_count=0_8
    do ipart=1,npart
       if(idp(ipart)<=0) cycle
       if(.not.ieee_is_finite(mp(ipart)) .or. mp(ipart)<=0.0d0 .or. &
            & .not.all(ieee_is_finite(xp(ipart,1:3))) .or. &
            & .not.all(ieee_is_finite(vp(ipart,1:3)))) &
            & call sidm_kinetic_diag_abort('invalid particle-reader state')
       write(state_unit,iostat=ios) int(idp(ipart),kind=8),mp(ipart), &
            & xp(ipart,1:3),vp(ipart,1:3)
       if(ios/=0) call sidm_kinetic_diag_abort('write particle-reader state')
       local_count=local_count+1_8
    end do
    flush(state_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('flush particle-reader state')
    close(state_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('close particle-reader state')
    inquire(file=trim(filename),size=file_bytes,iostat=ios)
    if(ios/=0 .or. file_bytes/=64_8*local_count) &
         & call sidm_kinetic_diag_abort('particle-reader state byte count mismatch')

    global_count=local_count
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_count,global_count,1,MPI_INTEGER8,MPI_SUM, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce particle-reader state count')
#endif
    if(global_count<=0_8) &
         & call sidm_kinetic_diag_abort('empty particle-reader state audit')

    endian_word=1_int32
    endian_bytes=transfer(endian_word,endian_bytes)
    if(endian_bytes(1)==1_int8) then
       byte_order='little-endian'
    else
       byte_order='big-endian'
    end if
    if(myid==1) then
       open(newunit=manifest_unit,file='sidm_reader_state.manifest', &
            & status='new',action='write',iostat=ios)
       if(ios==0) write(manifest_unit,'(A)') 'format_version=1'
       if(ios==0) write(manifest_unit,'(A,A)') 'byte_order=',trim(byte_order)
       if(ios==0) write(manifest_unit,'(A,I0)') 'ranks=',ncpu
       if(ios==0) write(manifest_unit,'(A,I0)') 'global_count=',global_count
       if(ios==0) write(manifest_unit,'(A,I0)') 'record_bytes=',64
       if(ios==0) write(manifest_unit,'(A)') &
            & 'fields=id:int64,mass:float64,x:float64,y:float64,z:float64,'// &
            & 'vx:float64,vy:float64,vz:float64'
       if(ios==0) flush(manifest_unit,iostat=ios)
       if(ios==0) close(manifest_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('write particle-reader manifest')
       write(*,'(A,I0,A,I0)') &
            & 'SIDM reader-state audit: count=',global_count,' ranks=',ncpu
    end if
#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'barrier after particle-reader audit')
#endif
  end subroutine sidm_kinetic_diag_export_input_state

  subroutine sidm_kinetic_diag_initial()
    real(dp) :: state_energy(nslab),mode_state(2),scale_mode,tol_mode
    real(dp) :: momentum_state(3),momentum_abs_state(3)
    integer(kind=8) :: state_count(nslab)
    integer :: ios,info,pending_unit

    call sidm_kinetic_diag_init()
    if(.not.target_control_enabled .or. initial_written) return
    if(abs(t)>1.0d-12*max(1.0d0,tau0_code)) &
         & call sidm_kinetic_diag_abort('initial state is not at t=0')
    call sidm_kinetic_diag_census(state_energy,state_count,mode_state, &
         & momentum_state,momentum_abs_state)
    if(nerror/=0) &
         & call sidm_kinetic_diag_abort('invalid initial particle census')
    if(sum(state_count)/=configured_particle_count .or. &
         & current_particle_mass_min/=current_particle_mass_max) &
         & call sidm_kinetic_diag_abort('unexpected count or unequal particle masses')
    scale_mode=2.0d0*sum(state_energy)/boxlen**3
    tol_mode=fourier_rtol*max(scale_mode,tiny(1.0d0))
    energy_previous=state_energy
    mode_previous=mode_state
    momentum_previous=momentum_state
    momentum_abs_previous=momentum_abs_state
    expected_particle_count=sum(state_count)
    expected_total_mass=current_total_mass
    expected_particle_mass=current_particle_mass_min
    if(myid==1 .and. enabled) then
       open(newunit=pending_unit,file='sidm_kinetic_transport.pending', &
            & status='new',action='write',iostat=ios)
       if(ios==0) write(pending_unit,'(A)',iostat=ios) &
            & 'observer_status=PENDING'
       if(ios==0) close(pending_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('create pending marker')
       open(newunit=output_unit,file='sidm_kinetic_transport.csv', &
            & status='new',action='write',iostat=ios)
       if(ios==0) write(output_unit,'(A)',iostat=ios) &
            & 'target,step,time,tau0_fraction,slab,N,events_interval,'// &
            & 'events_cumulative,mass,ux,uy,uz,'// &
            & 'Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,flux_advective,'// &
            & 'flux_pressure_work,qx,flux_total_moment,flux_decomp_residual,'// &
            & 'K_bin,deltaE,source_interval,'// &
            & 'flux_left,flux_right,residual,source_cumulative,'// &
            & 'flux_cumulative'
       if(ios/=0) call sidm_kinetic_diag_abort('create output file')
       open(newunit=fourier_unit,file='sidm_kinetic_fourier.csv', &
            & status='new',action='write',iostat=ios)
       if(ios==0) write(fourier_unit,'(A)',iostat=ios) &
            & 'target,step,time,Ehat0_re,Ehat0_im,Ehat_re,Ehat_im,'// &
            & 'source_re,source_im,flux_re,flux_im,residual_re,residual_im,'// &
            & 'source_total_re,source_total_im,flux_total_re,flux_total_im,'// &
            & 'deltaK_total,source_k0,source_abs,flux_path_abs,'// &
            & 'collision_prepost_energy_scale,'// &
            & 'deltaPx,deltaPy,deltaPz,momentum_abs_x,momentum_abs_y,'// &
            & 'momentum_abs_z,event_energy_residual,event_momentum_residual,'// &
            & 'closure_scale,tolerance'
       if(ios/=0) call sidm_kinetic_diag_abort('create Fourier output')
       open(newunit=step_unit,file='sidm_kinetic_step_checks.csv', &
            & status='new',action='write',iostat=ios)
       if(ios==0) write(step_unit,'(A)',iostat=ios) &
            & 'step,time,candidates,events,pmax,probability_bad,'// &
            & 'force_cells,force_nonfinite,max_cell_force,kick_count_1,'// &
            & 'kick_count_2,max_force_1,max_force_2,max_kick_1,max_kick_2,'// &
            & 'max_velocity_delta_1,max_velocity_delta_2,drift_count,'// &
         & 'drift_exposure,scatter_calls,scatter_dt_min,scatter_dt_max,'// &
         & 'preceding_drift_dt_min,preceding_drift_dt_max,'// &
         & 'current_drift_dt_min,current_drift_dt_max,step_status,'// &
         & 'step_errors,scatter_errors,kick_errors'
       if(ios/=0) call sidm_kinetic_diag_abort('create step summary')
    end if
    if(state_export_enabled) call sidm_kinetic_diag_export_state(0,0,t)
    if(enabled) then
       call sidm_kinetic_diag_write_rows(0,0,t,state_energy,state_count, &
            & state_energy,source_interval,flux_interval,source_total,flux_total, &
            & 0_8,0_8)
       call sidm_kinetic_diag_write_fourier(0,0,t,mode_state,mode_state, &
            & source_mode_interval,flux_mode_interval, &
            & source_mode_total,flux_mode_total,0.0d0,0.0d0, &
            & source_mode_abs_interval,flux_mode_abs_interval,0.0d0, &
            & (/0.0d0,0.0d0,0.0d0/),momentum_abs_state,0.0d0,0.0d0, &
            & scale_mode,tol_mode)
    end if
    initial_written=.true.
    target_index=2
  end subroutine sidm_kinetic_diag_initial

  subroutine sidm_kinetic_diag_clip_dt(dt_candidate)
    real(dp), intent(inout) :: dt_candidate
    real(dp) :: target_time,gap,tolerance

    call sidm_kinetic_diag_init()
    if(.not.target_control_enabled .or. target_index>ntarget) return
    dt_candidate=min(dt_candidate,dt_target_code)
    target_time=target_fraction(target_index)*tau0_code
    tolerance=1.0d-11*max(1.0d0,abs(target_time),abs(t))
    if(t>target_time+tolerance) &
         & call sidm_kinetic_diag_abort('missed observation target')
    gap=target_time-t
    if(gap>0.0d0 .and. gap<dt_candidate-tolerance) &
         & dt_candidate=gap
  end subroutine sidm_kinetic_diag_clip_dt

  subroutine sidm_kinetic_diag_target(ilevel,step,time_now)
    integer, intent(in) :: ilevel,step
    real(dp), intent(in) :: time_now
    real(dp) :: target_time,tolerance,state_energy(nslab),mode_state(2)
    real(dp) :: source_global(nslab),flux_global(nslab)
    real(dp) :: source_total_global(nslab),flux_total_global(nslab)
    real(dp) :: source_mode_global(2),source_mode_total_global(2)
    real(dp) :: flux_mode_global(2),flux_mode_total_global(2)
    real(dp) :: residual_mode(2),scale_mode,tol_mode
    real(dp) :: source_abs_global,flux_abs_global
    real(dp) :: collision_energy_scale_global
    real(dp) :: momentum_state(3),momentum_abs_state(3),delta_momentum(3)
    real(dp) :: delta_energy,source_k0,energy_conservation_scale
    real(dp) :: collision_e_residual_global,collision_p_residual_global
    real(dp) :: residual_slab,scale_slab,tol_slab
    real(dp) :: time_min,time_max
    integer(kind=8) :: state_count(nslab)
    integer(kind=8) :: event_interval_global,event_total_global
    integer :: info,error_global,ios,level_min,level_max
    integer :: target_min,target_max
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    call sidm_kinetic_diag_init()
    if(.not.target_control_enabled) return
    time_min=time_now
    time_max=time_now
    level_min=ilevel
    level_max=ilevel
    target_min=target_index
    target_max=target_index
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,time_min,1,MPI_DOUBLE_PRECISION, &
         & MPI_MIN,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce target time minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,time_max,1,MPI_DOUBLE_PRECISION, &
         & MPI_MAX,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce target time maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,level_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce target level minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,level_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce target level maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,target_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce target index minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,target_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce target index maximum')
#endif
    if(time_min/=time_max .or. level_min/=level_max .or. &
         & target_min/=target_max) &
         & call sidm_kinetic_diag_abort('rank time or target mismatch')
    if(ilevel/=levelmin .or. target_index>ntarget) return
    target_time=target_fraction(target_index)*tau0_code
    tolerance=1.0d-11*max(1.0d0,abs(target_time),abs(time_now))
    if(time_now<target_time-tolerance) return
    if(abs(time_now-target_time)>tolerance) &
         & call sidm_kinetic_diag_abort('observation target was crossed')

    if(state_export_enabled) &
         & call sidm_kinetic_diag_export_state(target_index-1,step,time_now)
    if(.not.enabled) then
       target_index=target_index+1
#ifndef WITHOUTMPI
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'synchronize control-only target state')
#endif
       if(target_index>ntarget) &
            & call sidm_kinetic_diag_write_control_marker(step,time_now)
#ifndef WITHOUTMPI
       if(target_index>ntarget) then
          call MPI_BARRIER(MPI_COMM_WORLD,info)
          call sidm_kinetic_diag_mpierr(info,'synchronize control-only completion')
       end if
#endif
       return
    end if

    call sidm_kinetic_diag_census(state_energy,state_count,mode_state, &
         & momentum_state,momentum_abs_state)
    source_global=source_interval
    flux_global=flux_interval
    source_total_global=source_total
    flux_total_global=flux_total
    source_mode_global=source_mode_interval
    source_mode_total_global=source_mode_total
    flux_mode_global=flux_mode_interval
    flux_mode_total_global=flux_mode_total
    source_abs_global=source_mode_abs_interval
    flux_abs_global=flux_mode_abs_interval
    collision_energy_scale_global=collision_energy_scale_interval
    collision_e_residual_global=collision_energy_residual_max
    collision_p_residual_global=collision_momentum_residual_max
    event_interval_global=events_interval
    event_total_global=events_total
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_global,nslab, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce interval collision source')
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux_global,nslab, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce interval streaming flux')
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_total_global,nslab, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce cumulative collision source')
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux_total_global,nslab, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce cumulative streaming flux')
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_mode_global,2, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce Fourier collision source')
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_mode_total_global,2, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce cumulative Fourier source')
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux_mode_global,2, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce Fourier streaming flux')
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux_mode_total_global,2, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce cumulative Fourier flux')
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_abs_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce absolute Fourier source scale')
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux_abs_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce absolute Fourier flux scale')
    call MPI_ALLREDUCE(MPI_IN_PLACE,collision_energy_scale_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce collision energy scale')
    call MPI_ALLREDUCE(MPI_IN_PLACE,collision_e_residual_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce collision energy residual')
    call MPI_ALLREDUCE(MPI_IN_PLACE,collision_p_residual_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce collision momentum residual')
    call MPI_ALLREDUCE(MPI_IN_PLACE,event_interval_global,1,MPI_INTEGER8, &
         & MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce interval collision count')
    call MPI_ALLREDUCE(MPI_IN_PLACE,event_total_global,1,MPI_INTEGER8, &
         & MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce cumulative collision count')
    error_global=nerror
    call MPI_ALLREDUCE(MPI_IN_PLACE,error_global,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce observer error count')
#else
    error_global=nerror
#endif
    if(.not.all(ieee_is_finite(state_energy)) .or. &
         & .not.all(ieee_is_finite(source_global)) .or. &
         & .not.all(ieee_is_finite(flux_global)) .or. &
         & .not.all(ieee_is_finite(source_total_global)) .or. &
         & .not.all(ieee_is_finite(flux_total_global)) .or. &
         & .not.all(ieee_is_finite(mode_state)) .or. &
         & .not.all(ieee_is_finite(source_mode_global)) .or. &
         & .not.all(ieee_is_finite(source_mode_total_global)) .or. &
         & .not.all(ieee_is_finite(flux_mode_global)) .or. &
         & .not.all(ieee_is_finite(flux_mode_total_global)) .or. &
         & .not.ieee_is_finite(source_abs_global) .or. &
         & .not.ieee_is_finite(flux_abs_global) .or. &
         & .not.ieee_is_finite(collision_energy_scale_global) .or. &
         & .not.all(ieee_is_finite(momentum_state)) .or. &
         & .not.all(ieee_is_finite(momentum_abs_state))) &
         & call sidm_kinetic_diag_abort('non-finite reduced observer fields')
    if(error_global/=0) &
         & call sidm_kinetic_diag_abort('invalid particle or path record')
    if(.not.sidm .and. (event_interval_global/=0_8 .or. &
         & any(source_global/=0.0d0) .or. &
         & any(source_mode_global/=0.0d0))) &
         & call sidm_kinetic_diag_abort('SIDM-off control has collision source')
    if(sum(state_count)/=expected_particle_count .or. &
         & abs(current_total_mass-expected_total_mass)> &
         & closure_rtol*expected_total_mass .or. &
         & abs(current_particle_mass_min-expected_particle_mass)> &
         & closure_rtol*expected_particle_mass .or. &
         & abs(current_particle_mass_max-expected_particle_mass)> &
         & closure_rtol*expected_particle_mass) &
         & call sidm_kinetic_diag_abort('particle inventory changed')
    do info=1,nslab
       residual_slab=state_energy(info)-energy_previous(info)- &
            & source_global(info)-flux_global(modulo(info-1,nslab)+1)+ &
            & flux_global(modulo(info,nslab)+1)
       scale_slab=max(abs(state_energy(info)),abs(energy_previous(info)), &
            & abs(source_global(info)),abs(flux_global(modulo(info-1,nslab)+1)), &
            & abs(flux_global(modulo(info,nslab)+1)),tiny(1.0d0))
       tol_slab=closure_rtol*scale_slab
       if(abs(residual_slab)>tol_slab) &
            & call sidm_kinetic_diag_abort('slab energy closure failed')
    end do
    residual_mode(1)=mode_state(1)-mode_previous(1)- &
         & fourier_k*flux_mode_global(2)-source_mode_global(1)
    residual_mode(2)=mode_state(2)-mode_previous(2)+ &
         & fourier_k*flux_mode_global(1)-source_mode_global(2)
    scale_mode=(sum(state_energy)+sum(energy_previous))/boxlen**3+ &
         & source_abs_global+fourier_k*flux_abs_global
    tol_mode=fourier_rtol*max(scale_mode,tiny(1.0d0))
    if(any(abs(residual_mode)>tol_mode)) &
         & call sidm_kinetic_diag_abort('Fourier energy closure failed')
    delta_energy=sum(state_energy)-sum(energy_previous)
    source_k0=sum(source_global)
    energy_conservation_scale=sum(state_energy)+sum(energy_previous)+ &
         & collision_energy_scale_global
    if(abs(source_k0)>1.0d-12*max(collision_energy_scale_global, &
         & tiny(1.0d0)) .or. abs(delta_energy)>1.0d-12* &
         & max(energy_conservation_scale,tiny(1.0d0))) &
         & call sidm_kinetic_diag_abort('elastic k=0 energy conservation failed')
    delta_momentum=momentum_state-momentum_previous
    if(any(abs(delta_momentum)>1.0d-12*max(momentum_abs_state, &
         & momentum_abs_previous,tiny(1.0d0)))) &
         & call sidm_kinetic_diag_abort('global momentum conservation failed')
    if(collision_e_residual_global>1.0d-12 .or. &
         & collision_p_residual_global>1.0d-12) &
         & call sidm_kinetic_diag_abort('elastic event conservation failed')
    call sidm_kinetic_diag_write_rows(target_index-1,step,time_now, &
         & state_energy,state_count,energy_previous,source_global, &
         & flux_global,source_total_global,flux_total_global, &
         & event_interval_global,event_total_global)
    call sidm_kinetic_diag_write_fourier(target_index-1,step,time_now, &
         & mode_previous,mode_state,source_mode_global,flux_mode_global, &
         & source_mode_total_global,flux_mode_total_global, &
         & delta_energy,source_k0,source_abs_global,flux_abs_global, &
         & collision_energy_scale_global, &
         & delta_momentum,momentum_abs_state,collision_e_residual_global, &
         & collision_p_residual_global,scale_mode,tol_mode)
    energy_previous=state_energy
    mode_previous=mode_state
    momentum_previous=momentum_state
    momentum_abs_previous=momentum_abs_state
    source_interval=0.0d0
    flux_interval=0.0d0
    events_interval=0_8
    source_mode_interval=0.0d0
    flux_mode_interval=0.0d0
    source_mode_abs_interval=0.0d0
    flux_mode_abs_interval=0.0d0
    collision_energy_scale_interval=0.0d0
    collision_energy_residual_max=0.0d0
    collision_momentum_residual_max=0.0d0
    target_index=target_index+1
    if(target_index>ntarget .and. myid==1) then
       close(output_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('close output file')
       output_unit=-1
       close(fourier_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('close Fourier output')
       fourier_unit=-1
    end if
#ifndef WITHOUTMPI
    if(target_index>ntarget) then
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'synchronize final observer output')
    end if
#endif
  end subroutine sidm_kinetic_diag_target

  subroutine sidm_kinetic_diag_begin_scatter(dt)
    real(dp), intent(in) :: dt
    logical :: accepted
    call sidm_kinetic_diag_init()
    if(.not.enabled) return
    if(.not.ieee_is_finite(dt) .or. dt<0.0d0) then
!$omp atomic update
       step_errors=step_errors+1_8
!$omp atomic update
       step_scatter_errors=step_scatter_errors+1_8
       return
    end if
!$omp critical(sidm_kinetic_scatter_interval)
    accepted=.true.
    if(dt==0.0d0) then
       accepted=sidm_kinetic_bootstrap_allowed(dt,t,initial_written, &
            & previous_drift_valid,bootstrap_zero_dt_seen,step_scatter_calls)
       if(accepted) then
          bootstrap_zero_dt_seen=.true.
          step_bootstrap_zero_dt=.true.
       else
          step_errors=step_errors+1_8
          step_scatter_errors=step_scatter_errors+1_8
       end if
    else if(.not.previous_drift_valid .or. .not.bootstrap_zero_dt_seen) then
       accepted=.false.
       step_errors=step_errors+1_8
       step_scatter_errors=step_scatter_errors+1_8
    end if
    if(accepted) then
       step_scatter_calls=step_scatter_calls+1_8
       step_scatter_dt_min=min(step_scatter_dt_min,dt)
       step_scatter_dt_max=max(step_scatter_dt_max,dt)
       if(previous_drift_valid) then
          step_preceding_drift_dt_min=previous_drift_dt_min
          step_preceding_drift_dt_max=previous_drift_dt_max
       end if
    end if
!$omp end critical(sidm_kinetic_scatter_interval)
  end subroutine sidm_kinetic_diag_begin_scatter

  subroutine sidm_kinetic_diag_step_complete(step,time_now,terminal_partial)
    integer, intent(in) :: step
    real(dp), intent(in) :: time_now
    logical, intent(in), optional :: terminal_partial
    integer(kind=8) :: counts(12)
    real(dp) :: maxima(8),exposure,dt_min(3),dt_max(3),expected_exposure
    integer :: partial_flag,partial_min,partial_max,info,ios,marker_unit
    integer :: dt_index
    integer :: pending_unit
    integer(kind=8) :: event_total_global
    integer :: step_min,step_max,previous_valid_flag,previous_valid_min
    integer :: previous_valid_max,bootstrap_flag,bootstrap_min,bootstrap_max
    integer :: bootstrap_history_flag,bootstrap_history_min,bootstrap_history_max
    integer(c_int) :: rename_status
    real(dp) :: time_min,time_max
    character(len=10) :: status
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    call sidm_kinetic_diag_init()
    if(.not.enabled) return
    counts=(/step_candidates,step_events,step_probability_bad, &
         & step_force_cells,step_force_bad,step_kick_count(1), &
         & step_kick_count(2),step_drift_count,step_scatter_calls,step_errors, &
         & step_scatter_errors,step_kick_errors/)
    maxima(1)=step_probability_max
    maxima(2)=step_force_cell_max
    maxima(3:4)=step_force_max
    maxima(5:6)=step_kick_max
    maxima(7:8)=step_velocity_delta_max
    exposure=step_drift_exposure
    dt_min=(/step_scatter_dt_min,step_preceding_drift_dt_min, &
         & step_drift_dt_min/)
    dt_max=(/step_scatter_dt_max,step_preceding_drift_dt_max, &
         & step_drift_dt_max/)
    partial_flag=0
    if(present(terminal_partial)) then
       if(terminal_partial) partial_flag=1
    end if
    partial_min=partial_flag
    partial_max=partial_flag
    step_min=step
    step_max=step
    time_min=time_now
    time_max=time_now
    previous_valid_flag=merge(1,0,previous_drift_valid)
    previous_valid_min=previous_valid_flag
    previous_valid_max=previous_valid_flag
    bootstrap_flag=merge(1,0,step_bootstrap_zero_dt)
    bootstrap_min=bootstrap_flag
    bootstrap_max=bootstrap_flag
    bootstrap_history_flag=merge(1,0,bootstrap_zero_dt_seen)
    bootstrap_history_min=bootstrap_history_flag
    bootstrap_history_max=bootstrap_history_flag
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,counts,12,MPI_INTEGER8,MPI_SUM, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step counts')
    call MPI_ALLREDUCE(MPI_IN_PLACE,maxima,8,MPI_DOUBLE_PRECISION,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step maxima')
    call MPI_ALLREDUCE(MPI_IN_PLACE,exposure,1,MPI_DOUBLE_PRECISION, &
         & MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step drift exposure')
    call MPI_ALLREDUCE(MPI_IN_PLACE,dt_min,3,MPI_DOUBLE_PRECISION,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce timestep minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,dt_max,3,MPI_DOUBLE_PRECISION,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce timestep maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,partial_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce partial-step flag minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,partial_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce partial-step flag maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,step_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step number minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,step_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step number maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,time_min,1,MPI_DOUBLE_PRECISION, &
         & MPI_MIN,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step time minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,time_max,1,MPI_DOUBLE_PRECISION, &
         & MPI_MAX,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce coarse-step time maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,previous_valid_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce prior drift coverage minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,previous_valid_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce prior drift coverage maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,bootstrap_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce zero-dt bootstrap minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,bootstrap_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce zero-dt bootstrap maximum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,bootstrap_history_min,1,MPI_INTEGER, &
         & MPI_MIN,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce bootstrap history minimum')
    call MPI_ALLREDUCE(MPI_IN_PLACE,bootstrap_history_max,1,MPI_INTEGER, &
         & MPI_MAX,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce bootstrap history maximum')
#endif
    if(partial_min/=partial_max) &
         & call sidm_kinetic_diag_abort('rank terminal-step mismatch')
    if(step_min/=step_max .or. time_min/=time_max) &
         & call sidm_kinetic_diag_abort('rank coarse-step mismatch')
    if(previous_valid_min/=previous_valid_max) &
         & call sidm_kinetic_diag_abort('rank prior drift state mismatch')
    if(bootstrap_min/=bootstrap_max) &
         & call sidm_kinetic_diag_abort('rank zero-dt bootstrap mismatch')
    if(bootstrap_history_min/=bootstrap_history_max) &
         & call sidm_kinetic_diag_abort('rank bootstrap history mismatch')
    if(bootstrap_history_min/=1) &
         & call sidm_kinetic_diag_abort('fresh-start zero-dt bootstrap was not recorded')
    if(counts(10)/=0_8) then
       if(myid==1) write(*,*) 'observer step errors total/scatter/kick=', &
            & counts(10),counts(11),counts(12),' step=',step, &
            & ' scatter_calls=',counts(9),' dt_min/max=',dt_min(1),dt_max(1), &
            & ' bootstrap/previous_valid=',bootstrap_max,previous_valid_min
       call sidm_kinetic_diag_abort('coarse-step observer checks failed')
    end if
    if(counts(3)/=0_8 .or. counts(5)/=0_8) then
       if(myid==1) write(*,*) 'observer probability/force checks bad/nonfinite=', &
            & counts(3),counts(5),' probability_max=',maxima(1), &
            & ' force_cells=',counts(4)
       call sidm_kinetic_diag_abort('coarse-step observer checks failed')
    end if
    if(counts(4)<=0_8 .or. any(maxima(2:8)/=0.0d0)) &
         & call sidm_kinetic_diag_abort('force-free pilot has a nonzero force or kick')
    if(sidm .and. counts(1)<=0_8 .and. bootstrap_max==0 .and. partial_max==0) &
         & call sidm_kinetic_diag_abort('active SIDM step has no candidate pairs')
    if(bootstrap_max==1 .and. &
         & (step/=1 .or. previous_valid_min/=0 .or. &
         & dt_min(1)/=0.0d0 .or. dt_max(1)/=0.0d0 .or. &
         & counts(2)/=0_8 .or. maxima(1)/=0.0d0)) &
         & call sidm_kinetic_diag_abort('invalid zero-dt bootstrap step')
    if(.not.sidm .and. (counts(1)/=0_8 .or. counts(2)/=0_8)) &
         & call sidm_kinetic_diag_abort('SIDM-off control recorded a collision')
    if(counts(9)/=int(ncpu,kind=8) .or. dt_min(1)/=dt_max(1)) &
         & call sidm_kinetic_diag_abort('missing synchronized scatter interval')
    if(bootstrap_max==0 .and. dt_min(1)<=0.0d0) &
         & call sidm_kinetic_diag_abort('nonpositive physical scatter interval')
    if(bootstrap_max==0 .and. partial_max==0) then
       do dt_index=1,3
          if(abs(dt_min(dt_index)-dt_target_code)> &
               & 1.0d-12*max(dt_target_code,tiny(1.0d0)) .or. &
               & abs(dt_max(dt_index)-dt_target_code)> &
               & 1.0d-12*max(dt_target_code,tiny(1.0d0))) &
               & call sidm_kinetic_diag_abort( &
               & 'measured interval differs from requested fixed timestep')
       end do
    end if
    if(previous_valid_min==1 .and. bootstrap_max==0) then
       if(.not.sidm_kinetic_interval_matches(dt_min(1),dt_max(1), &
            & dt_min(2),dt_max(2),1.0d-12)) &
            & call sidm_kinetic_diag_abort('scatter interval differs from preceding drift')
    end if
    if(partial_max==0) then
       if(counts(6)/=expected_particle_count .or. &
            & counts(7)/=expected_particle_count .or. &
            & counts(8)/=expected_particle_count) &
            & call sidm_kinetic_diag_abort('incomplete KDK or drift coverage')
       expected_exposure=dble(expected_particle_count)*dt_max(3)
       if(abs(exposure-expected_exposure)> &
            & 1.0d-12*max(expected_exposure,tiny(1.0d0)) .or. &
            & dt_min(3)<=0.0d0 .or. dt_min(3)/=dt_max(3)) &
            & call sidm_kinetic_diag_abort('drift exposure disagrees with particle count')
    else
       if(counts(6)/=expected_particle_count .or. counts(7)/=0_8 .or. &
            & counts(8)/=0_8 .or. exposure/=0.0d0) &
            & call sidm_kinetic_diag_abort('terminal sample is not sync-only')
    end if
    event_total_global=events_total
#ifndef WITHOUTMPI
    if(partial_max==1) then
       call MPI_ALLREDUCE(MPI_IN_PLACE,event_total_global,1,MPI_INTEGER8, &
            & MPI_SUM,MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce terminal collision count')
    end if
#endif
    status='PASS'
    if(bootstrap_max==1) status='BOOTSTRAP'
    if(partial_max==1) status='PARTIAL'
    if(myid==1) then
       write(step_unit,'(I0,",",ES24.16,",",I0,",",I0,",", &
            & ES24.16,",",I0,",",I0,",",I0,",",ES24.16,",", &
            & I0,",",I0,",",6(ES24.16,","),I0,",",ES24.16,",", &
            & I0,",",6(ES24.16,","),A,",",I0,",",I0,",",I0)', &
            & iostat=ios) step,time_now,counts(1),counts(2),maxima(1), &
            & counts(3),counts(4),counts(5),maxima(2),counts(6),counts(7), &
            & maxima(3:8),counts(8),exposure,counts(9),dt_min(1),dt_max(1), &
            & dt_min(2),dt_max(2),dt_min(3),dt_max(3),trim(status), &
            & counts(10),counts(11),counts(12)
       if(ios/=0) call sidm_kinetic_diag_abort('write coarse-step summary')
       flush(step_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('flush coarse-step summary')
       if(partial_max==1) then
          close(step_unit,iostat=ios)
          if(ios/=0) call sidm_kinetic_diag_abort('close coarse-step summary')
       end if
    end if
#ifndef WITHOUTMPI
    if(partial_max==1) then
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'synchronize final step summary')
    end if
#endif
    if(partial_max==1 .and. myid==1) then
       open(newunit=marker_unit, &
            & file='sidm_kinetic_transport.observer_validated.tmp', &
            & status='new',action='write',iostat=ios)
       if(ios==0) write(marker_unit,'(A)',iostat=ios) &
            & 'observer_status=VALIDATED;run_status=STOP_REQUESTED'
       if(ios==0) write(marker_unit,'(A,ES24.16,A,I0)',iostat=ios) &
            & 'terminal_time=',time_now,',terminal_step=',step
       if(ios==0) then
          if(sidm .and. event_total_global==0_8) then
             write(marker_unit,'(A)',iostat=ios) 'verdict=INCONCLUSIVE_ZERO_EVENTS'
          else if(sidm) then
             write(marker_unit,'(A,I0)',iostat=ios) 'verdict=ACTIVE_EVENTS,count=', &
                  & event_total_global
          else
             write(marker_unit,'(A)',iostat=ios) 'verdict=SIDM_OFF_CONTROL'
          end if
       end if
       if(ios==0) flush(marker_unit,iostat=ios)
       if(ios==0) close(marker_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('write observer validation marker')
       rename_status=sidm_kinetic_diag_rename( &
            & 'sidm_kinetic_transport.observer_validated.tmp'//c_null_char, &
            & 'sidm_kinetic_transport.observer_validated'//c_null_char)
       if(rename_status/=0_c_int) &
            & call sidm_kinetic_diag_abort('publish observer validation marker')
       open(newunit=pending_unit,file='sidm_kinetic_transport.pending', &
            & status='old',iostat=ios)
       if(ios==0) close(pending_unit,status='delete',iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('remove pending marker')
    end if
#ifndef WITHOUTMPI
    if(partial_max==1) then
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'synchronize completion marker')
    end if
#endif
    if(partial_max==0) then
       previous_drift_dt_min=dt_min(3)
       previous_drift_dt_max=dt_max(3)
       previous_drift_valid=.true.
    end if
    step_candidates=0_8
    step_events=0_8
    step_probability_bad=0_8
    step_probability_max=0.0d0
    step_force_cells=0_8
    step_force_bad=0_8
    step_force_cell_max=0.0d0
    step_kick_count=0_8
    step_kick_max=0.0d0
    step_force_max=0.0d0
    step_velocity_delta_max=0.0d0
    step_drift_count=0_8
    step_drift_exposure=0.0d0
    step_scatter_calls=0_8
    step_scatter_dt_min=huge(1.0d0)
    step_scatter_dt_max=0.0d0
    step_bootstrap_zero_dt=.false.
    step_preceding_drift_dt_min=0.0d0
    step_preceding_drift_dt_max=0.0d0
    step_drift_dt_min=huge(1.0d0)
    step_drift_dt_max=0.0d0
    step_errors=0_8
    step_scatter_errors=0_8
    step_kick_errors=0_8
  end subroutine sidm_kinetic_diag_step_complete

  subroutine sidm_kinetic_diag_scatter(x1,x2,m1,m2,v1_old,v2_old, &
       & v1_new,v2_new)
    real(dp), intent(in) :: x1(3),x2(3),m1,m2
    real(dp), intent(in) :: v1_old(3),v2_old(3),v1_new(3),v2_new(3)
    real(dp) :: dk1,dk2,source_mode(2),source_abs,volume
    real(dp) :: energy_scale,momentum_scale,event_energy_residual
    real(dp) :: event_momentum_residual,delta_p(3)
    integer :: b1,b2

    call sidm_kinetic_diag_init()
    if(.not.enabled) return
    if(m1<=0.0d0 .or. m2<=0.0d0 .or. &
         & .not.ieee_is_finite(m1) .or. .not.ieee_is_finite(m2) .or. &
         & .not.all(ieee_is_finite(x1)) .or. &
         & .not.all(ieee_is_finite(x2)) .or. &
         & .not.all(ieee_is_finite(v1_old)) .or. &
         & .not.all(ieee_is_finite(v2_old)) .or. &
         & .not.all(ieee_is_finite(v1_new)) .or. &
         & .not.all(ieee_is_finite(v2_new))) then
!$omp atomic update
       nerror=nerror+1
       return
    end if
    dk1=0.5d0*m1*(sum(v1_new**2)-sum(v1_old**2))
    dk2=0.5d0*m2*(sum(v2_new**2)-sum(v2_old**2))
    if(.not.ieee_is_finite(dk1) .or. .not.ieee_is_finite(dk2)) then
!$omp atomic update
       nerror=nerror+1
       return
    end if
    delta_p=m1*(v1_new-v1_old)+m2*(v2_new-v2_old)
    energy_scale=max(0.5d0*m1*sum(v1_old**2)+0.5d0*m2*sum(v2_old**2), &
         & 0.5d0*m1*sum(v1_new**2)+0.5d0*m2*sum(v2_new**2))
    if(.not.ieee_is_finite(energy_scale) .or. energy_scale<0.0d0) then
!$omp atomic update
       nerror=nerror+1
       return
    end if
    momentum_scale=max(sum(abs(m1*v1_old))+sum(abs(m2*v2_old)), &
         & sum(abs(m1*v1_new))+sum(abs(m2*v2_new)))
    event_energy_residual=abs(dk1+dk2)/max(energy_scale,tiny(1.0d0))
    event_momentum_residual=maxval(abs(delta_p))/ &
         & max(momentum_scale,tiny(1.0d0))
    if(.not.all(ieee_is_finite(delta_p)) .or. &
         & .not.ieee_is_finite(event_energy_residual) .or. &
         & .not.ieee_is_finite(event_momentum_residual)) then
!$omp atomic update
       nerror=nerror+1
       return
    end if
    if(event_energy_residual>1.0d-12 .or. &
         & event_momentum_residual>1.0d-12) then
!$omp atomic update
       nerror=nerror+1
    end if
    b1=sidm_kinetic_diag_bin(x1(1))
    b2=sidm_kinetic_diag_bin(x2(1))
    volume=boxlen**3
    call sidm_kinetic_event_source(x1(1),x2(1),dk1,dk2,volume, &
         & fourier_k,source_mode,source_abs)
!$omp critical(sidm_kinetic_source_update)
    source_interval(b1)=source_interval(b1)+dk1
    source_interval(b2)=source_interval(b2)+dk2
    source_total(b1)=source_total(b1)+dk1
    source_total(b2)=source_total(b2)+dk2
    events_interval=events_interval+1_8
    events_total=events_total+1_8
    source_mode_interval=source_mode_interval+source_mode
    source_mode_total=source_mode_total+source_mode
    source_mode_abs_interval=source_mode_abs_interval+source_abs
    collision_energy_residual_max=max(collision_energy_residual_max, &
         & event_energy_residual)
    collision_momentum_residual_max=max(collision_momentum_residual_max, &
         & event_momentum_residual)
    collision_energy_scale_interval=collision_energy_scale_interval+energy_scale
    step_events=step_events+1_8
!$omp end critical(sidm_kinetic_source_update)
  end subroutine sidm_kinetic_diag_scatter

  subroutine sidm_kinetic_diag_probability(probability)
    real(dp), intent(in) :: probability
    call sidm_kinetic_diag_init()
    if(.not.enabled) return
!$omp critical(sidm_kinetic_probability_update)
    step_candidates=step_candidates+1_8
    if(.not.ieee_is_finite(probability)) then
       step_probability_bad=step_probability_bad+1_8
    else
       step_probability_max=max(step_probability_max,probability)
       if(probability<0.0d0 .or. probability>0.05d0) &
            & step_probability_bad=step_probability_bad+1_8
    end if
!$omp end critical(sidm_kinetic_probability_update)
  end subroutine sidm_kinetic_diag_probability

  subroutine sidm_kinetic_diag_forcecell(nleaf,force_max,nonfinite_count)
    integer(kind=8), intent(in) :: nleaf,nonfinite_count
    real(dp), intent(in) :: force_max
    call sidm_kinetic_diag_init()
    if(.not.enabled) return
!$omp critical(sidm_kinetic_forcecell_update)
    step_force_cells=step_force_cells+nleaf
    step_force_bad=step_force_bad+nonfinite_count
    if(ieee_is_finite(force_max)) then
       step_force_cell_max=max(step_force_cell_max,force_max)
    else
       step_force_bad=step_force_bad+1_8
    end if
!$omp end critical(sidm_kinetic_forcecell_update)
  end subroutine sidm_kinetic_diag_forcecell

  subroutine sidm_kinetic_diag_gravitykick(site,force,kick,actual_delta)
    integer, intent(in) :: site
    real(dp), intent(in) :: force(3),kick(3),actual_delta(3)
    call sidm_kinetic_diag_init()
    if(.not.enabled) return
    if(site<1 .or. site>2) then
!$omp atomic update
       step_errors=step_errors+1_8
!$omp atomic update
       step_kick_errors=step_kick_errors+1_8
       return
    end if
!$omp critical(sidm_kinetic_kick_update)
    step_kick_count(site)=step_kick_count(site)+1_8
    if(.not.all(ieee_is_finite(force)) .or. &
         & .not.all(ieee_is_finite(kick)) .or. &
         & .not.all(ieee_is_finite(actual_delta))) then
       step_errors=step_errors+1_8
       step_kick_errors=step_kick_errors+1_8
    else
       step_force_max(site)=max(step_force_max(site),maxval(abs(force)))
       step_kick_max(site)=max(step_kick_max(site),maxval(abs(kick)))
       step_velocity_delta_max(site)=max(step_velocity_delta_max(site), &
            & maxval(abs(actual_delta)))
    end if
!$omp end critical(sidm_kinetic_kick_update)
  end subroutine sidm_kinetic_diag_gravitykick

  subroutine sidm_kinetic_diag_drift(mass,x0,x1,velocity,dt)
    real(dp), intent(in) :: mass,x0(3),x1(3),velocity(3),dt
    real(dp) :: energy,dx
    real(dp) :: mode_flux(2),mode_flux_abs
    integer(kind=8) :: kfirst,klast,k
    integer :: face,direction

    call sidm_kinetic_diag_init()
    if(.not.enabled) return
    if(mass<=0.0d0 .or. .not.ieee_is_finite(mass) .or. &
         & .not.all(ieee_is_finite(x0)) .or. &
         & .not.all(ieee_is_finite(x1)) .or. &
         & .not.all(ieee_is_finite(velocity)) .or. &
         & .not.ieee_is_finite(dt) .or. dt<=0.0d0) then
!$omp atomic update
       nerror=nerror+1
       return
    end if
    dx=boxlen/dble(nslab)
    energy=0.5d0*mass*sum(velocity**2)
    if(.not.ieee_is_finite(energy)) then
!$omp atomic update
       nerror=nerror+1
       return
    end if
    call sidm_kinetic_mode_flux(x0(1),x1(1),velocity(1),energy,dt, &
         & boxlen,fourier_k,mode_flux,mode_flux_abs)
    call sidm_kinetic_crossing_range(x0(1),x1(1),dx,kfirst,klast,direction)
    if(direction>0) then
       if(klast-kfirst>10000_8) then
!$omp atomic update
          nerror=nerror+1
          return
       end if
       do k=kfirst,klast
          face=1+int(modulo(k,int(nslab,kind=8)))
!$omp critical(sidm_kinetic_flux_update)
          flux_interval(face)=flux_interval(face)+energy
          flux_total(face)=flux_total(face)+energy
!$omp end critical(sidm_kinetic_flux_update)
       end do
    else if(direction<0) then
       if(kfirst-klast>10000_8) then
!$omp atomic update
          nerror=nerror+1
          return
       end if
       do k=kfirst,klast,-1_8
          face=1+int(modulo(k,int(nslab,kind=8)))
!$omp critical(sidm_kinetic_flux_update)
          flux_interval(face)=flux_interval(face)-energy
          flux_total(face)=flux_total(face)-energy
!$omp end critical(sidm_kinetic_flux_update)
       end do
    end if
!$omp critical(sidm_kinetic_flux_update)
    flux_mode_interval=flux_mode_interval+mode_flux
    flux_mode_total=flux_mode_total+mode_flux
    flux_mode_abs_interval=flux_mode_abs_interval+mode_flux_abs
    step_drift_count=step_drift_count+1_8
    step_drift_exposure=step_drift_exposure+dt
    step_drift_dt_min=min(step_drift_dt_min,dt)
    step_drift_dt_max=max(step_drift_dt_max,dt)
!$omp end critical(sidm_kinetic_flux_update)
  end subroutine sidm_kinetic_diag_drift

  subroutine sidm_kinetic_diag_census(energy,nparticle,mode,momentum, &
       & momentum_abs)
    real(dp), intent(out) :: energy(nslab)
    integer(kind=8), intent(out) :: nparticle(nslab)
    real(dp), intent(out) :: mode(2)
    real(dp), intent(out) :: momentum(3),momentum_abs(3)
    real(dp) :: raw(15,nslab),moment_global(15,nslab)
    integer(kind=8) :: count_local(nslab)
    real(dp) :: mode_local(2),mode_global(2),mass,kinetic,phase
    real(dp) :: mass_min_local,mass_max_local,mass_min_global,mass_max_global
    integer :: ip,b,info,igrid,jgrid,jpart,ipart
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    raw=0.0d0
    count_local=0_8
    mode_local=0.0d0
    mass_min_local=huge(1.0d0)
    mass_max_local=0.0d0
    if(numbl(myid,levelmin)>0) then
       igrid=headl(myid,levelmin)
       do jgrid=1,numbl(myid,levelmin)
          if(igrid<=0 .or. igrid>ngridmax) then
             nerror=nerror+1
             exit
          end if
          ipart=headp(igrid)
          do jpart=1,numbp(igrid)
             if(ipart<=0 .or. ipart>npartmax) then
                nerror=nerror+1
                exit
             end if
             if(idp(ipart)>0 .and. ptypep(ipart)==PTYPE_DM) then
                mass=mp(ipart)
                if(.not.ieee_is_finite(mass) .or. mass<=0.0d0 .or. &
                     & .not.all(ieee_is_finite(xp(ipart,1:3))) .or. &
                     & .not.all(ieee_is_finite(vp(ipart,1:3)))) then
                   nerror=nerror+1
                   ipart=nextp(ipart)
                   cycle
                end if
                b=sidm_kinetic_diag_bin(xp(ipart,1))
                raw(1,b)=raw(1,b)+mass
                raw(2:4,b)=raw(2:4,b)+mass*vp(ipart,1:3)
                raw(13:15,b)=raw(13:15,b)+abs(mass*vp(ipart,1:3))
                raw(5,b)=raw(5,b)+mass*vp(ipart,1)**2
                raw(6,b)=raw(6,b)+mass*vp(ipart,2)**2
                raw(7,b)=raw(7,b)+mass*vp(ipart,3)**2
                raw(8,b)=raw(8,b)+mass*vp(ipart,1)*vp(ipart,2)
                raw(9,b)=raw(9,b)+mass*vp(ipart,1)*vp(ipart,3)
                raw(10,b)=raw(10,b)+mass*vp(ipart,2)*vp(ipart,3)
                raw(11,b)=raw(11,b)+mass*sum(vp(ipart,1:3)**2)*vp(ipart,1)
                kinetic=0.5d0*mass*sum(vp(ipart,1:3)**2)
                raw(12,b)=raw(12,b)+kinetic
                mass_min_local=min(mass_min_local,mass)
                mass_max_local=max(mass_max_local,mass)
                phase=fourier_k*xp(ipart,1)
                mode_local(1)=mode_local(1)+kinetic*cos(phase)/boxlen**3
                mode_local(2)=mode_local(2)-kinetic*sin(phase)/boxlen**3
                count_local(b)=count_local(b)+1_8
             end if
             ipart=nextp(ipart)
          end do
          if(ipart/=0) nerror=nerror+1
          igrid=next(igrid)
       end do
    end if
    moment_global=raw
    nparticle=count_local
    mode_global=mode_local
    mass_min_global=mass_min_local
    mass_max_global=mass_max_local
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,moment_global,15*nslab, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce raw slab moments')
    call MPI_ALLREDUCE(MPI_IN_PLACE,nparticle,nslab,MPI_INTEGER8, &
         & MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce slab particle counts')
    call MPI_ALLREDUCE(MPI_IN_PLACE,mode_global,2,MPI_DOUBLE_PRECISION, &
         & MPI_SUM,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce actual-position Fourier mode')
    call MPI_ALLREDUCE(MPI_IN_PLACE,mass_min_global,1,MPI_DOUBLE_PRECISION, &
         & MPI_MIN,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce minimum particle mass')
    call MPI_ALLREDUCE(MPI_IN_PLACE,mass_max_global,1,MPI_DOUBLE_PRECISION, &
         & MPI_MAX,MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce maximum particle mass')
#endif
    if(.not.all(ieee_is_finite(moment_global))) then
       nerror=nerror+1
       moment_global=0.0d0
    end if
    energy=0.5d0*(moment_global(5,:)+moment_global(6,:)+ &
         & moment_global(7,:))
    mode=mode_global
    momentum=sum(moment_global(2:4,:),dim=2)
    momentum_abs=sum(moment_global(13:15,:),dim=2)
    momentum_current=momentum
    momentum_abs_current=momentum_abs
    current_total_mass=sum(moment_global(1,:))
    current_particle_mass_min=mass_min_global
    current_particle_mass_max=mass_max_global
    if(sum(nparticle)<=0_8) nerror=nerror+1
    call sidm_kinetic_diag_cache_moments(moment_global(1:12,:))
  end subroutine sidm_kinetic_diag_census

  subroutine sidm_kinetic_diag_cache_moments(moment_global)
    real(dp), intent(in) :: moment_global(12,nslab)
    integer :: i
    real(dp) :: mass,u(3),centered(6),qx,volume

    if(myid/=1) return
    volume=boxlen**3/dble(nslab)
    do i=1,nslab
       mass=moment_global(1,i)
       raw_cache_mass(i)=mass
       call sidm_kinetic_flux_decomposition(moment_global(:,i),volume,u, &
            & centered,qx,cache_flux_advective(i), &
            & cache_flux_pressure_work(i),cache_flux_total(i))
       cache_velocity(:,i)=u
       cache_pressure(:,i)=centered/volume
       cache_heat_flux(i)=qx
       if(.not.all(ieee_is_finite(u)) .or. &
            & .not.all(ieee_is_finite(centered)) .or. &
            & .not.ieee_is_finite(qx) .or. &
            & .not.ieee_is_finite(cache_flux_advective(i)) .or. &
            & .not.ieee_is_finite(cache_flux_pressure_work(i)) .or. &
            & .not.ieee_is_finite(cache_flux_total(i))) &
            & nerror=nerror+1
    end do
  end subroutine sidm_kinetic_diag_cache_moments

  subroutine sidm_kinetic_diag_write_rows(target,step,time_now,energy, &
       & nparticle,previous,source,flux,source_cum,flux_cum, &
       & nevent_interval,nevent_total)
    integer, intent(in) :: target,step
    real(dp), intent(in) :: time_now,energy(nslab),previous(nslab)
    real(dp), intent(in) :: source(nslab),flux(nslab)
    real(dp), intent(in) :: source_cum(nslab),flux_cum(nslab)
    integer(kind=8), intent(in) :: nparticle(nslab)
    integer(kind=8), intent(in) :: nevent_interval,nevent_total
    integer :: i,left,right,ios
    real(dp) :: volume,delta,residual,decomp_residual

    if(myid/=1) return
    volume=boxlen**3/dble(nslab)
    do i=1,nslab
       left=modulo(i-1,nslab)+1
       right=modulo(i,nslab)+1
       delta=energy(i)-previous(i)
       residual=delta-source(i)-flux(left)+flux(right)
       decomp_residual=cache_flux_total(i)-cache_flux_advective(i)- &
            & cache_flux_pressure_work(i)-cache_heat_flux(i)
       if(abs(decomp_residual)>closure_rtol*max( &
            & abs(cache_flux_total(i)),abs(cache_flux_advective(i)), &
            & abs(cache_flux_pressure_work(i)),abs(cache_heat_flux(i)), &
            & tiny(1.0d0))) &
            & call sidm_kinetic_diag_abort('kinetic flux decomposition failed')
       write(output_unit,'(I0,",",I0,",",ES24.16,",",ES24.16,",", &
            & 4(I0,","),23(ES24.16,:,","))',iostat=ios) &
            & target,step,time_now,time_now/tau0_code,i,nparticle(i), &
            & nevent_interval,nevent_total, &
            & raw_cache_mass(i),cache_velocity(:,i),cache_pressure(:,i), &
            & cache_flux_advective(i),cache_flux_pressure_work(i), &
            & cache_heat_flux(i),cache_flux_total(i),decomp_residual,energy(i), &
            & delta,source(i),flux(left),flux(right), &
            & residual,source_cum(i),flux_cum(right)
       if(ios/=0) call sidm_kinetic_diag_abort('write output row')
    end do
    flush(output_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('flush slab output')
  end subroutine sidm_kinetic_diag_write_rows

  subroutine sidm_kinetic_diag_write_fourier(target,step,time_now,mode0, &
       & mode1,source,flux,source_cum,flux_cum,delta_energy,source_k0, &
       & source_abs,flux_abs,collision_energy_scale,delta_momentum,momentum_abs, &
       & event_energy_residual,event_momentum_residual, &
       & closure_scale,closure_tolerance)
    integer, intent(in) :: target,step
    real(dp), intent(in) :: time_now,mode0(2),mode1(2)
    real(dp), intent(in) :: source(2),flux(2),source_cum(2),flux_cum(2)
    real(dp), intent(in) :: delta_energy,source_k0,source_abs,flux_abs
    real(dp), intent(in) :: collision_energy_scale
    real(dp), intent(in) :: delta_momentum(3),momentum_abs(3)
    real(dp), intent(in) :: event_energy_residual,event_momentum_residual
    real(dp), intent(in) :: closure_scale,closure_tolerance
    real(dp) :: residual(2)
    integer :: ios

    if(myid/=1) return
    residual(1)=mode1(1)-mode0(1)-fourier_k*flux(2)-source(1)
    residual(2)=mode1(2)-mode0(2)+fourier_k*flux(1)-source(2)
    write(fourier_unit,'(2(I0,","),30(ES24.16,:,","))',iostat=ios) &
         & target,step,time_now,mode0,mode1,source,flux,residual, &
         & source_cum,flux_cum,delta_energy,source_k0,source_abs,flux_abs, &
         & collision_energy_scale, &
         & delta_momentum,momentum_abs, &
         & event_energy_residual,event_momentum_residual,closure_scale, &
         & closure_tolerance
    if(ios/=0) call sidm_kinetic_diag_abort('write Fourier row')
    flush(fourier_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('flush Fourier output')
  end subroutine sidm_kinetic_diag_write_fourier

  integer function sidm_kinetic_diag_bin(x)
    real(dp), intent(in) :: x
    real(dp) :: wrapped,dx
    integer(kind=8) :: slab
    wrapped=modulo(x,boxlen)
    dx=boxlen/dble(nslab)
    slab=floor(wrapped/dx,kind=8)
    slab=max(0_8,min(int(nslab-1,kind=8),slab))
    sidm_kinetic_diag_bin=1+int(slab)
  end function sidm_kinetic_diag_bin

  subroutine sidm_kinetic_diag_export_state(target,step,time_now)
    integer, intent(in) :: target,step
    real(dp), intent(in) :: time_now
    integer :: state_unit,meta_unit,sidecar_unit,ios,info,igrid,jpart,ipart
    integer(kind=8) :: count_local,count_global
    integer(int32) :: endian_word
    integer(int8) :: endian_bytes(4)
    real(dp) :: interval_exposure
    character(len=16) :: byte_order
    character(len=80) :: filename
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    if(.not.state_export_enabled) return
    write(filename,'("sidm_kinetic_state_t",I2.2,"_rank",I4.4,".bin")') &
         & target,myid
    open(newunit=state_unit,file=trim(filename),status='new', &
         & action='write',access='stream',form='unformatted',iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('create ID-keyed state stream')
    count_local=0_8
    if(numbl(myid,levelmin)>0) then
       igrid=headl(myid,levelmin)
       do jpart=1,numbl(myid,levelmin)
          if(igrid<=0 .or. igrid>ngridmax) &
               & call sidm_kinetic_diag_abort('invalid grid in state export')
          ipart=headp(igrid)
          do info=1,numbp(igrid)
             if(ipart<=0 .or. ipart>npartmax) &
                  & call sidm_kinetic_diag_abort('invalid particle in state export')
             if(idp(ipart)>0 .and. ptypep(ipart)==PTYPE_DM) then
                if(.not.ieee_is_finite(mp(ipart)) .or. mp(ipart)<=0.0d0 .or. &
                     & .not.all(ieee_is_finite(xp(ipart,1:3))) .or. &
                     & .not.all(ieee_is_finite(vp(ipart,1:3)))) &
                     & call sidm_kinetic_diag_abort('non-finite ID-keyed state')
                write(state_unit,iostat=ios) int(idp(ipart),kind=8),mp(ipart), &
                     & xp(ipart,1:3),vp(ipart,1:3)
                if(ios/=0) call sidm_kinetic_diag_abort('write ID-keyed state')
                count_local=count_local+1_8
             end if
             ipart=nextp(ipart)
          end do
          if(ipart/=0) call sidm_kinetic_diag_abort('state particle-list tail mismatch')
          igrid=next(igrid)
       end do
    end if
    flush(state_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('flush ID-keyed state')
    close(state_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('close ID-keyed state')
    count_global=count_local
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,count_global,1,MPI_INTEGER8,MPI_SUM, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce ID-keyed state count')
#endif
    if(count_global/=expected_particle_count) &
         & call sidm_kinetic_diag_abort('ID-keyed state count mismatch')
    endian_word=1_int32
    endian_bytes=transfer(endian_word,endian_bytes)
    if(endian_bytes(1)==1_int8) then
       byte_order='little-endian'
    else
       byte_order='big-endian'
    end if
    interval_exposure=0.0d0
    if(target>0) interval_exposure=time_now- &
         & target_fraction(target)*tau0_code
    write(filename,'("sidm_kinetic_state_t",I2.2,"_rank",I4.4,".meta")') &
         & target,myid
    open(newunit=sidecar_unit,file=trim(filename),status='new', &
         & action='write',iostat=ios)
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) 'format_version=',1
    if(ios==0) write(sidecar_unit,'(A,A)',iostat=ios) &
         & 'byte_order=',trim(byte_order)
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'rank=',myid
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'target=',target
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'step=',step
    if(ios==0) write(sidecar_unit,'(A,ES24.16)',iostat=ios) &
         & 'time=',time_now
    if(ios==0) write(sidecar_unit,'(A,ES24.16)',iostat=ios) &
         & 'tau0_fraction=',time_now/tau0_code
    if(ios==0) write(sidecar_unit,'(A,ES24.16)',iostat=ios) &
         & 'interval_exposure=',interval_exposure
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'local_count=',count_local
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'id_bytes=',storage_size(0_8)/8
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'real_bytes=',storage_size(1.0d0)/8
    if(ios==0) write(sidecar_unit,'(A,I0)',iostat=ios) &
         & 'record_bytes=',(storage_size(0_8)+7*storage_size(1.0d0))/8
    if(ios==0) write(sidecar_unit,'(A)',iostat=ios) &
         & 'fields=int64_id,float64_mass,float64_x,float64_y,float64_z,'// &
         & 'float64_vx,float64_vy,float64_vz'
    if(ios==0) flush(sidecar_unit,iostat=ios)
    if(ios==0) close(sidecar_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('write state sidecar')
    if(myid==1) then
       if(target==0) then
          open(newunit=meta_unit,file='sidm_kinetic_state_targets.csv', &
               & status='new',action='write',iostat=ios)
          if(ios==0) write(meta_unit,'(A)',iostat=ios) &
               & 'target,step,time,tau0_fraction,N'
       else
          open(newunit=meta_unit,file='sidm_kinetic_state_targets.csv', &
               & status='old',position='append',action='write',iostat=ios)
       end if
       if(ios==0) write(meta_unit,'(I0,",",I0,",",ES24.16,",", &
            & ES24.16,",",I0)',iostat=ios) target,step,time_now, &
            & time_now/tau0_code,count_global
       if(ios==0) flush(meta_unit,iostat=ios)
       if(ios==0) close(meta_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('write state-target metadata')
    end if
  end subroutine sidm_kinetic_diag_export_state

  subroutine sidm_kinetic_diag_write_control_marker(step,time_now)
    integer, intent(in) :: step
    real(dp), intent(in) :: time_now
    integer :: marker_unit,ios

    if(.not.control_only .or. .not.state_export_enabled .or. &
         & target_index<=ntarget) &
         & call sidm_kinetic_diag_abort('invalid control-only completion state')
    if(abs(time_now-target_fraction(ntarget)*tau0_code)> &
         & 1.0d-11*max(1.0d0,abs(time_now))) &
         & call sidm_kinetic_diag_abort('control-only terminal time mismatch')
    if(myid/=1) return
    open(newunit=marker_unit,file='sidm_kinetic_control.observer_validated', &
         & status='new',action='write',iostat=ios)
    if(ios==0) write(marker_unit,'(A)',iostat=ios) &
         & 'status=CONTROL_COMPLETE;measurement=DISABLED;target_control=ACTIVE'
    if(ios==0) write(marker_unit,'(A,ES24.16,A,I0)',iostat=ios) &
         & 'terminal_time=',time_now,',terminal_step=',step
    if(ios==0) flush(marker_unit,iostat=ios)
    if(ios==0) close(marker_unit,iostat=ios)
    if(ios/=0) call sidm_kinetic_diag_abort('write control-only completion marker')
  end subroutine sidm_kinetic_diag_write_control_marker

  subroutine sidm_kinetic_diag_init()
    character(len=64) :: value
    integer :: status,length,ios,flags,flags_min,flags_max,info
    integer :: level_min,level_max
    integer(kind=8) :: count_min,count_max
    real(dp) :: tau0_min,tau0_max
    real(dp) :: dt_min,dt_max
    integer, external :: omp_get_max_threads
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    if(initialized) return
    initialized=.true.
    value=' '
    call get_environment_variable('SIDM_KINETIC_TRANSPORT',value, &
         & length=length,status=status)
    if(status==0 .and. trim(value)=='1') then
       enabled=.true.
    else if(status==1 .or. (status==0 .and. trim(value)=='0')) then
       enabled=.false.
    else
       call sidm_kinetic_diag_abort('invalid enable variable')
    end if
    value=' '
    call get_environment_variable('SIDM_KINETIC_CONTROL_ONLY',value, &
         & length=length,status=status)
    if(status==0 .and. trim(value)=='1') then
       control_only=.true.
    else if(status==1 .or. (status==0 .and. trim(value)=='0')) then
       control_only=.false.
    else
       call sidm_kinetic_diag_abort('invalid control-only variable')
    end if
    if(enabled .and. control_only) &
         & call sidm_kinetic_diag_abort('measurement and control-only modes conflict')
    target_control_enabled=enabled .or. control_only
    if(target_control_enabled) then
       value=' '
       call get_environment_variable('SIDM_KINETIC_TAU0_CODE',value, &
            & length=length,status=status)
       if(status/=0 .or. length==0) &
            & call sidm_kinetic_diag_abort('missing tau0 in code units')
       read(value,*,iostat=ios) tau0_code
       if(ios/=0 .or. .not.ieee_is_finite(tau0_code) .or. tau0_code<=0.0d0) &
            & call sidm_kinetic_diag_abort('invalid tau0 in code units')
       value=' '
       call get_environment_variable('SIDM_KINETIC_EXPECTED_LEVEL',value, &
            & length=length,status=status)
       if(status/=0 .or. length==0) &
            & call sidm_kinetic_diag_abort('missing expected uniform level')
       read(value,*,iostat=ios) configured_level
       if(ios/=0 .or. configured_level<1) &
            & call sidm_kinetic_diag_abort('invalid expected uniform level')
       value=' '
       call get_environment_variable('SIDM_KINETIC_EXPECTED_PARTICLES',value, &
            & length=length,status=status)
       if(status/=0 .or. length==0) &
            & call sidm_kinetic_diag_abort('missing expected particle count')
       read(value,*,iostat=ios) configured_particle_count
       if(ios/=0 .or. configured_particle_count<2_8) &
            & call sidm_kinetic_diag_abort('invalid expected particle count')
       value=' '
       call get_environment_variable('SIDM_KINETIC_DT_CODE',value, &
            & length=length,status=status)
       if(status/=0 .or. length==0) &
            & call sidm_kinetic_diag_abort('missing target timestep in code units')
       read(value,*,iostat=ios) dt_target_code
       if(ios/=0 .or. .not.ieee_is_finite(dt_target_code) .or. &
            & dt_target_code<=0.0d0) &
            & call sidm_kinetic_diag_abort('invalid target timestep in code units')
       if(ncpu/=2 .or. ndim/=3 .or. levelmin/=configured_level .or. &
            & nlevelmax/=configured_level .or. nrestart/=0 .or. static .or. tracer .or. &
            & omp_get_max_threads()/=1 .or. &
            & abs(boxlen-8.0d0*acos(-1.0d0))>1.0d-12) &
            & call sidm_kinetic_diag_abort( &
            & 'benchmark requires fresh 2-rank uniform 3D mesh and one OpenMP thread')
       if(.not.pic .or. .not.poisson .or. hydro .or. cosmo .or. &
            & star .or. sink .or. gpu_particle .or. use_adm .or. &
            & adm_hpm .or. use_fdm .or. use_mond .or. &
            & use_coupled_de .or. sidm_baryon .or. sidm_inelastic .or. &
            & sidm_fdiss/=0.0d0 .or. sidm_type/='constant' .or. &
            & sidm_a_type/='none' .or. sidm_angular/='isotropic' .or. &
            & gravity_type/=1 .or. any(gravity_params(1:3)/=0.0d0)) &
            & call sidm_kinetic_diag_abort('invalid force-free pilot physics')
       value=' '
       call get_environment_variable('SIDM_KINETIC_STATE_EXPORT',value, &
            & length=length,status=status)
       if(status==0 .and. trim(value)=='1') then
          state_export_enabled=.true.
       else if(status==1 .or. (status==0 .and. trim(value)=='0')) then
          state_export_enabled=.false.
       else
          call sidm_kinetic_diag_abort('invalid state-export variable')
       end if
       if(control_only .and. .not.state_export_enabled) &
            & call sidm_kinetic_diag_abort('control-only mode requires state export')
    end if
    flags=0
    if(enabled) flags=flags+1
    if(control_only) flags=flags+2
    if(state_export_enabled) flags=flags+4
    flags_min=flags
    flags_max=flags
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(flags,flags_min,1,MPI_INTEGER,MPI_MIN, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce enable flag minimum')
    call MPI_ALLREDUCE(flags,flags_max,1,MPI_INTEGER,MPI_MAX, &
         & MPI_COMM_WORLD,info)
    call sidm_kinetic_diag_mpierr(info,'reduce enable flag maximum')
    if(flags_min/=flags_max) &
         & call sidm_kinetic_diag_abort('rank configuration mismatch')
    if(target_control_enabled) then
       tau0_min=tau0_code
       tau0_max=tau0_code
       call MPI_ALLREDUCE(MPI_IN_PLACE,tau0_min,1,MPI_DOUBLE_PRECISION, &
            & MPI_MIN,MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce tau0 minimum')
       call MPI_ALLREDUCE(MPI_IN_PLACE,tau0_max,1,MPI_DOUBLE_PRECISION, &
            & MPI_MAX,MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce tau0 maximum')
       if(tau0_min/=tau0_max) &
            & call sidm_kinetic_diag_abort('rank tau0 mismatch')
       level_min=configured_level
       level_max=configured_level
       call MPI_ALLREDUCE(MPI_IN_PLACE,level_min,1,MPI_INTEGER,MPI_MIN, &
            & MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce expected level minimum')
       call MPI_ALLREDUCE(MPI_IN_PLACE,level_max,1,MPI_INTEGER,MPI_MAX, &
            & MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce expected level maximum')
       if(level_min/=level_max) &
            & call sidm_kinetic_diag_abort('rank expected-level mismatch')
       count_min=configured_particle_count
       count_max=configured_particle_count
       call MPI_ALLREDUCE(MPI_IN_PLACE,count_min,1,MPI_INTEGER8,MPI_MIN, &
            & MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce expected count minimum')
       call MPI_ALLREDUCE(MPI_IN_PLACE,count_max,1,MPI_INTEGER8,MPI_MAX, &
            & MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce expected count maximum')
       if(count_min/=count_max) &
            & call sidm_kinetic_diag_abort('rank expected-count mismatch')
       dt_min=dt_target_code
       dt_max=dt_target_code
       call MPI_ALLREDUCE(MPI_IN_PLACE,dt_min,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
            & MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce target timestep minimum')
       call MPI_ALLREDUCE(MPI_IN_PLACE,dt_max,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
            & MPI_COMM_WORLD,info)
       call sidm_kinetic_diag_mpierr(info,'reduce target timestep maximum')
       if(dt_min/=dt_max) &
            & call sidm_kinetic_diag_abort('rank target-timestep mismatch')
    end if
#endif
  end subroutine sidm_kinetic_diag_init

  subroutine sidm_kinetic_diag_abort(context)
    character(len=*), intent(in) :: context
    integer :: info
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    if(myid==1) write(*,'(A,A)') 'ERROR: SIDM kinetic observer: ',trim(context)
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop 1
#endif
  end subroutine sidm_kinetic_diag_abort

  subroutine sidm_kinetic_diag_mpierr(info,context)
    integer, intent(in) :: info
    character(len=*), intent(in) :: context
    if(info/=0) call sidm_kinetic_diag_abort(context)
  end subroutine sidm_kinetic_diag_mpierr

  subroutine sidm_kinetic_diag_orderly_stop()
    integer :: info,ios
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    if(output_unit/=-1) then
       flush(output_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('flush final transport output')
       close(output_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('close final transport output')
       output_unit=-1
    end if
    if(fourier_unit/=-1) then
       flush(fourier_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('flush final Fourier output')
       close(fourier_unit,iostat=ios)
       if(ios/=0) call sidm_kinetic_diag_abort('close final Fourier output')
       fourier_unit=-1
    end if
#ifndef WITHOUTMPI
    call MPI_FINALIZE(info)
    if(info/=0) then
       if(myid==1) write(*,'(A,I0)') &
            & 'ERROR: SIDM kinetic observer: MPI_FINALIZE returned ',info
       stop 1
    end if
#endif
    stop
  end subroutine sidm_kinetic_diag_orderly_stop

end module sidm_kinetic_transport_diag
