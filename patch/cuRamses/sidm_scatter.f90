!################################################################
!################################################################
!################################################################
! SIDM (Self-Interacting Dark Matter) Monte Carlo Scattering
!
! Features:
! - Velocity-dependent cross-section (constant/yukawa/power_law)
! - Isotropic or anisotropic (Rutherford-like) scattering angle
! - Inelastic scattering with mass splitting (iSIDM)
! - Timestep constraint via P_max tracking
!################################################################
subroutine sidm_scatter(ilevel)
  use pm_commons
  use pm_parameters, only: iseed
  use amr_parameters, only: sidm_estimator_diagnostics, &
       & sidm_probability_guard,sidm_rate_bound
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_commons
  use random
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  integer::icpu,info,sidm_n_substeps
  integer(int64)::sidm_particles_local,sidm_particles_global
  integer,parameter::sidm_max_subcycles=4096
  real(dp)::sidm_dt_interval,sidm_dt_remaining,sidm_dt_sub,sidm_dt_elapsed
  real(dp)::sidm_dt_previous,sidm_rate_current,sidm_p_remaining,sidm_p_substep
  real(dp)::sidm_scale_l,sidm_scale_t,sidm_scale_d,sidm_scale_v,sidm_scale_nH,sidm_scale_T2
  integer::sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  integer::sidm_n_vrel_reject
  integer::sidm_n_leaf_occupied,sidm_n_leaf_active,sidm_n_dm_leaf
  real(dp)::sidm_dp(3),sidm_dEk,sidm_dp_max,sidm_dEk_max
  real(dp)::sidm_Pmax_level
  real(dp)::sidm_dEdiss  ! Cumulative dissipated energy (dSIDM)
  real(dp)::sidm_expected_rate
  real(dp)::R_dummy
  integer,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer,dimension(IRandNumSize,0:0)::sidm_seeds
  logical,save::sidm_excited_init = .false.
  common /sidm_diag/ sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down,sidm_n_vrel_reject
  common /sidm_cons/ sidm_dp,sidm_dEk,sidm_dp_max,sidm_dEk_max
  common /sidm_pmax_c/ sidm_Pmax_level
  common /sidm_diss/ sidm_dEdiss
  common /sidm_estimator_counts/ sidm_n_leaf_occupied, &
       & sidm_n_leaf_active,sidm_n_dm_leaf
  common /sidm_estimator_rate/ sidm_expected_rate

  if(.not.sidm) return
  sidm_dt_interval=dtnew(ilevel)
  sidm_dt_elapsed=0.0d0
  sidm_n_substeps=0
  if(sidm_probability_guard) then
     ! All ranks participate even if a rank has no local grids.
     sidm_rate_bound(ilevel)=0.0d0
     call sidm_preflight_rate_gate(ilevel,sidm_dt_interval,sidm_dt_sub, &
          & sidm_rate_current,sidm_p_remaining,sidm_p_substep)
     if(npart<0) call sidm_probability_guard_abort('negative local particle count')
     sidm_particles_local=int(npart,kind=int64)
     call MPI_ALLREDUCE(sidm_particles_local,sidm_particles_global, &
          & 1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     if(info/=0 .or. sidm_particles_global<0_int64) &
          call sidm_probability_guard_abort('invalid global particle count')
     if(sidm_particles_global==0_int64) return
     call units(sidm_scale_l,sidm_scale_t,sidm_scale_d,sidm_scale_v, &
          & sidm_scale_nH,sidm_scale_T2)
  else
     if(numbtot(1,ilevel)==0) return
  end if
  if(verbose) write(*,111) ilevel

  ! Initialize excited DM fraction (once, for iSIDM)
  if(sidm_inelastic .and. .not.sidm_excited_init) then
     call sidm_init_excited()
     sidm_excited_init = .true.
  end if

  ! Initialize localseed if first call
  if(localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  ! Only thread 0's RNG state is retained between calls.  Reconstructing
  ! other thread streams from it repeats draws when their work is unequal.
  ! Until every worker stream is persistent across calls and restarts, run
  ! the local SIDM scatter serially, independently of other OpenMP regions.
  sidm_seeds(:,0) = localseed

  ! Reset level-wide counters
  sidm_n_scatter  = 0
  sidm_n_pairs    = 0
  sidm_n_up       = 0
  sidm_n_down     = 0
  sidm_n_vrel_reject = 0
  sidm_n_leaf_occupied = 0
  sidm_n_leaf_active = 0
  sidm_n_dm_leaf = 0
  sidm_dp(:)      = 0.0d0
  sidm_dEk        = 0.0d0
  sidm_dp_max     = 0.0d0
  sidm_dEk_max    = 0.0d0
  sidm_Pmax_level = 0.0d0
  sidm_dEdiss     = 0.0d0
  sidm_expected_rate = 0.0d0

#if NDIM==3
  ! Only local grids: after virtual_tree_fine all particles sit on
  ! grids owned by myid (virtual grids carry no particles here)
  icpu = myid
  do
     if(sidm_probability_guard) then
        sidm_dt_remaining=sidm_dt_interval-sidm_dt_elapsed
        if(sidm_n_substeps>0) then
           if(sidm_dt_remaining<=0.0d0) exit
           call sidm_preflight_rate_gate(ilevel,sidm_dt_remaining,sidm_dt_sub, &
                & sidm_rate_current,sidm_p_remaining,sidm_p_substep)
        end if
        if(sidm_n_substeps>=sidm_max_subcycles) &
             call sidm_probability_guard_abort('collision subcycle limit exceeded')
        sidm_n_substeps=sidm_n_substeps+1
        if(sidm_particles_global>int(huge(sidm_n_dm_leaf),kind=int64)/ &
             int(sidm_n_substeps,kind=int64)) &
             call sidm_probability_guard_abort('diagnostic visit counter may overflow')
        if(myid==1) write(*,'(A,I0,A,I0,A,I0,A,5ES25.16E3)') &
             ' SIDM_SUBCYCLE_PREFLIGHT: coarse=',nstep_coarse,' level=',ilevel, &
             ' index=',sidm_n_substeps, &
             ' rate_s,remaining_s,dt_s,P_remaining,P_substep=', &
             sidm_rate_current,sidm_dt_remaining*sidm_scale_t, &
             sidm_dt_sub*sidm_scale_t,sidm_p_remaining,sidm_p_substep
     else
        sidm_n_substeps=1
        sidm_dt_remaining=sidm_dt_interval
        sidm_dt_sub=sidm_dt_interval
     end if
     if(numbl(icpu,ilevel)>0) then
        call sub_sidm_scatter(ilevel,icpu,headl(icpu,ilevel),numbl(icpu,ilevel), &
             & sidm_seeds,0,1,sidm_dt_sub,sidm_dt_interval)
     end if
     if(.not.sidm_probability_guard) exit
     if(sidm_dt_interval==0.0d0) exit
     sidm_dt_previous=sidm_dt_elapsed
     if(sidm_dt_sub>=sidm_dt_remaining) then
        sidm_dt_elapsed=sidm_dt_interval
     else
        sidm_dt_elapsed=sidm_dt_elapsed+sidm_dt_sub
     end if
     if(.not.(sidm_dt_elapsed>sidm_dt_previous) .or. &
          sidm_dt_elapsed>sidm_dt_interval) &
          call sidm_probability_guard_abort('collision subcycle made invalid progress')
     if(sidm_dt_elapsed>=sidm_dt_interval) exit
  end do
  if(sidm_probability_guard .and. myid==1) then
     write(*,'(A,I0,A,I0,A,I0,A,2ES25.16E3)') &
          ' SIDM_SUBCYCLE_SUMMARY: coarse=',nstep_coarse,' level=',ilevel, &
          ' nsub=',sidm_n_substeps,' interval_s,covered_s=', &
          sidm_dt_interval*sidm_scale_t,sidm_dt_elapsed*sidm_scale_t
  end if
#endif

  ! Continue the rank stream from thread 0's final state (advance once
  ! so the seed always moves even if no pairs were processed)
  localseed = sidm_seeds(:,0)
  call ranf(localseed, R_dummy)

  ! MPI reduce all diagnostics
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_scatter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_pairs,  1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_up,     1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_down,   1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_vrel_reject,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_leaf_occupied,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_leaf_active,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_dm_leaf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dp,   3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dEk,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dp_max, 1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dEk_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_Pmax_level,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dEdiss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_expected_rate,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  if(sidm_probability_guard .and. myid==1) then
     write(*,'(A,I0,A,I0,A,ES25.16E3)') &
          ' SIDM_SUBCYCLE_ACTUAL: coarse=',nstep_coarse,' level=',ilevel, &
          ' Pmax=',sidm_Pmax_level
  end if

  ! Store P_max for timestep constraint in newdt_fine
  sidm_Pmax(ilevel) = sidm_Pmax_level

  ! Report statistics
  if(myid==1 .and. sidm_n_pairs>0) then
     write(*,'(A,I2,A,I8,A,I10,A,ES9.2)') &
          ' SIDM level ',ilevel,': scattered=',sidm_n_scatter, &
          ' pairs=',sidm_n_pairs,' Pmax=',sidm_Pmax_level
     if(sidm_n_scatter>0) then
        write(*,'(A,3ES12.4)') '   dp(x,y,z)=', sidm_dp(1:3)
        write(*,'(A,ES12.4,A,ES12.4)') '   dEk_total=', sidm_dEk, &
             '  dp_max=', sidm_dp_max
        write(*,'(A,ES12.4)') '   dEk_max  =', sidm_dEk_max
     end if
     if(sidm_inelastic .and. (sidm_n_up>0 .or. sidm_n_down>0)) then
        write(*,'(A,I6,A,I6,A,I6)') '   iSIDM: up=', &
             sidm_n_up,' down=',sidm_n_down, &
             ' net=',sidm_n_up-sidm_n_down
     end if
     if(sidm_fdiss > 0.0d0 .and. sidm_dEdiss /= 0.0d0) then
        write(*,'(A,ES12.4)') '   dEdiss_total=', sidm_dEdiss
     end if
     if(sidm_n_vrel_reject > 0) then
        write(*,'(A,I10,A,ES10.3,A)') &
             '   vrel_rejected=', sidm_n_vrel_reject, &
             ' (cap=', sidm_vrel_max, ' cm/s)'
     end if
  end if
  if(myid==1 .and. sidm_estimator_diagnostics) then
     if(sidm_probability_guard) then
        write(*,'(A,I0,A,I0)') ' SIDM estimator: counts are visits; substeps=',sidm_n_substeps, &
             ' level=',ilevel
     end if
     write(*,'(A,I2,A,I10,A,I10,A,I10,A,I10,A,ES12.4,A)') &
          ' SIDM estimator level ',ilevel,': dm=',sidm_n_dm_leaf, &
          ' occupied=',sidm_n_leaf_occupied, &
          ' active=',sidm_n_leaf_active,' sampled=',sidm_n_pairs, &
          ' rate=',sidm_expected_rate,' s^-1'
  end if

111 format('   Entering sidm_scatter for level ',I2)
end subroutine sidm_scatter
!################################################################
!################################################################
subroutine sidm_preflight_rate_gate(ilevel,dt_remaining_code,dt_sub_code, &
     & rate_out,p_remaining_out,p_substep_out)
  use pm_commons
  use amr_commons
  use amr_parameters, only: sidm_estimator_diagnostics,sidm_rate_bound, &
       & sidm_courant,sidm_cross_section,sidm_npart_min,sidm_type,sidm_a_type, &
       & sidm_inelastic,sidm_fdiss
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_remaining_code
  real(dp),intent(out)::dt_sub_code,rate_out,p_remaining_out,p_substep_out
  integer::igrid,jgrid,ind,icell,ipart,jp,npart1,ndm_cell,nx_loc
  integer::info,invalid_local,invalid_global,invalid_cell
  integer::invalid_reason_local,invalid_reason_global
  integer::gate_failed,gate_failed_global,sidm_dt_invalid
  integer::schedule_consistent
  integer::bad_ipart,bad_grid,bad_cell,bad_npart,bad_helper_status
  integer::bad_jp,bad_expected
  logical::geometry_valid
  integer,external::cell_index_from_part
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::dx,dx_loc,scale,vol_phys,dt_phys
  real(dp)::mass_max,mass_phys,rate_cell,rate_local,rate_global
  real(dp)::p_bound,p_sub_bound,p_limit,p_limit_global
  real(dp)::dt_code_min,dt_code_max,scale_t_min,scale_t_max
  real(dp)::courant_min,courant_max
  real(dp)::v_min(1:3),v_upper(1:3),v_phys(1:3)
  real(dp)::bad_mass,bad_rate,bad_velocity(1:3)
  real(dp)::bad_v_min(1:3),bad_v_upper(1:3)
  external sidm_probability_rate_bound
  external sidm_choose_collision_substep
  external sidm_collision_schedule_agrees

  dt_sub_code=0.0d0
  rate_out=0.0d0
  p_remaining_out=0.0d0
  p_substep_out=0.0d0

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  invalid_reason_local=0
  invalid_reason_global=0
  bad_ipart=0
  bad_grid=0
  bad_cell=0
  bad_npart=0
  bad_helper_status=0
  bad_jp=0
  bad_expected=0
  bad_mass=0.0d0
  bad_rate=0.0d0
  bad_velocity=0.0d0
  bad_v_min=0.0d0
  bad_v_upper=0.0d0
  geometry_valid=.false.
  invalid_local=0
  ! Check finiteness before ordered comparisons: Fortran does not promise
  ! short-circuit evaluation, and comparisons with NaN may trap.
  if(any(.not.ieee_is_finite((/scale_t,scale_l,scale_d,scale_v, &
       sidm_courant,sidm_cross_section,sidm_fdiss,dt_remaining_code,boxlen/)))) then
     invalid_local=1
     invalid_reason_local=1 ! nonfinite unit, parameter, or interval
  end if
  if(invalid_local==0) then
     if(scale_t<=0.0d0 .or. scale_l<=0.0d0 .or. scale_d<=0.0d0 .or. &
          scale_v<=0.0d0 .or. boxlen<=0.0d0 .or. dt_remaining_code<0.0d0 .or. &
          sidm_cross_section<=0.0d0 .or. sidm_courant<=0.0d0 .or. &
          sidm_courant>1.0d0 .or. sidm_fdiss/=0.0d0) then
        invalid_local=1
        invalid_reason_local=2 ! invalid unit, parameter, or interval
     end if
  end if
  if(sidm_npart_min<2 .or. trim(sidm_type)/='constant' .or. &
       trim(sidm_a_type)/='none' .or. sidm_inelastic) then
     invalid_local=1
     if(invalid_reason_local==0) invalid_reason_local=2 ! unsupported mode/parameters
  end if
#if NDIM != 3
  invalid_local=1
  if(invalid_reason_local==0) invalid_reason_local=3 ! unsupported dimension
#endif
  dt_phys=0.0d0
  if(invalid_local==0) then
     dt_phys=dt_remaining_code*scale_t
     if(.not.ieee_is_finite(dt_phys)) then
        invalid_local=1
        invalid_reason_local=4 ! invalid physical timestep
     end if
  end if

  rate_local=0.0d0
  dx=0.5d0**ilevel
  nx_loc=icoarse_max-icoarse_min+1
  if(nx_loc<=0 .or. invalid_local/=0) then
     invalid_local=1
     if(invalid_reason_local==0) invalid_reason_local=5 ! invalid grid geometry
     vol_phys=0.0d0
  else
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_phys=(dx_loc*scale_l)**3
     if(.not.ieee_is_finite(vol_phys) .or. vol_phys<=0.0d0) then
        invalid_local=1
        if(invalid_reason_local==0) invalid_reason_local=6 ! invalid cell volume
     else
        geometry_valid=.true.
     end if
  end if

#if NDIM == 3
  ! The particle classifier divides by boxlen/nx_loc.  Do not call it when
  ! the geometry preflight above has already found invalid inputs.
  if(geometry_valid) then
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     if(igrid<lbound(numbp,1) .or. igrid>ubound(numbp,1) .or. &
          igrid<lbound(next,1) .or. igrid>ubound(next,1) .or. &
          igrid<lbound(headp,1) .or. igrid>ubound(headp,1) .or. &
          igrid<lbound(xg,1) .or. igrid>ubound(xg,1)) then
        invalid_local=1
        if(invalid_reason_local==0) invalid_reason_local=12 ! bad grid index
        exit
     end if
     if(any(.not.ieee_is_finite(xg(igrid,1:3)))) then
        invalid_local=1
        if(invalid_reason_local==0) invalid_reason_local=13 ! bad grid position
        exit
     end if
     npart1=numbp(igrid)
     if(npart1<0) then
        invalid_local=1
        if(invalid_reason_local==0) invalid_reason_local=14 ! bad occupancy
        exit
     end if
     if(npart1==0) then
        igrid=next(igrid)
        cycle
     end if
     do ind=1,twotondim
        icell=ICELL_OF(igrid,ind)
        if(icell<lbound(son,1) .or. icell>ubound(son,1)) then
           invalid_local=1
           if(invalid_reason_local==0) invalid_reason_local=15 ! bad leaf index
           exit
        end if
        if(son(icell)/=0) cycle
        ndm_cell=0
        mass_max=0.0d0
        v_min=huge(1.0d0)
        v_upper=-huge(1.0d0)
        ipart=headp(igrid)
        jp=0
        ! The linked-list zero sentinel must never be dereferenced, even if
        ! numbp and the chain length disagree.
        do while(jp<npart1 .and. ipart>0)
           jp=jp+1
           if(ipart<lbound(idp,1) .or. ipart>ubound(idp,1) .or. &
                ipart<lbound(ptypep,1) .or. ipart>ubound(ptypep,1) .or. &
                ipart<lbound(nextp,1) .or. ipart>ubound(nextp,1) .or. &
                ipart<lbound(xp,1) .or. ipart>ubound(xp,1) .or. &
                ipart<lbound(mp,1) .or. ipart>ubound(mp,1) .or. &
                ipart<lbound(vp,1) .or. ipart>ubound(vp,1)) then
              invalid_local=1
              if(invalid_reason_local==0) then
                 invalid_reason_local=11 ! particle index outside array bounds
                 bad_ipart=ipart
                 bad_grid=igrid
                 bad_cell=ind
                 bad_npart=ndm_cell
                 bad_jp=jp
                 bad_expected=npart1
              end if
              exit
           end if
           if(any(.not.ieee_is_finite(xp(ipart,1:3)))) then
              invalid_local=1
              if(invalid_reason_local==0) invalid_reason_local=16 ! bad particle position
              exit
           end if
           if(cell_index_from_part(ipart,igrid,ilevel)==ind) then
              if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
                   ptypep(ipart)/=PTYPE_SINK) then
                 ndm_cell=ndm_cell+1
                 mass_phys=mp(ipart)*scale_d*scale_l**3
                 v_phys(1:3)=vp(ipart,1:3)*scale_v
                 if(.not.ieee_is_finite(mass_phys) .or. any(.not.ieee_is_finite(v_phys))) then
                    invalid_local=1
                    if(invalid_reason_local==0) then
                       invalid_reason_local=7 ! invalid particle mass or velocity
                       bad_ipart=ipart
                       bad_grid=igrid
                       bad_cell=ind
                       bad_npart=ndm_cell
                       bad_jp=jp
                       bad_expected=npart1
                       bad_mass=mass_phys
                       bad_velocity=v_phys
                       bad_v_min=v_phys
                       bad_v_upper=v_phys
                    end if
                 else if(mass_phys<=0.0d0) then
                    invalid_local=1
                    if(invalid_reason_local==0) invalid_reason_local=7
                 else
                    mass_max=max(mass_max,mass_phys)
                    v_min=min(v_min,v_phys)
                    v_upper=max(v_upper,v_phys)
                 end if
              end if
           end if
           ipart=nextp(ipart)
        end do
        if(jp/=npart1 .or. ipart>0) then
           invalid_local=1
           if(invalid_reason_local==0) then
              invalid_reason_local=10 ! numbp disagrees with linked-list length
              bad_grid=igrid
              bad_cell=ind
              bad_ipart=ipart
              bad_npart=ndm_cell
              bad_jp=jp
              bad_expected=npart1
           end if
        end if
        if(ndm_cell>=sidm_npart_min .and. mass_max>0.0d0) then
           call sidm_probability_rate_bound(sidm_cross_section,mass_max, &
                & v_min,v_upper,vol_phys,ndm_cell,rate_cell,invalid_cell)
           if(invalid_cell/=0 .or. .not.ieee_is_finite(rate_cell)) then
              invalid_local=1
              if(invalid_reason_local==0) then
                 invalid_reason_local=8 ! rate-helper rejection
                 bad_grid=igrid
                 bad_cell=ind
                 bad_npart=ndm_cell
                 bad_jp=jp
                 bad_expected=npart1
                 bad_helper_status=invalid_cell
                 bad_mass=mass_max
                 bad_rate=rate_cell
                 bad_v_min=v_min
                 bad_v_upper=v_upper
                 bad_velocity=v_upper-v_min
              end if
           else
              rate_local=max(rate_local,rate_cell)
           end if
        end if
     end do
     igrid=next(igrid)
  end do
  end if
#endif

  call MPI_ALLREDUCE(rate_local,rate_global,1,MPI_DOUBLE_PRECISION, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(invalid_local,invalid_global,1,MPI_INTEGER, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  if(.not.ieee_is_finite(rate_global)) then
     invalid_global=1
     if(invalid_reason_local==0) invalid_reason_local=9 ! invalid reduced rate
  else if(rate_global<0.0d0) then
     invalid_global=1
     if(invalid_reason_local==0) invalid_reason_local=9
  end if
  if(invalid_local/=0 .or. invalid_global/=0) then
     write(*,*) 'SIDM_PREFLIGHT_INVALID rank/reason=',myid,invalid_reason_local, &
          ' scales(t,l,d,v)=',scale_t,scale_l,scale_d,scale_v, &
          ' dt,nx,boxlen,volume,rate,sigma=',dt_phys,nx_loc,boxlen, &
          vol_phys,rate_local,sidm_cross_section, &
          ' grid,cell_index,particle,npart=', &
          bad_grid,bad_cell,bad_ipart,bad_npart, &
          ' linked_nodes,expected=',bad_jp,bad_expected, &
          ' helper_status,rate_cell,mass_max=',bad_helper_status,bad_rate,bad_mass, &
          ' v_min,v_upper=',bad_v_min,bad_v_upper, &
          ' mode,courant,npart_min=',trim(sidm_type),trim(sidm_a_type), &
          sidm_inelastic,sidm_fdiss,sidm_courant,sidm_npart_min
  end if
  call MPI_ALLREDUCE(invalid_reason_local,invalid_reason_global,1,MPI_INTEGER, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  if(invalid_global/=0) then
     if(myid==1) write(*,'(A,I0)') &
          'ERROR: invalid input or unsupported mode in SIDM probability preflight; reason=', &
          invalid_reason_global
     call sidm_probability_guard_abort('invalid preflight input or unsupported mode')
  end if

  call MPI_ALLREDUCE(dt_remaining_code,dt_code_min,1,MPI_DOUBLE_PRECISION, &
       & MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_remaining_code,dt_code_max,1,MPI_DOUBLE_PRECISION, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(scale_t,scale_t_min,1,MPI_DOUBLE_PRECISION, &
       & MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(scale_t,scale_t_max,1,MPI_DOUBLE_PRECISION, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sidm_courant,courant_min,1,MPI_DOUBLE_PRECISION, &
       & MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sidm_courant,courant_max,1,MPI_DOUBLE_PRECISION, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  ! newdt_fine broadcasts an identical interval.  Reject even a one-bit
  ! disagreement rather than allowing rank-local loop termination to split
  ! the collective substep schedule.
  call sidm_collision_schedule_agrees(dt_code_min,dt_code_max, &
       & scale_t_min,scale_t_max,courant_min,courant_max,schedule_consistent)
  gate_failed=1-schedule_consistent
  call MPI_ALLREDUCE(gate_failed,gate_failed_global,1,MPI_INTEGER, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  if(gate_failed_global/=0) then
     call sidm_probability_guard_abort('MPI ranks disagree on collision interval or guard limit')
  end if

  call sidm_choose_collision_substep(dt_code_min,rate_global,scale_t_min, &
       & courant_min,dt_sub_code,p_bound,p_sub_bound,sidm_dt_invalid)
  call MPI_ALLREDUCE(sidm_dt_invalid,gate_failed_global,1,MPI_INTEGER, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  if(gate_failed_global/=0) then
     call sidm_probability_guard_abort('invalid or unsafe collision substep')
  end if

  sidm_rate_bound(ilevel)=max(sidm_rate_bound(ilevel),rate_global)
  rate_out=rate_global
  p_remaining_out=p_bound
  p_substep_out=p_sub_bound
  p_limit=min(0.1d0,courant_min)
  p_limit_global=p_limit
  if(myid==1 .and. sidm_estimator_diagnostics) then
     write(*,'(A,I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
          ' SIDM_RATE_PREFLIGHT: coarse=',nstep_coarse,' level=',ilevel, &
          ' rate_s=',rate_global,' remaining_dt_s=',dt_code_min*scale_t_min, &
          ' P_remaining=',p_bound,' sub_dt_s=',dt_sub_code*scale_t_min, &
          ' P_substep=',p_sub_bound,' P_limit=',p_limit_global
  end if
end subroutine sidm_preflight_rate_gate

! Fail the complete MPI job with a nonzero status for P-safe violations.
subroutine sidm_probability_guard_abort(message)
  use amr_commons, only: myid
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=*),intent(in)::message
  integer::info

  if(myid==1) write(*,'(A,A)') 'ERROR: SIDM probability guard: ',trim(message)
  flush(6)
#ifndef WITHOUTMPI
  call MPI_ABORT(MPI_COMM_WORLD,1,info)
#endif
  stop 1
end subroutine sidm_probability_guard_abort
!################################################################
subroutine sub_sidm_scatter(ilevel,icpu,kgrid,subnump,thread_seeds, &
     & mythread_loc,nthreads_in,dt_step_code,dt_interval_code)
  use pm_commons
  use amr_commons
  use random
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel,icpu,kgrid,subnump
  integer,intent(in)::mythread_loc,nthreads_in
  real(dp),intent(in)::dt_step_code,dt_interval_code
  integer,dimension(IRandNumSize,0:nthreads_in-1)::thread_seeds
  external sidm_sample_rutherford_cosine,sidm_rotate_scattered_direction
  external sidm_pair_probability

  ! External function
  integer,external::cell_index_from_part

  ! Shared counters (common with sidm_scatter)
  integer::sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  integer::sidm_n_vrel_reject
  integer::sidm_n_leaf_occupied,sidm_n_leaf_active,sidm_n_dm_leaf
  real(dp)::sidm_dp(3),sidm_dEk,sidm_dp_max,sidm_dEk_max
  real(dp)::sidm_Pmax_level
  real(dp)::sidm_dEdiss
  real(dp)::sidm_expected_rate
  common /sidm_diag/ sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down,sidm_n_vrel_reject
  common /sidm_cons/ sidm_dp,sidm_dEk,sidm_dp_max,sidm_dEk_max
  common /sidm_pmax_c/ sidm_Pmax_level
  common /sidm_diss/ sidm_dEdiss
  common /sidm_estimator_counts/ sidm_n_leaf_occupied, &
       & sidm_n_leaf_active,sidm_n_dm_leaf
  common /sidm_estimator_rate/ sidm_expected_rate

  ! Local variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,dt_phys,dt_interval_phys,sigma_over_m
  real(dp)::P_scatter,twopi
  real(dp)::cos_theta,sin_theta,phi_rand,v_rel_mag,v_rel_kms
  real(dp),dimension(1:3)::v1,v2,v_cm,v_rel_vec,nhat,v_rel_new
  real(dp),dimension(1:3)::p_before,p_after,v1_new,v2_new
  real(dp)::Ek_before,Ek_after,dp_mag,dEk_evt
  real(dp)::m1,m2,mtot,R1,R2
  integer::igrid,jgrid,ind,iskip,icell,nx_loc
  integer::ipart,jpart,npart1,ndm_cell,ip,jp,npairs
  ! Per-cell particle index buffer (grown to numbp(igrid) as needed;
  ! a fixed nvector-sized buffer silently dropped particles beyond
  ! nvector and clamped the (N-1) rate factor in dense cells)
  integer,dimension(:),allocatable::ind_dm
  integer::itemp,ipair
  integer,dimension(IRandNumSize)::seed_loc
  logical::invalid_probability
#ifndef WITHOUTMPI
  integer::mpi_abort_ierr
#endif

  ! Scattering counters (thread-local)
  integer::n_scatter_loc,n_pairs_loc,n_up_loc,n_down_loc
  integer::n_vrel_reject_loc
  integer::n_leaf_occupied_loc,n_leaf_active_loc,n_dm_leaf_loc
  ! Conservation diagnostics (thread-local)
  real(dp)::dp_loc(3),dEk_loc,dp_max_loc,dEk_max_loc
  real(dp)::P_max_loc
  real(dp)::dEdiss_loc,KE_cm_before_diss
  real(dp)::expected_rate_loc

  ! Anisotropic scattering variables
  real(dp)::eps2
  logical::yukawa_angular  ! velocity-dependent angular distribution

  ! Inelastic scattering variables
  real(dp)::mu,KE_cm,KE_cm_new,delta_KE,v_rel_new_mag
  real(dp)::mchi_g,mu_chi
  logical::do_transition
  integer::istate_1,istate_2,new_state_1,new_state_2
  integer::n_channels,ichan,j1,j2
  real(dp)::R_inel,Ei_1,Ei_2,dE
  ! Multi-state channel arrays (max N*(N+1)/2 = 55 for N=10)
  integer::chan_j1(55),chan_j2(55)
  real(dp)::chan_dE(55)

  if(subnump<=0) return

  ! Unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  twopi = 2.0d0*ACOS(-1.0d0)

  ! Cell size at this level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale

  ! Physical (proper) cell volume [cm^3]: scale_l is already proper
  vol_phys = (dx_loc*scale_l)**3

  ! Physical timestep [s]
  dt_phys = dt_step_code*scale_t
  dt_interval_phys = dt_interval_code*scale_t

  ! (Multi-state energies are accessed directly as sidm_energy(:) [keV]
  !  and converted to erg inline: E [erg] = E [keV] * 1.602d-9)

  ! Precompute Rutherford parameters
  yukawa_angular = (trim(sidm_angular) == 'yukawa')
  if(trim(sidm_angular) == 'rutherford') then
     eps2 = 2.0d0*sidm_epsilon
  end if
  ! For 'yukawa' angular, eps2 is computed per pair below.

  ! Use per-thread seed from thread_seeds array
  seed_loc = thread_seeds(:,mythread_loc)

  n_scatter_loc = 0
  n_pairs_loc   = 0
  n_up_loc      = 0
  n_down_loc    = 0
  n_vrel_reject_loc = 0
  n_leaf_occupied_loc = 0
  n_leaf_active_loc = 0
  n_dm_leaf_loc = 0
  dp_loc(:)     = 0.0d0
  dEk_loc       = 0.0d0
  dp_max_loc    = 0.0d0
  dEk_max_loc   = 0.0d0
  dEdiss_loc    = 0.0d0
  P_max_loc     = 0.0d0
  expected_rate_loc = 0.0d0

  allocate(ind_dm(1024))

  ! Loop over grids assigned to this thread
  igrid = kgrid
  do jgrid=1,subnump

     npart1 = numbp(igrid)
     if(npart1<sidm_npart_min) then
        igrid = next(igrid)
        cycle
     end if
     if(npart1>size(ind_dm)) then
        deallocate(ind_dm); allocate(ind_dm(npart1))
     end if

     ! Loop over 8 subcells (twotondim=8 in 3D)
     do ind=1,twotondim
        icell = ICELL_OF(igrid,ind)

        ! Only process leaf cells
        if(son(icell)/=0) cycle

        ! Collect DM particles in this subcell
        ! DM-like: idp>0 and ptypep is DM or excited iSIDM (not star, not sink)
        ndm_cell = 0
        ipart = headp(igrid)
        do jp=1,npart1
           ! Check if this particle belongs to subcell 'ind'
           if(cell_index_from_part(ipart,igrid,ilevel)==ind) then
              if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. ptypep(ipart)/=PTYPE_SINK) then
                 ndm_cell = ndm_cell+1
                 ind_dm(ndm_cell) = ipart
              end if
           end if
           ipart = nextp(ipart)
        end do

        if(ndm_cell>0) then
           n_leaf_occupied_loc = n_leaf_occupied_loc + 1
           n_dm_leaf_loc = n_dm_leaf_loc + ndm_cell
        end if

        ! Need at least sidm_npart_min DM particles
        if(ndm_cell<sidm_npart_min) cycle
        n_leaf_active_loc = n_leaf_active_loc + 1

        ! Fisher-Yates shuffle
        do ip=ndm_cell,2,-1
           call ranf(seed_loc, R1)
           jp = 1 + int(R1*dble(ip))
           if(jp>ip) jp = ip  ! safety clamp
           if(jp/=ip) then
              itemp = ind_dm(ip)
              ind_dm(ip) = ind_dm(jp)
              ind_dm(jp) = itemp
           end if
        end do

        ! Sequential pairing: (1,2), (3,4), ...
        npairs = ndm_cell/2
        n_pairs_loc = n_pairs_loc + npairs

        do ip=1,npairs
           ! Pair index into ind_dm array
           ipair = 2*ip-1
           ! Particle masses (physical) [g]
           m1 = mp(ind_dm(ipair))   * scale_d * scale_l**3
           m2 = mp(ind_dm(ipair+1)) * scale_d * scale_l**3
           mtot = m1 + m2

           ! Particle velocities (code units -> proper peculiar [cm/s]);
           ! scale_v already carries the full supercomoving conversion
           v1(1:3) = vp(ind_dm(ipair),   1:3) * scale_v
           v2(1:3) = vp(ind_dm(ipair+1), 1:3) * scale_v

           ! [DEBUG] PRE-scatter NaN trap
           if(v1(1)/=v1(1).or.v1(2)/=v1(2).or.v1(3)/=v1(3).or. &
              v2(1)/=v2(1).or.v2(2)/=v2(2).or.v2(3)/=v2(3)) then
              write(*,'(A,I6,A,I3,A,I8,A,2I12)') &
                   ' SIDM PRE-SCATTER NaN: rank=',myid,' ilevel=',ilevel, &
                   ' icell=',icell,' part=',ind_dm(ipair),ind_dm(ipair+1)
              write(*,'(A,3ES14.6)') '   vp1 (code units)=',vp(ind_dm(ipair),  1:3)
              write(*,'(A,3ES14.6)') '   vp2 (code units)=',vp(ind_dm(ipair+1),1:3)
              write(*,'(A,2ES14.6)') '   mp1, mp2=',mp(ind_dm(ipair)),mp(ind_dm(ipair+1))
              write(*,'(A,2I4)')     '   ptypep=',ptypep(ind_dm(ipair)),ptypep(ind_dm(ipair+1))
              write(*,'(A,3ES14.6)') '   aexp,scale_v,scale_d=',aexp,scale_v,scale_d
              stop 'SIDM_PRE_NAN'
           end if

           ! Relative velocity
           v_rel_vec(1:3) = v1(1:3) - v2(1:3)
           v_rel_mag = sqrt(v_rel_vec(1)**2 + v_rel_vec(2)**2 &
                          + v_rel_vec(3)**2)

           if(v_rel_mag==0.0d0) cycle

           ! Pathological-pair rejection: skip if relative velocity exceeds cap
           if(sidm_vrel_max > 0.0d0 .and. v_rel_mag > sidm_vrel_max) then
              n_vrel_reject_loc = n_vrel_reject_loc + 1
              cycle
           end if

           ! --- Velocity-dependent cross-section ---
           v_rel_kms = v_rel_mag * 1.0d-5  ! cm/s -> km/s
           select case(trim(sidm_type))
           case('yukawa')
              ! sigma(v) = sigma_0 / (1 + (v/v0)^2)^2
              sigma_over_m = sidm_cross_section &
                   / (1.0d0 + (v_rel_kms/sidm_v0)**2)**2
           case('power_law')
              ! sigma(v) = sigma_0 * (v/v0)^n, capped
              if(v_rel_kms > 1.0d0) then
                 sigma_over_m = sidm_cross_section &
                      * (v_rel_kms/sidm_v0)**sidm_power
              else
                 sigma_over_m = sidm_cross_section &
                      * (1.0d0/sidm_v0)**sidm_power
              end if
              sigma_over_m = min(sigma_over_m, &
                   1.0d4*sidm_cross_section)
           case default  ! 'constant'
              sigma_over_m = sidm_cross_section
           end select

           ! --- Dark phase transition: sigma(v,a) = sigma(v) * g(a) ---
           select case(trim(sidm_a_type))
           case('step')
              if(aexp > sidm_a_transition) then
                 sigma_over_m = sigma_over_m * sidm_sigma_ratio
              end if
           case('sigmoid')
              sigma_over_m = sigma_over_m * (1.0d0 + &
                   (sidm_sigma_ratio - 1.0d0) * 0.5d0 * &
                   (1.0d0 + tanh((aexp - sidm_a_transition)/sidm_a_width)))
           case default  ! 'none'
              ! no modification
           end select

           ! Scattering probability
           ! P = (sigma/m) * m_p * v_rel * dt / V_cell * (N_dm - 1),
           ! times N/(2*floor(N/2)) so that floor(N/2) sampled pairs
           ! reproduce the exact N(N-1)/2 pair rate for odd N too
           call sidm_pair_probability(sigma_over_m,m1,m2,v_rel_mag, &
                & dt_phys,vol_phys,ndm_cell,npairs,P_scatter)
           ! A Bernoulli draw cannot represent an invalid probability.  The
           ! next-step Pmax constraint cannot repair the current sweep.
           invalid_probability=.not.ieee_is_finite(P_scatter)
           if(.not.invalid_probability) then
              invalid_probability=(P_scatter<0.0d0.or.P_scatter>1.0d0)
           end if
           if(invalid_probability) then
              write(error_unit,'(A,ES24.16,A,I0,A,I0,A,I0)') &
                   'ERROR: invalid SIDM probability ',P_scatter, &
                   ' rank=',myid,' level=',ilevel,' cell=',icell
              flush(error_unit)
#ifndef WITHOUTMPI
              ! Any rank may encounter the invalid pair; do not wait for a
              ! collective finalization while the other ranks continue.
              call MPI_ABORT(MPI_COMM_WORLD,1,mpi_abort_ierr)
#endif
              error stop 1
           end if
           if(sidm_probability_guard) then
              if(P_scatter>min(0.1d0,sidm_courant)) &
                   call sidm_probability_guard_abort('sampled P exceeds guarded limit')
           end if
           if(dt_interval_phys>0.0d0) then
              expected_rate_loc = expected_rate_loc + P_scatter/dt_interval_phys
           end if

           ! Track P_max for timestep constraint
           if(P_scatter > P_max_loc) P_max_loc = P_scatter

           ! Monte Carlo: scatter if random < P
           call ranf(seed_loc, R1)
           if(R1 >= P_scatter) cycle

           ! === SCATTER EVENT ===
           n_scatter_loc = n_scatter_loc + 1

           ! Before: momentum and kinetic energy (physical)
           p_before(1:3) = m1*v1(1:3) + m2*v2(1:3)
           Ek_before = 0.5d0*(m1*(v1(1)**2+v1(2)**2+v1(3)**2) &
                             +m2*(v2(1)**2+v2(2)**2+v2(3)**2))

           ! Center-of-mass velocity (physical)
           v_cm(1:3) = (m1*v1(1:3) + m2*v2(1:3)) / mtot

           ! --- Inelastic energy change (multi-state iSIDM) ---
           delta_KE = 0.0d0
           do_transition = .false.
           new_state_1 = 0
           new_state_2 = 0
           if(sidm_inelastic) then
              ! Current states: derived from ptypep (PTYPE_DM=0, ISIDM_EX1=1, ISIDM_EX2=2)
              istate_1 = ptype_to_isidm_state(ptypep(ind_dm(ipair)))
              istate_2 = ptype_to_isidm_state(ptypep(ind_dm(ipair+1)))
              istate_1 = max(0, min(istate_1, sidm_nstates-1))
              istate_2 = max(0, min(istate_2, sidm_nstates-1))
              Ei_1 = sidm_energy(istate_1) * 1.602d-9  ! keV -> erg
              Ei_2 = sidm_energy(istate_2) * 1.602d-9
              ! Kinematics of the MICROPHYSICAL chi pair: an N-body
              ! macro-particle scatters as if a single chi of mass
              ! sidm_mchi [GeV] does (velocities are shared), so the
              ! splitting compares against the per-chi-pair CM energy,
              ! not the astronomical macro-particle KE
              mchi_g = sidm_mchi * 1.783d-24  ! GeV -> g
              mu_chi = 0.5d0*mchi_g
              KE_cm = 0.5d0 * mu_chi * v_rel_mag**2

              ! Enumerate energetically accessible inelastic channels.
              ! 2-state model: the off-diagonal coupling flips both
              ! particles (gg->ee at threshold 2*delta, ee->gg
              ! exothermic); mixed pairs have no open channel and
              ! scatter elastically. N>2 states: phenomenological
              ! uniform branching over all open pair channels.
              ! The elastic channel always competes with equal weight.
              n_channels = 0
              do j1=0,sidm_nstates-1
                 do j2=j1,sidm_nstates-1
                    if(j1==istate_1 .and. j2==istate_2) cycle  ! elastic
                    if(j1==istate_2 .and. j2==istate_1) cycle  ! trivial swap
                    if(sidm_nstates==2 .and. &
                       .not.(j1==j2 .and. istate_1==istate_2)) cycle
                    ! Energy change: dE = (E_j1+E_j2) - (E_i1+E_i2)
                    ! dE>0 = endothermic (costs KE), dE<0 = exothermic (releases KE)
                    dE = (sidm_energy(j1)+sidm_energy(j2) &
                         -sidm_energy(istate_1)-sidm_energy(istate_2)) * 1.602d-9
                    if(dE > 0.0d0 .and. KE_cm < dE) cycle     ! not enough energy
                    n_channels = n_channels + 1
                    chan_j1(n_channels) = j1
                    chan_j2(n_channels) = j2
                    chan_dE(n_channels) = dE
                 end do
              end do

              if(n_channels > 0) then
                 ! ichan = 0 selects the elastic channel
                 call ranf(seed_loc, R_inel)
                 ichan = min(int(R_inel*dble(n_channels+1)), n_channels)
                 if(ichan >= 1) then
                    new_state_1 = chan_j1(ichan)
                    new_state_2 = chan_j2(ichan)
                    delta_KE = -chan_dE(ichan)  ! negative dE = KE gained
                    do_transition = .true.
                    ! Count up/down transitions
                    if(new_state_1 > istate_1) n_up_loc = n_up_loc + 1
                    if(new_state_1 < istate_1) n_down_loc = n_down_loc + 1
                    if(new_state_2 > istate_2) n_up_loc = n_up_loc + 1
                    if(new_state_2 < istate_2) n_down_loc = n_down_loc + 1
                 end if
              end if
           end if

           ! New relative velocity magnitude (per-chi kinematics:
           ! v'^2 = v^2 + 2*delta_KE/mu_chi)
           mu = m1*m2/mtot
           if(delta_KE == 0.0d0) then
              v_rel_new_mag = v_rel_mag  ! elastic
           else
              KE_cm_new = 0.5d0*mu_chi*v_rel_mag**2 + delta_KE
              if(KE_cm_new < 0.0d0) cycle  ! safety
              v_rel_new_mag = sqrt(2.0d0*KE_cm_new/mu_chi)
           end if

           ! Dissipative SIDM: remove fraction f_diss of CM kinetic energy
           if(sidm_fdiss > 0.0d0) then
              KE_cm_before_diss = 0.5d0*mu*v_rel_new_mag**2
              v_rel_new_mag = v_rel_new_mag * sqrt(1.0d0 - sidm_fdiss)
              dEdiss_loc = dEdiss_loc + sidm_fdiss*KE_cm_before_diss
           end if

           ! --- Scattering angle ---
           if(yukawa_angular) then
              ! Yukawa Born: dsigma/dOmega ~ 1/(1+R^2*sin^2(theta/2))^2
              ! equivalent to Rutherford with eps = v0^2/v_rel^2
              ! R = v_rel/v0 (physical km/s)
              eps2 = 2.0d0*sidm_v0**2/max(v_rel_kms**2, 1.0d0)
           end if

           select case(trim(sidm_angular))
           case('rutherford','yukawa')
              ! Rutherford-like: dsigma/dOmega ~ 1/(1-cos_theta+2*eps)^2
              ! For 'yukawa': eps = v0^2/v_rel^2 (velocity-dependent)
              ! The CDF inversion is tested independently in sidm_angular.
              call ranf(seed_loc, R1)
              call ranf(seed_loc, R2)
              call sidm_sample_rutherford_cosine(R1,eps2,cos_theta)
              sin_theta = sqrt(max(0.0d0, 1.0d0 - cos_theta**2))
              phi_rand = twopi*R2

              ! Rotate scattering angle from v_rel frame to lab frame
              call sidm_rotate_scattered_direction(v_rel_vec,v_rel_mag, &
                   & cos_theta,sin_theta,phi_rand,nhat)

           case default  ! 'isotropic'
              ! Random unit vector (isotropic)
              call ranf(seed_loc, R1)
              call ranf(seed_loc, R2)
              cos_theta = 2.0d0*R1 - 1.0d0
              sin_theta = sqrt(1.0d0 - cos_theta**2)
              phi_rand  = twopi*R2
              nhat(1) = sin_theta*cos(phi_rand)
              nhat(2) = sin_theta*sin(phi_rand)
              nhat(3) = cos_theta
           end select

           ! New relative velocity vector
           v_rel_new(1:3) = v_rel_new_mag * nhat(1:3)

           ! New velocities in COM frame
           v1_new(1:3) = v_cm(1:3) + (m2/mtot)*v_rel_new(1:3)
           v2_new(1:3) = v_cm(1:3) - (m1/mtot)*v_rel_new(1:3)

           ! [DEBUG] POST-scatter NaN trap
           if(v1_new(1)/=v1_new(1).or.v1_new(2)/=v1_new(2).or.v1_new(3)/=v1_new(3).or. &
              v2_new(1)/=v2_new(1).or.v2_new(2)/=v2_new(2).or.v2_new(3)/=v2_new(3)) then
              write(*,'(A,I6,A,I3,A,I8)') &
                   ' SIDM POST-SCATTER NaN: rank=',myid,' ilevel=',ilevel,' icell=',icell
              write(*,'(A,2I12)')    '   ind_dm pair=', ind_dm(ipair), ind_dm(ipair+1)
              write(*,'(A,2ES14.6)') '   m1, m2 [g]=', m1, m2
              write(*,'(A,3ES14.6)') '   v1=', v1(1:3)
              write(*,'(A,3ES14.6)') '   v2=', v2(1:3)
              write(*,'(A,3ES14.6)') '   v_cm=', v_cm(1:3)
              write(*,'(A,3ES14.6)') '   v_rel_vec=', v_rel_vec(1:3)
              write(*,'(A,2ES14.6)') '   v_rel_mag, v_rel_kms=', v_rel_mag, v_rel_kms
              write(*,'(A,3ES14.6)') '   cos_theta,sin_theta,phi_rand=', cos_theta,sin_theta,phi_rand
              write(*,'(A,3ES14.6)') '   nhat=', nhat(1:3)
              write(*,'(A,4ES14.6)') '   v_rel_new=', v_rel_new(1:3), v_rel_new_mag
              write(*,'(A,3ES14.6)') '   v1_new=', v1_new(1:3)
              write(*,'(A,3ES14.6)') '   v2_new=', v2_new(1:3)
              write(*,'(A,3ES14.6)') '   aexp,scale_v,sigma_over_m=', aexp,scale_v,sigma_over_m
              stop 'SIDM_POST_NAN'
           end if

           ! Apply state transitions (multi-state iSIDM)
           if(do_transition) then
              ptypep(ind_dm(ipair))   = isidm_state_to_ptype(new_state_1)
              ptypep(ind_dm(ipair+1)) = isidm_state_to_ptype(new_state_2)
           end if

           ! Conservation diagnostics (momentum: exact, energy: exact
           ! for elastic; intentional change for inelastic)
           p_after(1:3) = m1*v1_new(1:3) + m2*v2_new(1:3)
           Ek_after = 0.5d0*(m1*(v1_new(1)**2+v1_new(2)**2 &
                                +v1_new(3)**2) &
                            +m2*(v2_new(1)**2+v2_new(2)**2 &
                                +v2_new(3)**2))
           dp_loc(1:3) = dp_loc(1:3) + (p_after(1:3) - p_before(1:3))
           dEk_evt = Ek_after - Ek_before
           dEk_loc = dEk_loc + dEk_evt
           dp_mag = sqrt((p_after(1)-p_before(1))**2 &
                        +(p_after(2)-p_before(2))**2 &
                        +(p_after(3)-p_before(3))**2)
           if(dp_mag > dp_max_loc) dp_max_loc = dp_mag
           if(abs(dEk_evt) > dEk_max_loc) dEk_max_loc = abs(dEk_evt)

           ! Write new velocities (physical -> code units)
           vp(ind_dm(ipair),   1:3) = v1_new(1:3) / scale_v
           vp(ind_dm(ipair+1), 1:3) = v2_new(1:3) / scale_v
        end do  ! pairs

     end do  ! subcells (ind)

     igrid = next(igrid)
  end do  ! grids

  deallocate(ind_dm)

  ! Save seed back to per-thread array (no race)
  thread_seeds(:,mythread_loc) = seed_loc

  ! Accumulate thread-local counters to shared counters
  !$omp critical
  sidm_n_scatter = sidm_n_scatter + n_scatter_loc
  sidm_n_pairs   = sidm_n_pairs + n_pairs_loc
  sidm_n_up      = sidm_n_up + n_up_loc
  sidm_n_down    = sidm_n_down + n_down_loc
  sidm_n_vrel_reject = sidm_n_vrel_reject + n_vrel_reject_loc
  sidm_n_leaf_occupied = sidm_n_leaf_occupied + n_leaf_occupied_loc
  sidm_n_leaf_active = sidm_n_leaf_active + n_leaf_active_loc
  sidm_n_dm_leaf = sidm_n_dm_leaf + n_dm_leaf_loc
  sidm_dp(1:3)   = sidm_dp(1:3) + dp_loc(1:3)
  sidm_dEk       = sidm_dEk + dEk_loc
  if(dp_max_loc > sidm_dp_max) sidm_dp_max = dp_max_loc
  if(dEk_max_loc > sidm_dEk_max) sidm_dEk_max = dEk_max_loc
  if(P_max_loc > sidm_Pmax_level) sidm_Pmax_level = P_max_loc
  sidm_dEdiss = sidm_dEdiss + dEdiss_loc
  sidm_expected_rate = sidm_expected_rate + expected_rate_loc
  !$omp end critical

end subroutine sub_sidm_scatter
!################################################################
!################################################################
! Initialize excited DM fraction for iSIDM
!################################################################
subroutine sidm_init_excited()
  use pm_commons
  use pm_parameters, only: iseed
  use amr_commons
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icpu,igrid,jgrid,ipart,jpart,npart1,is
  integer::ndm_total,ndm_existing,info
  integer,dimension(0:9)::ndm_state
  integer,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer,dimension(1:IRandNumSize)::seed_init
  real(dp)::R1,cum_frac
  real(dp),dimension(0:9)::cum_prob  ! Cumulative probability for N states

  ! Check if any non-ground population is requested
  if(sum(sidm_frac_init(1:sidm_nstates-1)) <= 0.0d0) return

  ! Build cumulative probability
  cum_prob(0) = sidm_frac_init(0)
  do is=1,sidm_nstates-1
     cum_prob(is) = cum_prob(is-1) + sidm_frac_init(is)
  end do
  ! Normalize (in case fractions don't sum to 1)
  if(cum_prob(sidm_nstates-1) > 0.0d0) then
     cum_prob = cum_prob / cum_prob(sidm_nstates-1)
  end if

  ! Check if already initialized (restart case: some DM have excited ptype)
  ndm_existing = 0
  do ilevel=levelmin,nlevelmax
     igrid = headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        ipart = headp(igrid)
        do jpart=1,numbp(igrid)
           if(idp(ipart)>0 .and. (ptypep(ipart)==PTYPE_ISIDM_EX1 .or. ptypep(ipart)==PTYPE_ISIDM_EX2)) then
              ndm_existing = ndm_existing + 1
           end if
           ipart = nextp(ipart)
        end do
        igrid = next(igrid)
     end do
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_existing,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  if(ndm_existing > 0) then
     if(myid==1) write(*,'(A,I10,A)') &
          ' iSIDM: found ', ndm_existing, &
          ' pre-existing excited DM (skipping init)'
     return
  end if

  ! Initialize seed (offset to avoid correlation with main RNG)
  call rans(ncpu,iseed+12345,allseed)
  seed_init = allseed(myid,1:IRandNumSize)

  ndm_total = 0
  ndm_state = 0

  ! Loop over own grids and randomly assign states
  do ilevel=levelmin,nlevelmax
     igrid = headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        ipart = headp(igrid)
        do jpart=1,numbp(igrid)
           ! DM: idp>0, ptypep==PTYPE_DM (ground state, not yet assigned)
           if(idp(ipart)>0 .and. ptypep(ipart)==PTYPE_DM) then
              ndm_total = ndm_total + 1
              call ranf(seed_init, R1)
              ! Assign state by cumulative probability
              do is=0,sidm_nstates-1
                 if(R1 < cum_prob(is)) then
                    ptypep(ipart) = isidm_state_to_ptype(is)
                    ndm_state(is) = ndm_state(is) + 1
                    exit
                 end if
              end do
           end if
           ipart = nextp(ipart)
        end do
        igrid = next(igrid)
     end do
  end do

  ! Report
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_total,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_state,sidm_nstates,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  if(myid==1) then
     write(*,'(A,I10,A,I3,A)') &
          ' iSIDM: initialized ', ndm_total, ' DM into ', sidm_nstates, ' states:'
     do is=0,sidm_nstates-1
        write(*,'(A,I2,A,I10,A,F7.4)') &
             '   state ', is, ': N=', ndm_state(is), &
             '  frac=', dble(ndm_state(is))/dble(max(1,ndm_total))
     end do
  end if

end subroutine sidm_init_excited
!################################################################
!################################################################
! Report excited DM fraction (called from adaptive_loop)
!################################################################
subroutine sidm_report_excited_fraction()
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,igrid,jgrid,ipart,jpart,info,is,istate
  integer::ndm_total
  integer,dimension(0:9)::ndm_state

  if(.not.sidm_inelastic) return

  ndm_total = 0
  ndm_state = 0
  do ilevel=levelmin,nlevelmax
     igrid = headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        ipart = headp(igrid)
        do jpart=1,numbp(igrid)
           if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. ptypep(ipart)/=PTYPE_SINK) then
              ndm_total = ndm_total + 1
              istate = ptype_to_isidm_state(ptypep(ipart))
              istate = max(0, min(istate, 9))
              ndm_state(istate) = ndm_state(istate) + 1
           end if
           ipart = nextp(ipart)
        end do
        igrid = next(igrid)
     end do
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_total,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_state,sidm_nstates,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)

  if(myid==1 .and. ndm_total>0) then
     if(sidm_nstates==2) then
        ! Compact format for 2-state (backward compatible output)
        write(*,'(A,I10,A,I10,A,F8.5)') &
             ' iSIDM excited: ', ndm_state(1), ' / ', ndm_total, &
             '  f_exc=', dble(ndm_state(1))/dble(ndm_total)
     else
        write(*,'(A,I10,A,I3,A)') &
             ' iSIDM populations: N_DM=', ndm_total, '  (', sidm_nstates, ' states)'
        do is=0,sidm_nstates-1
           write(*,'(A,I2,A,I10,A,F7.4)') &
                '   state ', is, ': ', ndm_state(is), &
                '  f=', dble(ndm_state(is))/dble(ndm_total)
        end do
     end if
  end if

end subroutine sidm_report_excited_fraction
!################################################################
!################################################################
! DM-baryon drag force: continuous momentum exchange
!
! Physics: F/m_DM = -(sigma_db/m_DM) * rho_b * |v_rel| * v_rel
! Velocity-dependent: sigma(v) = sigma_0 * (v/100 km/s)^n
! Backreaction on gas conserves total momentum.
!################################################################
subroutine sidm_baryon_drag(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons, only: uold,nvar,gamma
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,vol_code,dt_phys
  real(dp)::rho_b_phys,sigma_eff
  real(dp),dimension(1:3)::v_gas_phys,v_dm_phys,v_rel,dv_dm,dv_gas
  real(dp)::v_rel_mag,v_rel_kms,drag_coeff
  real(dp)::m_dm_phys,m_gas_phys,dEk_dm,dv_max,dv_mag
  integer::igrid,jgrid,ind,iskip,icell,nx_loc
  integer::ipart,jpart,npart1,info
  integer,external::cell_index_from_part
  ! Diagnostics
  integer::n_drag_loc,n_drag_tot
  real(dp)::dp_drag_loc(3),dp_drag_tot(3),dv_max_loc,dv_max_tot

  if(.not.sidm_baryon) return
  if(.not.hydro) return
  if(numbtot(1,ilevel)==0) return

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale
  vol_code = dx_loc**3
  vol_phys = (dx_loc*scale_l)**3  ! proper volume (scale_l is proper)
  dt_phys  = dtnew(ilevel)*scale_t

  n_drag_loc  = 0
  dp_drag_loc = 0.0d0
  dv_max_loc  = 0.0d0

  ! Loop over local grids
  igrid = headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)

     npart1 = numbp(igrid)
     if(npart1==0) then
        igrid = next(igrid)
        cycle
     end if

     ! Loop over 8 subcells
     do ind=1,twotondim
        icell = ICELL_OF(igrid,ind)

        ! Only leaf cells with gas
        if(son(icell)/=0) cycle

        ! Gas density and velocity (code units -> proper cgs)
        if(uold(icell,1) <= 0.0d0) cycle
        rho_b_phys = uold(icell,1) * scale_d
        m_gas_phys = rho_b_phys * vol_phys

        v_gas_phys(1) = uold(icell,2)/uold(icell,1) * scale_v
        v_gas_phys(2) = uold(icell,3)/uold(icell,1) * scale_v
        v_gas_phys(3) = uold(icell,4)/uold(icell,1) * scale_v

        ! Loop over particles in this grid
        ipart = headp(igrid)
        do jpart=1,npart1
           ! DM only (idp>0, ptypep is DM or excited iSIDM)
           if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. ptypep(ipart)/=PTYPE_SINK) then
              ! Check subcell match
              if(cell_index_from_part(ipart,igrid,ilevel)==ind) then
                 ! DM velocity (proper peculiar)
                 v_dm_phys(1) = vp(ipart,1) * scale_v
                 v_dm_phys(2) = vp(ipart,2) * scale_v
                 v_dm_phys(3) = vp(ipart,3) * scale_v

                 ! Relative velocity
                 v_rel(1:3) = v_dm_phys(1:3) - v_gas_phys(1:3)
                 v_rel_mag = sqrt(v_rel(1)**2+v_rel(2)**2+v_rel(3)**2)
                 if(v_rel_mag == 0.0d0) then
                    ipart = nextp(ipart)
                    cycle
                 end if

                 ! Velocity-dependent cross-section
                 v_rel_kms = v_rel_mag * 1.0d-5
                 if(sidm_baryon_power == 0.0d0) then
                    sigma_eff = sidm_baryon_sigma
                 else
                    sigma_eff = sidm_baryon_sigma &
                         * (v_rel_kms/100.0d0)**sidm_baryon_power
                    sigma_eff = min(sigma_eff, 1.0d4*sidm_baryon_sigma)
                 end if

                 ! Drag: dv = -(sigma/m) * rho_b * |v_rel| * v_rel * dt
                 drag_coeff = sigma_eff * rho_b_phys * v_rel_mag * dt_phys
                 dv_dm(1:3) = -drag_coeff * v_rel(1:3)

                 ! Limit kick: <=50% of v_rel for the particle AND for
                 ! the implied gas kick (gas inertia matters in
                 ! DM-dominated cells); together the caps keep the
                 ! frictional heating strictly positive
                 m_dm_phys = mp(ipart) * scale_d * scale_l**3
                 dv_mag = sqrt(dv_dm(1)**2+dv_dm(2)**2+dv_dm(3)**2)
                 dv_max = 0.5d0 * v_rel_mag * min(1.0d0, m_gas_phys/m_dm_phys)
                 if(dv_mag > dv_max) then
                    dv_dm(1:3) = dv_dm(1:3) * (dv_max/dv_mag)
                    dv_mag = dv_max
                 end if

                 ! Apply to DM particle (physical -> code)
                 vp(ipart,1) = vp(ipart,1) + dv_dm(1)/scale_v
                 vp(ipart,2) = vp(ipart,2) + dv_dm(2)/scale_v
                 vp(ipart,3) = vp(ipart,3) + dv_dm(3)/scale_v

                 ! Backreaction on gas: momentum -m_DM*dv and total
                 ! energy -dKE_DM, so total DM+gas momentum and energy
                 ! are conserved (bulk + frictional heat)
                 dEk_dm = m_dm_phys*( v_dm_phys(1)*dv_dm(1) &
                        + v_dm_phys(2)*dv_dm(2) + v_dm_phys(3)*dv_dm(3) &
                        + 0.5d0*(dv_dm(1)**2+dv_dm(2)**2+dv_dm(3)**2) )
                 uold(icell,2) = uold(icell,2) &
                      - m_dm_phys*dv_dm(1)/vol_phys/(scale_d*scale_v)
                 uold(icell,3) = uold(icell,3) &
                      - m_dm_phys*dv_dm(2)/vol_phys/(scale_d*scale_v)
                 uold(icell,4) = uold(icell,4) &
                      - m_dm_phys*dv_dm(3)/vol_phys/(scale_d*scale_v)
                 uold(icell,5) = uold(icell,5) &
                      - dEk_dm/vol_phys/(scale_d*scale_v**2)

                 ! Keep the local gas velocity current so subsequent
                 ! particles in this cell see the updated bulk flow
                 dv_gas(1:3) = -m_dm_phys*dv_dm(1:3)/m_gas_phys
                 v_gas_phys(1:3) = v_gas_phys(1:3) + dv_gas(1:3)

                 ! Diagnostics
                 n_drag_loc = n_drag_loc + 1
                 dp_drag_loc(1:3) = dp_drag_loc(1:3) + m_dm_phys*dv_dm(1:3)
                 if(dv_mag > dv_max_loc) dv_max_loc = dv_mag
              end if
           end if
           ipart = nextp(ipart)
        end do  ! particles

     end do  ! subcells

     igrid = next(igrid)
  end do  ! grids

  ! MPI reduce diagnostics
  call MPI_ALLREDUCE(n_drag_loc, n_drag_tot, 1, MPI_INTEGER, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(dp_drag_loc, dp_drag_tot, 3, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(dv_max_loc, dv_max_tot, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, info)

  if(myid==1 .and. n_drag_tot>0) then
     write(*,'(A,I2,A,I10,A,ES9.2,A)') &
          ' IDM level ',ilevel,': n_drag=',n_drag_tot, &
          ' dv_max=',dv_max_tot*1d-5,' km/s'
  end if

end subroutine sidm_baryon_drag
