subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  use dust_mass_physics, only: dust_chimes_enabled
  use ksection, only: ksection_trim_done
  use morton_hash, only: morton_hash_epoch
#ifdef RT
  use rt_hydro_commons
#endif
#ifdef HYDRO_CUDA
  use poisson_cuda_interface, only: cuda_fft_print_timers_c
  use particle_cuda_interface, only: cuda_pm_report_c
  use cuda_commons, only: cuda_pool_init_f, cuda_available
  use iso_c_binding, only: c_int
#endif
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer(kind=8)::n_step
  integer::ilevel,idim,ivar,info,tot_pt
!jhshin1
  real(kind=8)::tt1,tt2,muspt,muspt_this_step,wallsec,dumpsec
!jhshin2
  real(kind=4)::real_mem,real_mem_tot,real_mem_min,real_mem_max
!jhshin1
  real(kind=8),save::tstart=0.0
!jhshin2
  ! SFR tracking
  real(dp)::scale_l_s,scale_t_s,scale_d_s,scale_v_s,scale_nH_s,scale_T2_s
  real(dp)::mstar_tot_glob,scale_m_sfr,V_box_sfr,sfr_inst
  real(dp),save::mstar_glob_old=0.0d0, t_sfr_old=-1.0d30
  ! SGS diagnostics
  real(dp)::sgs_loc(4),sgs_glob(4)
  real(dp)::esgs_tot,pturb_max,pth_at_max,csgs_ratio_max
  integer::icell,igrid_sgs,ind_sgs,ncell_sgs
  ! Sink/AGN coarse-step log
  integer::isink_log

#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
!jhshin1
  ! for calculating total run time
  if (tstart.eq.0.0) then
     tstart = MPI_WTIME()
  end if
!jhshin2
#endif

  call init_amr                      ! Initialize AMR variables

  ! Diagnostic: test MPI health after init_amr (build_comm may stress IB)
#ifndef WITHOUTMPI
  if(myid==1) write(*,*) 'Diagnostic: testing MPI after init_amr...'
  call flush(6)
  call MPI_BARRIER(MPI_COMM_WORLD, info)
  tot_pt = myid
  call MPI_ALLREDUCE(MPI_IN_PLACE, tot_pt, 1, MPI_INTEGER, &
       & MPI_SUM, MPI_COMM_WORLD, info)
  if(myid==1) write(*,'(A,I6)') &
       & ' Diagnostic: ALLREDUCE OK, sum=', tot_pt
  call flush(6)
#endif

  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       & call rt_init_hydro          ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
  if(use_fdm)call fdm_compute_hbar   ! Compute effective hbar in code units
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(cosmo)then
     ! Compute cooling table at current aexp
  endif
#else  
  ! CHIMES owns the thermal rates in this admitted mode. Retain init_time's
  ! initial gas temperature, but do not build the unused RAMSES CIE table.
  if(cooling.and..not.neq_chem.and..not.dust_chimes_enabled()) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again
#ifdef SNRT_CHIMES
  call initialize_chimes_live_state
#endif

  ! Initialize FDM wavefunction (after AMR grid and density are set)
  if(use_fdm .and. nrestart==0) call fdm_init_psi

  ! Rebuild Morton hash tables after full AMR grid is available
  if(nrestart==0) call morton_hash_build_and_verify()

#ifndef WITHOUTMPI
  muspt=0.
  tot_pt=-1
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif


  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  ! Initialize polytropic floor coefficients for eEOS in hydro solver
  call init_eeos_poly_coeff

#ifdef HYDRO_CUDA
  ! Early CUDA pool init: must happen before first multigrid_fine call so the
  ! cuFFT direct-solve gate (cuda_pool_is_initialized_c()/=0) passes on every
  ! rank. Previously init was lazy in godunov_fine, missing the first force solve.
  if(gpu_hydro .or. gpu_poisson .or. gpu_fft .or. gpu_sink .or. gpu_scalar &
       & .or. gpu_particle .or. mhd_gpu_faces) then
     call cuda_pool_init_f()
     if(myid==1) write(*,'(A,L1)') ' Adaptive loop: CUDA pool early-init, available=', cuda_available
  end if
#endif

  if(myid==1)write(*,*)'Starting time integration'
  call sidm_binary_orbit_monitor_init()

  do ! Main time loop
                               call timer('coarse levels','start')

#ifndef WITHOUTMPI
     tt1=MPI_WTIME()
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
              call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              call make_virtual_fine_dp_bulk(rtuold,nrtvar,ilevel)
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              call make_virtual_fine_dp_bulk(f,ndim,ilevel)
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
        ! Clear trim flag after exchange routines have processed it
        call ksection_trim_done()
     endif

     ! Call base level
     call amr_step(levelmin,1)

     ! make_grid/kill_grid and build_comm update Morton tables incrementally.
     ! A full rebuild is needed only after defrag has renumbered every grid.
     if(morton_hash_epoch /= amr_mesh_epoch) call morton_hash_rebuild
#ifdef FDMDEBUG
     ! Full self-verification is diagnostic-only and scans every local grid.
     ! Restrict the billion-cell FDM replay to the historical failure window;
     ! the hash rebuild itself remains unconditional because refinement can
     ! add or remove grids on every coarse step.
     if(.not.use_fdm .or. levelmin<10 .or. nstep_coarse>=269) &
          call morton_hash_verify('step')
#endif

                               call timer('coarse levels','start')

     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then

        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
              call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt)then
              call rt_upload_fine(ilevel)
              call make_virtual_fine_dp_bulk(rtuold,nrtvar,ilevel)
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              call make_virtual_fine_dp_bulk(f,ndim,ilevel)
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
!		if(lightcone) then
!			call output_cone_hydro(1)
!			call output_cone_part(1)
!			call output_cone_sink(1)
!			call output_cone_hydro(2)
!			call output_cone_part(2)
!			call output_cone_sink(2)
!		endif

 !       aexp_old2=aexp_old
!		if(myid .eq.1) write(*,*)aexp_old
!		if(myid .eq.1) write(*,*)aexp
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

     call check_jobcontrol

     ! Sink/AGN log every coarse step (append to sink_log.csv)
     if(sink .and. nsink>0 .and. myid==1) then
        call units(scale_l_s,scale_t_s,scale_d_s,scale_v_s,scale_nH_s,scale_T2_s)
        open(unit=667, file='sink_log.csv', position='append', &
             form='formatted', status='unknown')
        do isink_log=1, nsink
           write(667,'(I8,",",ES14.6,",",F8.4,",",I10,6(",",ES14.6))') &
                nstep_coarse, aexp, 1d0/aexp-1d0, &
                idsink(isink_log), &
                msink(isink_log)*scale_d_s*scale_l_s**3/1.98892d33, &
                xsink(isink_log,1), xsink(isink_log,2), &
                xsink(isink_log,3), &
                dMBH_coarse(isink_log), spinmag(isink_log)
        end do
        close(667)
     end if

!jhshin1 -walltime setting
#ifndef WITHOUTMPI
     tt2=MPI_WTIME()
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_max,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_min,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
        real_mem_tot=real_mem_max
        ! Global star mass for SFR diagnostic
        ! mstar_tot is already global (identical on all ranks from star_formation ALLREDUCE)
        ! No need for MPI_ALLREDUCE here — just use it directly
        mstar_tot_glob=mstar_tot
        if(myid==1)then
           if (tot_pt==0) muspt=0. ! dont count first timestep
           n_step = int(numbtot(1,levelmin),kind=8)*twotondim
           do ilevel=levelmin+1,nlevelmax
             n_step = n_step + int(numbtot(1,ilevel),kind=8)*product(nsubcycle(levelmin:ilevel-1))*(twotondim-1)
           enddo
           muspt_this_step = (tt2-tt1)*1e6/n_step*ncpu
           muspt = muspt + muspt_this_step
           tot_pt = tot_pt + 1
           write(*,'(a,f8.2,a,f12.2,a,f12.2,a)')' Time elapsed since last coarse step:',tt2-tt1 &
          ,' s',muspt_this_step,' mus/pt'  &
          ,muspt / max(tot_pt,1), ' mus/pt (av)'
           call writemem(real_mem_tot)
           call writemem_minmax(real_mem_min,real_mem_max)
           write(*,*)'Total running time:', NINT((tt2-tstart)*100.0)*0.01,'s'
#ifdef HYDRO_CUDA
           if(gpu_fft) call cuda_fft_print_timers_c(int(myid, c_int))
           if(gpu_particle .and. myid==1) call cuda_pm_report_c()
#endif
           ! SFR diagnostic
           if(nstar_tot>0)then
              call units(scale_l_s,scale_t_s,scale_d_s,scale_v_s,scale_nH_s,scale_T2_s)
              scale_m_sfr=scale_d_s*scale_l_s**3
              V_box_sfr=(boxlen_ini*3.08d24/(h0/100d0)/3.0857d24)**3
              sfr_inst=0.0d0
              if(t_sfr_old>-1.0d29 .and. t>t_sfr_old)then
                 sfr_inst=(mstar_tot_glob-mstar_glob_old)*scale_m_sfr/1.98892d33 &
                         /((t-t_sfr_old)*scale_t_s/3.15576d7)
              endif
              write(*,'(" SFR: N*=",I10," M*=",1PE10.3," Msun",' // &
                 '" SFR=",1PE10.3," Msun/yr",' // &
                 '" SFRD=",1PE10.3," Msun/yr/Mpc3  z=",0PF7.2)') &
                 nstar_tot, &
                 mstar_tot_glob*scale_m_sfr/1.98892d33, &
                 sfr_inst, sfr_inst/V_box_sfr, &
                 1.0d0/aexp-1.0d0
              mstar_glob_old=mstar_tot_glob
              t_sfr_old=t
           endif
           ! FPR diagnostic: per-level dx_phys and m_refine_eff/m_refine ratio
           if(dr_proper > 0.0d0 .and. cosmo) then
              call units(scale_l_s,scale_t_s,scale_d_s,scale_v_s,scale_nH_s,scale_T2_s)
              write(*,'(A)') ' FPR (Gnedin 2016):'
              do ilevel=levelmin,nlevelmax
                 call compute_fpr_m_refine_eff(ilevel)
                 if(m_refine(ilevel)>0.0d0) then
                    write(*,'(A,I2,A,F8.3,A,F8.2)') &
                         '   lv',ilevel,' dx_phys=', &
                         (0.5d0**ilevel)*boxlen/dble(icoarse_max-icoarse_min+1) &
                         *aexp*boxlen_ini*1000d0/(h0/100d0), &
                         ' kpc  m_eff/m=', m_refine_eff(ilevel)/m_refine(ilevel)
                 end if
              end do
           end if
        endif
        ! iSIDM excited fraction diagnostic (must be outside myid==1 for MPI_ALLREDUCE)
        if(sidm .and. sidm_inelastic) then
           call sidm_report_excited_fraction()
        end if
        ! FDM diagnostics (must be outside myid==1 for MPI_ALLREDUCE)
        if(use_fdm) call fdm_diagnostics()
        ! SGS turbulence diagnostics
        if(use_sgs .and. isgs>0) then
           ! Collect local SGS statistics over leaf cells at levelmin
           sgs_loc=0d0  ! (1)=sum(rho*esgs*vol), (2)=max(Pturb/Pth), (3)=ncells, (4)=sum(esgs)
           do ilevel=levelmin,nlevelmax
              do igrid_sgs=1,active(ilevel)%ngrid
                 do ind_sgs=1,twotondim
                    icell=ICELL_OF(active(ilevel)%igrid(igrid_sgs),ind_sgs)
                    if(son(icell)/=0) cycle
                    sgs_loc(3)=sgs_loc(3)+1d0
                    sgs_loc(1)=sgs_loc(1)+uold(icell,isgs)  ! sum(rho*esgs)
                    ! P_turb/P_th ratio
                    if(uold(icell,1)>0d0) then
                       esgs_tot=max(uold(icell,isgs)/uold(icell,1),0d0)
                       pturb_max=(2d0/3d0)*uold(icell,1)*esgs_tot
                       ! thermal pressure
                       pth_at_max=(gamma-1d0)*max(uold(icell,ndim+2) &
                            -0.5d0*(uold(icell,2)**2+uold(icell,3)**2+uold(icell,4)**2)/uold(icell,1),smallr*smallc**2)
                       if(pth_at_max>0d0) sgs_loc(2)=max(sgs_loc(2),pturb_max/pth_at_max)
                       sgs_loc(4)=sgs_loc(4)+esgs_tot
                    end if
                 end do
              end do
           end do
           call MPI_ALLREDUCE(sgs_loc,sgs_glob,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           ! Fix max: use MPI_MAX for element 2
           call MPI_ALLREDUCE(sgs_loc(2),sgs_glob(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
           if(myid==1 .and. sgs_glob(3)>0d0) then
              write(*,'(" SGS: E_sgs=",1PE10.3," <e_sgs>=",1PE10.3," max(Pt/Pth)=",1PE10.3)') &
                   sgs_glob(1), sgs_glob(4)/sgs_glob(3), sgs_glob(2)
           end if
        end if
        if(walltime_hrs.gt.0d0) then
           wallsec = walltime_hrs*3600.     ! Convert from hours to seconds
           dumpsec = minutes_dump*60.       ! Convert minutes before end to seconds
           if(wallsec-dumpsec.lt.tt2-tstart) then
              output_now=.true.
              if(myid==1) write(*,*) 'Dumping snapshot before walltime runs out'
              ! Now set walltime to a negative number so we don't keep printing
              ! outputs
              walltime_hrs = -1d0
           endif
        endif
     endif
#endif
!jhshin2

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop

! Optional target-binary monitor.  Checkpoints are emitted only from the
! same-time sink synchronization hook in synchro_fine; this initializer only
! validates the run and creates the ledger/configuration header.
subroutine sidm_binary_orbit_monitor_init
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use pm_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::local_error,ios,info,ipart,isink,ilevel,igrid,jgrid,jpart
  integer::sink_count(2),global_sink_count(2)
  logical::terminal_exists
  real(dp)::dx_min

  if(.not.binary_orbit_monitor)return
  local_error=0
  if(binary_orbit_monitor_initialized/=0) local_error=1
  if(ndim/=3.or.cosmo.or.hydro.or.sidm.or..not.sink.or..not.poisson.or.bondi.or. &
       .not.pic.or.static) &
       local_error=1
  if(nrestart/=0.or.create_sinks.or.nsink/=2) local_error=1
  if(binary_orbit_stop_time<=t.or..not.ieee_is_finite(binary_orbit_stop_time)) &
       local_error=1
  if(binary_orbit_sample_interval<=0d0.or. &
       .not.ieee_is_finite(binary_orbit_sample_interval)) local_error=1
  if(collisionless_courant_factor<=0d0.or. &
       collisionless_courant_factor>1d0.or. &
       .not.ieee_is_finite(collisionless_courant_factor)) local_error=1
  if(nstepmax<=0.or.foutput<=nstepmax.or.fbackup<=nstepmax) local_error=1
  if(.not.ieee_is_finite(walltime_hrs).or.walltime_hrs>0d0) local_error=1
  if(noutput<1.or.noutput>MAXOUT)then
     local_error=1
  else
     if(any(.not.ieee_is_finite(tout(1:noutput))).or. &
          any(.not.ieee_is_finite(aout(1:noutput)))) local_error=1
     if(any(tout(1:noutput)<=binary_orbit_stop_time).or. &
          any(aout(1:noutput)<=aexp)) local_error=1
  endif
  if(nlevelmax>levelmin)then
     if(any(nsubcycle(levelmin:nlevelmax-1)/=1)) local_error=1
  endif
  if(.not.allocated(idsink).or..not.allocated(msink).or. &
       .not.allocated(xsink).or..not.allocated(vsink).or. &
       .not.allocated(oksink_all)) local_error=1
  if(local_error==0)then
     if(size(idsink)<2.or.size(msink)<2.or.size(xsink,1)<2.or. &
          size(vsink,1)<2.or.size(oksink_all)<2.or. &
          size(xsink,2)<3.or.size(vsink,2)<3)then
        local_error=1
     else
        if(idsink(1)==idsink(2)) local_error=1
        if(any(.not.ieee_is_finite(msink(1:2))).or.any(msink(1:2)<=0d0)) &
             local_error=1
        if(any(.not.ieee_is_finite(xsink(1:2,1:3))).or. &
             any(.not.ieee_is_finite(vsink(1:2,1:3)))) local_error=1
        if(any(xsink(1:2,1:3)<0d0).or.any(xsink(1:2,1:3)>=boxlen)) &
             local_error=1
     endif
  endif
  if(boxlen<=0d0.or..not.ieee_is_finite(boxlen))local_error=1

  ! Count the two canonical seed particles in owned grid lists.  levelp is a
  ! particle time-integration tag, not its actual AMR mesh level.
  sink_count=0
  if(.not.allocated(idp).or..not.allocated(ptypep).or. &
       .not.allocated(headp).or..not.allocated(nextp).or. &
       .not.allocated(numbp).or..not.allocated(headl).or. &
       .not.allocated(numbl).or..not.allocated(next))then
     local_error=1
  else if(npart<0.or.npart>size(idp).or.npart>size(ptypep).or. &
       size(nextp)<npartmax.or.size(headp)<ngridmax.or. &
       size(numbp)<ngridmax.or.size(headl,1)<myid.or. &
       size(headl,2)<nlevelmax.or.size(numbl,1)<myid.or. &
       size(numbl,2)<nlevelmax.or.size(next)<ngridmax)then
     local_error=1
  else
     do ilevel=1,nlevelmax
        igrid=headl(myid,ilevel)
        do jgrid=1,numbl(myid,ilevel)
           ipart=headp(igrid)
           do jpart=1,numbp(igrid)
              if(ipart<=0.or.ipart>npartmax)then
                 local_error=1
                 exit
              endif
              if(ptypep(ipart)==PTYPE_SINK)then
                 if(idp(ipart)>=0_i8b.or.idp(ipart)<-2_i8b)then
                    local_error=1
                 else
                    isink=int(-idp(ipart))
                    sink_count(isink)=sink_count(isink)+1
                 endif
              endif
              ipart=nextp(ipart)
           enddo
           igrid=next(igrid)
        enddo
     enddo
  endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(sink_count,global_sink_count,2,MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD,info)
  if(info/=MPI_SUCCESS)local_error=1
#else
  global_sink_count=sink_count
#endif
  if(any(global_sink_count/=1))local_error=1
  call sidm_binary_orbit_monitor_check(local_error, &
       'binary monitor requires a fresh fixed-mass 3-D PIC run with two ' // &
       'uniquely owned seed sinks, nsubcycle=1, future output clocks and ' // &
       'a dump-free step cap')

  binary_orbit_monitor_sink_ids=idsink(1:2)
  binary_orbit_monitor_sink_masses=msink(1:2)
  binary_orbit_monitor_sink_level=0
  if(icoarse_max<icoarse_min)call sidm_binary_orbit_monitor_check(1, &
       'invalid coarse-grid bounds for the finest force-stencil margin')
  dx_min=boxlen/dble(icoarse_max-icoarse_min+1)*0.5d0**nlevelmax
  ! The initialized sink cloud reaches 4 dx_min; keep another 2 dx_min
  ! between its edge and the periodic boundary for the force stencil.
  binary_orbit_monitor_boundary_margin=6d0*dx_min
  if(.not.ieee_is_finite(dx_min).or. &
       dx_min<=0d0.or.binary_orbit_monitor_boundary_margin>=0.5d0*boxlen) &
       call sidm_binary_orbit_monitor_check(1, &
       'could not establish an interior margin for the finest force stencil')
  if(any(xsink(1:2,1:3)<binary_orbit_monitor_boundary_margin).or. &
       any(xsink(1:2,1:3)>=boxlen-binary_orbit_monitor_boundary_margin)) &
       call sidm_binary_orbit_monitor_check(1, &
       'initial sink lies inside the finest force-stencil boundary margin')
  binary_orbit_monitor_has_initial=.false.
  binary_orbit_monitor_last_step=-1
  ios=0
  if(myid==1)then
     inquire(file='sidm_binary_orbit_state.terminal.csv',exist=terminal_exists)
     ios=0
     if(terminal_exists)ios=1
     if(ios==0)open(newunit=binary_orbit_monitor_unit, &
          file='sidm_binary_orbit_state.csv',status='new', &
          action='write',form='formatted',iostat=ios)
     if(ios==0)write(binary_orbit_monitor_unit,'(A)',iostat=ios) &
          'record_type,step,time_code,dt_since_previous_checkpoint,' // &
          'dt_completed_coarse_step,aexp,requested_stop_time,level,phase,sink_id,' // &
          'mass,x,y,z,vx,vy,vz,overshoot'
     if(ios==0)write(binary_orbit_monitor_unit,'(A,9(1X,ES25.16E3),6(1X,I0))',iostat=ios) &
          '# CONFIG',binary_orbit_stop_time,binary_orbit_sample_interval, &
          collisionless_courant_factor,minval(tout(1:noutput)), &
          maxval(tout(1:noutput)),minval(aout(1:noutput)), &
          maxval(aout(1:noutput)),boxlen,binary_orbit_monitor_boundary_margin, &
          levelmin,nlevelmax,nstepmax,foutput,fbackup,noutput
     if(ios==0)flush(binary_orbit_monitor_unit,iostat=ios)
  endif
  call sidm_binary_orbit_monitor_check(merge(0,1,ios==0), &
       'could not create or initialize sidm_binary_orbit_state.csv')
  binary_orbit_monitor_initialized=1
  binary_orbit_monitor_next_time=0d0
  binary_orbit_monitor_last_time=t
  binary_orbit_monitor_last_sync_time=t
  if(myid==1)write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') &
       'SIDM_BINARY_ORBIT_MONITOR enabled; requested_stop=', &
       binary_orbit_stop_time,'sample_interval=',binary_orbit_sample_interval
end subroutine sidm_binary_orbit_monitor_init

subroutine sidm_binary_orbit_monitor_sync_step(ilevel,finished)
  use amr_commons
  use hydro_parameters
  use pm_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,intent(in)::ilevel
  logical,intent(out)::finished
  real(dp)::checkpoint_interval,completed_step_dt,overshoot,monitor_tol
  integer::completed_step

  finished=.false.
  if(.not.binary_orbit_monitor.or.binary_orbit_monitor_initialized/=1.or. &
       binary_orbit_monitor_terminal/=0)return
  if(ilevel==levelmin.and.numbtot(1,nlevelmax)<=0) &
       call sidm_binary_orbit_monitor_abort( &
       'finest-level mesh disappeared before synchronized checkpoint', &
       'finest_mesh_guard')
  if(binary_orbit_monitor_sink_level==0)then
     if(ilevel<nlevelmax)return
     if(ilevel/=nlevelmax.or.nsink/=2.or.any(oksink_all(1:2)/=1d0)) &
          call sidm_binary_orbit_monitor_abort( &
          'both canonical sinks are not present at the common finest sync', &
          'finest_sink_level_sync')
     binary_orbit_monitor_sink_level=ilevel
  endif
  if(ilevel/=binary_orbit_monitor_sink_level)return
  if(nsink/=2.or.any(oksink_all(1:2)/=1d0)) &
       call sidm_binary_orbit_monitor_abort( &
       'both sinks were not synchronized at their common finest level', &
       'sink_level_sync')
  if(any(xsink(1:2,1:3)<binary_orbit_monitor_boundary_margin).or. &
       any(xsink(1:2,1:3)>=boxlen-binary_orbit_monitor_boundary_margin)) &
       call sidm_binary_orbit_monitor_abort( &
       'sink cloud or force stencil reached the periodic boundary margin', &
       'sink_boundary_guard')

  if(.not.binary_orbit_monitor_has_initial)then
     call sidm_binary_orbit_monitor_write_rows('initial',0d0,0d0,0d0,ilevel)
     binary_orbit_monitor_has_initial=.true.
     binary_orbit_monitor_last_step=nstep_coarse
     binary_orbit_monitor_last_time=t
     binary_orbit_monitor_last_sync_time=t
     binary_orbit_monitor_next_time=t+binary_orbit_sample_interval
     return
  endif

  completed_step=nstep_coarse
  if(completed_step<=binary_orbit_monitor_last_step.or. &
       .not.ieee_is_finite(t).or.t<=binary_orbit_monitor_last_time) &
       call sidm_binary_orbit_monitor_abort( &
       'synchronized step/time did not advance monotonically', &
       'sink_level_sync')

  checkpoint_interval=t-binary_orbit_monitor_last_time
  completed_step_dt=dtnew(ilevel)
  if(.not.ieee_is_finite(completed_step_dt).or.completed_step_dt<=0d0) &
       call sidm_binary_orbit_monitor_abort( &
       'completed coarse timestep is not finite and positive', &
       'sink_level_sync')
  monitor_tol=5d-12*max(1d-14,max(abs(checkpoint_interval), &
       abs(completed_step_dt),abs(binary_orbit_sample_interval)))
  if(.not.ieee_is_finite(binary_orbit_monitor_last_sync_time).or. &
       t<=binary_orbit_monitor_last_sync_time.or. &
       abs((t-binary_orbit_monitor_last_sync_time)-completed_step_dt)>monitor_tol.or. &
       checkpoint_interval<=0d0.or. &
       completed_step_dt>binary_orbit_sample_interval+monitor_tol.or. &
       checkpoint_interval>2d0*binary_orbit_sample_interval+monitor_tol) &
       call sidm_binary_orbit_monitor_abort( &
       'completed timestep or checkpoint cadence exceeded its validated bound', &
       'sink_level_sync')

  if(t>=binary_orbit_stop_time)then
     overshoot=t-binary_orbit_stop_time
     if(overshoot<0d0.or.overshoot>completed_step_dt+monitor_tol.or. &
          overshoot>binary_orbit_sample_interval+monitor_tol) &
          call sidm_binary_orbit_monitor_abort( &
          'physical-stop overshoot exceeded one completed timestep', &
          'sink_level_sync')
     call sidm_binary_orbit_monitor_write_rows('terminal', &
          checkpoint_interval,completed_step_dt,overshoot,ilevel)
     call sidm_binary_orbit_monitor_finish(.true.,'physical_time', &
          'centered_sink_level_sync')
     finished=.true.
     return
  endif

  monitor_tol=5d-12*max(1d-14,max(abs(t), &
       abs(binary_orbit_monitor_next_time)))
  if(t+monitor_tol>=binary_orbit_monitor_next_time)then
     call sidm_binary_orbit_monitor_write_rows('sample', &
          checkpoint_interval,completed_step_dt,0d0,ilevel)
     do while(binary_orbit_monitor_next_time<=t+monitor_tol)
        binary_orbit_monitor_next_time= &
             binary_orbit_monitor_next_time+binary_orbit_sample_interval
     enddo
  endif
  if(completed_step>=nstepmax) &
       call sidm_binary_orbit_monitor_abort( &
       'safety_step_cap_reached_before_physical_stop','sink_level_sync')
  binary_orbit_monitor_last_sync_time=t
end subroutine sidm_binary_orbit_monitor_sync_step

subroutine sidm_binary_orbit_monitor_write_rows(record_type,checkpoint_interval, &
     completed_step_dt,overshoot,ilevel)
  use amr_commons
  use hydro_parameters
  use pm_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=*),intent(in)::record_type
  real(dp),intent(in)::checkpoint_interval,completed_step_dt,overshoot
  integer,intent(in)::ilevel
  integer::isink,ios,local_error
  integer::stored_index
  integer::completed_step

  local_error=0
  if(nsink/=2.or..not.allocated(idsink).or..not.allocated(msink).or. &
       .not.allocated(xsink).or..not.allocated(vsink)) local_error=1
  if(local_error==0)then
     if(idsink(1)==idsink(2))local_error=1
     if(any(.not.ieee_is_finite(msink(1:2))).or. &
          any(.not.ieee_is_finite(xsink(1:2,1:3))).or. &
          any(.not.ieee_is_finite(vsink(1:2,1:3)))) local_error=1
     if(any(xsink(1:2,1:3)<binary_orbit_monitor_boundary_margin).or. &
          any(xsink(1:2,1:3)>= &
          boxlen-binary_orbit_monitor_boundary_margin)) &
          local_error=1
     do isink=1,2
        stored_index=0
        if(idsink(isink)==binary_orbit_monitor_sink_ids(1))stored_index=1
        if(idsink(isink)==binary_orbit_monitor_sink_ids(2))stored_index=2
        if(stored_index==0)then
           local_error=1
        else if(msink(isink)<=0d0.or. &
             abs(msink(isink)-binary_orbit_monitor_sink_masses(stored_index))> &
             5d-12*abs(binary_orbit_monitor_sink_masses(stored_index)))then
           local_error=1
        endif
     enddo
  endif
  call sidm_binary_orbit_monitor_check(local_error, &
       'sink identities, masses, or canonical phase became invalid during orbit monitor')

  ios=0
  completed_step=nstep_coarse
  if(myid==1)then
     if(binary_orbit_monitor_unit==0)then
        ios=1
     else
        do isink=1,2
           write(binary_orbit_monitor_unit, &
                '(A,",",I0,5(",",ES25.16E3),",",I0,",",A,",",I0,8(",",ES25.16E3))', &
                iostat=ios) trim(record_type),completed_step,t,checkpoint_interval, &
                completed_step_dt,aexp,binary_orbit_stop_time, &
                ilevel,'centered_sink_level_sync',idsink(isink),msink(isink), &
                xsink(isink,1),xsink(isink,2),xsink(isink,3), &
                vsink(isink,1),vsink(isink,2),vsink(isink,3),overshoot
           if(ios/=0)exit
        enddo
        if(ios==0)flush(binary_orbit_monitor_unit,iostat=ios)
     endif
  endif
  call sidm_binary_orbit_monitor_check(merge(0,1,ios==0), &
       'failed writing canonical sink rows to sidm_binary_orbit_state.csv')
  binary_orbit_monitor_last_step=completed_step
  binary_orbit_monitor_last_time=t
end subroutine sidm_binary_orbit_monitor_write_rows

subroutine sidm_binary_orbit_monitor_abort(reason,phase)
  use hydro_parameters
  implicit none
  character(len=*),intent(in)::reason,phase
  call sidm_binary_orbit_monitor_finish(.false.,reason,phase)
  call sidm_binary_orbit_monitor_check(1,reason)
end subroutine sidm_binary_orbit_monitor_abort

subroutine sidm_binary_orbit_monitor_finish(success,reason,phase)
  use amr_commons
  use hydro_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical,intent(in)::success
  character(len=*),intent(in)::reason,phase
  integer::ios,local_error,terminal_unit,completed_step
  character(len=8)::terminal_state

  if(binary_orbit_monitor_initialized/=1.or. &
       binary_orbit_monitor_terminal/=0)return
  terminal_state='failed'
  if(success)terminal_state='complete'
  local_error=0
  completed_step=nstep_coarse
  if(myid==1)then
     if(binary_orbit_monitor_unit/=0)then
        flush(binary_orbit_monitor_unit,iostat=ios)
        if(ios==0)close(binary_orbit_monitor_unit,iostat=ios)
        if(ios/=0)local_error=1
        binary_orbit_monitor_unit=0
     else
        local_error=1
     endif
     if(local_error==0)then
        open(newunit=terminal_unit,file='sidm_binary_orbit_state.terminal.csv', &
             status='new',action='write',form='formatted',iostat=ios)
        if(ios==0)write(terminal_unit,'(A)',iostat=ios) &
             'terminal_state,reason,step_at_record,time_code,requested_stop_time,overshoot,phase'
        if(ios==0)write(terminal_unit,'(A,",",A,",",I0,3(",",ES25.16E3),",",A)', &
             iostat=ios) trim(terminal_state),trim(reason),completed_step,t, &
             binary_orbit_stop_time,t-binary_orbit_stop_time,trim(phase)
        if(ios==0)flush(terminal_unit,iostat=ios)
        if(ios==0)close(terminal_unit,iostat=ios)
        if(ios/=0)local_error=1
     endif
  endif
  call sidm_binary_orbit_monitor_check(local_error, &
       'failed to close the orbit ledger or write its terminal commit')
  binary_orbit_monitor_terminal=1
end subroutine sidm_binary_orbit_monitor_finish

subroutine sidm_binary_orbit_monitor_check(local_error,message)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::local_error
  character(len=*),intent(in)::message
  integer::global_error,info

  global_error=local_error
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(local_error,global_error,1,MPI_INTEGER,MPI_MAX, &
       MPI_COMM_WORLD,info)
  if(info/=MPI_SUCCESS)global_error=1
#endif
  if(global_error/=0)then
     if(myid==1)write(*,'(A,1X,A)') &
          'ERROR: SIDM binary orbit monitor:',trim(message)
#ifndef WITHOUTMPI
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
     stop 1
#endif
  endif
end subroutine sidm_binary_orbit_monitor_check

subroutine check_jobcontrol
  use amr_commons
  use amr_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::mpi_err
#endif
  integer::istep,iaction,ios,action_todo,lb_mode_todo,lb_force_todo
  integer::lb_mode_old,ipos
  logical::fexist
  character(len=256)::line
  character(len=8)::old_mode,new_mode

  if(len_trim(jobcontrolfile)==0) return

  action_todo = 0  ! default: do nothing
  lb_mode_todo = -1 ! 0=memory, 1=work, 2=timed work
  lb_force_todo = 0

  if(myid==1) then
    inquire(file=trim(jobcontrolfile), exist=fexist)
    if(fexist) then
      open(unit=99, file=trim(jobcontrolfile), form='formatted', status='old', iostat=ios)
      if(ios==0) then
        do  ! read all lines; blank lines, comments, and malformed lines are skipped
          read(99, '(A)', iostat=ios) line
          if(ios < 0) exit
          if(ios > 0) cycle
          ipos=index(line,'#')
          if(ipos==1)then
            line=''
          else if(ipos>1)then
            line=line(:ipos-1)
          endif
          if(len_trim(line)==0) cycle
          read(line, *, iostat=ios) istep, iaction
          if(ios /= 0) then
            write(*,'(A,A)') ' Job control: ignoring malformed line: ',trim(line)
            cycle
          end if
          if(istep == 0 .or. istep == nstep_coarse) then
            select case(iaction)
            case(0)
              continue
            case(-1)
              action_todo = -1  ! stop always wins
            case(1)
              if(action_todo /= -1) action_todo = 1
            case(2:4)
              lb_mode_todo=iaction-2
            case(5)
              ! A wildcard force-remap would repeat forever because the file is
              ! intentionally persistent.  Require an explicit coarse step.
              if(istep>0) then
                lb_force_todo=1
              else
                write(*,'(A)') ' Job control: ignoring wildcard (step 0) force-remap action'
              end if
            case default
              write(*,'(A,I0)') ' Job control: ignoring unknown action ',iaction
            end select
          endif
        end do
        close(99)
      endif
    endif
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(action_todo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
  call MPI_BCAST(lb_mode_todo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
  call MPI_BCAST(lb_force_todo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#endif

  if(memory_balance)then
    lb_mode_old=0
    old_mode='memory'
  else if(time_balance_alpha>0d0)then
    lb_mode_old=2
    old_mode='timed'
  else
    lb_mode_old=1
    old_mode='work'
  endif

  if(lb_mode_todo>=0 .and. lb_mode_todo/=lb_mode_old)then
    select case(lb_mode_todo)
    case(0)
      memory_balance=.true.
      time_balance_alpha=0d0
      new_mode='memory'
    case(1)
      memory_balance=.false.
      cost_weighting=.true.
      time_balance_alpha=0d0
      new_mode='work'
    case(2)
      memory_balance=.false.
      cost_weighting=.true.
      time_balance_alpha=0.3d0
      new_mode='timed'
    end select
    ! Samples gathered under another cost model are not comparable.  Reset
    ! only the work/imbalance history; the measured remap cost remains useful.
    level_time_loc=0d0
    level_ncells_loc=0_8
    level_cell_time_ema=0d0
    level_mesh_scale_ema=1d0
    level_rank_scale_ema=1d0
    lb_step_time_ema=0d0
    lb_imbalance_ema=0d0
    lb_imbalance_ema_valid=.false.
    lb_force_remap=.true.
    if(myid==1) write(*,'(A,A,A,A,A,I0)') &
         ' Job control: load-balance mode ',trim(old_mode),' -> ', &
         trim(new_mode),' at step ',nstep_coarse
  endif
  if(lb_force_todo==1)then
    lb_force_remap=.true.
    if(myid==1) write(*,'(A,I0)') &
         ' Job control: forced load balance requested at step ',nstep_coarse
  endif

  if(action_todo == 1) then
    output_now = .true.
    if(myid==1) write(*,'(A,I7)') ' Job control: extra output requested at step', nstep_coarse
  else if(action_todo == -1) then
    output_now = .true.
    nstepmax = nstep_coarse  ! triggers clean_stop in update_time
    if(myid==1) write(*,'(A,I7)') ' Job control: stop requested at step', nstep_coarse
  endif
end subroutine check_jobcontrol
