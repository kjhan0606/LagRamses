subroutine newdt_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY: gravity_type
  use fdm_commons, ONLY: hbar_code
#ifdef RT
  use rt_parameters, ONLY: rt_advect, rt_nsubcycle
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 3 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp 
  ! 4- maximum step time for ATON
  ! This routine also compute the particle kinetic energy.
  !-----------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,nx_loc
  integer::npart1,ip,info,ilev
#ifndef _OPENMP
  integer,dimension(1:nvector),save::ind_part
#else
  integer, dimension(:), allocatable:: nparticles, ptrhead
  integer mythread, nthreads,nwork,icount,jcount,npart3,subnump
  common /openmpthreads/ mythread, nthreads
!$omp threadprivate(/openmpthreads/)
#endif

  real(kind=8)::dt_loc,dt_all,ekin_loc,ekin_all,dt_acc_min
  real(kind=8)::idt_loc,iekin_loc
  real(kind=8)::jdt_loc,jekin_loc
  real(dp)::tff,fourpi,threepi2
  real(dp)::aton_time_step,dt_aton,dt_rt
  real(dp)::dx_min,dx,scale,dt_fact,limiting_dt_fact
  real(dp)::dx_fdm,k2_nyq,k2_eff,k2_dB,dt_acc,dt_nyq,dt_adv
  real(dp)::dtheta_max,c1_max,eps_c1,dt_qp
  logical::highest_level

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifdef _OPENMP
!$omp parallel
  mythread = omp_get_thread_num()
  nthreads = omp_get_num_threads()
!$omp end parallel
  allocate(ptrhead(0:nthreads-1), nparticles(0:nthreads-1))
#endif


  threepi2=3.0d0*ACOS(-1.0d0)**2

  ! Save old time step
  dtold(ilevel)=dtnew(ilevel)

  ! Maximum time step
  dtnew(ilevel)=boxlen/smallc
  if(poisson.and.gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     if(cosmo)fourpi=1.5d0*omega_m*aexp
     tff=sqrt(threepi2/8./fourpi/rho_max(ilevel))
     dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*tff)
  end if
  if(cosmo)then
     dtnew(ilevel)=MIN(dtnew(ilevel),0.1/hexp)
  end if
  ! HJM CFL debug: after gravity+cosmo constraints
  if(myid==1 .and. nstep_coarse_old < 36 .and. ilevel==levelmin &
     .and. fdm_use_hjm) then
     write(*,'(" HJM_CFL[grav]: dt=",1PE12.5," tff=",1PE12.5," hexp=",1PE12.5)') &
          dtnew(ilevel), courant_factor*tff, 0.1d0/hexp
  end if
  ! FDM kinetic timestep. The classic CFL dt < fdm_courant*dx^2*a^2/(6*hbar)
  ! is the worst case: it equals fdm_courant*a^2/(hbar*k2) with k2 = 6/dx^2,
  ! the 7-point stencil's max eigenvalue (~Nyquist). Two knobs relax it:
  !   fdm_kinetic=1 (Crank-Nicolson): fine levels are unconditionally STABLE,
  !       so they need no dx^2 limit (only the base, spectral, sets accuracy).
  !   fdm_hybrid=.true. (Madelung/de Broglie): replace the Nyquist k2 by the
  !       actual max |grad theta|^2 in the field (capped at Nyquist). Smooth
  !       (fluid-like) regions have small k -> large dt; dt tightens adaptively
  !       only where quantum interference (high phase gradient) appears.
  ! Default (fdm_kinetic=0, fdm_hybrid=.false.) reproduces the original CFL
  ! at every level, byte-for-byte.
  if(use_fdm .and. hbar_code > 0.0d0)then
     dx_fdm  = 0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1)
     k2_nyq  = 6.0d0/dx_fdm**2

     if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
        ! HJM fluid level: advection CFL from Madelung velocity
        call fdm_vmax_level(ilevel, dtheta_max)
        dt_adv = huge(1.0d0)
        if(dtheta_max > 0.0d0) then
           dt_adv = fdm_courant * dx_fdm * aexp**2 &
                / (hbar_code * dtheta_max)
           dtnew(ilevel) = MIN(dtnew(ilevel), dt_adv)
        end if
        if(myid==1 .and. nstep_coarse_old < 36 .and. ilevel==levelmin) then
           write(*,'(" HJM_CFL[fdm]: dtheta=",1PE10.3," dt_adv=",1PE10.3,&
                &" dt=",1PE10.3)') dtheta_max, dt_adv, dtnew(ilevel)
        end if
     else
        ! Wave levels: existing Schrodinger CFL
        if(fdm_hybrid)then
           call fdm_kdb_max_level(ilevel, k2_dB)
           k2_eff = max(1.0d-30, min(k2_dB, k2_nyq))
        else
           k2_eff = k2_nyq
        end if
        dt_acc = fdm_courant * aexp**2 / (hbar_code * k2_eff)
        dt_nyq = fdm_courant * aexp**2 / (hbar_code * k2_nyq)
        if(ilevel==levelmin)then
           dtnew(ilevel)=MIN(dtnew(ilevel), dt_acc)
        else
           if(fdm_kinetic==1)then
              if(fdm_hybrid) dtnew(ilevel)=MIN(dtnew(ilevel), dt_acc)
           else
              dtnew(ilevel)=MIN(dtnew(ilevel), dt_nyq)
           end if
        end if
     end if
  end if

#ifdef ATON
  ! Maximum time step for ATON
  if(aton)then
     dt_aton = aton_time_step()
     if(dt_aton>0d0)then
        dtnew(ilevel)=MIN(dtnew(ilevel),dt_aton)
     end if
  end if
#endif

#ifdef RT
  ! Maximum time step for radiative transfer
  if(rt_advect)then
     call get_rt_courant_coarse(dt_rt)
     dtnew(ilevel) = 0.99999 * &
          MIN(dtnew(ilevel), dt_rt/2.0**(ilevel-levelmin) * rt_nsubcycle)
     if(static) RETURN
  endif
#endif

  if(pic) then

     dt_all=dtnew(ilevel); dt_loc=dt_all
     ekin_all=0.0; ekin_loc=0.0
     
!jhshin1
#ifdef _OPENMP
     ! Compute maximum time step on active region
     if(numbl(myid,ilevel)>0)then
        call pthreadLinkedList(headl(myid,ilevel),numbl(myid,ilevel),nthreads,nparticles,ptrhead,next)
        iekin_loc = 0
        idt_loc = dt_all
!$omp parallel private(subnump,igrid,jdt_loc, jekin_loc) reduction(+:iekin_loc), reduction(MIN: idt_loc)
        subnump = nparticles(mythread)
        igrid = ptrhead(mythread)
	    call sub_newdt_fine(ilevel, igrid,subnump, jdt_loc, jekin_loc, dtnew(ilevel))
        iekin_loc = iekin_loc + jekin_loc
        idt_loc = MIN(idt_loc, jdt_loc)
!$omp end parallel
        ekin_loc = iekin_loc
        dt_loc = idt_loc
     end if
#else 
     ! Compute maximum time step on active region
     if(numbl(myid,ilevel)>0)then
        ! Loop over grids
        ip=0
        igrid=headl(myid,ilevel)
        do jgrid=1,numbl(myid,ilevel)
           npart1=numbp(igrid)   ! Number of particles in the grid
           if(npart1>0)then
              ! Loop over particles
              ipart=headp(igrid)
              do jpart=1,npart1
                 ip=ip+1
                 ind_part(ip)=ipart
                 if(ip==nvector)then
                    call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
                    ip=0
                 end if
                 ipart=nextp(ipart)    ! Go to next particle
              end do
              ! End loop over particles
           end if
           igrid=next(igrid)   ! Go to next grid
        end do
        ! End loop over grids
        if(ip>0)call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
     end if
#endif
!jhshin2

     ! Minimize time step over all cpus
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
          & MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     dt_all=dt_loc
     ekin_all=ekin_loc
#endif
     ekin_tot=ekin_tot+ekin_all
     dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)
     if(myid==1 .and. nstep_coarse_old < 36 .and. ilevel==levelmin &
        .and. fdm_use_hjm) then
        write(*,'(" HJM_CFL[pic]: dt_all=",1PE12.5," dt=",1PE12.5)') &
             dt_all, dtnew(ilevel)
     end if

  end if
!jhshin1
#ifdef _OPENMP
  deallocate(ptrhead, nparticles)
#endif
!jhshin2

  if(hydro)call courant_fine(ilevel)
  if(myid==1 .and. nstep_coarse_old < 36 .and. ilevel==levelmin &
     .and. fdm_use_hjm) then
     write(*,'(" HJM_CFL[hydro]: dt=",1PE12.5)') dtnew(ilevel)
  end if

  ! SIDM timestep constraint: keep P_scatter < sidm_courant
  if(sidm .and. sidm_Pmax(ilevel) > 0.0d0) then
     dtnew(ilevel) = MIN(dtnew(ilevel), &
          sidm_courant / sidm_Pmax(ilevel) * dtold(ilevel))
  end if

  ! SGS dt diagnostic: print per-level dt for first 20 coarse steps
  if(myid==1 .and. nstep_coarse_old < 80) then
     write(*,'(" SGS_DT: level=",I2," dt=",1PE10.3)') ilevel, dtnew(ilevel)
  end if

  ! HJM CFL debug: print breakdown for first 5 steps at levelmin
  if(myid==1 .and. nstep_coarse_old < 36 .and. ilevel==levelmin &
     .and. fdm_use_hjm) then
     write(*,'(" HJM_CFL_DBG: dtold=",1PE10.3," dt_final=",1PE10.3)') &
          dtold(ilevel), dtnew(ilevel)
  end if

#ifdef _OPENMP
111 format('  +Entering newdt_fine for level ',I2)
#else
111 format('   Entering newdt_fine for level ',I2)
#endif

end subroutine newdt_fine
!#####################################################################
!#####################################################################
subroutine sub_newdt_fine(ilevel, igrid, subnump, dt_loc, ekin_loc,dt_all)
  use pm_commons
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY: gravity_type
#ifdef RT
  use rt_parameters, ONLY: rt_advect, rt_nsubcycle
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 3 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp 
  ! 4- maximum step time for ATON
  ! This routine also compute the particle kinetic energy.
  !-----------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,nx_loc
  integer::npart1,ip,info,ilev
  integer,dimension(1:nvector)::ind_part
  integer subnump

  real(kind=8)::dt_loc,dt_all,ekin_loc,ekin_all,dt_acc_min
  real(kind=8)::idt_loc,iekin_loc
  real(dp)::tff,fourpi,threepi2
  real(dp)::aton_time_step,dt_aton,dt_rt
  real(dp)::dx_min,dx,scale,dt_fact,limiting_dt_fact
  logical::highest_level
  ip=0
  dt_loc = dt_all
  ekin_loc = 0
  do jgrid=1,subnump
     npart1=numbp(igrid)   ! Number of particles in the grid
     if(npart1>0)then
        ! Loop over particles
        ipart=headp(igrid)
        do jpart=1,npart1
           ip=ip+1
           ind_part(ip)=ipart
           if(ip==nvector)then
              call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
              ip=0
           end if
           ipart=nextp(ipart)    ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0) then
     call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
  endif

end subroutine sub_newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt2(ind_part,dt_loc,ekin_loc,nn,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(kind=8)::dt_loc,ekin_loc
  integer::nn,ilevel
  integer,dimension(1:nvector)::ind_part

  integer::i,idim,nx_loc
  real(dp)::dx,dx_loc,scale,dtpart
  real(dp),dimension(1:nvector)::v2

  ! Compute time step
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  v2(1:nn)=0.0D0
  do idim=1,ndim
     do i=1,nn
        v2(i)=MAX(v2(i),vp(ind_part(i),idim)**2)
!        v2(i)=v2(i)+vp(ind_part(i),idim)**2
     end do
  end do
  do i=1,nn
     if(v2(i)>0.0D0)then
        dtpart=courant_factor*dx_loc/sqrt(v2(i))
        dt_loc=MIN(dt_loc,dtpart)
     end if
  end do

  ! Compute kinetic energy
  do idim=1,ndim
     do i=1,nn
        ekin_loc=ekin_loc+0.5D0*mp(ind_part(i))*vp(ind_part(i),idim)**2
     end do
  end do
    
end subroutine newdt2




