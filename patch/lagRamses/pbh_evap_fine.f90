!=======================================================================
! pbh_evap_fine: evaporating-PBH dark matter, per-level update.
!
! Called once per level step from amr_step, in the unew window right
! after thermal_feedback (after set_unew, before godunov_fine), so the
! existing make_virtual_reverse_dp(unew) collects cross-domain deposits.
!
! Per step (paper appendix A):
!   heat : dE = mp*efac/vol  ->  unew(:,ndim+2) (+enew if pressure_fix)
!          with efac = f_PBH*boost*dQtilde / (w0*scale_v**2)
!   mass : mp = mp * w(a1)/w(a0)          (exact mixed-mass table ratio)
! Deposit uses the step-start mp (m_initial = mp/w0), then the mass is
! updated, so the bookkeeping matches the cumulative table exactly.
!
! Thread safety WITHOUT atomics (owner-writes + deferred remote buffer):
! after the per-step particle-tree sort a particle's NGP containing cell
! belongs to its own grid, and one grid is processed by one thread, so
! the common-case deposit is single-writer. The rare non-owned targets
! (parent-cell fallback, boundary-drift neighbours) are queued in
! thread-private buffers and applied serially after the parallel loop,
! sorted by (cell, particle id). The injection tally is accumulated per
! grid chunk and summed in chunk order. The result is therefore bitwise
! independent of the OpenMP thread count and schedule.
!=======================================================================
subroutine pbh_evap_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use pbh_commons
!$ use omp_lib
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::ichunk,ncache,info,nthr,nch,k,i,j,e,ntot
  real(dp)::dt_phys,ratio,dQ,dQcr,wbeg,efac,efaccr,conv
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::einj_all,ecr_all,ecr_mesh,ecr_mesh_all
  integer(kind=8)::nfall_all
  logical::do_heat,do_cr
  integer,parameter::NGRID_CHUNK=64
  real(dp),allocatable::einj_ch(:),ecr_ch(:)
  integer(kind=8),allocatable::nfall_ch(:)
  ! gathered remote deposits (sorted, applied serially)
  integer,allocatable::gi(:)
  integer(kind=8),allocatable::gp(:)
  real(dp),allocatable::gd(:),gc(:)
  integer::tmpi
  integer(kind=8)::tmpp
  real(dp)::tmpd,tmpc

  if(.not.use_pbh)return
  if(.not.pic)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  call pbh_lazy_init(aexp,aexp_ini,nrestart,myid,cosmo)

  ! step-global factors (aexp is the END of this level step here,
  ! update_time having already run)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  dt_phys=dtnew(ilevel)*scale_t
  call pbh_step(ilevel,nlevelmax,aexp,dt_phys,ratio,dQ,wbeg,dQcr)
  do_heat=hydro .and. trim(pbh_energy_sink)=='local_heat' .and. dQ>0.0d0
  do_cr=hydro .and. pbh_cr_ivar>ndim+2 .and. pbh_cr_ivar<=nvar .and. dQcr>0.0d0
  ! f_PBH is the PBH fraction of DARK MATTER. In hydro cosmological runs the
  ! grafic IC already normalises the particle mass to the dark-matter share
  ! (init_part: mp = 0.5**(3*ilevel)*(1-omega_b/omega_m)), so mp IS the DM
  ! mass and f_PBH*mp*dQ is the PBH energy with no further share factor.
  ! Adding (1-Omega_b/Omega_m) here double-counts the share and under-heats
  ! by 15 percent (verified against f*M_dm*dQtilde on the P1 uniform box).
  efac=pbh_fraction*dQ/(wbeg*scale_v**2)   ! per unit current code mass
  efaccr=pbh_fraction*dQcr/(wbeg*scale_v**2)
  conv=scale_v**2*scale_d*scale_l**3        ! code energy -> erg

  ! thread-private remote-deposit buffers (one per possible thread)
  nthr=1
!$ nthr=omp_get_max_threads()
  if(.not.allocated(pbh_rbuf))allocate(pbh_rbuf(0:nthr-1))

  ncache=active(ilevel)%ngrid
  nch=max(1,(ncache+NGRID_CHUNK-1)/NGRID_CHUNK)
  allocate(einj_ch(nch),ecr_ch(nch),nfall_ch(nch))
  einj_ch=0.0d0
  ecr_ch=0.0d0
  nfall_ch=0

!$omp parallel do private(ichunk,k) schedule(dynamic)
  do ichunk=1,ncache,NGRID_CHUNK
     k=(ichunk-1)/NGRID_CHUNK+1
     call sub_pbh_evap_fine(ilevel,ichunk,min(NGRID_CHUNK,ncache-ichunk+1), &
          & ratio,efac,efaccr,do_heat,do_cr,einj_ch(k),ecr_ch(k),nfall_ch(k))
  end do
!$omp end parallel do

  ! apply buffered remote deposits, sorted by (cell, particle id) so the
  ! per-cell addition order never depends on threads or schedule
  ntot=0
  do i=0,ubound(pbh_rbuf,1)
     ntot=ntot+pbh_rbuf(i)%n
  end do
  if(ntot>0)then
     allocate(gi(ntot),gp(ntot),gd(ntot),gc(ntot))
     e=0
     do i=0,ubound(pbh_rbuf,1)
        do j=1,pbh_rbuf(i)%n
           e=e+1
           gi(e)=pbh_rbuf(i)%icell(j)
           gp(e)=pbh_rbuf(i)%pid(j)
           gd(e)=pbh_rbuf(i)%de(j)
           gc(e)=pbh_rbuf(i)%dcr(j)
        end do
        pbh_rbuf(i)%n=0
     end do
     ! insertion sort (remote deposits are rare, ntot is tiny)
     do i=2,ntot
        tmpi=gi(i); tmpp=gp(i); tmpd=gd(i); tmpc=gc(i)
        j=i-1
        do while(j>=1)
           if(gi(j)>tmpi .or. (gi(j)==tmpi .and. gp(j)>tmpp))then
              gi(j+1)=gi(j); gp(j+1)=gp(j); gd(j+1)=gd(j); gc(j+1)=gc(j)
              j=j-1
           else
              exit
           end if
        end do
        gi(j+1)=tmpi; gp(j+1)=tmpp; gd(j+1)=tmpd; gc(j+1)=tmpc
     end do
     do e=1,ntot
        if(gd(e)/=0.0d0)then
           unew(gi(e),ndim+2)=unew(gi(e),ndim+2)+gd(e)
           if(pressure_fix)then
              enew(gi(e))=enew(gi(e))+gd(e)
           end if
        end if
        if(gc(e)/=0.0d0)then
           unew(gi(e),pbh_cr_ivar)=unew(gi(e),pbh_cr_ivar)+gc(e)
        end if
     end do
     deallocate(gi,gp,gd,gc)
  end if

  ! tallies: chunk partials summed in fixed chunk order (deterministic)
  do k=1,nch
     pbh_einj_loc=pbh_einj_loc+einj_ch(k)*conv
     pbh_ecr_loc=pbh_ecr_loc+ecr_ch(k)*conv
     pbh_nfallback_loc=pbh_nfallback_loc+nfall_ch(k)
  end do
  deallocate(einj_ch,ecr_ch,nfall_ch)

  ! coarse-step diagnostics (all ranks reach this when numbtot>0)
  if(ilevel==levelmin)then
     ! in-situ mesh content of the CR reservoir (code units, this level)
     ecr_mesh=0.0d0
     if(do_cr)then
        do i=1,active(ilevel)%ngrid
           do j=1,twotondim
              e=ncoarse+(j-1)*ngridmax+active(ilevel)%igrid(i)
              ecr_mesh=ecr_mesh+uold(e,pbh_cr_ivar)
           end do
        end do
        ecr_mesh=ecr_mesh*(0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1))**ndim
     end if
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(pbh_einj_loc,einj_all,1,MPI_DOUBLE_PRECISION, &
          & MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(pbh_ecr_loc,ecr_all,1,MPI_DOUBLE_PRECISION, &
          & MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(pbh_nfallback_loc,nfall_all,1,MPI_INTEGER8, &
          & MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ecr_mesh,ecr_mesh_all,1,MPI_DOUBLE_PRECISION, &
          & MPI_SUM,MPI_COMM_WORLD,info)
#else
     einj_all=pbh_einj_loc
     ecr_all=pbh_ecr_loc
     nfall_all=pbh_nfallback_loc
     ecr_mesh_all=ecr_mesh
#endif
     pbh_einj_tot=einj_all
     if(myid==1)write(*,112)aexp,wbeg*ratio,rho_tot,einj_all,nfall_all,eint_tot, &
          & ecr_all,ecr_mesh_all
  end if

111 format('   Entering pbh_evap_fine for level ',I2)
112 format(' PBHDIAG aexp=',ES14.7,' w=',ES16.9,' rho_tot=',ES16.9, &
         & ' einj[erg]=',ES14.7,' nfallback=',I12,' eint=',ES16.9, &
         & ' ecr[erg]=',ES14.7,' ecrmesh=',ES14.7)
end subroutine pbh_evap_fine

!=======================================================================
subroutine sub_pbh_evap_fine(ilevel,istart,nchunk,ratio,efac,efaccr, &
     & do_heat,do_cr,einj_c,ecr_c,nfall_c)
  ! Gather DM particles of one chunk of active grids into nvector
  ! batches (classic feedback pattern) and hand them to pbh_dump.
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel,istart,nchunk
  real(dp)::ratio,efac,efaccr,einj_c,ecr_c
  integer(kind=8)::nfall_c
  logical::do_heat,do_cr
  integer::i,jgrid,igrid,jpart,ipart,next_part,npart1,npart2,ig,ip
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part

  ig=0
  ip=0
  do jgrid=1,nchunk
     igrid=active(ilevel)%igrid(istart+jgrid-1)
     npart1=numbp(igrid)
     npart2=0
     ! count DM particles
     if(npart1>0)then
        ipart=headp(igrid)
        do jpart=1,npart1
           next_part=nextp(ipart)
           if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
                & ptypep(ipart)/=PTYPE_SINK)then
              npart2=npart2+1
           end if
           ipart=next_part
        end do
     end if
     ! gather DM particles
     if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        do jpart=1,npart1
           next_part=nextp(ipart)
           if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
                & ptypep(ipart)/=PTYPE_SINK)then
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
           end if
           if(ip==nvector)then
              call pbh_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel, &
                   & ratio,efac,efaccr,do_heat,do_cr,einj_c,ecr_c,nfall_c)
              ip=0
              ig=0
           end if
           ipart=next_part
        end do
     end if
  end do
  if(ip>0)call pbh_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel, &
       & ratio,efac,efaccr,do_heat,do_cr,einj_c,ecr_c,nfall_c)
end subroutine sub_pbh_evap_fine

!=======================================================================
subroutine pbh_dump(ind_grid,ind_part,ind_grid_part,ng,np,ilevel, &
     & ratio,efac,efaccr,do_heat,do_cr,einj_c,ecr_c,nfall_c)
  ! NGP energy deposit + exact mass update for one particle batch.
  ! Owned targets (cell inside the particle's own grid) are written
  ! directly; everything else goes to the thread-private remote buffer.
  use amr_commons
  use pm_commons
  use hydro_commons
  use pbh_commons, only: pbh_rbuf,pbh_rbuf_push,pbh_cr_ivar
!$ use omp_lib
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::ratio,efac,efaccr,einj_c,ecr_c
  integer(kind=8)::nfall_c
  logical::do_heat,do_cr
  integer::i,j,idim,nx_loc,tid
  real(dp)::dx,dx_loc,scale,vol_loc,dE,dEcr,wv
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::x0
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector)::ok,own
  real(dp),dimension(1:nvector,1:ndim)::x
  integer ,dimension(1:nvector,1:ndim)::id,igd,icd
  integer ,dimension(1:nvector)::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

#if NDIM==3
  if(do_heat.or.do_cr)then
     tid=0
!$   tid=omp_get_thread_num()
     ! Lower left corner of 3x3x3 grid-cube
     do idim=1,ndim
        do i=1,ng
           x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
        end do
     end do
     ! Gather 27 neighboring father cells (should be present anytime !)
     do i=1,ng
        ind_cell(i)=father(ind_grid(i))
     end do
     call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)
     ! Rescale position at level ilevel
     do idim=1,ndim
        do j=1,np
           x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
           x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
           x(j,idim)=x(j,idim)/dx
        end do
     end do
     ! NGP at level ilevel
     do idim=1,ndim
        do j=1,np
           id(j,idim)=x(j,idim)
        end do
     end do
     ! Compute parent grids
     do idim=1,ndim
        do j=1,np
           igd(j,idim)=id(j,idim)/2
        end do
     end do
     do j=1,np
        kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     end do
     do j=1,np
        igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
     end do
     ! Particles whose child grid exists deposit at ilevel
     do j=1,np
        ok(j)=igrid(j)>0
     end do
     ! Owned = the target cell lies in the particle's own grid, which this
     ! thread is processing exclusively (single-writer, no protection needed)
     do j=1,np
        own(j)=ok(j).and.(igrid(j)==ind_grid(ind_grid_part(j)))
     end do
     ! Compute parent cell position and address
     do idim=1,ndim
        do j=1,np
           if(ok(j))icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end do
     end do
     do j=1,np
        if(ok(j))then
           icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
           indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
        else
           ! fallback: deposit into the coarser parent cell (ilevel-1)
           ! so that no particle's heat is ever silently dropped
           indp(j)=nbors_father_cells(ind_grid_part(j),kg(j))
        end if
     end do

     ! Deposit thermal energy (total energy index ndim+2; the auxiliary
     ! internal-energy field enew must receive the same increment when
     ! the pressure fix is active, or set_uold discards the heat in
     ! truncation-dominated cells)
     do j=1,np
        if(ok(j))then
           wv=vol_loc
        else
           wv=vol_loc*2**ndim
           nfall_c=nfall_c+1
        end if
        dE=0.0d0
        dEcr=0.0d0
        if(do_heat)dE=mp(ind_part(j))*efac/wv
        if(do_cr)dEcr=mp(ind_part(j))*efaccr/wv
        if(own(j))then
           if(do_heat)then
              unew(indp(j),ndim+2)=unew(indp(j),ndim+2)+dE
              if(pressure_fix)then
                 enew(indp(j))=enew(indp(j))+dE
              end if
           end if
           if(do_cr)then
              unew(indp(j),pbh_cr_ivar)=unew(indp(j),pbh_cr_ivar)+dEcr
           end if
        else
           call pbh_rbuf_push(pbh_rbuf(tid),indp(j), &
                & int(idp(ind_part(j)),kind=8),dE,dEcr)
        end if
        if(do_heat)einj_c=einj_c+mp(ind_part(j))*efac
        if(do_cr)ecr_c=ecr_c+mp(ind_part(j))*efaccr
     end do
  end if
#else
  if(do_heat.or.do_cr)then
     write(*,*)'PBH ERROR: pbh_evap_fine requires NDIM=3'
     call clean_stop
  end if
#endif

  ! Exact mixed-mass update (each particle belongs to exactly one grid,
  ! hence one thread; no protection needed)
  do j=1,np
     mp(ind_part(j))=mp(ind_part(j))*ratio
  end do

end subroutine pbh_dump
