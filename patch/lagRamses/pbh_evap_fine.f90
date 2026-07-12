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
! NGP deposit follows feedbk (feedback.kjhan3.f90) with two changes:
! particles whose child grid is absent fall back to the coarser parent
! cell instead of being skipped, and the deposit lines are atomic.
!=======================================================================
subroutine pbh_evap_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use pbh_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::ichunk,ncache,info
  real(dp)::dt_phys,ratio,dQ,wbeg,efac,conv
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::einj_all
  integer(kind=8)::nfall_all
  logical::do_heat
  integer,parameter::NGRID_CHUNK=64

  if(.not.use_pbh)return
  if(.not.pic)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  call pbh_lazy_init(aexp,aexp_ini,nrestart,myid,cosmo)

  ! step-global factors (aexp is the END of this level step here,
  ! update_time having already run)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  dt_phys=dtnew(ilevel)*scale_t
  call pbh_step(ilevel,nlevelmax,aexp,dt_phys,ratio,dQ,wbeg)
  do_heat=hydro .and. trim(pbh_energy_sink)=='local_heat' .and. dQ>0.0d0
  efac=pbh_fraction*dQ/(wbeg*scale_v**2)      ! per unit current code mass
  conv=scale_v**2*scale_d*scale_l**3        ! code energy -> erg

  ncache=active(ilevel)%ngrid
!$omp parallel do private(ichunk) schedule(dynamic)
  do ichunk=1,ncache,NGRID_CHUNK
     call sub_pbh_evap_fine(ilevel,ichunk,min(NGRID_CHUNK,ncache-ichunk+1), &
          & ratio,efac,conv,do_heat)
  end do
!$omp end parallel do

  ! coarse-step diagnostics (all ranks reach this when numbtot>0)
  if(ilevel==levelmin)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(pbh_einj_loc,einj_all,1,MPI_DOUBLE_PRECISION, &
          & MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(pbh_nfallback_loc,nfall_all,1,MPI_INTEGER8, &
          & MPI_SUM,MPI_COMM_WORLD,info)
#else
     einj_all=pbh_einj_loc
     nfall_all=pbh_nfallback_loc
#endif
     pbh_einj_tot=einj_all
     if(myid==1)write(*,112)aexp,wbeg*ratio,rho_tot,einj_all,nfall_all,eint_tot
  end if

111 format('   Entering pbh_evap_fine for level ',I2)
112 format(' PBHDIAG aexp=',ES14.7,' w=',ES16.9,' rho_tot=',ES16.9, &
         & ' einj[erg]=',ES14.7,' nfallback=',I12,' eint=',ES16.9)
end subroutine pbh_evap_fine

!=======================================================================
subroutine sub_pbh_evap_fine(ilevel,istart,nchunk,ratio,efac,conv,do_heat)
  ! Gather DM particles of one chunk of active grids into nvector
  ! batches (classic feedback pattern) and hand them to pbh_dump.
  use amr_commons
  use pm_commons
  use pbh_commons, only: pbh_energy_sink
  implicit none
  integer::ilevel,istart,nchunk
  real(dp)::ratio,efac,conv
  logical::do_heat
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
                   & ratio,efac,conv,do_heat)
              ip=0
              ig=0
           end if
           ipart=next_part
        end do
     end if
  end do
  if(ip>0)call pbh_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel, &
       & ratio,efac,conv,do_heat)
end subroutine sub_pbh_evap_fine

!=======================================================================
subroutine pbh_dump(ind_grid,ind_part,ind_grid_part,ng,np,ilevel, &
     & ratio,efac,conv,do_heat)
  ! NGP energy deposit + exact mass update for one particle batch.
  use amr_commons
  use pm_commons
  use hydro_commons
  use pbh_commons, only: pbh_einj_loc,pbh_nfallback_loc
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::ratio,efac,conv
  logical::do_heat
  integer::i,j,idim,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dE,einj,wv
  integer(kind=8)::nfall
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::x0
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector)::ok
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
  if(do_heat)then
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
     einj=0.0d0
     nfall=0
     do j=1,np
        if(ok(j))then
           wv=vol_loc
        else
           wv=vol_loc*2**ndim
           nfall=nfall+1
        end if
        dE=mp(ind_part(j))*efac/wv
!$omp atomic
        unew(indp(j),ndim+2)=unew(indp(j),ndim+2)+dE
        if(pressure_fix)then
!$omp atomic
           enew(indp(j))=enew(indp(j))+dE
        end if
        einj=einj+mp(ind_part(j))*efac
     end do
!$omp atomic
     pbh_einj_loc=pbh_einj_loc+einj*conv
!$omp atomic
     pbh_nfallback_loc=pbh_nfallback_loc+nfall
  end if
#else
  if(do_heat)then
     write(*,*)'PBH ERROR: pbh_evap_fine requires NDIM=3'
     call clean_stop
  end if
#endif

  ! Exact mixed-mass update (each particle belongs to exactly one grid,
  ! hence one thread; no atomic needed)
  do j=1,np
     mp(ind_part(j))=mp(ind_part(j))*ratio
  end do

end subroutine pbh_dump
