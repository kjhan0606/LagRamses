!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag
  use amr_commons
  implicit none
  integer::ilevel

  if(verbose)write(*,*)'Entering flag'
  do ilevel=nlevelmax-1,1,-1
     call flag_fine(ilevel,2)
  end do
  call flag_coarse
  if(verbose)write(*,*)'Complete flag'

end subroutine flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag_coarse
  use amr_commons
  implicit none
  !--------------------------------------------------------------
  ! This routine compute the refinement map at the coarse level.
  !--------------------------------------------------------------
  integer::ind,nxny,ix,iy,iz

  if(verbose)write(*,*)'  Entering flag_coarse'
  ! Constants
  nxny=nx*ny
  ! Reset flag1 array at coarse level
  flag1(0:ncoarse)=0
  ! Set flag1 to 1 at coarse level for inner cells only
  nflag=0
  do iz=kcoarse_min,kcoarse_max
     do iy=jcoarse_min,jcoarse_max
        do ix=icoarse_min,icoarse_max
           ind=1+ix+iy*nx+iz*nxny
           flag1(ind)=1
           nflag=nflag+1
        end do
     end do
  end do
  if(verbose)write(*,112)nflag  
  call make_virtual_coarse_int(flag1(1))
  if(simple_boundary)call make_boundary_coarse

112 format('   ==> Flag ',i6,' cells')

end subroutine flag_coarse
!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag_fine(ilevel,icount)
  use amr_commons
  implicit none
  integer::ilevel,icount
  !--------------------------------------------------------
  ! This routine builds the refinement map at level ilevel.
  !--------------------------------------------------------
  integer::iexpand

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Step 1: initialize refinement map to minimal refinement rules
  call init_flag(ilevel)
  if(verbose)write(*,*) ' +==> end step 1',nflag

  ! If ilevel < levelmin, exit routine
  if(ilevel<levelmin)return
  if(balance)return

  ! Precompute neighbor grids for this level
  call precompute_nbor_active(ilevel)

  ! Step 2: make one cubic buffer around flagged cells,
  ! in order to enforce numerical rule.
  call smooth_fine(ilevel)
  if(verbose)write(*,*) ' +==> end step 2',nflag

  ! Step 3: if cell satisfies user-defined physical citeria,
  ! then flag cell for refinement.
  call userflag_fine(ilevel)    
  if(verbose)write(*,*) ' +==> end step 3',nflag

  ! Step 4: make nexpand cubic buffers around flagged cells.
  do iexpand=1,nexpand(ilevel)
     call smooth_fine(ilevel)
  end do
  if(verbose)write(*,*) ' +==> end step 4',nflag

  ! Cleanup neighbor cache
  call cleanup_nbor_active()

  if(verbose)write(*,112)nflag

  ! In case of adaptive time step ONLY, check for refinement rules.
  if(ilevel>levelmin)then
     if(icount<nsubcycle(ilevel-1))then
        call ensure_ref_rules(ilevel)
     end if
  end if

111 format('  +Entering flag_fine for level ',I2)
112 format('   ==> Flag ',i6,' cells')

end subroutine flag_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flag(ilevel)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ilevel
  !-------------------------------------------
  ! This routine initialize the refinement map
  ! to a minimal state in order to satisfy the
  ! refinement rules.
  !-------------------------------------------
  integer::i,ind,mflag

  ! Initialize flag1 to 0
  nflag=0
! do ind=1,twotondim
!    iskip=ncoarse+(ind-1)*ngridmax
!    do i=1,active(ilevel)%ngrid
!       flag1(active(ilevel)%igrid(i)+iskip)=0
!    end do
! end do
! ! If load balancing operations, flag only refined cells
! if(balance)then
!!$omp parallel do private(ind,iskip,i) reduction(+: nflag)
!    do ind=1,twotondim
!       iskip=ncoarse+(ind-1)*ngridmax
!       do i=1,active(ilevel)%ngrid
!          if(son(active(ilevel)%igrid(i)+iskip)>0)then
!             flag1(active(ilevel)%igrid(i)+iskip)=1
!             nflag=nflag+1
!          end if
!       end do
!    end do
! else
!    ! If cell is refined and contains a flagged son
!    ! or a refined son, then flag cell for refinement.
!    if(ilevel>=levelmin)then
!       call test_flag(ilevel)
!    else
!       ! If ilevel < levelmin, set flag to 1 for all cells
!!$omp parallel do private(ind,iskip,i) reduction(+: nflag)
!       do ind=1,twotondim
!          iskip=ncoarse+(ind-1)*ngridmax
!          do i=1,active(ilevel)%ngrid
!             flag1(active(ilevel)%igrid(i)+iskip)=1
!          end do
!          nflag=nflag+active(ilevel)%ngrid
!       end do
!    end if
! end if
!$omp parallel do private(i,ind)
  do i=1,active(ilevel)%ngrid
     do ind=1,twotondim
        flag1(ICELL_OF(active(ilevel)%igrid(i),ind))=0
     end do
  end do

  ! If load balancing operations, flag only refined cells
  if(balance)then
     mflag = 0
!$omp parallel do private(ind,i) reduction(+: mflag)
     do i=1,active(ilevel)%ngrid
        do ind=1,twotondim
           if(son(ICELL_OF(active(ilevel)%igrid(i),ind))>0)then
              flag1(ICELL_OF(active(ilevel)%igrid(i),ind))=1
              mflag=mflag+1
           end if
        end do
     end do
     nflag = nflag + mflag
  else
     ! If cell is refined and contains a flagged son
     ! or a refined son, then flag cell for refinement.
     if(ilevel>=levelmin)then
        call test_flag(ilevel)
     else
        ! If ilevel < levelmin, set flag to 1 for all cells
        mflag = 0
!$omp parallel do private(ind,i) reduction(+: mflag)
        do i=1,active(ilevel)%ngrid
           do ind=1,twotondim
              flag1(ICELL_OF(active(ilevel)%igrid(i),ind))=1
           end do
           mflag=mflag+twotondim
        end do
        nflag = nflag + mflag
     end if
  end if
  
  ! Update boundaries
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

end subroutine init_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine test_flag(ilevel)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets flag1 to 1 if cell is refined and 
  ! contains a flagged son or a refined son.
  ! This ensures that refinement rules are satisfied.
  !---------------------------------------------------------
  integer::i,ind_son,ind
  integer::ind_grid_son,ind_cell_son
  logical::ok

  ! Loop over cells
  do ind=1,twotondim
     ! Test all refined cells
     do i=1,active(ilevel)%ngrid
        ! Gather child grid number
        ind_grid_son=son(ICELL_OF(active(ilevel)%igrid(i),ind))
        ! Test child if it exists
        ok=.false.
        if(ind_grid_son>0)then
           ! Loop over children cells
           do ind_son=1,twotondim
              ind_cell_son=ICELL_OF(ind_grid_son,ind_son)
              ok=(ok.or.(son  (ind_cell_son)> 0))
              ok=(ok.or.(flag1(ind_cell_son)==1))
           end do
        end if
        ! If ok, then flag1 cells.
        if(ok)then
           flag1(ICELL_OF(active(ilevel)%igrid(i),ind))=1
           nflag=nflag+1
        end if
     end do
  end do
  ! End loop over cells

end subroutine test_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine ensure_ref_rules(ilevel)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------
  ! This routine determines if all grids at level ilevel are 
  ! surrounded by 26 neighboring grids, in order to enforce the 
  ! strict refinement rule. 
  ! Used in case of adaptive time steps only.
  !-----------------------------------------------------------------
  integer::i,ind,igrid,ngrid,ncache
  integer,dimension(1:nvector)::ind_cell,ind_grid
  integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim)::nbors_father_grids
  logical,dimension(1:nvector)::ok

  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Gather neighboring father cells (should be present anytime !)
     do i=1,ngrid
        ind_cell(i)=father(ind_grid(i))
     end do
     call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids &
          & ,ngrid,ilevel)
     
     do i=1,ngrid
        ok(i)=.true.
     end do

     do ind=1,threetondim
        do i=1,ngrid
           ind_cell(i)=nbors_father_cells(i,ind)
           if(ind_cell(i)==0)ok(i)=.false.
        end do
        do i=1,ngrid
           if(ind_cell(i)>0)then
              if(son(ind_cell(i))==0)ok(i)=.false.
           endif
        end do
     end do
     
     do ind=1,twotondim
        do i=1,ngrid
           ind_cell(i)=ICELL_OF(ind_grid(i),ind)
        end do
        do i=1,ngrid
           if(.not.ok(i))flag1(ind_cell(i))=0
        end do
     end do

  end do

  ! Update boundaries
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

end subroutine ensure_ref_rules 
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine sub_userflag_fine(ilevel,skip_loc,scale, igrid,ngrid,iflag)
  use amr_commons
  use hydro_commons
  use pm_commons    ! headp, nextp, idp, tp, xp for sink particle check
  use cooling_module
#include "amr_index.h"
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel.
  ! -------------------------------------------------------------------
  integer::i,j,ncache,nok,ix,iy,iz,iflag,jflag
  integer::igrid,ind,idim,ngrid,ivar
  integer::nx_loc
  integer::ipart,ind_part
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector)::sink_mask

  logical,dimension(1:nvector)::ok

  real(dp)::dx,dx_loc,scale
  real(dp)::d0,dx_min,vol_min,mstar,msnk,nISM,nCOM
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx

  iflag = 0

  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do

  ! Pre-compute sink particle child cell bitmask per grid
  sink_mask(1:ngrid) = 0
  if(sink_refine .and. sink .and. pic)then
     do i=1,ngrid
        ipart = headp(ind_grid(i))
        do while(ipart > 0)
           if(ptypep(ipart) == PTYPE_SINK)then
              ! Determine child cell index (1..8) from particle position
              ind_part = 1
              if(xp(ipart,1) > xg(ind_grid(i),1)) ind_part = ind_part + 1
              if(xp(ipart,2) > xg(ind_grid(i),2)) ind_part = ind_part + 2
              if(xp(ipart,3) > xg(ind_grid(i),3)) ind_part = ind_part + 4
              sink_mask(i) = ior(sink_mask(i), ishft(1, ind_part-1))
           end if
           ipart = nextp(ipart)
        end do
     end do
  end if

  ! Loop over cells
  do ind=1,twotondim

     do i=1,ngrid
        ind_cell(i)=ICELL_OF(ind_grid(i),ind)
     end do

     ! Initialize refinement to false
     do i=1,ngrid
        ok(i)=.false.
     end do

     ! Apply purely local Lagrangian refinement criteria
     if(m_refine_eff(ilevel)>-1.0d0)then
        call poisson_refine(ind_cell,ok,ngrid,ilevel)
        ! Apply sink particle refinement from pre-computed bitmask
        if(sink_refine .and. sink .and. pic)then
           do i=1,ngrid
              if(btest(sink_mask(i), ind-1)) ok(i) = .true.
           end do
        end if
        ! Apply geometry-based refinement criteria
        if(r_refine(ilevel)>-1.0)then
           ! Compute cell center in code units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do
           call geometry_refine(xx,ind_cell,ok,ngrid,ilevel)
        end if
     end if

     ! Count newly flagged cells
     nok=0
     do i=1,ngrid
        if(flag1(ind_cell(i))==0.and.ok(i))then
           nok=nok+1
        end if
     end do
     
     do i=1,ngrid
        if(ok(i))flag1(ind_cell(i))=1
     end do

     iflag=iflag+nok
  end do
     ! End loop over cells
end subroutine sub_userflag_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine userflag_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  integer::i,j,ncache,nok,ix,iy,iz,iskip,iflag,jflag
  integer::igrid,ind,idim,ngrid,ivar
  integer::nx_loc
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  integer,dimension(1:nvector,1:twondim)::indn

  logical,dimension(1:nvector)::ok

  real(dp)::dx,dx_loc,scale
  real(dp)::d0,dx_min,vol_min,mstar,msnk,nISM,nCOM
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx

  logical::prevent_refine

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return

  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set positions of cell centers relative to each grid center.  The void
  ! floor is evaluated before the cooling holdback check below.
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Do we prevent the whole level from refining ?
  prevent_refine=.false.
  
  ! Prevent over refinement due to gas cooling
  ! This translates into :
  ! - a constant physical resolution at low redshift (ilevel<=nlevelmax_part+nlevel_collapse)
  ! - a constant comobile resolution at high redshift (ilevel>nlevelmax_part+nlevel_collapse)
  if(cosmo.and.cooling.and.q_refine_holdback)then
     ! Binary holdback: prevent refinement at high levels at early times
     dx_min=(0.5D0**nlevelmax)*scale
     if(ilevel.gt.nlevelmax_part+nlevel_collapse)then
         if(dx_loc<(4d0**(1d0/ndim))*(dx_min/aexp)) prevent_refine=.true.
     endif
  endif

  ! The void floor is independent of density and remains active during the
  ! cooling holdback.  It is evaluated first so the early return below only
  ! suppresses the ordinary high-level criteria.
  if(void_refine .and. ilevel<void_refine_min_level)then
     call voidflag_fine(ilevel,skip_loc,scale,xc)
  end if

  ! This is a separate, opt-in path for the void simulations.  Its mesh
  ! floors remain active during the ordinary cooling holdback, but it is a
  ! complete no-op (including no allocation or communication) when disabled.
  if(void_web_refine .and. ilevel<void_web_wall_level)then
     call void_webflag_fine(ilevel)
  end if

  if(prevent_refine) then
     call make_virtual_fine_int(flag1(1),ilevel)
     if(simple_boundary)call make_boundary_flag(ilevel)
     return
  end if

  ! FDM refinement: de Broglie (wave levels) + Madelung (fluid levels)
  ! Placed after prevent_refine so gas holdback gates FDM refinement too
  if(use_fdm .and. .not.fdm_refine_matched) call fdm_refine_flag(ilevel)
  if(use_fdm .and. .not.fdm_refine_matched) call fdm_madelung_refine_flag(ilevel)

  ! Compute FPR-adjusted effective m_refine for this level
  call compute_fpr_m_refine_eff(ilevel)

  ! Loop over active grids
  ncache=active(ilevel)%ngrid
  iflag = 0
!$omp parallel do private(igrid,ngrid,jflag) reduction(+:iflag)
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
	 call sub_userflag_fine(ilevel, skip_loc,scale, igrid,ngrid,jflag)
	 iflag = iflag + jflag
  enddo
  nflag = nflag + iflag



  ! Do the same for hydro solver
  if(hydro)call hydro_flag(ilevel)

#ifdef RT
  ! Do the same for RT solver
  if(rt)call rt_hydro_flag(ilevel)
#endif

  ! Update boundaries
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

end subroutine userflag_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine voidflag_fine(ilevel,skip_loc,scale,xc)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::scale
  real(dp),dimension(1:3),intent(in)::skip_loc
  real(dp),dimension(1:twotondim,1:3),intent(in)::xc
  integer::i,ind,idim,ind_grid,ind_cell,nnew
  real(dp),dimension(1:ndim)::xx
  logical::refine_region_contains

  nnew=0
!$omp parallel do private(i,ind,idim,ind_grid,ind_cell,xx) reduction(+:nnew)
  do i=1,active(ilevel)%ngrid
     ind_grid=active(ilevel)%igrid(i)
     do ind=1,twotondim
        ind_cell=ICELL_OF(ind_grid,ind)
        do idim=1,ndim
           xx(idim)=(xg(ind_grid,idim)+xc(ind,idim)-skip_loc(idim))*scale
        end do
        if(refine_region_contains(xx,ilevel))then
           if(flag1(ind_cell)==0)then
              flag1(ind_cell)=1
              nnew=nnew+1
           end if
        end if
     end do
  end do
!$omp end parallel do
  nflag=nflag+nnew

end subroutine voidflag_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine void_webflag_fine(ilevel)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  integer::i,ind,ind_grid,ind_cell,ancestor_cell
  integer::ancestor_level,ancestor_grid,state,nnew
  logical::refine_cell

  if(.not.void_web_refine)return
  if(ilevel>=void_web_wall_level)return

  call ensure_void_web_state
  if(.not.void_web_state_valid)return

  nnew=0
!$omp parallel do private(i,ind,ind_grid,ind_cell,ancestor_cell) &
!$omp& private(ancestor_level,ancestor_grid,state,refine_cell) reduction(+:nnew)
  do i=1,active(ilevel)%ngrid
     ind_grid=active(ilevel)%igrid(i)
     do ind=1,twotondim
        ind_cell=ICELL_OF(ind_grid,ind)

        ! Follow the AMR tree to the fixed environmental scale.  The state
        ! therefore does not depend on the current leaf level.
        ancestor_cell=ind_cell
        ancestor_level=ilevel
        do while(ancestor_level>void_web_env_level .and. ancestor_cell>ncoarse)
           ancestor_grid=IGRID_OF(ancestor_cell)
           ancestor_cell=father(ancestor_grid)
           ancestor_level=ancestor_level-1
        end do

        state=0
        if(ancestor_level==void_web_env_level .and. ancestor_cell>0) &
             & state=void_web_state(ancestor_cell)
        refine_cell=.false.
        if(ilevel<void_web_base_level .and. iand(state,1)/=0)refine_cell=.true.
        if(ilevel<void_web_wall_level .and. iand(state,2)/=0)refine_cell=.true.
        if(refine_cell .and. flag1(ind_cell)==0)then
           flag1(ind_cell)=1
           nnew=nnew+1
        end if
     end do
  end do
!$omp end parallel do
  nflag=nflag+nnew

end subroutine void_webflag_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine ensure_void_web_state
  use amr_commons
  implicit none
  logical::update_due

  if(.not.void_web_refine)return

  update_due=.not.void_web_state_valid
  update_due=update_due.or.(void_web_state_epoch/=amr_mesh_epoch)
  if(void_web_state_valid)then
     update_due=update_due.or. &
          & (nstep_coarse-void_web_state_step>=void_web_update_interval)
  end if
  if(update_due)call update_void_web_state

end subroutine ensure_void_web_state
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine update_void_web_state
  use amr_commons
  use hydro_commons
  use omp_lib, only: omp_get_wtime
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,ngrid,ncache,array_size,mpi_err
  integer::nscope_chunk,nwall_chunk,nvalid_chunk
  integer::nscope_loc,nwall_loc,nvalid_loc
  integer::nscope_all,nwall_all,nvalid_all
  integer,dimension(1:3)::count_loc,count_all
  real(dp)::dx_loc,rate_norm,lambda_min_chunk,lambda_max_chunk
  real(dp)::lambda_min_loc,lambda_max_loc,lambda_min_all,lambda_max_all
  real(dp)::update_time_start,update_time_loc,update_time_all
  real(dp),dimension(1:3)::diagnostic_loc,diagnostic_all
  logical::reset_state

  if(.not.void_web_refine)return
  update_time_start=omp_get_wtime()

  array_size=ncoarse+twotondim*ngridmax
  reset_state=.not.allocated(void_web_state)
  if(allocated(void_web_state))then
     if(size(void_web_state)/=array_size)then
        deallocate(void_web_state)
        reset_state=.true.
     end if
  end if
  if(.not.allocated(void_web_state))allocate(void_web_state(1:array_size))
  if(void_web_state_epoch/=amr_mesh_epoch)reset_state=.true.
  if(reset_state)void_web_state=0

  ! In RAMSES cosmological units v_code/x_code divided by a*hexp is
  ! (dv_phys/dx_comoving)/(a H/h).  This makes the eigenvalue threshold
  ! dimensionless and identical at all epochs.
  dx_loc=0.5d0**void_web_env_level*boxlen/dble(icoarse_max-icoarse_min+1)
  rate_norm=aexp*hexp
  if(rate_norm<=tiny(1.0d0))then
     if(myid==1)write(*,*)'Invalid aexp*hexp in void V-web calculation'
     call clean_stop
  end if

  nscope_loc=0
  nwall_loc=0
  nvalid_loc=0
  lambda_min_loc=huge(1.0d0)
  lambda_max_loc=-huge(1.0d0)
  ncache=active(void_web_env_level)%ngrid
!$omp parallel do private(igrid,ngrid,nscope_chunk,nwall_chunk,nvalid_chunk) &
!$omp& private(lambda_min_chunk,lambda_max_chunk) &
!$omp& reduction(+:nscope_loc,nwall_loc,nvalid_loc) &
!$omp& reduction(min:lambda_min_loc) reduction(max:lambda_max_loc)
  do igrid=1,ncache,nvector
     ngrid=min(nvector,ncache-igrid+1)
     call sub_update_void_web_state(igrid,ngrid,dx_loc,rate_norm, &
          & nscope_chunk,nwall_chunk,nvalid_chunk, &
          & lambda_min_chunk,lambda_max_chunk)
     nscope_loc=nscope_loc+nscope_chunk
     nwall_loc=nwall_loc+nwall_chunk
     nvalid_loc=nvalid_loc+nvalid_chunk
     lambda_min_loc=min(lambda_min_loc,lambda_min_chunk)
     lambda_max_loc=max(lambda_max_loc,lambda_max_chunk)
  end do
!$omp end parallel do

  call make_virtual_fine_int(void_web_state(1),void_web_env_level)
  if(simple_boundary)then
     if(myid==1)write(*,*)'void_web_refine currently expects periodic boundaries'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  update_time_loc=omp_get_wtime()-update_time_start
  ! Pack diagnostics to keep the feature at two scalar collectives per update:
  ! one sum and one minimum (maxima are represented by their negatives).
  count_loc=(/nscope_loc,nwall_loc,nvalid_loc/)
  call MPI_ALLREDUCE(count_loc,count_all,3,MPI_INTEGER,MPI_SUM, &
       & MPI_COMM_WORLD,mpi_err)
  nscope_all=count_all(1)
  nwall_all=count_all(2)
  nvalid_all=count_all(3)
  diagnostic_loc=(/lambda_min_loc,-lambda_max_loc,-update_time_loc/)
  call MPI_ALLREDUCE(diagnostic_loc,diagnostic_all,3,MPI_DOUBLE_PRECISION, &
       & MPI_MIN,MPI_COMM_WORLD,mpi_err)
  lambda_min_all=diagnostic_all(1)
  lambda_max_all=-diagnostic_all(2)
  update_time_all=-diagnostic_all(3)
#else
  nscope_all=nscope_loc
  nwall_all=nwall_loc
  nvalid_all=nvalid_loc
  lambda_min_all=lambda_min_loc
  lambda_max_all=lambda_max_loc
  update_time_all=omp_get_wtime()-update_time_start
#endif

  void_web_state_step=nstep_coarse
  void_web_state_epoch=amr_mesh_epoch
  void_web_state_valid=.true.
  if(myid==1)then
     write(*,'(A,I8,A,I3,A,I12,A,I12)')' Void V-web update: step=', &
          & nstep_coarse,' level=',void_web_env_level,' scope=',nscope_all, &
          & ' wall=',nwall_all
     if(nvalid_all>0)write(*,'(A,2ES12.4)')'   lambda_max range = ', &
          & lambda_min_all,lambda_max_all
     write(*,'(A,F10.4,A)')'   V-web update wall time = ',update_time_all,' s'
     if(nscope_all==0)write(*,'(A)') &
          & '   WARNING: void V-web scope is empty; no mesh floor will be applied'
  end if

end subroutine update_void_web_state
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine sub_update_void_web_state(igrid_start,ngrid,dx_loc,rate_norm, &
     & nscope,nwall,nvalid,lambda_min,lambda_max)
  use amr_commons
  use hydro_commons
  use hydro_parameters, only: smallr
#include "amr_index.h"
  implicit none
  integer,intent(in)::igrid_start,ngrid
  real(dp),intent(in)::dx_loc,rate_norm
  integer,intent(out)::nscope,nwall,nvalid
  real(dp),intent(out)::lambda_min,lambda_max
  integer::i,ind,idim,jdim,ind_cell,old_state,new_state
  integer::cell_minus,cell_plus
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector,0:twondim)::igridn
  integer,dimension(1:nvector,1:twondim)::indn
  real(dp)::rho_minus,rho_plus,lambda
  real(dp),dimension(1:ndim,1:ndim)::grad,sigma
  logical::in_scope,neighbors_ok,wall_on
  real(dp)::void_web_lambda_max_3d

  do i=1,ngrid
     ind_grid(i)=active(void_web_env_level)%igrid(igrid_start+i-1)
  end do
  call getnborgrids(ind_grid,igridn,ngrid)

  nscope=0
  nwall=0
  nvalid=0
  lambda_min=huge(1.0d0)
  lambda_max=-huge(1.0d0)
  do ind=1,twotondim
     call getnborcells(igridn,ind,indn,ngrid)
     do i=1,ngrid
        ind_cell=ICELL_OF(ind_grid(i),ind)
        select case(void_web_scope_ivar)
        case(-1)
           in_scope=.true.
        case(0)
           in_scope=(cpu_map2(ind_cell)==1)
        case default
           in_scope=(uold(ind_cell,void_web_scope_ivar) / &
                & max(uold(ind_cell,1),smallr)>void_web_scope_cut)
        end select

        old_state=void_web_state(ind_cell)
        new_state=0
        if(in_scope)then
           new_state=1
           nscope=nscope+1
           neighbors_ok=.true.
           grad=0.0d0
           do jdim=1,ndim
              cell_minus=indn(i,2*jdim-1)
              cell_plus =indn(i,2*jdim)
              if(cell_minus<=0 .or. cell_plus<=0)then
                 neighbors_ok=.false.
              else
                 rho_minus=max(uold(cell_minus,1),smallr)
                 rho_plus =max(uold(cell_plus ,1),smallr)
                 do idim=1,ndim
                    grad(idim,jdim)= &
                         & (uold(cell_plus,idim+1)/rho_plus- &
                         &  uold(cell_minus,idim+1)/rho_minus)/(2.0d0*dx_loc)
                 end do
              end if
           end do

           if(neighbors_ok)then
              do jdim=1,ndim
                 do idim=1,ndim
                    sigma(idim,jdim)=-0.5d0*(grad(idim,jdim)+ &
                         & grad(jdim,idim))/rate_norm
                 end do
              end do
              lambda=void_web_lambda_max_3d(sigma)
              nvalid=nvalid+1
              lambda_min=min(lambda_min,lambda)
              lambda_max=max(lambda_max,lambda)
              wall_on=(iand(old_state,2)/=0)
              if(wall_on)then
                 wall_on=(lambda>=void_web_lambda_off)
              else
                 wall_on=(lambda>=void_web_lambda_on)
              end if
              if(wall_on)then
                 new_state=ior(new_state,2)
                 nwall=nwall+1
              end if
           end if
        end if
        void_web_state(ind_cell)=new_state
     end do
  end do

end subroutine sub_update_void_web_state
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
function void_web_lambda_max_3d(a) result(lambda_max)
  use amr_parameters, only: dp
  implicit none
  real(dp),dimension(1:3,1:3),intent(in)::a
  real(dp)::lambda_max
  real(dp)::q,p2,p,r,phi,detb
  real(dp),dimension(1:3,1:3)::b
  ! Stable closed-form eigenvalue for a real symmetric 3x3 matrix.
  q=(a(1,1)+a(2,2)+a(3,3))/3.0d0
  p2=(a(1,1)-q)**2+(a(2,2)-q)**2+(a(3,3)-q)**2+ &
       & 2.0d0*(a(1,2)**2+a(1,3)**2+a(2,3)**2)
  if(p2<=tiny(1.0d0))then
     lambda_max=q
     return
  end if
  p=sqrt(p2/6.0d0)
  b=a
  b(1,1)=b(1,1)-q
  b(2,2)=b(2,2)-q
  b(3,3)=b(3,3)-q
  b=b/p
  detb=b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2))- &
       & b(1,2)*(b(2,1)*b(3,3)-b(2,3)*b(3,1))+ &
       & b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))
  r=max(-1.0d0,min(1.0d0,0.5d0*detb))
  phi=acos(r)/3.0d0
  lambda_max=q+2.0d0*p*cos(phi)

end function void_web_lambda_max_3d
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine poisson_refine(ind_cell,ok,ncell,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  logical,dimension(1:nvector)::ok
  !-------------------------------------------------
  ! This routine sets flag1 to 1 if cell statisfy
  ! user-defined physical criterion for refinement.
  !-------------------------------------------------
  integer::i,nx_loc
  real(dp)::d_scale,d_scale_fdm,scale,dx,dx_loc,vol_loc

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale
  vol_loc=dx_loc**3
  ! FDM quasi-Lagrangian reference: mean DM mass per levelmin cell.
  ! With hydro the gas carries Omega_b and |psi|^2 averages
  ! (1 - Omega_b/Omega_m); without hydro |psi|^2 averages 1.
  d_scale_fdm=0.5d0**(ndim*levelmin)/vol_loc
  if(hydro .and. omega_m > 0.0d0) d_scale_fdm = d_scale_fdm*(1.0d0-omega_b/omega_m)

  if(poisson)then

     if(.not. init) then
        if(ivar_refine < 0 .and. m_refine_eff(ilevel) > 0.0d0)then
           ! Zoom-in: use density criterion (quasi-Lagrangian)
           d_scale=mass_sph/vol_loc
           if(use_fdm)then
              ! FDM: density is |psi|^2; uold (hydro) is empty under hydro=.false.
              do i=1,ncell
                 if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
                    ok(i)=ok(i).or.(psi_re(ind_cell(i))>=m_refine_eff(ilevel)*d_scale_fdm)
                 else
                    ok(i)=ok(i).or.((psi_re(ind_cell(i))**2+psi_im(ind_cell(i))**2)>=m_refine_eff(ilevel)*d_scale_fdm)
                 end if
              end do
           else if(hydro)then
              do i=1,ncell
                 ok(i)=ok(i).or.(uold(ind_cell(i),1)>=m_refine_eff(ilevel)*d_scale)
              end do
           else
              ! dmonly: uold is NOT allocated (hydro=.false.) — use the
              ! particle density rho (mean-density units, mean=1). The
              ! particle-count criterion is rho >= m_refine * d_scale_fdm
              ! (same reference mass per levelmin cell as the FDM branch).
              do i=1,ncell
                 ok(i)=ok(i).or.(rho(ind_cell(i))>=m_refine_eff(ilevel)*d_scale_fdm)
              end do
           end if
        else
           if(use_fdm)then
              do i=1,ncell
                 if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
                    ok(i)=ok(i).or.(psi_re(ind_cell(i))>=m_refine_eff(ilevel)*d_scale_fdm)
                 else
                    ok(i)=ok(i).or.((psi_re(ind_cell(i))**2+psi_im(ind_cell(i))**2)>=m_refine_eff(ilevel)*d_scale_fdm)
                 end if
              end do
           else
              do i=1,ncell
                 ok(i)=ok(i).or.(cpu_map2(ind_cell(i))==1)
              end do
           end if
        end if
     else
        if(ivar_refine==0)then
           if(use_fdm)then
              do i=1,ncell
                 if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
                    ok(i)=ok(i).or.(psi_re(ind_cell(i))>=m_refine_eff(ilevel)*d_scale_fdm)
                 else
                    ok(i)=ok(i).or.((psi_re(ind_cell(i))**2+psi_im(ind_cell(i))**2)>=m_refine_eff(ilevel)*d_scale_fdm)
                 end if
              end do
           else
              do i=1,ncell
                 ok(i)=ok(i).or.(cpu_map2(ind_cell(i))==1)
              end do
           end if
        else if(ivar_refine>0)then
           do i=1,ncell
              ok(i)=ok(i).or. &
                   & (uold(ind_cell(i),ivar_refine)/uold(ind_cell(i),1) &
                   & > var_cut_refine)
           end do
        else if(m_refine_eff(ilevel)==0.0)then
           do i=1,ncell
              ok(i)=.true.
           end do
        else
           ! ivar_refine<0, m_refine>0: density-based (zoom-in init)
           d_scale=mass_sph/vol_loc
           if(use_fdm)then
              ! FDM: density is |psi|^2; uold (hydro) is empty under hydro=.false.
              do i=1,ncell
                 if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
                    ok(i)=ok(i).or.(psi_re(ind_cell(i))>=m_refine_eff(ilevel)*d_scale_fdm)
                 else
                    ok(i)=ok(i).or.((psi_re(ind_cell(i))**2+psi_im(ind_cell(i))**2)>=m_refine_eff(ilevel)*d_scale_fdm)
                 end if
              end do
           else if(hydro)then
              do i=1,ncell
                 ok(i)=ok(i).or.(uold(ind_cell(i),1)>=m_refine_eff(ilevel)*d_scale)
              end do
           else
              ! dmonly: uold is NOT allocated (hydro=.false.) — use the
              ! particle density rho (mean-density units, mean=1). The
              ! particle-count criterion is rho >= m_refine * d_scale_fdm
              ! (same reference mass per levelmin cell as the FDM branch).
              do i=1,ncell
                 ok(i)=ok(i).or.(rho(ind_cell(i))>=m_refine_eff(ilevel)*d_scale_fdm)
              end do
           end if
        endif
     endif

  else

     if(hydro)then
        d_scale=mass_sph/vol_loc
        do i=1,ncell
           ok(i)=ok(i).or.(uold(ind_cell(i),1)>=m_refine_eff(ilevel)*d_scale)
        end do
     endif

  end if

end subroutine poisson_refine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine geometry_refine(xx,ind_cell,ok,ncell,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:ndim)::xx
  logical ,dimension(1:nvector)::ok
  !-------------------------------------------------
  ! This routine sets flag1 to 1 if cell statisfy
  ! user-defined physical criterion for refinement.
  !-------------------------------------------------
  integer ::i
  logical::refine_region_contains

  ! Authorize refinement if cell lies within region,
  ! otherwise unmark cell (no refinement outside region)
  if(r_refine(ilevel)>-1.0)then
     do i=1,ncell
        ok(i)=ok(i).and.refine_region_contains(xx(i,1:ndim),ilevel)
     end do
  endif

end subroutine geometry_refine
!############################################################
!############################################################
!############################################################
!############################################################
logical function refine_region_contains(xx,ilevel)
  use amr_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),dimension(1:ndim),intent(in)::xx
  real(dp)::er,xr,yr,zr,rr,xn,yn,zn,r,aa,bb

  er=exp_refine(ilevel)
  xr=x_refine(ilevel)
  yr=y_refine(ilevel)
  zr=z_refine(ilevel)
  rr=r_refine(ilevel)
  aa=a_refine(ilevel)
  bb=b_refine(ilevel)

  if(rr<=0.0d0 .or. er<=0.0d0 .or. aa<=0.0d0 .or. bb<=0.0d0)then
     refine_region_contains=.false.
     return
  end if

  xn=abs(xx(1)-xr)
  if(cosmo .and. xn>0.5d0)xn=1.0d0-xn
  xn=2.0d0*xn/rr
  yn=0.0d0
  zn=0.0d0
#if NDIM > 1
  yn=abs(xx(2)-yr)
  if(cosmo .and. yn>0.5d0)yn=1.0d0-yn
  yn=2.0d0*yn/(aa*rr)
#endif
#if NDIM > 2
  zn=abs(xx(3)-zr)
  if(cosmo .and. zn>0.5d0)zn=1.0d0-zn
  zn=2.0d0*zn/(bb*rr)
#endif
  if(er<10.0d0)then
     r=(xn**er+yn**er+zn**er)**(1.0d0/er)
  else
     r=max(xn,yn,zn)
  end if
  refine_region_contains=(r<1.0d0)

end function refine_region_contains
!############################################################
!############################################################
!############################################################
!############################################################
subroutine sub1_smooth_fine(ilevel, igrid,ngrid)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth
  integer::i,ncache,ngrid
  integer::igrid,ind
  integer,dimension(1:nvector)::ind_grid,ind_cell
  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do
  do ind=1,twotondim
     do i=1,ngrid
        ind_cell(i)=ICELL_OF(ind_grid(i),ind)
     end do
     do i=1,ngrid
        flag2(ind_cell(i))=0
     end do
  end do

end subroutine sub1_smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine sub2_smooth_fine(ilevel, igrid,ngrid,iflag)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ilevel,iflag
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth
  integer::i,ncache,ngrid
  integer::igrid,ind
  integer,dimension(1:nvector)::ind_grid,ind_cell
  iflag = 0
  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do
  do ind=1,twotondim
     do i=1,ngrid
        ind_cell(i)=ICELL_OF(ind_grid(i),ind)
     end do
     do i=1,ngrid
        if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=0
     end do
     do i=1,ngrid
        if(flag2(ind_cell(i))==1)then
           flag1(ind_cell(i))=1
           iflag=iflag+1
        end if
     end do
  end do
end subroutine sub2_smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine smooth_fine(ilevel)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth,mflag
  integer::i,j,ncache,ngrid,iflag,jflag
  integer::igrid,ind
  integer,dimension(1:3)::n_nbor
! integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  integer,target, allocatable, dimension(:,:)::Pind_grid,Pind_cell
  integer,pointer, dimension(:)::ind_grid,ind_cell
  integer mythread, nthreads,nwork,icount,jcount
  common /omp_threads/ mythread, nthreads
!$omp threadprivate(/omp_threads/)
  common /omp_smooth_fine/ ind_grid,ind_cell
!$omp threadprivate(/omp_smooth_fine/)


  
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return

!$omp parallel
  mythread = omp_get_thread_num()
  nthreads = omp_get_num_threads()
!$omp end parallel
  allocate(Pind_grid(1:nvector,0:nthreads-1), Pind_cell(1:nvector,0:nthreads-1))
!$omp parallel
  ind_grid => Pind_grid(:, mythread)
  ind_cell => Pind_cell(:, mythread)
!$omp end parallel


  n_nbor(1:3)=(/1,2,2/)
  flag1(0)=0
  ncache=active(ilevel)%ngrid
  ! Loop over steps
  do ismooth=1,ndim
     ! Initialize flag2 to 0
!$omp parallel do private(igrid,ngrid,i,ind)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           do i=1,ngrid
              ind_cell(i)=ICELL_OF(ind_grid(i),ind)
           end do
           do i=1,ngrid
              flag2(ind_cell(i))=0
           end do
        end do
     end do
     ! Count neighbors and set flag2 accordingly (cache + OpenMP)
!$omp parallel do private(igrid,ngrid,i,j,ind,igridn)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           do j=0,twondim
              igridn(i,j) = nbor_active_cache(j, igrid+i-1)
           end do
        end do
        do ind=1,twotondim
           call count_nbors(igridn,ind,n_nbor(ismooth),ngrid)
        end do
     end do
     ! Set flag1=1 for cells with flag2=1
     mflag = 0
!$omp parallel do private(igrid,ngrid,i,ind) reduction(+:mflag)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           do i=1,ngrid
              ind_cell(i)=ICELL_OF(ind_grid(i),ind)
           end do
           do i=1,ngrid
              if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=0
           end do
           do i=1,ngrid
              if(flag2(ind_cell(i))==1)then
                 flag1(ind_cell(i))=1
                 mflag=mflag+1
              end if
           end do
        end do
     end do
     nflag = nflag + mflag
     ! Update boundaries
     call make_virtual_fine_int(flag1(1),ilevel)
     if(simple_boundary)call make_boundary_flag(ilevel)
  end do
  ! End loop over steps
  deallocate(Pind_grid,Pind_cell)

end subroutine smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine count_nbors(igridn,ind,n_nbor,nn)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ind,nn,n_nbor
  integer,dimension(1:nvector,0:twondim)::igridn
  !----------------------------------------------------
  ! This routine computes the number of neighbors 
  ! for cell ind in grid igridn(:,0) for which flag1=1.
  ! The user must provide the neighboring grids index
  ! stored in igridn(:,:) and the threshold n_nbor
  ! If the number of flag1 neighbors exceeds n_nbor, 
  ! then cell is marked with flag2=1
  !----------------------------------------------------
  integer::i,in
  integer,dimension(1:nvector)::ind_cell,i_nbor
  integer,dimension(1:nvector,1:twondim)::indn
  ! Compute cell number
  do i=1,nn
     ind_cell(i)=ICELL_OF(igridn(i,0),ind)
  end do
  ! Gather neighbors
  call getnborcells(igridn,ind,indn,nn)
  ! Check if neighboring cell exists and count it as a flagged neighbor
  i_nbor(1:nn)=0
  do in=1,twondim
     do i=1,nn
        i_nbor(i)=i_nbor(i)+flag1(indn(i,in))
     end do
  end do
  ! flag2 cell if necessary
  do i=1,nn
     if(i_nbor(i)>=n_nbor)flag2(ind_cell(i))=1
  end do
end subroutine count_nbors
!############################################################
!############################################################
!############################################################
!############################################################
subroutine count_nbors2(igridn,ind,n_nbor,nn)
  use amr_commons
#include "amr_index.h"
  implicit none
  integer::ind,nn,n_nbor
  integer,dimension(1:nvector,0:twondim)::igridn
  !----------------------------------------------------
  ! This routine computes the number of neighbors 
  ! for cell ind in grid igridn(:,0) for which flag2=1.
  ! The user must provide the neighboring grids index
  ! stored in igridn(:,:) and the threshold n_nbor
  ! If the number of flag2 neighbors exceeds n_nbor, 
  ! then cell is marked with flag1=1
  !----------------------------------------------------
  integer::i,in
  integer,dimension(1:nvector)::ind_cell,i_nbor
  integer,dimension(1:nvector,1:twondim)::indn
  ! Compute cell number
  do i=1,nn
     ind_cell(i)=ICELL_OF(igridn(i,0),ind)
  end do
  ! Gather neighbors
  call getnborcells(igridn,ind,indn,nn)
  ! Check if neighboring cell exists and count it as a flagged neighbor
  i_nbor(1:nn)=0
  do in=1,twondim
     do i=1,nn
        i_nbor(i)=i_nbor(i)+flag2(indn(i,in))
     end do
  end do
  ! flag2 cell if necessary
  do i=1,nn
     if(i_nbor(i)>=n_nbor)flag1(ind_cell(i))=1
  end do
end subroutine count_nbors2
!############################################################
!############################################################
!############################################################
subroutine init_refmap
  use amr_commons
  implicit none

  integer::ilevel,ivar

  if(verbose)write(*,*)'Entering init_refmap'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call init_refmap_fine(ilevel)
     call make_virtual_fine_int(cpu_map2(1),ilevel)
  end do
  if(verbose)write(*,*)'Complete init_refmap'

end subroutine init_refmap
!############################################################
!############################################################
!############################################################
!############################################################
subroutine init_refmap_fine(ilevel)
  use amr_commons
  use hydro_commons
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::i,icell,igrid,ncache,ngrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::i1_lo,i1_hi,i2_lo,i2_hi
  integer::buf_count,info,nvar_in
  integer ,dimension(1:nvector)::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ::vv
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:nvar)::uu

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file1,ok_file2,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar
  integer,parameter::tag=1103
  integer::dummy_io,info2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  if(multiple)then
     filename=TRIM(initfile(ilevel))//'/dir_refmap/ic_refmap.00001'
     INQUIRE(file=filename,exist=ok_file2)
  else
     filename=TRIM(initfile(ilevel))//'/ic_refmap'
     INQUIRE(file=filename,exist=ok_file2)
  endif
  if (.not. ok_file2)then
     if(myid==1)write(*,*)'File '//TRIM(filename)//' not found'
     call clean_stop
  endif
     
  !-------------------------------------------------------------------------
  ! First step: compute level boundaries in terms of initial condition array
  !-------------------------------------------------------------------------
  if(ncache>0)then
     i1_min=n1(ilevel)+1; i1_max=0
     i2_min=n2(ilevel)+1; i2_max=0
     i3_min=n3(ilevel)+1; i3_max=0
     do ind=1,twotondim           
        do i=1,ncache
           igrid=active(ilevel)%igrid(i)
           xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
           xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
           xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
           xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
           xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
           xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
           i1_min=MIN(i1_min,int(xx1)+1)
           i1_max=MAX(i1_max,int(xx1)+1)
           i2_min=MIN(i2_min,int(xx2)+1)
           i2_max=MAX(i2_max,int(xx2)+1)
           i3_min=MIN(i3_min,int(xx3)+1)
           i3_max=MAX(i3_max,int(xx3)+1)
        end do
     end do
     error=.false.
     if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
     if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
     if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
     if(error) then
        write(*,*)'Some grid are outside initial conditions sub-volume'
        write(*,*)'for ilevel=',ilevel
        write(*,*)i1_min,i1_max
        write(*,*)i2_min,i2_max
        write(*,*)i3_min,i3_max
        write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
!jhshin1
!        call clean_stop
!jhshin2
     end if
  endif

  !-----------------------------------------
  ! Second step: read initial condition file
  !-----------------------------------------
  ! Allocate initial conditions array
  if(ncache>0)then
     allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
     ! A rank may own coarse background cells outside a zoom IC sub-volume.
     ! Those cells must retain the default (unrefined) map value; copying the
     ! full rank bounding box from init_plane used to read outside its bounds.
     init_array=0.0_dp
     i1_lo=max(1,i1_min); i1_hi=min(n1(ilevel),i1_max)
     i2_lo=max(1,i2_min); i2_hi=min(n2(ilevel),i2_max)
  endif
  allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))

  if(myid==1)write(*,*)'Reading file '//TRIM(filename)
  if(multiple)then
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif
     
     ilun=ncpu+myid+10
     open(ilun,file=filename,form='unformatted')
     rewind ilun
     read(ilun) ! skip first line
     do i3=1,n3(ilevel)
        read(ilun) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
        if(i3.ge.i3_min.and.i3.le.i3_max.and. &
             & i1_lo.le.i1_hi.and.i2_lo.le.i2_hi)then
           init_array(i1_lo:i1_hi,i2_lo:i2_hi,i3) = &
                & init_plane(i1_lo:i1_hi,i2_lo:i2_hi)
        end if
     end do
     close(ilun)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif

  else
     if(myid==1)then
        open(10,file=filename,form='unformatted')
        rewind 10
        read(10) ! skip first line
     endif
     do i3=1,n3(ilevel)
        if(myid==1)then
           read(10) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
        else
           init_plane=0.0
        endif
        buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
        call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
        if(ncache>0)then
           if(i3.ge.i3_min.and.i3.le.i3_max.and. &
                & i1_lo.le.i1_hi.and.i2_lo.le.i2_hi)then
              init_array(i1_lo:i1_hi,i2_lo:i2_hi,i3) = &
                   & init_plane(i1_lo:i1_hi,i2_lo:i2_hi)
           end if
        endif
     end do
     if(myid==1)close(10)
  endif
  
  if(ncache>0)then
     
     ! Loop over cells
!$omp parallel do private(ind,i,igrid,icell,xx1,xx2,xx3,i1,i2,i3)
     do ind=1,twotondim
        do i=1,ncache
           igrid=active(ilevel)%igrid(i)
           icell=ICELL_OF(igrid,ind)
           xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
           xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
           xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
           xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
           xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
           xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
           i1=int(xx1)+1
           i1=int(xx1)+1
           i2=int(xx2)+1
           i2=int(xx2)+1
           i3=int(xx3)+1
           i3=int(xx3)+1
           ! Scatter to corresponding primitive variable
           cpu_map2(icell)=int(init_array(i1,i2,i3))
        end do
     end do
     ! End loop over cells
  endif

  ! Deallocate initial conditions array
  if(ncache>0)deallocate(init_array)
  deallocate(init_plane) 
  
111 format('   Entering init_refmap_fine ',I2)

end subroutine init_refmap_fine
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine precompute_nbor_active(ilevel)
  use amr_commons
  use morton_keys
  use morton_hash
  implicit none
  integer, intent(in) :: ilevel
  integer :: ncache, igrid, igrid_amr, j
  integer :: l
  integer(8) :: nmax_x, nmax_y, nmax_z
  type(mkey_t) :: mkey, nkey

  ncache = active(ilevel)%ngrid
  if(ncache == 0) return

  if(allocated(nbor_active_cache)) deallocate(nbor_active_cache)
  allocate(nbor_active_cache(0:twondim, 1:ncache))

  l = ilevel
  nmax_x = int(nx, 8) * 2_8**(l-1)
  nmax_y = int(ny, 8) * 2_8**(l-1)
  nmax_z = int(nz, 8) * 2_8**(l-1)

!$omp parallel do private(igrid, igrid_amr, j, mkey, nkey)
  do igrid = 1, ncache
     igrid_amr = active(ilevel)%igrid(igrid)
     nbor_active_cache(0, igrid) = igrid_amr
     mkey = grid_to_morton(igrid_amr, l)
     do j = 1, twondim
        nkey = morton_neighbor(mkey, j, nmax_x, nmax_y, nmax_z)
        nbor_active_cache(j, igrid) = morton_hash_lookup(mort_table(l), nkey)
     end do
  end do

end subroutine precompute_nbor_active
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine cleanup_nbor_active()
  use amr_commons
  implicit none
  if(allocated(nbor_active_cache)) deallocate(nbor_active_cache)
end subroutine cleanup_nbor_active
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine compute_fpr_m_refine_eff(ilevel)
  !-----------------------------------------------------------------
  ! Compute FPR-adjusted effective m_refine for level ilevel.
  ! If dx_phys < dr_proper, increase threshold by (dr_proper/dx_phys)^3
  ! so that refinement is naturally suppressed when physical resolution
  ! exceeds the target (Gnedin 2016, ApJ 821, 50, Appendix A).
  ! When FPR is disabled (dr_proper=0), m_refine_eff = m_refine.
  !-----------------------------------------------------------------
  use amr_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp)::dx,dx_loc,scale,dx_phys_kpc,fpr_factor
  integer::nx_loc

  ! Default: no modification
  m_refine_eff(ilevel) = m_refine(ilevel)

  ! FPR only active for cosmo runs with dr_proper > 0
  if(.not. cosmo .or. dr_proper <= 0.0d0) return
  if(m_refine(ilevel) <= 0.0d0) return

  ! Compute physical cell size in kpc
  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale
  ! dx_phys [kpc] = dx_code * aexp * boxlen_ini [h^-1 Mpc] * 1000 / (h0/100)
  dx_phys_kpc = dx_loc * aexp * boxlen_ini * 1000.0d0 / (h0 / 100.0d0)

  if(dx_phys_kpc < dr_proper) then
     fpr_factor = (dr_proper / dx_phys_kpc)**3
     m_refine_eff(ilevel) = m_refine(ilevel) * fpr_factor
  end if

end subroutine compute_fpr_m_refine_eff
