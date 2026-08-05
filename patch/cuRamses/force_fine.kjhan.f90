!#########################################################
!#########################################################
! Dynamic OMP/CUDA hybrid dispatch for force_fine
! Overrides patch/Horizon5-master-2/force_fine.kjhan.f90
!#########################################################
!#########################################################

#ifdef HYDRO_CUDA
!=========================================================
! Module: GPU state for gradient_phi superbatch
!=========================================================
module force_hybrid_commons
  use amr_parameters, only: dp, nvector, ndim, twotondim, twondim
  implicit none

  integer, parameter :: FORCE_SUPER_SIZE = 4096

  type force_gpu_state_t
     ! phi stencil buffer: phi_buf(4, cap*24)
     ! Layout: phi_buf(1:4, slot) where slot=(g-1)*24+(ind-1)*3+idim
     real(dp), allocatable :: phi_buf(:,:)
     ! Force output buffer: f_buf(cap*24)
     real(dp), allocatable :: f_buf(:)
     ! Cell index buffer for scatter: cell_buf(cap, 8)
     integer, allocatable :: cell_buf(:,:)
     integer :: off = 0     ! number of grids accumulated
     integer :: cap = 0     ! current capacity
  end type

  type(force_gpu_state_t), allocatable, save, target :: force_gstates(:)
  integer, save :: force_hybrid_inited = 0

contains

  subroutine force_gstate_ensure(gs, needed)
    type(force_gpu_state_t), intent(inout) :: gs
    integer, intent(in) :: needed
    if (needed <= gs%cap) return
    gs%cap = max(needed, FORCE_SUPER_SIZE)
    if (allocated(gs%phi_buf)) deallocate(gs%phi_buf)
    if (allocated(gs%f_buf)) deallocate(gs%f_buf)
    if (allocated(gs%cell_buf)) deallocate(gs%cell_buf)
    allocate(gs%phi_buf(4, gs%cap * 24))
    allocate(gs%f_buf(gs%cap * 24))
    allocate(gs%cell_buf(gs%cap, 8))
  end subroutine

end module force_hybrid_commons
#endif

!#########################################################
!#########################################################
subroutine force_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use dark_energy_commons, only: cosmo_poisson_fourpi
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  ! Work arrays (thread-private via OMP private clause)
  integer ,dimension(1:nvector)::ind_grid_w,ind_cell_w
  real(dp),dimension(1:nvector,1:ndim)::xx_w,ff_w

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then

     ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_grid_w,ind_cell_w,xx_w,ff_w) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do
           do idim=1,ndim
              do i=1,ngrid
                 xx_w(i,idim)=xg(ind_grid_w(i),idim)+xc(ind,idim)
              end do
           end do
           do idim=1,ndim
              do i=1,ngrid
                 xx_w(i,idim)=(xx_w(i,idim)-skip_loc(idim))*scale
              end do
           end do
           call gravana(xx_w,ff_w,dx_loc,ngrid)
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell_w(i),idim)=ff_w(i,idim)
              end do
           end do
        end do
     end do

     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  !------------------------------
  ! Compute gradient of potential
  !------------------------------
  else
     call make_boundary_phi(ilevel)

     ncache=active(ilevel)%ngrid
#ifdef HYDRO_CUDA
     call force_gradient_hybrid(ilevel, icount, ncache)
#else
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        call gradient_phi(ilevel,igrid,ngrid,icount)
     end do
#endif

     ! Apply MOND (QUMOND) correction to Newtonian force
     if(use_mond .and. mond_type == 0) then
        call apply_mond_force(ilevel)
     end if

     ! Apply Coupled Dark Energy force enhancement
     if(use_coupled_de) then
        call apply_coupled_de_force(ilevel)
     end if

     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  rho_loc =0.0; rho_all =0.0
  epot_loc=0.0; epot_all=0.0
  fourpi=cosmo_poisson_fourpi(aexp,1.0d0)
  fact=-dx_loc**ndim/fourpi/2.0D0

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_grid_w,ind_cell_w) &
!$omp& reduction(+:epot_loc) reduction(max:rho_loc)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              if(son(ind_cell_w(i))==0)then
                 epot_loc=epot_loc+fact*f(ind_cell_w(i),idim)**2
              end if
           end do
        end do
        do i=1,ngrid
           rho_loc=MAX(rho_loc,dble(abs(rho(ind_cell_w(i)))))
        end do
     end do
  end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc =rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

111 format('   Entering force_fine for level ',I2)
end subroutine force_fine

#ifdef HYDRO_CUDA
!#########################################################
!#########################################################
! Hybrid dispatch for gradient_phi
!#########################################################
!#########################################################
subroutine force_gradient_hybrid(ilevel, icount, ncache)
  use amr_commons
  use poisson_commons
  use force_hybrid_commons
  use cuda_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: ilevel, icount, ncache

  integer :: igrid, ngrid, stream_slot
  type(force_gpu_state_t), pointer :: gs

  ! First-call: allocate GPU states
  if (force_hybrid_inited == 0) then
     allocate(force_gstates(0:7))
     force_hybrid_inited = 1
  end if

  !$omp parallel private(igrid, ngrid, stream_slot, gs)
  stream_slot = cuda_acquire_stream_c()
  if (stream_slot >= 0) then
     gs => force_gstates(stream_slot)
     call force_gstate_ensure(gs, FORCE_SUPER_SIZE)
     gs%off = 0
  end if

  !$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     if (stream_slot >= 0) then
        call force_gpu_gather_batch(gs, ilevel, icount, igrid, ngrid, stream_slot)
     else
        call gradient_phi(ilevel, igrid, ngrid, icount)
     end if
  end do
  !$omp end do nowait

  if (stream_slot >= 0) then
     if (gs%off > 0) call force_gpu_flush_scatter(gs, stream_slot, ilevel)
     call cuda_release_stream_c(stream_slot)
  end if
  !$omp end parallel

end subroutine force_gradient_hybrid

!#########################################################
! GPU gather: fill phi stencil buffer for superbatch
!#########################################################
subroutine force_gpu_gather_batch(gs, ilevel, icount, igrid_start, ngrid, stream_slot)
  use amr_commons
  use poisson_commons
  use morton_hash
  use force_hybrid_commons
  implicit none

  type(force_gpu_state_t), intent(inout) :: gs
  integer, intent(in) :: ilevel, icount, igrid_start, ngrid, stream_slot

  integer :: i, ind, idim, nx_loc, iskip
  integer :: id1, id2, id3, id4
  integer :: ig1, ig2, ig3, ig4
  integer :: ih1, ih2, ih3, ih4
  real(dp) :: dx, scale, dx_loc
  integer, dimension(1:3,1:4,1:8) :: ggg, hhh

  integer, dimension(1:nvector) :: ind_grid
  integer, dimension(1:nvector, 0:twondim) :: igridn
  integer, dimension(1:nvector, 1:ndim) :: ind_left, ind_right
  real(dp), dimension(1:nvector, 1:twotondim, 1:ndim) :: phi_left, phi_right

  integer :: g, slot, off

  ! Load grid indices
  do i = 1, ngrid
     ind_grid(i) = active(ilevel)%igrid(igrid_start + i - 1)
  end do

  ! Mesh size
  dx = 0.5D0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  ! Stencil lookup tables
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i = 1, ngrid
     igridn(i, 0) = ind_grid(i)
  end do
  do idim = 1, ndim
     do i = 1, ngrid
        igridn(i, 2*idim-1) = morton_nbor_grid(ind_grid(i), ilevel, 2*idim-1)
        igridn(i, 2*idim  ) = morton_nbor_grid(ind_grid(i), ilevel, 2*idim  )
        ind_left (i, idim) = morton_nbor_cell(ind_grid(i), ilevel, 2*idim-1)
        ind_right(i, idim) = morton_nbor_cell(ind_grid(i), ilevel, 2*idim  )
     end do
  end do

  ! Interpolate phi from upper level
  if (ilevel > levelmin) then
     do idim = 1, ndim
        call interpol_phi(ind_left(1,idim), phi_left(1,1,idim), ngrid, ilevel, icount)
        call interpol_phi(ind_right(1,idim), phi_right(1,1,idim), ngrid, ilevel, icount)
     end do
  end if

  ! Fill phi_buf and cell_buf
  off = gs%off
  do ind = 1, twotondim
     iskip = ncoarse + (ind - 1) * ngridmax
     do i = 1, ngrid
        gs%cell_buf(off + i, ind) = iskip + ind_grid(i)
     end do

     do idim = 1, ndim
        id1 = hhh(idim,1,ind); ig1 = ggg(idim,1,ind); ih1 = ncoarse+(id1-1)*ngridmax
        id2 = hhh(idim,2,ind); ig2 = ggg(idim,2,ind); ih2 = ncoarse+(id2-1)*ngridmax
        id3 = hhh(idim,3,ind); ig3 = ggg(idim,3,ind); ih3 = ncoarse+(id3-1)*ngridmax
        id4 = hhh(idim,4,ind); ig4 = ggg(idim,4,ind); ih4 = ncoarse+(id4-1)*ngridmax

        do i = 1, ngrid
           g = off + i
           slot = (g - 1) * 24 + (ind - 1) * 3 + idim

           ! phi1
           if (igridn(i, ig1) > 0) then
              gs%phi_buf(1, slot) = phi(igridn(i, ig1) + ih1)
           else
              gs%phi_buf(1, slot) = phi_left(i, id1, idim)
           end if
           ! phi2
           if (igridn(i, ig2) > 0) then
              gs%phi_buf(2, slot) = phi(igridn(i, ig2) + ih2)
           else
              gs%phi_buf(2, slot) = phi_right(i, id2, idim)
           end if
           ! phi3
           if (igridn(i, ig3) > 0) then
              gs%phi_buf(3, slot) = phi(igridn(i, ig3) + ih3)
           else
              gs%phi_buf(3, slot) = phi_left(i, id3, idim)
           end if
           ! phi4
           if (igridn(i, ig4) > 0) then
              gs%phi_buf(4, slot) = phi(igridn(i, ig4) + ih4)
           else
              gs%phi_buf(4, slot) = phi_right(i, id4, idim)
           end if
        end do
     end do
  end do

  gs%off = off + ngrid

  ! Flush if buffer is near full
  if (gs%off + nvector > gs%cap) then
     call force_gpu_flush_scatter(gs, stream_slot, ilevel)
  end if

end subroutine force_gpu_gather_batch

!#########################################################
! GPU flush: launch kernel, scatter results
!#########################################################
subroutine force_gpu_flush_scatter(gs, stream_slot, ilevel)
  use amr_commons
  use poisson_commons
  use force_hybrid_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none

  type(force_gpu_state_t), intent(inout) :: gs
  integer, intent(in) :: stream_slot, ilevel

  real(dp) :: dx, a, b
  integer :: g, ind, idim, slot

  if (gs%off == 0) return

  dx = 0.5D0**ilevel
  a = 0.50D0 * 4.0D0 / 3.0D0 / dx
  b = 0.25D0 * 1.0D0 / 3.0D0 / dx

  ! Launch GPU kernel
  call gradient_phi_cuda_async_f(gs%phi_buf, gs%f_buf, a, b, &
       int(gs%off, c_int), int(stream_slot, c_int))
  call gradient_phi_cuda_sync_f(int(stream_slot, c_int))

  ! Scatter force values to global f() array
  do g = 1, gs%off
     do ind = 1, twotondim
        do idim = 1, ndim
           slot = (g - 1) * 24 + (ind - 1) * 3 + idim
           f(gs%cell_buf(g, ind), idim) = gs%f_buf(slot)
        end do
     end do
  end do

  gs%off = 0

end subroutine force_gpu_flush_scatter
#endif

!#########################################################
!#########################################################
! gradient_phi (unchanged from base version)
!#########################################################
!#########################################################
subroutine gradient_phi(ilevel,igrid,ngrid,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc,igrid
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim)::igridn
  real(dp),dimension(1:nvector)::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::phi_left,phi_right

  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        igridn(i,2*idim-1)=morton_nbor_grid(ind_grid(i),ilevel,2*idim-1)
        igridn(i,2*idim  )=morton_nbor_grid(ind_grid(i),ilevel,2*idim  )
        ind_left (i,idim)=morton_nbor_cell(ind_grid(i),ilevel,2*idim-1)
        ind_right(i,idim)=morton_nbor_cell(ind_grid(i),ilevel,2*idim  )
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=phi(igridn(i,ig1)+ih1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=phi(igridn(i,ig2)+ih2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=phi(igridn(i,ig3)+ih3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=phi(igridn(i,ig4)+ih4)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           f(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                &             -b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine gradient_phi
!#########################################################
!#########################################################
subroutine apply_mond_force(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  !-------------------------------------------------------
  ! QUMOND: multiply Newtonian force by nu(|g_N+g_ext|/a0)
  ! With EFE: f_d = nu * f_d + (nu-1) * g_ext_d
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale_a,a0_code,gnorm,x,nu
  real(dp)::g_ext_code(1:3)
  integer,dimension(1:nvector)::ind_cell

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  ! Convert a0 and g_ext from CGS to code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_a = scale_l / scale_t**2
  a0_code = a0_mond / scale_a
  g_ext_code(1:3) = g_ext_mond(1:3) / scale_a

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_cell, &
!$omp&  gnorm,x,nu) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+active(ilevel)%igrid(igrid+i-1)
        end do

        do i=1,ngrid
           gnorm = sqrt((f(ind_cell(i),1)+g_ext_code(1))**2 &
                      + (f(ind_cell(i),2)+g_ext_code(2))**2 &
                      + (f(ind_cell(i),3)+g_ext_code(3))**2)

           if(gnorm > 1d-30*a0_code) then
              x = gnorm / a0_code
              ! Compute nu-function
              if(mond_mu_type == 1) then
                 nu = 0.5d0 + 0.5d0*sqrt(1d0 + 4d0/x)
              else
                 nu = sqrt(0.5d0 + 0.5d0*sqrt(1d0 + 4d0/(x*x)))
              end if
              ! f_d = nu * f_d + (nu-1) * g_ext_d
              do idim=1,ndim
                 f(ind_cell(i),idim) = nu * f(ind_cell(i),idim) &
                      + (nu - 1d0) * g_ext_code(idim)
              end do
           end if
        end do
     end do
  end do

end subroutine apply_mond_force
!#########################################################
!#########################################################
subroutine compute_mond_phantom_density(ilevel, is_aqual)
  use amr_commons
  use poisson_commons
  use morton_hash
  use dark_energy_commons, only: cosmo_poisson_fourpi
  implicit none
  integer,intent(in)::ilevel
  logical,intent(in)::is_aqual
  !-------------------------------------------------------
  ! Compute phantom density and add to rho() in-place.
  !
  ! QUMOND (is_aqual=.false.):
  !   rho_ph = -div[(nu-1)*(f+g_ext)] / fourpi  (nu-1 > 0)
  ! AQUAL (is_aqual=.true.):
  !   rho_ph = +div[(mu-1)*(f+g_ext)] / fourpi  (mu-1 < 0)
  !
  ! Uses f() (force, already boundary-exchanged)
  ! with 6-neighbor centered difference for divergence.
  ! Cells at coarse-fine boundaries (igridn=0) are skipped.
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale_a,a0_code
  real(dp)::g_ext_code(1:3)
  real(dp)::dx,scale,dx_loc,fourpi
  integer::nx_loc
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::gnorm_nb,x_nb,func_m1
  real(dp)::h_right,h_left,div_h
  logical::skip_cell

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  ! Unit conversion
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_a = scale_l / scale_t**2
  a0_code = a0_mond / scale_a
  g_ext_code(1:3) = g_ext_mond(1:3) / scale_a

  ! Mesh size
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! 4piG factor (same convention as force_fine)
  fourpi=cosmo_poisson_fourpi(aexp,1.0d0)

  ! Neighbor lookup tables (left=1, right=2 per dimension)
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  gnorm_nb,x_nb,func_m1, &
!$omp&  h_right,h_left,div_h,skip_cell) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids (6 faces)
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do

        do i=1,ngrid
           ! Check if all neighbors exist (skip at coarse-fine boundaries)
           skip_cell=.false.
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              if(igridn_w(i,ig_left)==0 .or. igridn_w(i,ig_right)==0) then
                 skip_cell=.true.
                 exit
              end if
           end do
           if(skip_cell) cycle

           ! Compute divergence of (func-1)*(f+g_ext)
           div_h = 0d0
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax

              ind_nb_left =igridn_w(i,ig_left )+ih_left
              ind_nb_right=igridn_w(i,ig_right)+ih_right

              ! Left neighbor
              gnorm_nb = sqrt((f(ind_nb_left,1)+g_ext_code(1))**2 &
                            + (f(ind_nb_left,2)+g_ext_code(2))**2 &
                            + (f(ind_nb_left,3)+g_ext_code(3))**2)
              x_nb = gnorm_nb / a0_code
              if(is_aqual) then
                 call get_mond_mu_minus1(x_nb, mond_mu_type, func_m1)
              else
                 call get_mond_nu_minus1(x_nb, mond_mu_type, func_m1)
              end if
              h_left = func_m1 * (f(ind_nb_left, idim) + g_ext_code(idim))

              ! Right neighbor
              gnorm_nb = sqrt((f(ind_nb_right,1)+g_ext_code(1))**2 &
                            + (f(ind_nb_right,2)+g_ext_code(2))**2 &
                            + (f(ind_nb_right,3)+g_ext_code(3))**2)
              x_nb = gnorm_nb / a0_code
              if(is_aqual) then
                 call get_mond_mu_minus1(x_nb, mond_mu_type, func_m1)
              else
                 call get_mond_nu_minus1(x_nb, mond_mu_type, func_m1)
              end if
              h_right = func_m1 * (f(ind_nb_right, idim) + g_ext_code(idim))

              div_h = div_h + (h_right - h_left)
           end do

           ! QUMOND: rho_ph = -div[(nu-1)*(f+g_ext)] / (2*dx) / fourpi
           ! AQUAL:  rho_ph = +div[(mu-1)*(f+g_ext)] / (2*dx) / fourpi
           if(is_aqual) then
              rho(ind_cell_w(i)) = rho(ind_cell_w(i)) + div_h / (2d0*dx_loc) / fourpi
           else
              rho(ind_cell_w(i)) = rho(ind_cell_w(i)) - div_h / (2d0*dx_loc) / fourpi
           end if

        end do  ! i
     end do  ! ind
  end do  ! igrid

end subroutine compute_mond_phantom_density
!#########################################################
!#########################################################
subroutine get_mond_nu_minus1(x, mu_type, num1)
  use amr_parameters, only: dp
  implicit none
  real(dp),intent(in)::x
  integer,intent(in)::mu_type
  real(dp),intent(out)::num1
  ! Returns nu(x) - 1 for MOND interpolation function
  if(x < 1d-30) then
     num1 = 0d0
     return
  end if
  if(mu_type == 1) then
     ! Simple: mu=x/(1+x) -> nu = 0.5*(1+sqrt(1+4/x))
     num1 = 0.5d0*(-1d0 + sqrt(1d0 + 4d0/x))
  else
     ! Standard: mu=x/sqrt(1+x^2) -> nu = sqrt(0.5+0.5*sqrt(1+4/x^2))
     num1 = sqrt(0.5d0 + 0.5d0*sqrt(1d0 + 4d0/(x*x))) - 1d0
  end if
end subroutine get_mond_nu_minus1
!#########################################################
!#########################################################
subroutine get_mond_mu_minus1(x, mu_type, mum1)
  use amr_parameters, only: dp
  implicit none
  real(dp),intent(in)::x
  integer,intent(in)::mu_type
  real(dp),intent(out)::mum1
  ! Returns mu(x) - 1 for MOND interpolation function
  ! mu(x) - 1 < 0 always (mu < 1 in deep-MOND)
  if(x < 1d-30) then
     mum1 = -1d0   ! mu(0) = 0
     return
  end if
  if(mu_type == 1) then
     ! Simple: mu=x/(1+x) -> mu-1 = -1/(1+x)
     mum1 = -1d0 / (1d0 + x)
  else
     ! Standard: mu=x/sqrt(1+x^2) -> mu-1 = x/sqrt(1+x^2) - 1
     mum1 = x / sqrt(1d0 + x*x) - 1d0
  end if
end subroutine get_mond_mu_minus1
!#########################################################
!#########################################################
subroutine aqual_iterate(ilevel, icount)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel, icount
  !-------------------------------------------------------
  ! AQUAL Phase 2: iterative fixed-point solver
  ! Solves: div[mu(|grad(Phi)|/a0) * grad(Phi)] = 4piG*rho
  ! by iterating:
  !   1. Compute AQUAL phantom density from current forces
  !   2. Re-solve Poisson with rho_bary + rho_phantom
  !   3. Re-compute forces
  !   4. Check convergence: max|f_new - f_old| / max|f_new|
  !-------------------------------------------------------
  integer :: iter, ncache, ind, iskip, i, idx, idim, info
  real(dp) :: delta_f_local, f_norm_local
  real(dp) :: delta_f_global, f_norm_global, rel_change
  logical :: converged

  real(dp), allocatable :: rho_bary(:)
  real(dp), allocatable :: f_old(:,:)
  integer :: ncell_active
  integer, allocatable :: ind_cell_list(:)

  ncache = active(ilevel)%ngrid
  if(ncache == 0) return

  ncell_active = ncache * twotondim
  allocate(rho_bary(1:ncell_active))
  allocate(f_old(1:ncell_active, 1:ndim))
  allocate(ind_cell_list(1:ncell_active))

  ! Build cell list and save rho_bary
  idx = 0
  do ind = 1, twotondim
     iskip = ncoarse + (ind-1) * ngridmax
     do i = 1, ncache
        idx = idx + 1
        ind_cell_list(idx) = iskip + active(ilevel)%igrid(i)
        rho_bary(idx) = rho(ind_cell_list(idx))
     end do
  end do

  ! Save current forces (Newtonian) as f_old
  do idx = 1, ncell_active
     do idim = 1, ndim
        f_old(idx, idim) = f(ind_cell_list(idx), idim)
     end do
  end do

  converged = .false.
  do iter = 1, n_iter_mond

     ! (a) Restore rho to baryonic
     do idx = 1, ncell_active
        rho(ind_cell_list(idx)) = rho_bary(idx)
     end do

     ! (b) Compute AQUAL phantom density (mu-1, current forces)
     call compute_mond_phantom_density(ilevel, .true.)
     call make_virtual_fine_dp(rho(1), ilevel)

     ! (c) Re-solve Poisson
     if(ilevel > levelmin) then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel, icount)
        else
           call multigrid_fine(ilevel, icount)
        end if
     else
        call multigrid_fine(levelmin, icount)
     end if

     ! (d) Re-compute forces
     call force_fine(ilevel, icount)

     ! (e) Convergence check
     delta_f_local = 0d0
     f_norm_local = 0d0
!$omp parallel do private(idx, idim) &
!$omp& reduction(max:delta_f_local, f_norm_local)
     do idx = 1, ncell_active
        do idim = 1, ndim
           delta_f_local = max(delta_f_local, &
                abs(f(ind_cell_list(idx), idim) - f_old(idx, idim)))
           f_norm_local = max(f_norm_local, &
                abs(f(ind_cell_list(idx), idim)))
        end do
     end do

     ! Save forces for next convergence check
     do idx = 1, ncell_active
        do idim = 1, ndim
           f_old(idx, idim) = f(ind_cell_list(idx), idim)
        end do
     end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(delta_f_local, delta_f_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(f_norm_local, f_norm_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     delta_f_global = delta_f_local
     f_norm_global = f_norm_local
#endif

     if(f_norm_global > 0d0) then
        rel_change = delta_f_global / f_norm_global
     else
        rel_change = 0d0
     end if

     if(myid == 1) then
        write(*,'(A,I3,A,ES12.4)') &
             ' AQUAL iter ', iter, ': delta_f/f_norm = ', rel_change
     end if

     if(rel_change < mond_eps) then
        converged = .true.
        if(myid == 1) write(*,'(A,I3,A)') &
             ' AQUAL converged in ', iter, ' iterations'
        exit
     end if
  end do

  if(.not. converged .and. myid == 1) then
     write(*,'(A,I3,A,ES12.4)') &
          ' WARNING: AQUAL did NOT converge after ', n_iter_mond, &
          ' iters, rel_change=', rel_change
  end if

  deallocate(rho_bary, f_old, ind_cell_list)
end subroutine aqual_iterate
!#########################################################
!#########################################################
! f(R) Hu-Sawicki gravity + nDGP gravity solvers
!#########################################################
!#########################################################

!=========================================================
! compute_fifth_force: add gradient of scalar_gr to f()
! factor = -0.5 for f(R) and nDGP; the nDGP 1/beta
! coupling is already part of the solved field equation.
!=========================================================
subroutine compute_fifth_force(ilevel, factor)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::factor
  !-------------------------------------------------------
  ! Compute F5_d = factor * (scalar_gr(right)-scalar_gr(left))/(2*dx)
  ! and add to f(icell, idim).
  ! Uses simple 2-point centered difference (same as rho divergence).
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc
  integer::nx_loc
  real(dp)::grad_u,u_cen,u_left,u_right

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ! NOTE: no early return on ncache==0 — the final make_virtual_fine_dp
  ! must be entered by every rank (matched communication)
  ncache=active(ilevel)%ngrid

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Neighbor lookup tables (left=1, right=2)
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right,grad_u,u_cen,u_left,u_right) schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do

        do i=1,ngrid
           u_cen = scalar_gr(ind_cell_w(i))
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax

              ! Same-level neighbours are already cached above.  The old
              ! path decoded the cell Morton coordinate and performed a
              ! hash lookup six times per active cell even on a complete
              ! uniform level.  Retain parent-CIC only for a genuinely
              ! missing AMR neighbour.
              if(igridn_w(i,ig_left)>0) then
                 u_left=scalar_gr(igridn_w(i,ig_left)+ih_left)
              else
                 call scalar_sample_axis( &
                      & ind_grid_w(i),ind,ilevel,idim,-1,u_cen,u_left)
              end if
              if(igridn_w(i,ig_right)>0) then
                 u_right=scalar_gr(igridn_w(i,ig_right)+ih_right)
              else
                 call scalar_sample_axis( &
                      & ind_grid_w(i),ind,ilevel,idim,1,u_cen,u_right)
              end if
              grad_u=(u_right-u_left)/(2d0*dx_loc)
              f(ind_cell_w(i),idim) = f(ind_cell_w(i),idim) + factor * grad_u
           end do
        end do
     end do
  end do

  ! Update MPI virtual boundaries: force_fine synced f BEFORE the
  ! fifth force was added; without this, reception cells keep the
  ! Newtonian-only force and boundary particles miss F5
  do idim=1,ndim
     call make_virtual_fine_dp(f(1,idim),ilevel)
  end do

end subroutine compute_fifth_force

!=========================================================
! scalar_lookup_cell: look up a scalar cell by its integer
! coordinates at an AMR level.  Coordinates are periodic.
!=========================================================
subroutine scalar_lookup_cell(ilevel, ix_in, iy_in, iz_in, value, found)
  use amr_commons
  use poisson_commons
  use morton_hash
  use morton_keys, only: mkey_t, morton_encode
  implicit none
  integer,intent(in)::ilevel
  integer(8),intent(in)::ix_in,iy_in,iz_in
  real(dp),intent(out)::value
  logical,intent(out)::found
  integer(8)::ix,iy,iz,ncx,ncy,ncz,gx,gy,gz
  integer::igrid,ind,icell
  type(mkey_t)::key

  found=.false.
  value=0d0
  if(.not. allocated(mort_table)) return
  if(ilevel < 1 .or. ilevel > size(mort_table)) return

  ncx=int(nx,8)*2_8**ilevel
  ncy=int(ny,8)*2_8**ilevel
  ncz=int(nz,8)*2_8**ilevel
  ix=modulo(ix_in,ncx)
  iy=modulo(iy_in,ncy)
  iz=modulo(iz_in,ncz)
  gx=ix/2_8
  gy=iy/2_8
  gz=iz/2_8
  key=morton_encode(gx,gy,gz)
  igrid=morton_hash_lookup(mort_table(ilevel),key)
  if(igrid <= 0) return

  ind=1+int(modulo(ix,2_8))+2*int(modulo(iy,2_8))+4*int(modulo(iz,2_8))
  icell=ncoarse+(ind-1)*ngridmax+igrid
  value=scalar_gr(icell)
  found=.true.
end subroutine scalar_lookup_cell

!=========================================================
! scalar_sample_offset: sample a same-level scalar neighbor.
! If that fine cell does not exist, impose a CIC-interpolated
! Dirichlet value from the parent level instead of the old
! zero-gradient closure.  This is also used for diagonal Hessian
! samples in Vainshtein operators.
!=========================================================
subroutine scalar_sample_offset(igrid, ind, ilevel, ox, oy, oz, fallback, value)
  use amr_commons
  use morton_hash
  use morton_keys, only: mkey_t, grid_to_morton, morton_decode
  implicit none
  integer,intent(in)::igrid,ind,ilevel,ox,oy,oz
  real(dp),intent(in)::fallback
  real(dp),intent(out)::value
  integer(8)::gx,gy,gz,cx,cy,cz,tx,ty,tz,x0,y0,z0
  integer::bx,by,bz
  real(dp)::wx,wy,wz,w,val
  logical::found
  type(mkey_t)::key

  key=grid_to_morton(igrid,ilevel)
  call morton_decode(key,gx,gy,gz)
  cx=2_8*gx+int(mod(ind-1,2),8)
  cy=2_8*gy+int(mod((ind-1)/2,2),8)
  cz=2_8*gz+int((ind-1)/4,8)
  tx=cx+int(ox,8)
  ty=cy+int(oy,8)
  tz=cz+int(oz,8)

  call scalar_lookup_cell(ilevel,tx,ty,tz,value,found)
  if(found) return
  if(ilevel <= 1) then
     value=fallback
     return
  end if

  if(modulo(tx,2_8)==0_8) then
     x0=tx/2_8-1_8
     wx=0.75d0
  else
     x0=(tx-1_8)/2_8
     wx=0.25d0
  end if
  if(modulo(ty,2_8)==0_8) then
     y0=ty/2_8-1_8
     wy=0.75d0
  else
     y0=(ty-1_8)/2_8
     wy=0.25d0
  end if
  if(modulo(tz,2_8)==0_8) then
     z0=tz/2_8-1_8
     wz=0.75d0
  else
     z0=(tz-1_8)/2_8
     wz=0.25d0
  end if

  value=0d0
  do bz=0,1
     do by=0,1
        do bx=0,1
           call scalar_lookup_cell(ilevel-1,x0+bx,y0+by,z0+bz,val,found)
           if(.not. found) then
              value=fallback
              return
           end if
           w=merge(wx,1d0-wx,bx==1)*merge(wy,1d0-wy,by==1) &
                & *merge(wz,1d0-wz,bz==1)
           value=value+w*val
        end do
     end do
  end do
end subroutine scalar_sample_offset

subroutine scalar_sample_axis(igrid,ind,ilevel,idim,side,fallback,value)
  use amr_parameters, only: dp
  implicit none
  integer,intent(in)::igrid,ind,ilevel,idim,side
  real(dp),intent(in)::fallback
  real(dp),intent(out)::value
  integer::ox,oy,oz

  ox=0; oy=0; oz=0
  select case(idim)
  case(1)
     ox=side
  case(2)
     oy=side
  case(3)
     oz=side
  end select
  call scalar_sample_offset(igrid,ind,ilevel,ox,oy,oz,fallback,value)
end subroutine scalar_sample_axis

!=========================================================
! Stop a strict scalar solve with a non-zero scheduler status.
! The cuRamses VPATH clean_stop finalizes MPI and executes a
! status-zero STOP, which can make a failed production run look
! successful to Slurm.
!=========================================================
subroutine scalar_solver_abort
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  call MPI_ABORT(MPI_COMM_WORLD,92,info)
#else
  stop 92
#endif
end subroutine scalar_solver_abort

#ifdef USE_FFTW
!#########################################################
!  FFT ACCELERATION FOR THE SCALAR SOLVERS (domain level)
!
!  Plain level GS cannot converge wavelengths >> dx in
!  n_iter sweeps. On the fully refined periodic domain
!  level we solve the (linearized) equation spectrally:
!    - f(R)/symmetron/dilaton: Newton step
!        (lap - m2bar) du = -residual
!      with m2bar the level-averaged Jacobian mass term
!    - nDGP/galileon: operator splitting: solve the local
!      cell quadratic for L = lap(phi) (branch with L->S
!      as coeff->0; clamped at the extremum if the
!      discriminant is negative), then invert lap du = L - lap(phi)
!  The rhs is staged in scalar_gr_old (otherwise dead
!  between solves); the correction is ADDED to scalar_gr.
!  Discrete 7-point eigenvalues match the GS stencil.
!  Strategy: retain the single-rank local FFT reference path through
!  256^3.  Multi-rank runs reuse the production base-Poisson FFTW-MPI+OMP
!  engine, including its in-place transposed plan and sparse slab exchange;
!  the distributed path has no physics-motivated mesh-size ceiling.
!#########################################################
logical function level_fft_ok(ilevel)
  use amr_commons
  implicit none
  integer,intent(in)::ilevel
  integer(kind=8)::ntot
  level_fft_ok = .false.
  if(ilevel /= levelmin) return
  ntot = int(nx,8)*int(ny,8)*int(nz,8)*(int(2,8)**(3*ilevel))
  ! The local path replicates the complete transform and remains capped at
  ! 256^3.  Multi-rank jobs use the shared distributed Poisson FFT engine.
  if(ncpu <= 1) then
     if(ntot > 16777216_8) return
  end if
  level_fft_ok = .true.
end function level_fft_ok

!=========================================================
! level_fft_helmholtz: solve (lap - m2) x = b on the fully
! refined periodic level; b read from scalar_gr_old, x is
! ADDED into scalar_gr.  relax multiplies the correction;
! step_frac>0 additionally applies a pointwise
! trust region |du| <= step_frac*|u|.  This is required for
! fields with a one-sided physical branch (f_R<0, dilaton>0):
! a spatially averaged Newton mass is only a preconditioner and
! its unrestricted correction can cross the singular branch.
!=========================================================
subroutine level_fft_helmholtz(ilevel, m2, step_frac, relax)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(in)::m2,step_frac,relax

  integer::fft_Nx,fft_Ny,fft_Nz,igrid,igrid_amr,ind,iskip,icell,info
  integer::Kx,Ky,Kz,ix,iy,iz,i,j,k,ngrid_loc,nx_loc
  integer(kind=8)::ntot
  integer(kind=8),save::plan_f=0_8,plan_b=0_8
  integer,save::saved_Nx=0,saved_Ny=0,saved_Nz=0
  real(dp),save::saved_dx=-1d0
  real(dp)::dx_loc,scale,kd2,denom,twopi,du,uold
  real(dp),allocatable,save::b3(:,:,:)
  real(dp),allocatable,save::eigx(:),eigy(:),eigz(:)
  complex(kind=8),allocatable,save::bk(:,:,:)
  integer,parameter::FFTW_EST=64

#ifdef USE_FFTW
  ! Reuse the production Poisson FFTW-MPI+OMP engine on multi-rank jobs.
  ! This shares its in-place transposed plans, sparse slab communication and
  ! validated 1024^3 path.  Keep the single-rank local implementation below
  ! only as a compact reference path.
  if(ncpu > 1) then
     call fftw_scalar_solve_uniform(ilevel,m2,step_frac,relax)
     return
  end if
#endif

  fft_Nx = nx*2**ilevel; fft_Ny = ny*2**ilevel; fft_Nz = nz*2**ilevel
  ntot = int(fft_Nx,8)*int(fft_Ny,8)*int(fft_Nz,8)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=0.5d0**ilevel*scale
  twopi=2d0*acos(-1d0)

  ! Cache FFT work arrays, plans, and discrete one-dimensional
  ! eigenvalues.  The previous path allocated about 3*N^3 reals,
  ! recreated both FFTW plans, and recomputed millions of cosines on
  ! every nonlinear correction.
  if(fft_Nx/=saved_Nx .or. fft_Ny/=saved_Ny .or. fft_Nz/=saved_Nz &
       & .or. dx_loc/=saved_dx) then
     if(plan_f/=0_8) call dfftw_destroy_plan(plan_f)
     if(plan_b/=0_8) call dfftw_destroy_plan(plan_b)
     plan_f=0_8
     plan_b=0_8
     if(allocated(b3)) deallocate(b3,bk,eigx,eigy,eigz)
     allocate(b3(fft_Nx,fft_Ny,fft_Nz))
     allocate(bk(fft_Nx/2+1,fft_Ny,fft_Nz))
     allocate(eigx(fft_Nx/2+1),eigy(fft_Ny),eigz(fft_Nz))
     do i=1,fft_Nx/2+1
        eigx(i)=(2d0-2d0*cos(twopi*dble(i-1)/dble(fft_Nx)))/dx_loc**2
     end do
     do j=1,fft_Ny
        eigy(j)=(2d0-2d0*cos(twopi*dble(j-1)/dble(fft_Ny)))/dx_loc**2
     end do
     do k=1,fft_Nz
        eigz(k)=(2d0-2d0*cos(twopi*dble(k-1)/dble(fft_Nz)))/dx_loc**2
     end do
     call dfftw_plan_dft_r2c_3d(plan_f,fft_Nx,fft_Ny,fft_Nz,b3,bk,FFTW_EST)
     call dfftw_plan_dft_c2r_3d(plan_b,fft_Nx,fft_Ny,fft_Nz,bk,b3,FFTW_EST)
     saved_Nx=fft_Nx
     saved_Ny=fft_Ny
     saved_Nz=fft_Nz
     saved_dx=dx_loc
  end if

  ! Gather the staged rhs in place.  Every active cell has a unique
  ! owner, so zero-fill plus SUM is equivalent to the old bl -> b3
  ! out-of-place allreduce while eliminating one full-grid work array.
  b3=0d0
  ngrid_loc=active(ilevel)%ngrid
  do igrid=1,ngrid_loc
     igrid_amr=active(ilevel)%igrid(igrid)
     Kx=nint(xg(igrid_amr,1)*dble(fft_Nx))
     Ky=nint(xg(igrid_amr,2)*dble(fft_Ny))
     Kz=nint(xg(igrid_amr,3)*dble(fft_Nz))
     do ind=1,twotondim
        ix=modulo(Kx-1+mod(ind-1,2),  fft_Nx)
        iy=modulo(Ky-1+mod((ind-1)/2,2), fft_Ny)
        iz=modulo(Kz-1+(ind-1)/4,     fft_Nz)
        iskip=ncoarse+(ind-1)*ngridmax
        icell=iskip+igrid_amr
        b3(ix+1,iy+1,iz+1)=scalar_gr_old(icell)
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(MPI_IN_PLACE,b3,int(ntot),MPI_DOUBLE_PRECISION,MPI_SUM, &
       & MPI_COMM_WORLD,info)
#endif

  ! Spectral solve with the discrete 7-point eigenvalues
  call dfftw_execute_dft_r2c(plan_f,b3,bk)
!$omp parallel do private(i,j,k,kd2,denom) collapse(2)
  do k=1,fft_Nz
     do j=1,fft_Ny
        do i=1,fft_Nx/2+1
           kd2=eigx(i)+eigy(j)+eigz(k)
           denom=-(kd2+m2)
           if(abs(denom) > 1d-30) then
              bk(i,j,k)=bk(i,j,k)/denom/dble(ntot)
           else
              bk(i,j,k)=(0d0,0d0)   ! k=0 with m2=0: gauge mode
           end if
        end do
     end do
  end do
  call dfftw_execute_dft_c2r(plan_b,bk,b3)

  ! Add the correction to the local cells
  do igrid=1,ngrid_loc
     igrid_amr=active(ilevel)%igrid(igrid)
     Kx=nint(xg(igrid_amr,1)*dble(fft_Nx))
     Ky=nint(xg(igrid_amr,2)*dble(fft_Ny))
     Kz=nint(xg(igrid_amr,3)*dble(fft_Nz))
     do ind=1,twotondim
        ix=modulo(Kx-1+mod(ind-1,2),  fft_Nx)
        iy=modulo(Ky-1+mod((ind-1)/2,2), fft_Ny)
        iz=modulo(Kz-1+(ind-1)/4,     fft_Nz)
        iskip=ncoarse+(ind-1)*ngridmax
        icell=iskip+igrid_amr
        uold=scalar_gr(icell)
        du=relax*b3(ix+1,iy+1,iz+1)
        if(step_frac > 0d0 .and. abs(uold) > 0d0) then
           du=max(-step_frac*abs(uold),min(step_frac*abs(uold),du))
        end if
        scalar_gr(icell)=uold+du
     end do
  end do

end subroutine level_fft_helmholtz

#ifdef USE_FFTW
!=========================================================
! level_fft_helmholtz_mpi: distributed FFTW-MPI solve of
! (lap-m2) du = scalar_gr_old on a fully refined periodic
! level.  RAMSES cells are exchanged with FFTW x-slabs as
! (slab index,value) pairs and returned in the same packed
! order after the inverse transform.
!=========================================================
subroutine level_fft_helmholtz_mpi(ilevel,m2,step_frac,relax)
  use amr_commons
  use poisson_commons
  use iso_c_binding
  use omp_lib
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  include 'fftw3-mpi.f03'

  integer,intent(in)::ilevel
  real(dp),intent(in)::m2,step_frac,relax

  logical,save::initialized=.false.
  integer,save::saved_Nx=0,saved_Ny=0,saved_Nz=0
  type(C_PTR),save::plan_f=C_NULL_PTR,plan_b=C_NULL_PTR,p_data=C_NULL_PTR
  real(C_DOUBLE),pointer,save::rdata(:)=>null()
  complex(C_DOUBLE_COMPLEX),pointer,save::cdata(:)=>null()
  integer(C_INTPTR_T),save::local_nx=0,nx_start=0,alloc_local=0
  real(dp),allocatable,save::eigx(:),eigy(:),eigz(:)

  integer::fft_Nx,fft_Ny,fft_Nz,fftw_block,slab_real_size
  integer::igrid,igrid_amr,ngrid_loc,ind,iskip,icell
  integer::Kx,Ky,Kz,ix,iy,iz,x_local,dest,idx_3d
  integer::irank,info,n_send,n_recv
  integer::i,j,k
  integer(kind=8)::ntot,ncomplex,idx_c
  integer,allocatable::sendcounts(:),recvcounts(:),sdispls(:),rdispls(:)
  integer,allocatable::sendpos(:)
  real(dp),allocatable::sendbuf(:),recvbuf(:)
  real(dp)::dx_loc,scale,kd2,denom,twopi,du,uold
  integer::nx_loc
  integer(C_INT)::thread_ok

  fft_Nx=nx*2**ilevel
  fft_Ny=ny*2**ilevel
  fft_Nz=nz*2**ilevel
  ntot=int(fft_Nx,8)*int(fft_Ny,8)*int(fft_Nz,8)
  fftw_block=(fft_Nx+ncpu-1)/ncpu
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=0.5d0**ilevel*scale
  twopi=2d0*acos(-1d0)

  if(.not.initialized) then
     thread_ok=fftw_init_threads()
     call fftw_plan_with_nthreads(int(omp_get_max_threads(),C_INT))
     call fftw_mpi_init()
     initialized=.true.
  end if

  if(fft_Nx/=saved_Nx .or. fft_Ny/=saved_Ny .or. fft_Nz/=saved_Nz) then
     if(c_associated(plan_f)) call fftw_destroy_plan(plan_f)
     if(c_associated(plan_b)) call fftw_destroy_plan(plan_b)
     if(c_associated(p_data)) call fftw_free(p_data)
     plan_f=C_NULL_PTR
     plan_b=C_NULL_PTR
     p_data=C_NULL_PTR
     nullify(rdata,cdata)
     if(allocated(eigx)) deallocate(eigx,eigy,eigz)

     alloc_local=fftw_mpi_local_size_3d( &
          int(fft_Nx,C_INTPTR_T),int(fft_Ny,C_INTPTR_T), &
          int(fft_Nz/2+1,C_INTPTR_T),MPI_COMM_WORLD,local_nx,nx_start)
     p_data=fftw_alloc_complex(alloc_local)
     slab_real_size=2*int(alloc_local)
     call c_f_pointer(p_data,rdata,[slab_real_size])
     call c_f_pointer(p_data,cdata,[int(alloc_local)])

     plan_f=fftw_mpi_plan_dft_r2c_3d( &
          int(fft_Nx,C_INTPTR_T),int(fft_Ny,C_INTPTR_T), &
          int(fft_Nz,C_INTPTR_T),rdata,cdata,MPI_COMM_WORLD,FFTW_ESTIMATE)
     plan_b=fftw_mpi_plan_dft_c2r_3d( &
          int(fft_Nx,C_INTPTR_T),int(fft_Ny,C_INTPTR_T), &
          int(fft_Nz,C_INTPTR_T),cdata,rdata,MPI_COMM_WORLD,FFTW_ESTIMATE)

     allocate(eigx(0:int(local_nx)-1),eigy(0:fft_Ny-1),eigz(0:fft_Nz/2))
     do i=0,int(local_nx)-1
        eigx(i)=(2d0-2d0*cos(twopi*dble(int(nx_start)+i)/dble(fft_Nx))) &
             & /dx_loc**2
     end do
     do j=0,fft_Ny-1
        eigy(j)=(2d0-2d0*cos(twopi*dble(j)/dble(fft_Ny)))/dx_loc**2
     end do
     do k=0,fft_Nz/2
        eigz(k)=(2d0-2d0*cos(twopi*dble(k)/dble(fft_Nz)))/dx_loc**2
     end do
     saved_Nx=fft_Nx
     saved_Ny=fft_Ny
     saved_Nz=fft_Nz
     if(myid==1) write(*,'(A,I5,A,I5,A,I5,A)') &
          ' scalar FFTW-MPI plan: ',fft_Nx,'x',fft_Ny,'x',fft_Nz,' grid'
  end if

  allocate(sendcounts(0:ncpu-1),recvcounts(0:ncpu-1))
  allocate(sdispls(0:ncpu-1),rdispls(0:ncpu-1),sendpos(0:ncpu-1))
  sendcounts=0
  ngrid_loc=active(ilevel)%ngrid
  do igrid=1,ngrid_loc
     igrid_amr=active(ilevel)%igrid(igrid)
     Kx=nint(xg(igrid_amr,1)*dble(fft_Nx))
     do ind=1,twotondim
        ix=modulo(Kx-1+mod(ind-1,2),fft_Nx)
        dest=min(ix/fftw_block,ncpu-1)
        sendcounts(dest)=sendcounts(dest)+1
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLTOALL(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER, &
       & MPI_COMM_WORLD,info)
#else
  recvcounts=sendcounts
#endif
  sdispls(0)=0
  rdispls(0)=0
  do irank=1,ncpu-1
     sdispls(irank)=sdispls(irank-1)+sendcounts(irank-1)
     rdispls(irank)=rdispls(irank-1)+recvcounts(irank-1)
  end do
  n_send=sdispls(ncpu-1)+sendcounts(ncpu-1)
  n_recv=rdispls(ncpu-1)+recvcounts(ncpu-1)
  allocate(sendbuf(0:max(1,2*n_send)-1),recvbuf(0:max(1,2*n_recv)-1))
  sendpos=sdispls

  do igrid=1,ngrid_loc
     igrid_amr=active(ilevel)%igrid(igrid)
     Kx=nint(xg(igrid_amr,1)*dble(fft_Nx))
     Ky=nint(xg(igrid_amr,2)*dble(fft_Ny))
     Kz=nint(xg(igrid_amr,3)*dble(fft_Nz))
     do ind=1,twotondim
        ix=modulo(Kx-1+mod(ind-1,2),fft_Nx)
        iy=modulo(Ky-1+mod((ind-1)/2,2),fft_Ny)
        iz=modulo(Kz-1+(ind-1)/4,fft_Nz)
        dest=min(ix/fftw_block,ncpu-1)
        x_local=ix-dest*fftw_block
        idx_3d=x_local*fft_Ny*2*(fft_Nz/2+1) &
             & +iy*2*(fft_Nz/2+1)+iz
        iskip=ncoarse+(ind-1)*ngridmax
        icell=iskip+igrid_amr
        sendbuf(2*sendpos(dest))=dble(idx_3d)
        sendbuf(2*sendpos(dest)+1)=scalar_gr_old(icell)
        sendpos(dest)=sendpos(dest)+1
     end do
  end do

  sendcounts=2*sendcounts
  recvcounts=2*recvcounts
  sdispls=2*sdispls
  rdispls=2*rdispls
#ifndef WITHOUTMPI
  call MPI_ALLTOALLV(sendbuf,sendcounts,sdispls,MPI_DOUBLE_PRECISION, &
       & recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
#else
  recvbuf=sendbuf
#endif

  rdata=0d0
  do i=0,n_recv-1
     idx_3d=nint(recvbuf(2*i))
     rdata(idx_3d+1)=rdata(idx_3d+1)+recvbuf(2*i+1)
  end do
  call fftw_mpi_execute_dft_r2c(plan_f,rdata,cdata)

  ncomplex=int(local_nx,8)*int(fft_Ny,8)*int(fft_Nz/2+1,8)
!$omp parallel do private(idx_c,i,j,k,kd2,denom) schedule(static)
  do idx_c=0,ncomplex-1
     i=int(idx_c/(int(fft_Ny,8)*int(fft_Nz/2+1,8)))
     j=int(mod(idx_c/int(fft_Nz/2+1,8),int(fft_Ny,8)))
     k=int(mod(idx_c,int(fft_Nz/2+1,8)))
     kd2=eigx(i)+eigy(j)+eigz(k)
     denom=-(kd2+m2)
     if(abs(denom)>1d-30) then
        cdata(idx_c+1)=cdata(idx_c+1)/denom/dble(ntot)
     else
        cdata(idx_c+1)=(0d0,0d0)
     end if
  end do
  call fftw_mpi_execute_dft_c2r(plan_b,cdata,rdata)

  do i=0,n_recv-1
     idx_3d=nint(recvbuf(2*i))
     recvbuf(2*i+1)=rdata(idx_3d+1)
  end do
#ifndef WITHOUTMPI
  call MPI_ALLTOALLV(recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION, &
       & sendbuf,sendcounts,sdispls,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
#else
  sendbuf=recvbuf
#endif

  sendpos=sdispls/2
  do igrid=1,ngrid_loc
     igrid_amr=active(ilevel)%igrid(igrid)
     Kx=nint(xg(igrid_amr,1)*dble(fft_Nx))
     do ind=1,twotondim
        ix=modulo(Kx-1+mod(ind-1,2),fft_Nx)
        dest=min(ix/fftw_block,ncpu-1)
        iskip=ncoarse+(ind-1)*ngridmax
        icell=iskip+igrid_amr
        uold=scalar_gr(icell)
        du=relax*sendbuf(2*sendpos(dest)+1)
        if(step_frac>0d0 .and. abs(uold)>0d0) then
           du=max(-step_frac*abs(uold),min(step_frac*abs(uold),du))
        end if
        scalar_gr(icell)=uold+du
        sendpos(dest)=sendpos(dest)+1
     end do
  end do

  deallocate(sendcounts,recvcounts,sdispls,rdispls,sendpos,sendbuf,recvbuf)
end subroutine level_fft_helmholtz_mpi
#endif

!=========================================================
! fR_build_fft_rhs: stage b = -(lap u - S(u)) in scalar_gr_old
! and return the level-mean Newton mass term m2bar
!=========================================================
subroutine fR_build_fft_rhs(ilevel, R_bar, fR_bar, m2bar)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(in)::R_bar, fR_bar
  real(dp),intent(out)::m2bar

  integer::igrid,ngrid,ncache,i,ind,iskip,idim,info
  integer::ig_left,ig_right,ih_left,ih_right
  real(dp)::dx,scale,dx_loc,dx2_inv,nx_frac
  integer::nx_loc,np1
  real(dp)::u_c,u_abs,u_nb_l,u_nb_r,lapl,source,R_of_u,dR_du
  real(dp)::a2_over_3,rho_coeff,boxratio_sq,R_bar0
  real(dp)::fR0_abs,small_fR,inv_np1
  real(dp),dimension(2)::acc_loc,acc_glob
  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2_inv=1d0/dx_loc**2

  np1 = fR_n + 1
  fR0_abs=abs(fR0)
  small_fR=1d-30*fR0_abs
  inv_np1=1d0/dble(np1)
  boxratio_sq = (boxlen_ini / 2997.92458d0)**2
  a2_over_3 = aexp**2 * boxratio_sq / 3d0
  R_bar0 = 3d0 * (omega_m + 4d0 * omega_l)
  rho_coeff = omega_m * boxratio_sq / aexp

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  acc_loc=0d0

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w,ig_left,ig_right,ih_left,ih_right, &
!$omp&  u_c,u_abs,u_nb_l,u_nb_r,lapl,source,R_of_u,dR_du) &
!$omp& reduction(+:acc_loc) schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=vain_face_grid(igrid+i-1,2*idim-1)
           igridn_w(i,2*idim  )=vain_face_grid(igrid+i-1,2*idim  )
        end do
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do
        do i=1,ngrid
           u_c = scalar_gr(ind_cell_w(i))
           u_abs = abs(u_c)
           lapl = 0d0
           do idim=1,ndim
              ig_left =ggg(idim,1,ind); ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax
              u_nb_l = u_c; u_nb_r = u_c
              if(igridn_w(i,ig_left ) > 0) u_nb_l = scalar_gr(igridn_w(i,ig_left )+ih_left)
              if(igridn_w(i,ig_right) > 0) u_nb_r = scalar_gr(igridn_w(i,ig_right)+ih_right)
              lapl = lapl + (u_nb_l + u_nb_r - 2d0*u_c) * dx2_inv
           end do
           if(u_abs > small_fR) then
              select case(np1)
              case(1)
                 R_of_u = R_bar0 * fR0_abs/u_abs
              case(2)
                 R_of_u = R_bar0 * sqrt(fR0_abs/u_abs)
              case default
                 R_of_u = R_bar0 * (fR0_abs/u_abs)**inv_np1
              end select
              dR_du = -R_of_u / (dble(np1) * u_c)
           else
              R_of_u = R_bar
              dR_du = 0d0
           end if
           source = a2_over_3*(R_of_u - R_bar) &
                & - rho_coeff*(rho(ind_cell_w(i)) - rho_tot)
           scalar_gr_old(ind_cell_w(i))=-(lapl-source)
           acc_loc(1) = acc_loc(1) + a2_over_3*dR_du
           acc_loc(2) = acc_loc(2) + 1d0
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(acc_loc,acc_glob,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
       & MPI_COMM_WORLD,info)
#else
  acc_glob=acc_loc
#endif
  m2bar = acc_glob(1)/max(acc_glob(2),1d0)
  m2bar = max(m2bar, 0d0)

end subroutine fR_build_fft_rhs

!=========================================================
! sb_build_fft_rhs: symmetron/dilaton Newton rhs (shared);
! b = -(lap u - mass*u - cubic*u^3) in scalar_gr_old,
! m2bar = level-mean of (mass + 3*cubic*u^2), clipped >= 0
!=========================================================
subroutine sb_build_fft_rhs(ilevel, assb_in, L_in, m2bar)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(in)::assb_in, L_in
  real(dp),intent(out)::m2bar

  integer::igrid,ngrid,ncache,i,ind,iskip,idim,info
  integer::ig_left,ig_right,ih_left,ih_right
  real(dp)::dx,scale,dx_loc,dx2_inv
  integer::nx_loc
  real(dp)::u_c,u_nb_l,u_nb_r,lapl,mass_term,cubic_coeff,a2_over_2L2
  real(dp),dimension(2)::acc_loc,acc_glob
  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2_inv=1d0/dx_loc**2
  a2_over_2L2 = aexp**2 / (2d0*(L_in/boxlen_ini)**2)
  cubic_coeff = a2_over_2L2

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  acc_loc=0d0

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w,ig_left,ig_right,ih_left,ih_right, &
!$omp&  u_c,u_nb_l,u_nb_r,lapl,mass_term) &
!$omp& reduction(+:acc_loc) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=vain_face_grid(igrid+i-1,2*idim-1)
           igridn_w(i,2*idim  )=vain_face_grid(igrid+i-1,2*idim  )
        end do
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do
        do i=1,ngrid
           u_c = scalar_gr(ind_cell_w(i))
           lapl = 0d0
           do idim=1,ndim
              ig_left =ggg(idim,1,ind); ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax
              u_nb_l = u_c; u_nb_r = u_c
              if(igridn_w(i,ig_left ) > 0) u_nb_l = scalar_gr(igridn_w(i,ig_left )+ih_left)
              if(igridn_w(i,ig_right) > 0) u_nb_r = scalar_gr(igridn_w(i,ig_right)+ih_right)
              lapl = lapl + (u_nb_l + u_nb_r - 2d0*u_c) * dx2_inv
           end do
           mass_term = a2_over_2L2 * (rho(ind_cell_w(i))*(assb_in/aexp)**3 - 1d0)
           scalar_gr_old(ind_cell_w(i)) = &
                & -(lapl - mass_term*u_c - cubic_coeff*u_c**3)
           acc_loc(1) = acc_loc(1) + mass_term + 3d0*cubic_coeff*u_c**2
           acc_loc(2) = acc_loc(2) + 1d0
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(acc_loc,acc_glob,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
       & MPI_COMM_WORLD,info)
#else
  acc_glob=acc_loc
#endif
  m2bar = max(acc_glob(1)/max(acc_glob(2),1d0), 0d0)

end subroutine sb_build_fft_rhs

!=========================================================
! dil_build_fft_rhs: Brax+12 dilaton Newton rhs;
! b = -(lap chi - S(chi)) in scalar_gr_old,
! m2bar = level-mean of dS/dchi (>0)
!=========================================================
subroutine dil_build_fft_rhs(ilevel, A2_d, s_d, chibar_d, m2bar)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(in)::A2_d,s_d,chibar_d
  real(dp),intent(out)::m2bar

  integer::igrid,ngrid,ncache,i,ind,iskip,idim,info
  integer::ig_left,ig_right,ih_left,ih_right
  real(dp)::dx,scale,dx_loc,dx2_inv
  integer::nx_loc
  real(dp)::u_c,u_nb_l,u_nb_r,lapl
  real(dp)::boxratio_sq,cA,cV,pexp,wfac,vphi,dvphi,vbar,source
  real(dp),dimension(2)::acc_loc,acc_glob
  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2_inv=1d0/dx_loc**2

  boxratio_sq=(boxlen_ini/2997.92458d0)**2
  cA = 3d0*omega_m*A2_d*boxratio_sq/aexp
  cV = aexp**2*boxratio_sq
  pexp = 1d0 - 3d0/s_d
  vbar = -3d0*omega_m*beta_dilaton*(A2_d*chibar_d/beta_dilaton)**pexp

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  acc_loc=0d0

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w,ig_left,ig_right,ih_left,ih_right, &
!$omp&  u_c,u_nb_l,u_nb_r,lapl,wfac,vphi,dvphi,source) &
!$omp& reduction(+:acc_loc) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do
        do i=1,ngrid
           u_c = scalar_gr(ind_cell_w(i))
           lapl = 0d0
           do idim=1,ndim
              ig_left =ggg(idim,1,ind); ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax
              u_nb_l = u_c; u_nb_r = u_c
              if(igridn_w(i,ig_left ) > 0) u_nb_l = scalar_gr(igridn_w(i,ig_left )+ih_left)
              if(igridn_w(i,ig_right) > 0) u_nb_r = scalar_gr(igridn_w(i,ig_right)+ih_right)
              lapl = lapl + (u_nb_l + u_nb_r - 2d0*u_c) * dx2_inv
           end do
           wfac = A2_d*max(u_c,1d-30)/beta_dilaton
           vphi = -3d0*omega_m*beta_dilaton*wfac**pexp
           dvphi = -3d0*omega_m*A2_d*pexp*wfac**(pexp-1d0)
           source = cA*(rho(ind_cell_w(i))*u_c - chibar_d) + cV*(vphi - vbar)
           scalar_gr_old(ind_cell_w(i)) = -(lapl - source)
           acc_loc(1) = acc_loc(1) + cA*rho(ind_cell_w(i)) + cV*dvphi
           acc_loc(2) = acc_loc(2) + 1d0
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(acc_loc,acc_glob,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
       & MPI_COMM_WORLD,info)
#else
  acc_glob=acc_loc
#endif
  m2bar = max(acc_glob(1)/max(acc_glob(2),1d0), 0d0)

end subroutine dil_build_fft_rhs

!=========================================================
! vain_prepare_uniform_cache: build the face/edge grid topology once
! for a uniform spectral scalar solve.  The active-grid order and domain
! decomposition cannot change inside one nDGP/Galileon nonlinear loop.
! Rebuilding once at the next coarse step is intentionally conservative:
! it remains correct across a RAMSES load balance without signatures or
! invalidation hooks.
!=========================================================
subroutine vain_prepare_uniform_cache(ilevel)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  integer::ia,igrid,idim,sx,sy,sz,ixside,iyside,izside,iedge,igrid_diag
  integer::ncache,nalloc

  ncache=active(ilevel)%ngrid
  nalloc=max(ncache,1)
  vain_cache_ready=.false.

  if(allocated(vain_face_grid)) then
     if(size(vain_face_grid,1)/=nalloc) then
        deallocate(vain_face_grid,vain_xy_grid,vain_xz_grid,vain_yz_grid)
     end if
  end if
  if(.not.allocated(vain_face_grid)) then
     allocate(vain_face_grid(nalloc,6))
     allocate(vain_xy_grid(nalloc,4),vain_xz_grid(nalloc,4),vain_yz_grid(nalloc,4))
  end if

!$omp parallel do private(ia,igrid,idim,sx,sy,sz,ixside,iyside,izside, &
!$omp& iedge,igrid_diag) schedule(static)
  do ia=1,ncache
     igrid=active(ilevel)%igrid(ia)
     do idim=1,ndim
        vain_face_grid(ia,2*idim-1)=morton_nbor_grid(igrid,ilevel,2*idim-1)
        vain_face_grid(ia,2*idim  )=morton_nbor_grid(igrid,ilevel,2*idim  )
     end do
     do sx=-1,1,2
        ixside=merge(1,2,sx<0)
        do sy=-1,1,2
           iyside=merge(3,4,sy<0)
           iedge=((sx+1)/2)*2+(sy+1)/2+1
           igrid_diag=vain_face_grid(ia,ixside)
           if(igrid_diag>0) igrid_diag=morton_nbor_grid(igrid_diag,ilevel,iyside)
           vain_xy_grid(ia,iedge)=igrid_diag
        end do
        do sz=-1,1,2
           izside=merge(5,6,sz<0)
           iedge=((sx+1)/2)*2+(sz+1)/2+1
           igrid_diag=vain_face_grid(ia,ixside)
           if(igrid_diag>0) igrid_diag=morton_nbor_grid(igrid_diag,ilevel,izside)
           vain_xz_grid(ia,iedge)=igrid_diag
        end do
     end do
     do sy=-1,1,2
        iyside=merge(3,4,sy<0)
        do sz=-1,1,2
           izside=merge(5,6,sz<0)
           iedge=((sy+1)/2)*2+(sz+1)/2+1
           igrid_diag=vain_face_grid(ia,iyside)
           if(igrid_diag>0) igrid_diag=morton_nbor_grid(igrid_diag,ilevel,izside)
           vain_yz_grid(ia,iedge)=igrid_diag
        end do
     end do
  end do
!$omp end parallel do

  vain_cache_level=ilevel
  vain_cache_ngrid=ncache
  vain_cache_ready=.true.
end subroutine vain_prepare_uniform_cache

!=========================================================
! vain_build_fft_rhs: nDGP/galileon operator splitting;
! split H_ij into trace and traceless pieces,
!   H_ij H_ij = Hbar_ij Hbar_ij + A^2/3,
! and solve the local quadratic
!   (2 coeff/3)*A^2 + A = S + coeff*Hbar_ij Hbar_ij
! for A = lap(phi) (branch A->S as coeff->0; clamped at the
! extremum when the discriminant is negative) and stage
! b = A - lap(phi) in scalar_gr_old (then lap(dphi) = b).
!=========================================================
subroutine vain_build_fft_rhs(ilevel, srcfac, coeff, centered_mixed, rhs_rel)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(in)::srcfac, coeff
  logical,intent(in)::centered_mixed
  real(dp),intent(out)::rhs_rel

  integer::igrid,ngrid,ncache,i,ind,iskip,idim,icell,info
  integer::ig_left,ig_right
  integer::d,bx,by,bz,tx,ty,tz,sx,sy,sz,ixside,iyside,izside
  integer::iedge,igrid_diag,ind_diag
  real(dp)::dx,scale,dx_loc,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,s_src,t_ij,tbar_ij,qcoeff,disc,a_tgt
  real(dp)::rhs_sum_local,rhs_sum_global,ncell_local,ncell_global,rhs_mean
  real(dp)::rhs_max_local,rhs_max_global,src_max_local,src_max_global
  real(dp)::phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp
  real(dp)::phi_xx,phi_yy,phi_zz,mix_xy2,mix_xz2,mix_yz2
  real(dp)::phi_pp,phi_pm,phi_mp,phi_mm
  real(dp)::dpp,dpm,dmp,dmm
  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w
  integer,dimension(1:nvector,1:4)::igridxy_w,igridxz_w,igridyz_w
  integer,dimension(1:nvector,1:12)::ind_diag_w
  integer,dimension(12),parameter::odx=(/ 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0/)
  integer,dimension(12),parameter::ody=(/ 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1/)
  integer,dimension(12),parameter::odz=(/ 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)

  ncache=active(ilevel)%ngrid
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2_inv=1d0/dx_loc**2

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  if(.not.vain_cache_ready .or. vain_cache_level/=ilevel .or. &
       & vain_cache_ngrid/=ncache) call vain_prepare_uniform_cache(ilevel)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w,igridxy_w,igridxz_w,igridyz_w, &
!$omp&  ind_diag_w,ig_left,ig_right, &
!$omp&  d,bx,by,bz,tx,ty,tz,sx,sy,sz,ixside,iyside,izside, &
!$omp&  iedge,igrid_diag,ind_diag, &
!$omp&  u_c,lapl,s_src,t_ij,tbar_ij,qcoeff,disc,a_tgt, &
!$omp&  phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp,phi_xx,phi_yy,phi_zz, &
!$omp&  mix_xy2,mix_xz2,mix_yz2,phi_pp,phi_pm,phi_mp,phi_mm, &
!$omp&  dpp,dpm,dmp,dmm) &
!$omp& schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=vain_face_grid(igrid+i-1,2*idim-1)
           igridn_w(i,2*idim  )=vain_face_grid(igrid+i-1,2*idim  )
        end do
     end do
     do i=1,ngrid
        igridxy_w(i,1:4)=vain_xy_grid(igrid+i-1,1:4)
        igridxz_w(i,1:4)=vain_xz_grid(igrid+i-1,1:4)
        igridyz_w(i,1:4)=vain_yz_grid(igrid+i-1,1:4)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do
        bx=mod(ind-1,2)
        by=mod((ind-1)/2,2)
        bz=(ind-1)/4
        do d=1,12
           tx=bx+odx(d); ty=by+ody(d); tz=bz+odz(d)
           sx=0; sy=0; sz=0
           if(tx<0) sx=-1
           if(tx>1) sx= 1
           if(ty<0) sy=-1
           if(ty>1) sy= 1
           if(tz<0) sz=-1
           if(tz>1) sz= 1
           tx=modulo(tx,2); ty=modulo(ty,2); tz=modulo(tz,2)
           ind_diag=1+tx+2*ty+4*tz
           do i=1,ngrid
              if(sx/=0 .and. sy/=0) then
                 iedge=((sx+1)/2)*2+(sy+1)/2+1
                 igrid_diag=igridxy_w(i,iedge)
              else if(sx/=0 .and. sz/=0) then
                 iedge=((sx+1)/2)*2+(sz+1)/2+1
                 igrid_diag=igridxz_w(i,iedge)
              else if(sy/=0 .and. sz/=0) then
                 iedge=((sy+1)/2)*2+(sz+1)/2+1
                 igrid_diag=igridyz_w(i,iedge)
              else if(sx/=0) then
                 ixside=merge(1,2,sx<0)
                 igrid_diag=igridn_w(i,ixside)
              else if(sy/=0) then
                 iyside=merge(3,4,sy<0)
                 igrid_diag=igridn_w(i,iyside)
              else if(sz/=0) then
                 izside=merge(5,6,sz<0)
                 igrid_diag=igridn_w(i,izside)
              else
                 igrid_diag=ind_grid_w(i)
              end if
              if(igrid_diag>0) then
                 ind_diag_w(i,d)=ncoarse+(ind_diag-1)*ngridmax+igrid_diag
              else
                 ind_diag_w(i,d)=0
              end if
           end do
        end do
        do i=1,ngrid
           u_c = scalar_gr(ind_cell_w(i))

           ig_left =ggg(1,1,ind); ig_right=ggg(1,2,ind)
           phi_xm = u_c; phi_xp = u_c
           if(igridn_w(i,ig_left ) > 0) phi_xm = scalar_gr(igridn_w(i,ig_left) +ncoarse+(hhh(1,1,ind)-1)*ngridmax)
           if(igridn_w(i,ig_right) > 0) phi_xp = scalar_gr(igridn_w(i,ig_right)+ncoarse+(hhh(1,2,ind)-1)*ngridmax)
           ig_left =ggg(2,1,ind); ig_right=ggg(2,2,ind)
           phi_ym = u_c; phi_yp = u_c
           if(igridn_w(i,ig_left ) > 0) phi_ym = scalar_gr(igridn_w(i,ig_left) +ncoarse+(hhh(2,1,ind)-1)*ngridmax)
           if(igridn_w(i,ig_right) > 0) phi_yp = scalar_gr(igridn_w(i,ig_right)+ncoarse+(hhh(2,2,ind)-1)*ngridmax)
           ig_left =ggg(3,1,ind); ig_right=ggg(3,2,ind)
           phi_zm = u_c; phi_zp = u_c
           if(igridn_w(i,ig_left ) > 0) phi_zm = scalar_gr(igridn_w(i,ig_left) +ncoarse+(hhh(3,1,ind)-1)*ngridmax)
           if(igridn_w(i,ig_right) > 0) phi_zp = scalar_gr(igridn_w(i,ig_right)+ncoarse+(hhh(3,2,ind)-1)*ngridmax)

           lapl = (phi_xp+phi_xm+phi_yp+phi_ym+phi_zp+phi_zm-6d0*u_c)*dx2_inv
           phi_xx = (phi_xp+phi_xm-2d0*u_c)*dx2_inv
           phi_yy = (phi_yp+phi_ym-2d0*u_c)*dx2_inv
           phi_zz = (phi_zp+phi_zm-2d0*u_c)*dx2_inv
           if(ind_diag_w(i,1)>0) then
              phi_pp=scalar_gr(ind_diag_w(i,1))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 1, 0,u_c,phi_pp)
           end if
           if(ind_diag_w(i,2)>0) then
              phi_pm=scalar_gr(ind_diag_w(i,2))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1,-1, 0,u_c,phi_pm)
           end if
           if(ind_diag_w(i,3)>0) then
              phi_mp=scalar_gr(ind_diag_w(i,3))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 1, 0,u_c,phi_mp)
           end if
           if(ind_diag_w(i,4)>0) then
              phi_mm=scalar_gr(ind_diag_w(i,4))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1,-1, 0,u_c,phi_mm)
           end if
           if(centered_mixed) then
              mix_xy2=(0.25d0*(phi_pp-phi_pm-phi_mp+phi_mm)*dx2_inv)**2
           else
              dpp=(phi_pp-phi_xp-phi_yp+u_c)*dx2_inv
              dpm=(phi_xp-phi_pm-u_c+phi_ym)*dx2_inv
              dmp=(phi_yp-phi_mp-u_c+phi_xm)*dx2_inv
              dmm=(u_c-phi_ym-phi_xm+phi_mm)*dx2_inv
              mix_xy2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
           end if
           if(ind_diag_w(i,5)>0) then
              phi_pp=scalar_gr(ind_diag_w(i,5))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0, 1,u_c,phi_pp)
           end if
           if(ind_diag_w(i,6)>0) then
              phi_pm=scalar_gr(ind_diag_w(i,6))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0,-1,u_c,phi_pm)
           end if
           if(ind_diag_w(i,7)>0) then
              phi_mp=scalar_gr(ind_diag_w(i,7))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0, 1,u_c,phi_mp)
           end if
           if(ind_diag_w(i,8)>0) then
              phi_mm=scalar_gr(ind_diag_w(i,8))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0,-1,u_c,phi_mm)
           end if
           if(centered_mixed) then
              mix_xz2=(0.25d0*(phi_pp-phi_pm-phi_mp+phi_mm)*dx2_inv)**2
           else
              dpp=(phi_pp-phi_xp-phi_zp+u_c)*dx2_inv
              dpm=(phi_xp-phi_pm-u_c+phi_zm)*dx2_inv
              dmp=(phi_zp-phi_mp-u_c+phi_xm)*dx2_inv
              dmm=(u_c-phi_zm-phi_xm+phi_mm)*dx2_inv
              mix_xz2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
           end if
           if(ind_diag_w(i,9)>0) then
              phi_pp=scalar_gr(ind_diag_w(i,9))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1, 1,u_c,phi_pp)
           end if
           if(ind_diag_w(i,10)>0) then
              phi_pm=scalar_gr(ind_diag_w(i,10))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1,-1,u_c,phi_pm)
           end if
           if(ind_diag_w(i,11)>0) then
              phi_mp=scalar_gr(ind_diag_w(i,11))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1, 1,u_c,phi_mp)
           end if
           if(ind_diag_w(i,12)>0) then
              phi_mm=scalar_gr(ind_diag_w(i,12))
           else
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1,-1,u_c,phi_mm)
           end if
           if(centered_mixed) then
              mix_yz2=(0.25d0*(phi_pp-phi_pm-phi_mp+phi_mm)*dx2_inv)**2
           else
              dpp=(phi_pp-phi_yp-phi_zp+u_c)*dx2_inv
              dpm=(phi_yp-phi_pm-u_c+phi_zm)*dx2_inv
              dmp=(phi_zp-phi_mp-u_c+phi_ym)*dx2_inv
              dmm=(u_c-phi_zm-phi_ym+phi_mm)*dx2_inv
              mix_yz2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
           end if
           t_ij = phi_xx**2 + phi_yy**2 + phi_zz**2 &
                & + 2d0*(mix_xy2+mix_xz2+mix_yz2)
           tbar_ij=max(t_ij-lapl**2/3d0,0d0)

           s_src = srcfac*(rho(ind_cell_w(i)) - rho_tot)
           if(abs(coeff) > 1d-12) then
              qcoeff=2d0*coeff/3d0
              disc = 1d0 + 4d0*qcoeff*(s_src + coeff*tbar_ij)
              if(disc > 0d0) then
                 a_tgt = (-1d0 + sqrt(disc))/(2d0*qcoeff)
              else
                 a_tgt = -1d0/(2d0*qcoeff)
              end if
           else
              a_tgt = s_src + coeff*tbar_ij
           end if
           scalar_gr_old(ind_cell_w(i)) = a_tgt - lapl
        end do
     end do
  end do

  ! A periodic Laplacian has no k=0 mode.  The nonlinear local
  ! quadratic solve does not preserve this compatibility condition
  ! at an intermediate operator-split iterate, so explicitly project
  ! the staged correction onto its zero-mean subspace.  FFTW would
  ! otherwise discard the same mode silently.
  rhs_sum_local=0d0
  ncell_local=0d0
!$omp parallel do private(igrid,ngrid,i,ind,iskip,icell) &
!$omp& reduction(+:rhs_sum_local,ncell_local) schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           rhs_sum_local=rhs_sum_local+scalar_gr_old(icell)
           ncell_local=ncell_local+1d0
        end do
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_sum_local,rhs_sum_global,1,MPI_DOUBLE_PRECISION, &
       & MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ncell_local,ncell_global,1,MPI_DOUBLE_PRECISION, &
       & MPI_SUM,MPI_COMM_WORLD,info)
#else
  rhs_sum_global=rhs_sum_local
  ncell_global=ncell_local
#endif
  rhs_mean=rhs_sum_global/max(ncell_global,1d0)
  rhs_max_local=0d0
  src_max_local=0d0
!$omp parallel do private(igrid,ngrid,i,ind,iskip,icell,s_src) &
!$omp& reduction(max:rhs_max_local,src_max_local) schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell)=scalar_gr_old(icell)-rhs_mean
           rhs_max_local=max(rhs_max_local,abs(scalar_gr_old(icell)))
           s_src=srcfac*(rho(icell)-rho_tot)
           src_max_local=max(src_max_local,abs(s_src))
        end do
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_max_local,rhs_max_global,1,MPI_DOUBLE_PRECISION, &
       & MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(src_max_local,src_max_global,1,MPI_DOUBLE_PRECISION, &
       & MPI_MAX,MPI_COMM_WORLD,info)
#else
  rhs_max_global=rhs_max_local
  src_max_global=src_max_local
#endif
  rhs_rel=rhs_max_global/max(src_max_global,1d-30)

end subroutine vain_build_fft_rhs
#endif

!=========================================================
! fR_background: compute background Ricci scalar and f_R
!=========================================================
subroutine fR_background(aa, R_bar, fR_bar)
  use amr_parameters, only: dp, omega_m, omega_l, fR0, fR_n
  implicit none
  real(dp),intent(in) ::aa
  real(dp),intent(out)::R_bar, fR_bar
  !-------------------------------------------------------
  ! R_bar = 3*H0^2*(Omega_m/a^3 + 4*Omega_Lambda)
  ! In code units where H0=1:
  !   R_bar = 3*(omega_m/a^3 + 4*omega_l)
  ! f_R_bar = fR0 * (R_bar0/R_bar)^(n+1),  R_bar0 = R_bar(a=1)
  ! Standard convention (Hu & Sawicki 2007): fR0 IS the field value
  ! today, f_R_bar(a=1) = fR0 (no extra factor n).
  !-------------------------------------------------------
  real(dp)::R_bar0
  integer::np1

  np1 = fR_n + 1
  R_bar  = 3d0 * (omega_m / aa**3 + 4d0 * omega_l)
  R_bar0 = 3d0 * (omega_m + 4d0 * omega_l)
  fR_bar = fR0 * (R_bar0 / R_bar)**np1

end subroutine fR_background

!=========================================================
! fR_solve_level: top-level f(R) solver for one AMR level
!=========================================================
subroutine fR_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
#ifdef HYDRO_CUDA
  use scalar_cuda_interface, only: SCAL_MODEL_FR
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount
  !-------------------------------------------------------
  ! 1. Initialize scalar_gr to background fR_bar
  ! 2. Newton-GS iteration
  ! 3. Add fifth force F5 = -0.5 * grad(fR)
  !-------------------------------------------------------
  integer::iter,ncache,info
  real(dp)::R_bar,fR_bar
  logical::gscal_ok
  real(dp),dimension(1:12)::sparams
  real(dp)::gs_dx2
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res
  real(dp)::background_ratio
  logical::converged
  real(dp),save::fR_bar_previous(1:MAXLEVEL)=0d0
#ifdef USE_FFTW
  real(dp)::m2bar
  logical,external::level_fft_ok
#endif

  ! NOTE: do NOT return when ncache==0 — the relaxation loop and the
  ! FFT stage contain MPI collectives (ALLREDUCE) that every rank
  ! must enter; ranks without grids on this level contribute nothing.
  ncache=active(ilevel)%ngrid

  ! Get background values
  call fR_background(aexp, R_bar, fR_bar)

  ! Predict the new-time solution by evolving the previous solution with
  ! the homogeneous background.  In Hu-Sawicki gravity fR_bar changes
  ! rapidly at early times; reusing the absolute old field without this
  ! rescaling can leave Newton-GS thousands of sweeps away from the new
  ! solution.  A zero previous value marks the first solve on this rank.
  ! On restart scalar_gr itself is restored, while this level-local
  ! background predictor is rebuilt on the first post-restart solve.
  if(fR_bar_previous(ilevel) /= 0d0) then
     background_ratio=fR_bar/fR_bar_previous(ilevel)
     if(background_ratio > 0d0) &
          & call fR_rescale_warm_start(ilevel,background_ratio)
  end if
  fR_bar_previous(ilevel)=fR_bar

  ! Seed every cell still at exactly 0 with the background value.
  ! Covers first step, legacy restarts without scalar data, and cells
  ! created by refinement after step 0; converged f_R is
  ! strictly negative so 0 uniquely marks uninitialized cells.
  call fR_seed_scalar(ilevel, fR_bar)
  call vain_prepare_uniform_cache(ilevel)

#ifdef USE_FFTW
  ! Spectral Newton on the uniform domain level: converges the
  ! long-wavelength modes that the GS sweeps below cannot reach
  if(level_fft_ok(ilevel)) then
     call make_virtual_fine_dp(scalar_gr(1), ilevel)
     do iter=1,8
        call fR_build_fft_rhs(ilevel, R_bar, fR_bar, m2bar)
        call level_fft_helmholtz(ilevel, m2bar, 0.25d0, 1d0)
        call make_virtual_fine_dp(scalar_gr(1), ilevel)
     end do
  end if
#endif

  ! Newton-GS relaxation (GPU sweeps when gpu_scalar is active)
  gscal_ok=.false.
#ifdef HYDRO_CUDA
  call scalar_gpu_begin(ilevel, .false., gscal_ok)
  if(gscal_ok) then
     gs_dx2=(0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1))**2
     sparams=0d0
     sparams(1)=1d0/gs_dx2
     sparams(2)=aexp**2*(boxlen_ini/2997.92458d0)**2/3d0
     sparams(3)=omega_m*(boxlen_ini/2997.92458d0)**2/aexp
     sparams(4)=R_bar
     sparams(5)=3d0*(omega_m+4d0*omega_l)
     sparams(6)=abs(fR0)
     sparams(7)=1d-30*abs(fR0)
     sparams(8)=1d0/dble(fR_n+1)
     sparams(9)=dble(fR_n+1)
     sparams(10)=rho_tot
  end if
#endif
  converged = .false.
  do iter=1,n_iter_fR

#ifdef HYDRO_CUDA
     if(gscal_ok) then
        call scalar_gpu_sweep_halo(ilevel, SCAL_MODEL_FR, sparams, 0, &
             & res_max_local, src_max_local)
     else
#endif
     call fR_gauss_seidel(ilevel, R_bar, fR_bar, res_max_local, src_max_local)

     ! Exchange boundaries after each sweep
     call make_virtual_fine_dp(scalar_gr(1), ilevel)
#ifdef HYDRO_CUDA
     end if
#endif

     ! Global convergence check
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < fR_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I6,A,ES10.3)') &
             ' f(R) level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do
#ifdef HYDRO_CUDA
  if(gscal_ok) call scalar_gpu_end(ilevel)
#endif

  if(.not. converged) then
     if(myid==1) write(*,'(A,I2,A,I6,A,ES10.3)') &
          & ' WARNING: f(R) level ',ilevel,' NOT converged after ', &
          & n_iter_fR,' iters, res=',rel_res
     if(scalar_solver_strict) call scalar_solver_abort
  end if

  ! Save scalar_gr for warm start next step
  call fR_save_old(ilevel)

  ! Add fifth force. Proper-frame F5 = +(c^2/2) grad(delta f_R);
  ! with scalar_gr = f_R (dimensionless) and box-unit gradients this
  ! becomes F5_code = +(a^2/2) (c/(H0*L_box))^2 * grad(f_R)
  ! (background gradient vanishes spatially). Sub-Compton unscreened
  ! limit recovers F5 = F_N/3 exactly.
  call compute_fifth_force(ilevel, &
       & 0.5d0*aexp**2*(2997.92458d0/boxlen_ini)**2)

end subroutine fR_solve_level

!=========================================================
! fR_rescale_warm_start: advance an initialized f_R field
! with the homogeneous-background ratio before Newton solve.
!=========================================================
subroutine fR_rescale_warm_start(ilevel, ratio)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::ratio

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(static)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           if(scalar_gr(icell) /= 0d0) &
                & scalar_gr(icell)=scalar_gr(icell)*ratio
        end do
     end do
  end do

end subroutine fR_rescale_warm_start

!=========================================================
! fR_init_scalar: initialize scalar_gr at a level
!=========================================================
subroutine fR_init_scalar(ilevel, fR_bar)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::fR_bar

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = fR_bar
           scalar_gr_old(icell) = fR_bar
        end do
     end do
  end do

end subroutine fR_init_scalar

!=========================================================
! fR_seed_scalar: set cells still at exactly 0 to fR_bar
! (uninitialized: first step, restart, newly refined grids)
!=========================================================
subroutine fR_seed_scalar(ilevel, fR_bar)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::fR_bar

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           if(scalar_gr(icell) == 0d0) then
              scalar_gr(icell) = fR_bar
              scalar_gr_old(icell) = fR_bar
           end if
        end do
     end do
  end do

end subroutine fR_seed_scalar

!=========================================================
! fR_save_old: save scalar_gr → scalar_gr_old
!=========================================================
subroutine fR_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do

end subroutine fR_save_old

!=========================================================
! fR_gauss_seidel: one Newton-GS sweep for f(R) equation
!
! Physical quasi-static equation (Oyaizu 2008; ECOSMOG):
!   nabla_com^2 delta f_R = (a^2/3c^2) [delta R - 8 pi G delta rho]
!
! In supercomoving box units (x in L_box, R~ = R/H0^2, rho = 1+delta):
!   lap f_R = +(a^2/3)*(H0*L/c)^2*(R~(f_R) - R~_bar)
!             - omega_m*(H0*L/c)^2*(rho - rho_tot)/a
! with (H0*L/c)^2 = (boxlen_ini/2997.92458)^2.
! Linearized mass term is +Yukawa (dR/df_R > 0 for f_R < 0), so the
! Newton Jacobian -6/dx^2 - (a^2/3)(H0L/c)^2 dR/df_R is negative
! definite. Chameleon limit: R~(f_R) -> R~_bar + 3 omega_m delta/a^3.
!
! Hu-Sawicki inversion: R = R_bar0 * (fR0/f_R)^(1/(n+1))
!=========================================================
subroutine fR_gauss_seidel(ilevel, R_bar, fR_bar, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::R_bar, fR_bar
  real(dp),intent(out)::res_max, src_max
  !-------------------------------------------------------
  ! Newton-GS update: u_new = u - F/J
  !   F = laplacian(u) - source(u)
  !   J = -6/dx² - dS/du
  ! Red-black ordering
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,u_abs,lapl,source,R_of_u,dR_du,residual,jacobian,delta_u
  real(dp)::u_nb_l,u_nb_r
  real(dp)::a2_over_3,rho_coeff,boxratio_sq
  real(dp)::R_bar0,fR0_abs,small_fR,inv_np1
  integer::np1,icolor

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  np1 = fR_n + 1
  fR0_abs=abs(fR0)
  small_fR=1d-30*fR0_abs
  inv_np1=1d0/dble(np1)
  ! (H0 L_box / c)^2 converts H0-unit curvature/density terms to
  ! box-unit Laplacians (same pattern as dark_energy_commons)
  boxratio_sq = (boxlen_ini / 2997.92458d0)**2
  a2_over_3 = aexp**2 * boxratio_sq / 3d0
  R_bar0 = 3d0 * (omega_m + 4d0 * omega_l)
  ! Matter source coefficient: omega_m*(H0 L/c)^2/a
  rho_coeff = omega_m * boxratio_sq / aexp

  ! Neighbor lookup
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  ! Red-black sweep (icolor=0: red, icolor=1: black)
  do icolor=0,1

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  u_c,u_abs,lapl,source,R_of_u,dR_du,residual,jacobian,delta_u,u_nb_l,u_nb_r) &
!$omp& reduction(max:res_max,src_max) schedule(static)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Gather neighbors
        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=vain_face_grid(igrid+i-1,2*idim-1)
              igridn_w(i,2*idim  )=vain_face_grid(igrid+i-1,2*idim  )
           end do
        end do

        do ind=1,twotondim
           ! Red-black coloring
           ! True 3D red-black: color = parity of (i+j+k) of the cell,
           ! i.e. popcount of the oct-local index (oct origins are even)
           if(mod(popcnt(ind-1), 2) /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              u_c = scalar_gr(ind_cell_w(i))
              u_abs = abs(u_c)

              ! Same-level neighbors, or parent-CIC Dirichlet data at
              ! a coarse-fine interface.
              lapl = 0d0
              do idim=1,ndim
                 ig_left =ggg(idim,1,ind)
                 ig_right=ggg(idim,2,ind)
                 ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
                 ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax
                 if(igridn_w(i,ig_left) > 0) then
                    u_nb_l=scalar_gr(igridn_w(i,ig_left)+ih_left)
                 else
                    call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,-1,u_c,u_nb_l)
                 end if
                 if(igridn_w(i,ig_right) > 0) then
                    u_nb_r=scalar_gr(igridn_w(i,ig_right)+ih_right)
                 else
                    call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,1,u_c,u_nb_r)
                 end if
                 lapl = lapl + (u_nb_l + u_nb_r - 2d0*u_c) * dx2_inv
              end do

              ! R(f_R) from Hu-Sawicki inversion:
              ! f_R = fR0 * (R_bar0/R)^(n+1)
              ! → R = R_bar0 * (fR0/f_R)^(1/(n+1))
              ! Since fR0<0 and f_R<0, use absolute values
              if(u_abs > small_fR) then
                 select case(np1)
                 case(1)
                    R_of_u = R_bar0 * fR0_abs/u_abs
                 case(2)
                    R_of_u = R_bar0 * sqrt(fR0_abs/u_abs)
                 case default
                    R_of_u = R_bar0 * (fR0_abs/u_abs)**inv_np1
                 end select
              else
                 R_of_u = R_bar  ! fallback to background
              end if

              ! Source (see subroutine header):
              !   S = +(a²/3)(H0L/c)²*(R(fR)-R_bar) - Ωm(H0L/c)²*(rho-rho_tot)/a
              ! rho has box mean rho_tot ≈ 1 (the Poisson solver subtracts
              ! it too); a2_over_3 and rho_coeff carry the (H0L/c)² factor.
              source = a2_over_3*(R_of_u - R_bar) &
                   & - rho_coeff*(rho(ind_cell_w(i)) - rho_tot)

              ! Residual: F = laplacian - source
              residual = lapl - source

              ! Jacobian of source w.r.t. u_c:
              ! dS/du = +(a²/3)(H0L/c)² * dR/dfR
              ! dR/dfR = -R/((n+1)*fR) > 0 for fR < 0
              if(u_abs > small_fR) then
                 dR_du = -R_of_u / (dble(np1) * u_c)
              else
                 dR_du = 0d0
              end if

              ! Jacobian: J = -6/dx² - dS/du < 0 (negative definite)
              jacobian = -6d0*dx2_inv - a2_over_3*dR_du

              ! Newton update
              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 ! Clamp update to prevent overshoot
                 if(abs(delta_u) > 0.5d0*u_abs .and. u_abs > small_fR) then
                    delta_u = sign(0.5d0*u_abs, delta_u)
                 end if
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
                 ! Enforce f_R < 0 (physical constraint)
                 if(scalar_gr(ind_cell_w(i)) > 0d0) &
                      & scalar_gr(ind_cell_w(i)) = -0.5d0*max(u_abs,small_fR)
              end if

              ! Track convergence
              res_max = max(res_max, abs(residual))
              src_max = max(src_max, abs(source))

           end do  ! i
        end do  ! ind
     end do  ! igrid

  end do  ! icolor

end subroutine fR_gauss_seidel

!=========================================================
! nDGP_beta: compute β(a) parameter
!=========================================================
function nDGP_beta(aa, orc, branch) result(beta)
  use amr_parameters, only: dp, omega_m, omega_l
  implicit none
  real(dp),intent(in)::aa, orc
  integer,intent(in)::branch
  real(dp)::beta
  !-------------------------------------------------------
  ! H(a) = H0 * sqrt(Ω_m/a³ + Ω_Λ)  [flat ΛCDM]
  ! Ḣ/H² = -(3/2)*(Ω_m/a³)/(Ω_m/a³ + Ω_Λ)
  ! r_c = 1/(2*sqrt(omega_rc)) / H0
  ! β = 1 + branch * 2*H*r_c*(1 + Ḣ/(3H²))
  !   = 1 + branch / sqrt(omega_rc) * sqrt(Ω_m/a³+Ω_Λ)
  !         * (1 - (Ω_m/a³)/(2*(Ω_m/a³+Ω_Λ)))
  !-------------------------------------------------------
  real(dp)::Oma3, E2, Hdot_over_H2, HrC

  Oma3 = omega_m / aa**3
  E2   = Oma3 + omega_l          ! H²/H0²
  Hdot_over_H2 = -1.5d0 * Oma3 / E2
  ! H*r_c = sqrt(E2)/H0 * 1/(2*sqrt(omega_rc)*H0) = sqrt(E2)/(2*sqrt(omega_rc))
  HrC = sqrt(E2) / (2d0 * sqrt(orc))
  beta = 1d0 + dble(branch) * 2d0 * HrC * (1d0 + Hdot_over_H2 / 3d0)

end function nDGP_beta

!=========================================================
! nDGP_solve_level: top-level nDGP solver for one AMR level
!=========================================================
subroutine nDGP_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
#ifdef HYDRO_CUDA
  use scalar_cuda_interface, only: SCAL_MODEL_NDGP
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount

  integer::iter,ncache,info,gs_iter_max
  real(dp)::beta
  real(dp)::nDGP_beta  ! external function
  logical::gscal_ok
  real(dp),dimension(1:12)::sparams
  real(dp)::gs_dx2
  integer::trk
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res,fft_rel
  logical::converged
#ifdef USE_FFTW
  logical::fft_attempted
  logical,external::level_fft_ok
#endif

  ! NOTE: do NOT return when ncache==0 — the relaxation loop and the
  ! FFT stage contain MPI collectives (ALLREDUCE) that every rank
  ! must enter; ranks without grids on this level contribute nothing.
  ncache=active(ilevel)%ngrid

  ! Compute β(a)
  beta = nDGP_beta(aexp, omega_rc, nDGP_branch)

  if(myid==1 .and. nstep==0 .and. ilevel==levelmin) then
     write(*,'(A,F8.4,A,F8.4)') ' nDGP: beta(a)=', beta, ' at a=', aexp
  end if

  ! Initialize scalar_gr on first step
  if(nstep==0) then
     call nDGP_init_scalar(ilevel)
  end if

  converged = .false.
#ifdef USE_FFTW
  fft_attempted = .false.
  ! Operator-split spectral solve on the uniform domain level
  if(level_fft_ok(ilevel)) then
     fft_attempted = .true.
     call vain_prepare_uniform_cache(ilevel)
     call make_virtual_fine_dp(scalar_gr(1), ilevel)
     do iter=1,50
        call vain_build_fft_rhs(ilevel, omega_m*aexp/beta, &
             & 1d0/(12d0*omega_rc*beta*aexp**4), .false., fft_rel)
        if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
             & ' nDGP level ',ilevel,' FFT iteration ',iter-1, &
             & ' residual=',fft_rel
        if(fft_rel < nDGP_eps) then
           converged=.true.
           if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
                & ' nDGP level ',ilevel,' FFT converged in ',iter-1, &
                & ' iters, res=',fft_rel
           exit
        end if
        call level_fft_helmholtz(ilevel, 0d0, 0d0, 1d0)
        call make_virtual_fine_dp(scalar_gr(1), ilevel)
     end do
  end if

  ! On the fully refined periodic base level the distributed spectral
  ! solver is the production algorithm.  Silently starting thousands of
  ! global Newton-GS sweeps after a failed FFT iteration made 1024^3 runs
  ! look hung while consuming an entire node.  Strict production runs now
  ! stop with the last residual; non-strict diagnostic runs retain only a
  ! bounded GS safety path below.
  if(ilevel==levelmin .and. .not.fft_attempted .and. scalar_solver_strict) then
     if(myid==1) write(*,'(A,I2,A)') &
          & ' ERROR: nDGP level ',ilevel, &
          & ' scalar FFT unavailable for this MPI layout; refusing global GS fallback'
     call scalar_solver_abort
  end if
  if(fft_attempted .and. .not.converged .and. scalar_solver_strict) then
     if(myid==1) write(*,'(A,I2,A,ES10.3)') &
          & ' ERROR: nDGP level ',ilevel, &
          & ' scalar FFT failed; refusing global GS fallback, residual=',fft_rel
     call scalar_solver_abort
  end if
#endif

  ! Newton-GS relaxation
  gs_iter_max=n_iter_nDGP
#ifdef USE_FFTW
  if(fft_attempted) gs_iter_max=min(gs_iter_max,100)
#endif
  gscal_ok=.false.
  trk=0
  if(.not.converged) then
#ifdef HYDRO_CUDA
  call scalar_gpu_begin(ilevel, .true., gscal_ok)
  if(gscal_ok) then
     gs_dx2=(0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1))**2
     sparams=0d0
     sparams(1)=1d0/gs_dx2
     sparams(2)=1d0/(12d0*omega_rc*beta*aexp**4)
     sparams(3)=omega_m*aexp/beta
     sparams(4)=rho_tot
     sparams(5)=1d-2*omega_m*aexp/beta
     sparams(6)=gs_dx2
     if(galileon_tracker) trk=1
  end if
#endif
  do iter=1,gs_iter_max

     if(myid==1 .and. (iter==1 .or. mod(iter,10)==0)) &
          & write(*,'(A,I2,A,I5,A,I5)') ' nDGP level ',ilevel, &
          & ' starting GS iteration ',iter,' / ',gs_iter_max

#ifdef HYDRO_CUDA
     if(gscal_ok) then
        call scalar_gpu_sweep_halo(ilevel, SCAL_MODEL_NDGP, sparams, trk, &
             & res_max_local, src_max_local)
     else
#endif
     call nDGP_gauss_seidel(ilevel, beta, res_max_local, src_max_local)

     call make_virtual_fine_dp(scalar_gr(1), ilevel)
#ifdef HYDRO_CUDA
     end if
#endif

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < nDGP_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
             ' nDGP level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do
#ifdef HYDRO_CUDA
  if(gscal_ok) call scalar_gpu_end(ilevel)
#endif
  end if

  if(.not. converged) then
     if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
          & ' WARNING: nDGP level ',ilevel,' NOT converged after ', &
          & gs_iter_max,' iters, res=',rel_res
     if(scalar_solver_strict) call scalar_solver_abort
  end if

  ! Save scalar_gr for warm start
  call nDGP_save_old(ilevel)

  ! Fifth force: F5 = -(1/2) * grad(φ)
  ! (the 1/β coupling is already in the field-equation source;
  !  linear limit must give G_eff/G = 1 + 1/(3β), Winther+15 eq. 16-17)
  call compute_fifth_force(ilevel, -0.5d0)

end subroutine nDGP_solve_level

!=========================================================
! nDGP_init_scalar: initialize scalar_gr = 0 at a level
!=========================================================
subroutine nDGP_init_scalar(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = 0d0
           scalar_gr_old(icell) = 0d0
        end do
     end do
  end do

end subroutine nDGP_init_scalar

!=========================================================
! nDGP_save_old: save scalar_gr → scalar_gr_old
!=========================================================
subroutine nDGP_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do

end subroutine nDGP_save_old

!=========================================================
! nDGP_gauss_seidel: one Newton-GS sweep for nDGP equation
!
! ∇²φ + coeff*[(∇²φ)² - (∇ᵢ∇ⱼφ)²] = source
!
! coeff = r_c²/(3β a² c²) → in code units (c=1, H0=1):
!   coeff = 1/(12 * omega_rc * beta * a²)
!
! source = 8πG/(3β) * a² * δρ
!   = 2*fourpi/(3β) * δ where fourpi=1.5*Ωm*a
!   = Ωm*a / β * δ
!
! Mixed derivatives via diagonal neighbors (double morton_nbor_grid)
!=========================================================
subroutine nDGP_gauss_seidel(ilevel, beta, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::beta
  real(dp),intent(out)::res_max, src_max

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,source,residual,jacobian,delta_u,sclamp
  real(dp)::coeff
  real(dp)::phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp
  real(dp)::phi_xx,phi_yy,phi_zz,mix_xy2,mix_xz2,mix_yz2
  real(dp)::lapl2,trace_ij2,vain_term
  real(dp)::dpp,dpm,dmp,dmm,dmix_du
  integer::icolor

  real(dp)::phi_pp,phi_pm,phi_mp,phi_mm

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  ! Vainshtein coefficient: physical r_c²/(3β a²) with
  ! r_c² = 1/(4*omega_rc*H0²). In supercomoving units the quadratic
  ! term carries one extra (H0²/a²) from the second Laplacian, so
  ! coeff_code = [1/(12 omega_rc beta a²)]*(H0²/a²)/H0² = 1/(12 Ω_rc β a⁴)
  coeff = 1d0 / (12d0 * omega_rc * beta * aexp**4)

  ! Neighbor lookup
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  do icolor=0,1

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  u_c,lapl,source,residual,jacobian,delta_u, &
!$omp&  phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp, &
!$omp&  phi_xx,phi_yy,phi_zz,mix_xy2,mix_xz2,mix_yz2, &
!$omp&  lapl2,trace_ij2,vain_term,dpp,dpm,dmp,dmm,dmix_du, &
!$omp&  phi_pp,phi_pm,phi_mp,phi_mm) &
!$omp& reduction(max:res_max,src_max) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Gather face neighbors
        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
              igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
           end do
        end do

        do ind=1,twotondim
           ! Compatible mixed-derivative squares are relaxed with the
           ! same red-black ordering used by the original nDGP path.
           if(mod(popcnt(ind-1),2) /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              u_c = scalar_gr(ind_cell_w(i))

              ! Parent-CIC Dirichlet closure at coarse-fine interfaces.
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0, 0,u_c,phi_xm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0, 0,u_c,phi_xp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1, 0,u_c,phi_ym)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1, 0,u_c,phi_yp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 0,-1,u_c,phi_zm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 0, 1,u_c,phi_zp)

              ! Laplacian
              lapl = (phi_xp + phi_xm + phi_yp + phi_ym + phi_zp + phi_zm - 6d0*u_c) * dx2_inv

              ! Diagonal second derivatives
              phi_xx = (phi_xp + phi_xm - 2d0*u_c) * dx2_inv
              phi_yy = (phi_yp + phi_ym - 2d0*u_c) * dx2_inv
              phi_zz = (phi_zp + phi_zm - 2d0*u_c) * dx2_inv

              ! Full mixed Hessian from centered diagonal samples.
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 1, 0,u_c,phi_pp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1,-1, 0,u_c,phi_pm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 1, 0,u_c,phi_mp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1,-1, 0,u_c,phi_mm)
              dpp=(phi_pp-phi_xp-phi_yp+u_c)*dx2_inv
              dpm=(phi_xp-phi_pm-u_c+phi_ym)*dx2_inv
              dmp=(phi_yp-phi_mp-u_c+phi_xm)*dx2_inv
              dmm=(u_c-phi_ym-phi_xm+phi_mm)*dx2_inv
              mix_xy2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
              if(galileon_tracker) &
                   & mix_xy2=(0.25d0*(phi_pp-phi_pm-phi_mp+phi_mm)*dx2_inv)**2
              dmix_du=0.5d0*dx2_inv*(dpp-dpm-dmp+dmm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0, 1,u_c,phi_pp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0,-1,u_c,phi_pm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0, 1,u_c,phi_mp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0,-1,u_c,phi_mm)
              dpp=(phi_pp-phi_xp-phi_zp+u_c)*dx2_inv
              dpm=(phi_xp-phi_pm-u_c+phi_zm)*dx2_inv
              dmp=(phi_zp-phi_mp-u_c+phi_xm)*dx2_inv
              dmm=(u_c-phi_zm-phi_xm+phi_mm)*dx2_inv
              mix_xz2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
              if(galileon_tracker) &
                   & mix_xz2=(0.25d0*(phi_pp-phi_pm-phi_mp+phi_mm)*dx2_inv)**2
              dmix_du=dmix_du+0.5d0*dx2_inv*(dpp-dpm-dmp+dmm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1, 1,u_c,phi_pp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1,-1,u_c,phi_pm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1, 1,u_c,phi_mp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1,-1,u_c,phi_mm)
              dpp=(phi_pp-phi_yp-phi_zp+u_c)*dx2_inv
              dpm=(phi_yp-phi_pm-u_c+phi_zm)*dx2_inv
              dmp=(phi_zp-phi_mp-u_c+phi_ym)*dx2_inv
              dmm=(u_c-phi_zm-phi_ym+phi_mm)*dx2_inv
              mix_yz2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
              if(galileon_tracker) &
                   & mix_yz2=(0.25d0*(phi_pp-phi_pm-phi_mp+phi_mm)*dx2_inv)**2
              dmix_du=dmix_du+0.5d0*dx2_inv*(dpp-dpm-dmp+dmm)

              ! Vainshtein operator: (trace H)^2 - trace(H^2).
              lapl2 = lapl**2
              trace_ij2 = phi_xx**2 + phi_yy**2 + phi_zz**2 &
                   & + 2d0*(mix_xy2+mix_xz2+mix_yz2)
              vain_term = lapl2 - trace_ij2

              ! Source = Ωm*a/β * δ  (in code units)
              ! rho has box mean rho_tot ≈ 1; subtract it (δ = rho - rho_tot)
              source = omega_m * aexp / beta * (rho(ind_cell_w(i)) - rho_tot)

              ! Residual: F = lapl + coeff * vain_term - source
              residual = lapl + coeff * vain_term - source

              ! Jacobian: dF/du_c
              ! dlapl/du_c = -6/dx²
              ! d(vain_term)/du_c = d(lapl²-trace_ij²)/du_c
              !   d(lapl²)/du_c = 2*lapl*(-6*dx2_inv) ... but actually dlapl/du = -(2*ndim)/dx²
              ! Let's use: dlapl/du = -6/dx² (each of 6 face neighbors has +1, center has -6)
              ! d(lapl²)/du = 2*lapl*(-6/dx²)
              ! d(phi_xx²)/du = 2*phi_xx*(-2/dx²), similarly for yy, zz
              ! d(trace_ij2)/du = 2*(-2/dx²)*(phi_xx+phi_yy+phi_zz) = 2*(-2/dx²)*lapl
              ! The compatible mixed stencil contributes dmix_du.
              jacobian = -6d0*dx2_inv + coeff * &
                   & (-8d0*lapl*dx2_inv-2d0*dmix_du)

              ! Newton update
              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 ! Damped Newton; floor the clamp scale at 1% of the
                 ! delta=1 source so cells with delta~0 still relax
                 sclamp = 0.5d0*dx2*max(abs(source), 1d-2*omega_m*aexp/beta)
                 if(abs(delta_u) > sclamp) delta_u = sign(sclamp, delta_u)
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
              end if

              res_max = max(res_max, abs(residual))
              src_max = max(src_max, abs(source))

           end do  ! i
        end do  ! ind
     end do  ! igrid

  end do  ! icolor

end subroutine nDGP_gauss_seidel

!#########################################################
!#########################################################
!  SYMMETRON SCALAR FIELD SOLVER
!#########################################################
!#########################################################

!=========================================================
! symmetron_solve_level: top-level Symmetron solver
!=========================================================
subroutine symmetron_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
#ifdef HYDRO_CUDA
  use scalar_cuda_interface, only: SCAL_MODEL_SYMMETRON
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount

  integer::iter,ncache,info
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res
  logical::converged
  logical::gscal_ok
  real(dp),dimension(1:12)::sparams
  real(dp)::gs_dx2
#ifdef USE_FFTW
  real(dp)::m2bar
  logical,external::level_fft_ok
#endif

  ! NOTE: do NOT return when ncache==0 — the relaxation loop and the
  ! FFT stage contain MPI collectives (ALLREDUCE) that every rank
  ! must enter; ranks without grids on this level contribute nothing.
  ncache=active(ilevel)%ngrid

  ! Seed the broken-phase VEV chi_bar(a) = sqrt(1-(a_ssb/a)^3) into
  ! cells still at exactly 0. The field equation is homogeneous in
  ! chi, so chi=0 is an exact (unphysical) fixed point of Newton-GS:
  ! without seeding the fifth force stays identically zero. Covers
  ! first step, restarts and newly refined cells; pre-SSB (a<=a_ssb)
  ! chi_bar=0 and chi=0 is the true solution.
  call symmetron_seed_scalar(ilevel)
  call vain_prepare_uniform_cache(ilevel)

#ifdef USE_FFTW
  ! Spectral Newton on the uniform domain level (broken phase only)
  if(level_fft_ok(ilevel) .and. aexp > a_ssb) then
     call make_virtual_fine_dp(scalar_gr(1), ilevel)
     do iter=1,3
        call sb_build_fft_rhs(ilevel, a_ssb, L_symmetron, m2bar)
        call level_fft_helmholtz(ilevel, m2bar, 0d0, 1d0)
        call make_virtual_fine_dp(scalar_gr(1), ilevel)
     end do
  end if
#endif

  ! Newton-GS relaxation (GPU sweeps when gpu_scalar is active)
  gscal_ok=.false.
#ifdef HYDRO_CUDA
  call scalar_gpu_begin(ilevel, .false., gscal_ok)
  if(gscal_ok) then
     gs_dx2=(0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1))**2
     sparams=0d0
     sparams(1)=1d0/gs_dx2
     sparams(2)=aexp**2/(2d0*(L_symmetron/boxlen_ini)**2)
     sparams(3)=(a_ssb/aexp)**3
  end if
#endif
  converged = .false.
  do iter=1,n_iter_symmetron

#ifdef HYDRO_CUDA
     if(gscal_ok) then
        call scalar_gpu_sweep_halo(ilevel, SCAL_MODEL_SYMMETRON, sparams, 0, &
             & res_max_local, src_max_local)
     else
#endif
     call symmetron_gauss_seidel(ilevel, res_max_local, src_max_local)

     call make_virtual_fine_dp(scalar_gr(1), ilevel)
#ifdef HYDRO_CUDA
     end if
#endif

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < symmetron_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I6,A,ES10.3)') &
             ' Symmetron level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do
#ifdef HYDRO_CUDA
  if(gscal_ok) call scalar_gpu_end(ilevel)
#endif

  if(.not. converged) then
     if(myid==1) write(*,'(A,I2,A,I6,A,ES10.3)') &
          & ' WARNING: Symmetron level ',ilevel,' NOT converged after ', &
          & n_iter_symmetron,' iters, res=',rel_res
     if(scalar_solver_strict) call scalar_solver_abort
  end if

  ! Save for warm start
  call symmetron_save_old(ilevel)

  ! Fifth force: F5 = -6*Ωm*β²*(L/L_box)²*(a²/a_ssb³)*χ*∇χ
  call compute_fifth_force_symmetron(ilevel)

end subroutine symmetron_solve_level

!=========================================================
! symmetron_init_scalar: initialize scalar_gr = 0
!=========================================================
subroutine symmetron_init_scalar(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = 0d0
           scalar_gr_old(icell) = 0d0
        end do
     end do
  end do

end subroutine symmetron_init_scalar

!=========================================================
! symmetron_seed_scalar: seed cells still at exactly 0 with
! the broken-phase VEV chi_bar(a) (0 pre-SSB). chi=0 is an
! exact fixed point of the homogeneous equation, so seeding
! is required for the field to leave the symmetric phase.
!=========================================================
subroutine symmetron_seed_scalar(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell
  real(dp)::chibar

  chibar = 0d0
  if(aexp > a_ssb) chibar = sqrt(max(1d0 - (a_ssb/aexp)**3, 0d0))
  if(chibar == 0d0) return

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           if(scalar_gr(icell) == 0d0) then
              scalar_gr(icell) = chibar
              scalar_gr_old(icell) = chibar
           end if
        end do
     end do
  end do

end subroutine symmetron_seed_scalar

!=========================================================
! symmetron_save_old: save scalar_gr → scalar_gr_old
!=========================================================
subroutine symmetron_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do

end subroutine symmetron_save_old

!=========================================================
! symmetron_gauss_seidel: one Newton-GS sweep
! (Davis, Li, Mota, Winther 2012 / ISIS conventions)
!
! ∇²χ = (a²/2L²)[(ρ/ρ_ssb) - 1]·χ + (a²/2L²)·χ³
!
! In code units:
!   L_code = L_symmetron / boxlen_ini (Mpc/h → box units)
!   ρ_ssb = mean matter density at a_ssb
!   rho(cell) = ρ/ρ̄(a) = 1+δ (box mean rho_tot ≈ 1)
!   ρ/ρ_ssb = rho(cell)*(a_ssb/a)³
!
! Broken-phase VEV: χ̄² = 1-(a_ssb/a)³ for a > a_ssb.
!=========================================================
subroutine symmetron_gauss_seidel(ilevel, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(out)::res_max, src_max

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,residual,jacobian,delta_u,u_nb_l,u_nb_r
  real(dp)::L_code,L2_inv,a2_over_2L2
  real(dp)::rho_ratio,mass_term,cubic_coeff
  integer::icolor

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  ! Convert L_symmetron from Mpc/h to code units
  L_code = L_symmetron / boxlen_ini
  L2_inv = 1d0 / L_code**2
  a2_over_2L2 = aexp**2 * L2_inv / 2d0
  cubic_coeff = a2_over_2L2

  ! Neighbor lookup
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  do icolor=0,1

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  u_c,lapl,residual,jacobian,delta_u,u_nb_l,u_nb_r, &
!$omp&  rho_ratio,mass_term) &
!$omp& reduction(max:res_max,src_max) schedule(static)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=vain_face_grid(igrid+i-1,2*idim-1)
              igridn_w(i,2*idim  )=vain_face_grid(igrid+i-1,2*idim  )
           end do
        end do

        do ind=1,twotondim
           ! True 3D red-black: color = parity of (i+j+k) of the cell,
           ! i.e. popcount of the oct-local index (oct origins are even)
           if(mod(popcnt(ind-1), 2) /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              u_c = scalar_gr(ind_cell_w(i))

              ! Same-level neighbours are already gathered for this vector.
              ! Use the Morton/parent-CIC path only at a true AMR boundary.
              lapl = 0d0
              do idim=1,ndim
                 ig_left =ggg(idim,1,ind)
                 ig_right=ggg(idim,2,ind)
                 ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
                 ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax
                 if(igridn_w(i,ig_left) > 0) then
                    u_nb_l=scalar_gr(igridn_w(i,ig_left)+ih_left)
                 else
                    call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,-1,u_c,u_nb_l)
                 end if
                 if(igridn_w(i,ig_right) > 0) then
                    u_nb_r=scalar_gr(igridn_w(i,ig_right)+ih_right)
                 else
                    call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,1,u_c,u_nb_r)
                 end if
                 lapl=lapl+u_nb_l+u_nb_r
              end do

              lapl = (lapl - 6d0*u_c) * dx2_inv

              ! ρ/ρ_ssb = rho*(a_ssb/a)³  (rho = ρ/ρ̄(a) = 1+δ already)
              rho_ratio = rho(ind_cell_w(i)) * (a_ssb/aexp)**3

              ! mass_term = a2_over_2L2 * (rho_ratio - 1) (relative to VEV)
              mass_term = a2_over_2L2 * (rho_ratio - 1d0)

              ! F = lapl - mass_term*χ - cubic_coeff*χ³
              residual = lapl - mass_term * u_c - cubic_coeff * u_c**3

              ! dF/du = -6/dx² - mass_term - 3*cubic_coeff*χ²
              jacobian = -6d0*dx2_inv - mass_term - 3d0*cubic_coeff*u_c**2

              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
              end if

              src_max = max(src_max, abs(mass_term*u_c) + abs(cubic_coeff*u_c**3))
              res_max = max(res_max, abs(residual))

           end do  ! i
        end do  ! ind
     end do  ! igrid

  end do  ! icolor

end subroutine symmetron_gauss_seidel

!=========================================================
! compute_fifth_force_symmetron: field-dependent fifth force
! F₅_d = -6·Ωm·β²·(a_ssb/a)·χ·∂χ/∂x_d
! (field-dependent coupling → cannot use generic compute_fifth_force)
!=========================================================
subroutine compute_fifth_force_symmetron(ilevel)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc
  integer::nx_loc
  real(dp)::grad_u,chi_c,factor,u_left,u_right

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ! NOTE: no early return on ncache==0 — the final make_virtual_fine_dp
  ! must be entered by every rank (matched communication)
  ncache=active(ilevel)%ngrid

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Factor: -6*Ωm*β²*(L_symmetron/L_box)²*a²/a_ssb³
  ! Derived from A(φ)=1+φ²/2M², φ_ssb²/M² = 6Ωmβ²H0²λ0²/a_ssb³ and
  ! the supercomoving force conversion f_code = a³ g/(L_box H0²);
  ! unscreened small-scale limit then gives F5/F_N = 2β²χ̄² exactly.
  factor = -6d0 * omega_m * beta_symmetron**2 &
       & * (L_symmetron/boxlen_ini)**2 * aexp**2 / a_ssb**3

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right,grad_u,chi_c,u_left,u_right) schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=vain_face_grid(igrid+i-1,2*idim-1)
           igridn_w(i,2*idim  )=vain_face_grid(igrid+i-1,2*idim  )
        end do
     end do

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do

        do i=1,ngrid
           chi_c = scalar_gr(ind_cell_w(i))
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax

              if(igridn_w(i,ig_left) > 0) then
                 u_left=scalar_gr(igridn_w(i,ig_left)+ih_left)
              else
                 call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,-1,chi_c,u_left)
              end if
              if(igridn_w(i,ig_right) > 0) then
                 u_right=scalar_gr(igridn_w(i,ig_right)+ih_right)
              else
                 call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,1,chi_c,u_right)
              end if
              grad_u=(u_right-u_left)/(2d0*dx_loc)
              f(ind_cell_w(i),idim) = f(ind_cell_w(i),idim) + factor * chi_c * grad_u
           end do
        end do
     end do
  end do

  ! Update MPI virtual boundaries (f was synced before F5 was added)
  do idim=1,ndim
     call make_virtual_fine_dp(f(1,idim),ilevel)
  end do

end subroutine compute_fifth_force_symmetron

!#########################################################
!#########################################################
!  DILATON SCALAR FIELD SOLVER
!#########################################################
!#########################################################

!=========================================================
! dilaton_solve_level: top-level Dilaton solver
!=========================================================
subroutine dilaton_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
#ifdef HYDRO_CUDA
  use scalar_cuda_interface, only: SCAL_MODEL_DILATON
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount
  !-------------------------------------------------------
  ! Environment-dependent dilaton (Brax, van de Bruck, Davis,
  ! Shaw 2010; N-body form: Brax, Davis, Li, Winther, Zhao 2012,
  ! arXiv:1206.3568, original r=3/2 model):
  !   A(phi) = 1 + (A2/2) phi^2/Mpl^2,  V = V0 exp(-gamma phi/Mpl)
  !   m^2(a) = 3 A2 H^2(a)  =>  xi = H0/m0 = 1/sqrt(3 A2)
  !   beta(a) = beta0 * a^s,  s = 9 Omega_m A2 xi^2 = 3 Omega_m
  ! Working variable chi = phi/Mpl; background chibar = beta0 a^s/A2.
  ! Code-unit field equation (x in box units, rho mean-normalized):
  !   lap chi = (3 Om A2 B2/a) (rho*chi - chibar)
  !           + a^2 B2 [vphi(chi) - vphi(chibar)]
  !   vphi(chi) = -3 Om beta0 (A2 chi/beta0)^(1-3/s)
  ! with B2 = (boxlen_ini/2997.92458)^2 = 1/ctilde^2.
  ! Parameters: beta_dilaton = beta0 (cosmological coupling today),
  ! L_dilaton = 2998*xi = fifth-force range today [Mpc/h].
  ! (a0_dilaton is ignored — legacy of the old symmetron-clone.)
  !-------------------------------------------------------
  integer::iter,ncache,info
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res
  logical::converged
  real(dp)::xi_d,A2_d,s_d,chibar_d,fac5
  logical::gscal_ok
  real(dp),dimension(1:12)::sparams
  real(dp)::gs_dx2
#ifdef USE_FFTW
  real(dp)::m2bar
  logical,external::level_fft_ok
#endif

  ! NOTE: do NOT return when ncache==0 — the relaxation loop and the
  ! FFT stage contain MPI collectives (ALLREDUCE) that every rank
  ! must enter; ranks without grids on this level contribute nothing.
  ncache=active(ilevel)%ngrid

  ! Model constants from (beta0, range)
  xi_d = L_dilaton/2997.92458d0
  A2_d = 1d0/(3d0*xi_d**2)
  s_d  = 3d0*omega_m
  chibar_d = beta_dilaton*aexp**s_d/A2_d

  ! Seed cells still at exactly 0 with the background value
  ! (chi=0 is a singular point of vphi; also covers restarts and
  ! newly refined cells)
  call dilaton_seed_scalar(ilevel, chibar_d)

#ifdef USE_FFTW
  ! Spectral Newton on the uniform domain level
  if(level_fft_ok(ilevel)) then
     call make_virtual_fine_dp(scalar_gr(1), ilevel)
     do iter=1,3
        call dil_build_fft_rhs(ilevel, A2_d, s_d, chibar_d, m2bar)
        call level_fft_helmholtz(ilevel, m2bar, 0.25d0, 1d0)
        call make_virtual_fine_dp(scalar_gr(1), ilevel)
     end do
  end if
#endif

  gscal_ok=.false.
#ifdef HYDRO_CUDA
  call scalar_gpu_begin(ilevel, .false., gscal_ok)
  if(gscal_ok) then
     gs_dx2=(0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1))**2
     sparams=0d0
     sparams(1)=1d0/gs_dx2
     sparams(2)=3d0*omega_m*A2_d*(boxlen_ini/2997.92458d0)**2/aexp
     sparams(3)=aexp**2*(boxlen_ini/2997.92458d0)**2
     sparams(4)=1d0-3d0/s_d
     sparams(5)=-3d0*omega_m*beta_dilaton &
          & *(A2_d*chibar_d/beta_dilaton)**(1d0-3d0/s_d)
     sparams(6)=chibar_d
     sparams(7)=A2_d/beta_dilaton
     sparams(8)=-3d0*omega_m*beta_dilaton
     sparams(9)=-3d0*omega_m*A2_d*(1d0-3d0/s_d)
     sparams(10)=1d-30*chibar_d
  end if
#endif
  converged = .false.
  do iter=1,n_iter_dilaton

#ifdef HYDRO_CUDA
     if(gscal_ok) then
        call scalar_gpu_sweep_halo(ilevel, SCAL_MODEL_DILATON, sparams, 0, &
             & res_max_local, src_max_local)
     else
#endif
     call dilaton_gauss_seidel(ilevel, A2_d, s_d, chibar_d, &
          & res_max_local, src_max_local)

     call make_virtual_fine_dp(scalar_gr(1), ilevel)
#ifdef HYDRO_CUDA
     end if
#endif

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < dilaton_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
             ' Dilaton level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do
#ifdef HYDRO_CUDA
  if(gscal_ok) call scalar_gpu_end(ilevel)
#endif

  if(.not. converged) then
     if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
          & ' WARNING: Dilaton level ',ilevel,' NOT converged after ', &
          & n_iter_dilaton,' iters, res=',rel_res
     if(scalar_solver_strict) call scalar_solver_abort
  end if

  call dilaton_save_old(ilevel)

  ! Fifth force: F5 = -ctilde^2 a^2 A2 * chi * grad(chi)
  ! (unscreened linear limit gives F5/FN = 2 beta(a)^2 exactly)
  fac5 = -(2997.92458d0/boxlen_ini)**2 * aexp**2 * A2_d
  call compute_fifth_force_dilaton(ilevel, fac5)

end subroutine dilaton_solve_level

!=========================================================
! dilaton_init_scalar / dilaton_save_old
!=========================================================
subroutine dilaton_init_scalar(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = 0d0
           scalar_gr_old(icell) = 0d0
        end do
     end do
  end do
end subroutine dilaton_init_scalar

!=========================================================
! dilaton_seed_scalar: seed cells still at exactly 0 with the
! broken-phase VEV (see symmetron_seed_scalar)
!=========================================================
subroutine dilaton_seed_scalar(ilevel, chibar_in)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::chibar_in

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           if(scalar_gr(icell) == 0d0) then
              scalar_gr(icell) = chibar_in
              scalar_gr_old(icell) = chibar_in
           end if
        end do
     end do
  end do

end subroutine dilaton_seed_scalar

subroutine dilaton_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do
end subroutine dilaton_save_old

!=========================================================
! dilaton_gauss_seidel: one Newton-GS sweep for the Brax+12
! dilaton equation (see dilaton_solve_level header):
!   F = lap chi - cA*(rho*chi - chibar) - cV*[vphi(chi)-vphi(chibar)]
!   vphi(chi) = -3 Om beta0 wfac^pexp, wfac = A2 chi/beta0,
!   pexp = 1 - 3/s < 0  =>  dvphi/dchi > 0 and the Newton
!   Jacobian -6/h^2 - dS/dchi is negative definite.
! chi > 0 is enforced with the same halving guard as f(R).
!=========================================================
subroutine dilaton_gauss_seidel(ilevel, A2_d, s_d, chibar_d, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::A2_d,s_d,chibar_d
  real(dp),intent(out)::res_max, src_max

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,residual,jacobian,delta_u,u_nb_l,u_nb_r
  real(dp)::boxratio_sq,cA,cV,pexp,wfac,vphi,dvphi,vbar,source
  integer::icolor

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  boxratio_sq=(boxlen_ini/2997.92458d0)**2
  cA = 3d0*omega_m*A2_d*boxratio_sq/aexp
  cV = aexp**2*boxratio_sq
  pexp = 1d0 - 3d0/s_d
  vbar = -3d0*omega_m*beta_dilaton*(A2_d*chibar_d/beta_dilaton)**pexp

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  do icolor=0,1

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  u_c,lapl,residual,jacobian,delta_u,u_nb_l,u_nb_r, &
!$omp&  wfac,vphi,dvphi,source) &
!$omp& reduction(max:res_max,src_max) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
              igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
           end do
        end do

        do ind=1,twotondim
           ! True 3D red-black: color = parity of (i+j+k) of the cell
           if(mod(popcnt(ind-1), 2) /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              u_c = scalar_gr(ind_cell_w(i))

              ! Parent-CIC Dirichlet data at coarse-fine interfaces.
              lapl = 0d0
              do idim=1,ndim
                 call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,-1,u_c,u_nb_l)
                 call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim, 1,u_c,u_nb_r)
                 lapl=lapl+u_nb_l+u_nb_r
              end do
              lapl = (lapl - 6d0*u_c) * dx2_inv

              wfac = A2_d*max(u_c,1d-30)/beta_dilaton
              vphi = -3d0*omega_m*beta_dilaton*wfac**pexp
              dvphi = -3d0*omega_m*A2_d*pexp*wfac**(pexp-1d0)

              source = cA*(rho(ind_cell_w(i))*u_c - chibar_d) &
                   & + cV*(vphi - vbar)

              residual = lapl - source
              jacobian = -6d0*dx2_inv - cA*rho(ind_cell_w(i)) - cV*dvphi

              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 ! Clamp and keep chi > 0 (vphi is singular at 0)
                 if(abs(delta_u) > 0.5d0*abs(u_c)) &
                      & delta_u = sign(0.5d0*abs(u_c), delta_u)
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
                 if(scalar_gr(ind_cell_w(i)) <= 0d0) &
                      & scalar_gr(ind_cell_w(i)) = 0.5d0*max(abs(u_c),1d-30*chibar_d)
              end if

              src_max = max(src_max, abs(source))
              res_max = max(res_max, abs(residual))

           end do
        end do
     end do

  end do

end subroutine dilaton_gauss_seidel

!=========================================================
! compute_fifth_force_dilaton: F₅ = -2β²·φ_D·∇φ_D
!=========================================================
subroutine compute_fifth_force_dilaton(ilevel, factor_in)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::factor_in

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc
  integer::nx_loc
  real(dp)::grad_u,phi_c,factor,u_left,u_right

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ! NOTE: no early return on ncache==0 — the final make_virtual_fine_dp
  ! must be entered by every rank (matched communication)
  ncache=active(ilevel)%ngrid

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! F5 = -ctilde^2 a^2 A2 * chi * grad(chi)  (see dilaton_solve_level)
  factor = factor_in

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right,grad_u,phi_c,u_left,u_right) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do

        do i=1,ngrid
           phi_c = scalar_gr(ind_cell_w(i))
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax

              call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim,-1,phi_c,u_left)
              call scalar_sample_axis(ind_grid_w(i),ind,ilevel,idim, 1,phi_c,u_right)
              grad_u=(u_right-u_left)/(2d0*dx_loc)
              f(ind_cell_w(i),idim) = f(ind_cell_w(i),idim) + factor * phi_c * grad_u
           end do
        end do
     end do
  end do

  ! Update MPI virtual boundaries (f was synced before F5 was added)
  do idim=1,ndim
     call make_virtual_fine_dp(f(1,idim),ilevel)
  end do

end subroutine compute_fifth_force_dilaton

!#########################################################
!#########################################################
!  GALILEON (CUBIC) SCALAR FIELD SOLVER
!#########################################################
!#########################################################

!=========================================================
! galileon_beta: β_G(a) = c₂/(6·c₃·H)
! H(a) = H0*sqrt(Ωm/a³ + ΩΛ) in code units (H0=1)
!=========================================================
function galileon_beta(aa)
  use amr_parameters, only: dp, c2_galileon, c3_galileon, omega_m, omega_l
  implicit none
  real(dp)::galileon_beta
  real(dp),intent(in)::aa
  real(dp)::Ha

  Ha = sqrt(omega_m / aa**3 + omega_l)
  galileon_beta = c2_galileon / (6d0 * c3_galileon * Ha)

end function galileon_beta

!=========================================================
! galileon_solve_level: top-level cubic Galileon solver
! Reuses nDGP Vainshtein structure with different coeff/β
!=========================================================
subroutine galileon_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
#ifdef HYDRO_CUDA
  use scalar_cuda_interface, only: SCAL_MODEL_GALILEON
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount

  integer::iter,ncache,info
  real(dp)::beta_G,coeff_G
  real(dp)::galileon_beta  ! external function
  logical::gscal_ok
  real(dp),dimension(1:12)::sparams
  real(dp)::gs_dx2
  integer::trk
  real(dp)::Ha
  real(dp)::xi_t,E2_t,hd_t,beta1_t
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res,fft_rel
  logical::converged
#ifdef USE_FFTW
  logical,external::level_fft_ok
#endif

  ! NOTE: do NOT return when ncache==0 — the relaxation loop and the
  ! FFT stage contain MPI collectives (ALLREDUCE) that every rank
  ! must enter; ranks without grids on this level contribute nothing.
  ncache=active(ilevel)%ngrid

  ! Compute β_G(a) and Vainshtein coeff
  if(galileon_tracker) then
     ! Barreira+13 tracker (parameter-free): H*phidot = xi*H0^2*Mpl,
     ! c2 = -6*c3*xi, xi = sqrt(6(1-Om)), c3 = 1/(6 xi). Closed forms:
     !   E^2(a) = [Om a^-3 + sqrt(Om^2 a^-6 + 4(1-Om))]/2
     !   Hdot/H^2 = -(3/2) Om a^-3 / (2E^2 - Om a^-3)
     !   beta1 = (xi/3)[2 Hdot/H^2 - 1 + (1-Om)/E^4]   (Barreira eq. 14)
     !   beta2 = 2 E^2 beta1 / xi^2                    (Barreira eq. 15)
     ! Code-unit field eq (u = a^2 phi/(Mpl H0^2 L^2)):
     !   lap u + [1/(3 beta1 a^4)][(lap u)^2-(didj u)^2] = (Om a/beta2) delta
     xi_t   = sqrt(6d0*(1d0-omega_m))
     E2_t   = 0.5d0*(omega_m/aexp**3 &
          & + sqrt((omega_m/aexp**3)**2 + 4d0*(1d0-omega_m)))
     hd_t   = -1.5d0*(omega_m/aexp**3)/(2d0*E2_t - omega_m/aexp**3)
     beta1_t = (xi_t/3d0)*(2d0*hd_t - 1d0 + (1d0-omega_m)/E2_t**2)
     beta_G  = 2d0*E2_t*beta1_t/xi_t**2      ! beta2: source coupling
     coeff_G = 1d0/(3d0*beta1_t*aexp**4)     ! Vainshtein coefficient
     ! Skip the solve while the linear coupling is negligible: the
     ! unscreened G_eff/G-1 = -xi/(9 beta2 E^2) decays as 1/E^4 into
     ! the past (< 1e-3 for z > 2.5). This avoids both wasted work
     ! and the pathological |coeff| ~ a^-4 regime at high redshift.
     if(-xi_t/(9d0*beta_G*E2_t) < 1d-3) return
  else
     ! LEGACY simplified nDGP-template coefficients (experimental)
     beta_G = galileon_beta(aexp)
     coeff_G = c3_galileon / (c2_galileon * aexp**4)
     xi_t = 0d0; E2_t = 1d0
  end if

  if(myid==1 .and. nstep==0 .and. ilevel==levelmin) then
     write(*,'(A,F8.4,A,ES10.3,A,F8.4)') &
          ' Galileon: beta=', beta_G, ' coeff=', coeff_G, ' at a=', aexp
  end if

  if(nstep==0) then
     call galileon_init_scalar(ilevel)
  end if

  converged = .false.
#ifdef USE_FFTW
  ! Operator-split spectral solve on the uniform domain level
  if(level_fft_ok(ilevel)) then
     call vain_prepare_uniform_cache(ilevel)
     call make_virtual_fine_dp(scalar_gr(1), ilevel)
     ! The tracker becomes only weakly elliptic near a~0.8 and needs
     ! damped trace/traceless Picard steps.  Test the mean-projected
     ! residual staged for the periodic FFT: after the negative-
     ! discriminant prescription the unprojected local targets need
     ! not have zero mean and therefore cannot themselves be a
     ! periodic Laplacian.
     do iter=1,500
        call vain_build_fft_rhs(ilevel, omega_m*aexp/beta_G, coeff_G, &
             & .true., fft_rel)
        if(fft_rel < galileon_eps) then
           converged=.true.
           if(myid==1) write(*,'(A,I2,A,I5,A,ES10.3)') &
                & ' Galileon level ',ilevel,' FFT converged in ',iter-1, &
                & ' iters, res=',fft_rel
           exit
        end if
        call level_fft_helmholtz(ilevel, 0d0, 0d0, 0.25d0)
        call make_virtual_fine_dp(scalar_gr(1), ilevel)
     end do
  end if
#endif

  gscal_ok=.false.
  trk=0
  if(.not.converged) then
#ifdef HYDRO_CUDA
  call scalar_gpu_begin(ilevel, .true., gscal_ok)
  if(gscal_ok) then
     gs_dx2=(0.5d0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1))**2
     sparams=0d0
     sparams(1)=1d0/gs_dx2
     sparams(2)=coeff_G
     sparams(3)=omega_m*aexp/beta_G
     sparams(4)=rho_tot
     sparams(5)=1d-2*omega_m*aexp/abs(beta_G)
     sparams(6)=gs_dx2
     if(galileon_tracker) trk=1
  end if
#endif
  do iter=1,n_iter_galileon

#ifdef HYDRO_CUDA
     if(gscal_ok) then
        call scalar_gpu_sweep_halo(ilevel, SCAL_MODEL_GALILEON, sparams, trk, &
             & res_max_local, src_max_local)
     else
#endif
     call galileon_gauss_seidel(ilevel, beta_G, coeff_G, res_max_local, src_max_local)

     call make_virtual_fine_dp(scalar_gr(1), ilevel)
#ifdef HYDRO_CUDA
     end if
#endif

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < galileon_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I5,A,ES10.3)') &
             ' Galileon level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do
#ifdef HYDRO_CUDA
  if(gscal_ok) call scalar_gpu_end(ilevel)
#endif
  end if

  if(.not. converged) then
     if(myid==1) write(*,'(A,I2,A,I5,A,ES10.3)') &
          & ' WARNING: Galileon level ',ilevel,' NOT converged after ', &
          & n_iter_galileon,' iters, res=',rel_res
     if(scalar_solver_strict) call scalar_solver_abort
  end if

  call galileon_save_old(ilevel)

  if(galileon_tracker) then
     ! Poisson back-reaction (Barreira eq. 11): the extra term
     ! -(kappa c3/M^3) phidot^2 lap(phi) integrates to a potential
     ! -(c3 xi^2/E^2) u, so the WHOLE fifth force is a gradient of u:
     !   F5 = +(xi/(6 E^2)) grad(u)   [c3 xi^2 = xi/6]
     ! Unscreened linear limit: F5/FN = -xi/(9 beta2 E^2)
     ! (= +0.84 at a=1 for Om=0.3, decaying as 1/E^4 into the past).
     call compute_fifth_force(ilevel, xi_t/(6d0*E2_t))
  else
     ! LEGACY: F5 = -(1/2) grad(u)
     call compute_fifth_force(ilevel, -0.5d0)
  end if

end subroutine galileon_solve_level

!=========================================================
! galileon_init_scalar / galileon_save_old
!=========================================================
subroutine galileon_init_scalar(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = 0d0
           scalar_gr_old(icell) = 0d0
        end do
     end do
  end do
end subroutine galileon_init_scalar

subroutine galileon_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do
end subroutine galileon_save_old

!=========================================================
! galileon_gauss_seidel: one Newton-GS sweep
! Same Vainshtein operator as nDGP, different coeff/source
!
! ∇²φ + coeff_G*[(∇²φ)² - (∇ᵢ∇ⱼφ)²] = source
! source = Ωm*a/β_G * δ
!=========================================================
subroutine galileon_gauss_seidel(ilevel, beta_G, coeff_G, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::beta_G, coeff_G
  real(dp),intent(out)::res_max, src_max

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,source,residual,jacobian,delta_u,sclamp
  real(dp)::phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp
  real(dp)::phi_xx,phi_yy,phi_zz,mix_xy2,mix_xz2,mix_yz2
  real(dp)::phi_pp,phi_pm,phi_mp,phi_mm
  real(dp)::lapl2,trace_ij2,vain_term
  real(dp)::dpp,dpm,dmp,dmm,dmix_du
  real(dp)::tbar_ij,qcoeff,disc,a_tgt
  integer::icolor

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  ! Mixed Hessian terms couple checkerboard-equal edge diagonals.
  ! Eight cell-parity colours remove that same-sweep dependence.
  do icolor=0,7

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  u_c,lapl,source,residual,jacobian,delta_u, &
!$omp&  phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp, &
!$omp&  phi_xx,phi_yy,phi_zz,mix_xy2,mix_xz2,mix_yz2, &
!$omp&  phi_pp,phi_pm,phi_mp,phi_mm,dpp,dpm,dmp,dmm,dmix_du, &
!$omp&  lapl2,trace_ij2,vain_term,tbar_ij,qcoeff,disc,a_tgt,sclamp) &
!$omp& reduction(max:res_max,src_max) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
              igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
           end do
        end do

        do ind=1,twotondim
           if(ind-1 /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              u_c = scalar_gr(ind_cell_w(i))

              ! Parent-CIC Dirichlet closure at coarse-fine interfaces.
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0, 0,u_c,phi_xm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0, 0,u_c,phi_xp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1, 0,u_c,phi_ym)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1, 0,u_c,phi_yp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 0,-1,u_c,phi_zm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 0, 1,u_c,phi_zp)

              ! Laplacian
              lapl = (phi_xp + phi_xm + phi_yp + phi_ym + phi_zp + phi_zm - 6d0*u_c) * dx2_inv

              ! Diagonal second derivatives
              phi_xx = (phi_xp + phi_xm - 2d0*u_c) * dx2_inv
              phi_yy = (phi_yp + phi_ym - 2d0*u_c) * dx2_inv
              phi_zz = (phi_zp + phi_zm - 2d0*u_c) * dx2_inv

              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 1, 0,u_c,phi_pp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1,-1, 0,u_c,phi_pm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 1, 0,u_c,phi_mp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1,-1, 0,u_c,phi_mm)
              dpp=(phi_pp-phi_xp-phi_yp+u_c)*dx2_inv
              dpm=(phi_xp-phi_pm-u_c+phi_ym)*dx2_inv
              dmp=(phi_yp-phi_mp-u_c+phi_xm)*dx2_inv
              dmm=(u_c-phi_ym-phi_xm+phi_mm)*dx2_inv
              mix_xy2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
              dmix_du=0.5d0*dx2_inv*(dpp-dpm-dmp+dmm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0, 1,u_c,phi_pp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 1, 0,-1,u_c,phi_pm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0, 1,u_c,phi_mp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel,-1, 0,-1,u_c,phi_mm)
              dpp=(phi_pp-phi_xp-phi_zp+u_c)*dx2_inv
              dpm=(phi_xp-phi_pm-u_c+phi_zm)*dx2_inv
              dmp=(phi_zp-phi_mp-u_c+phi_xm)*dx2_inv
              dmm=(u_c-phi_zm-phi_xm+phi_mm)*dx2_inv
              mix_xz2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
              dmix_du=dmix_du+0.5d0*dx2_inv*(dpp-dpm-dmp+dmm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1, 1,u_c,phi_pp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0, 1,-1,u_c,phi_pm)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1, 1,u_c,phi_mp)
              call scalar_sample_offset(ind_grid_w(i),ind,ilevel, 0,-1,-1,u_c,phi_mm)
              dpp=(phi_pp-phi_yp-phi_zp+u_c)*dx2_inv
              dpm=(phi_yp-phi_pm-u_c+phi_zm)*dx2_inv
              dmp=(phi_zp-phi_mp-u_c+phi_ym)*dx2_inv
              dmm=(u_c-phi_zm-phi_ym+phi_mm)*dx2_inv
              mix_yz2=0.25d0*(dpp**2+dpm**2+dmp**2+dmm**2)
              dmix_du=dmix_du+0.5d0*dx2_inv*(dpp-dpm-dmp+dmm)

              ! Vainshtein: (trace H)^2 - trace(H^2)
              lapl2 = lapl**2
              trace_ij2 = phi_xx**2 + phi_yy**2 + phi_zz**2 &
                   & + 2d0*(mix_xy2+mix_xz2+mix_yz2)
              vain_term = lapl2 - trace_ij2

              ! Source = Ωm*a/β_G * δ
              ! rho has box mean rho_tot ≈ 1; subtract it (δ = rho - rho_tot)
              source = omega_m * aexp / beta_G * (rho(ind_cell_w(i)) - rho_tot)

              ! Cubic-Galileon QSA loses real roots in sufficiently
              ! underdense cells at late times (Barreira et al. 2013).
              ! Match their ECOSMOG prescription: when the local
              ! trace/traceless discriminant is negative, set its
              ! square root to zero and relax to the double root.
              ! The FFT preconditioner applies the identical projection.
              tbar_ij=max(trace_ij2-lapl2/3d0,0d0)
              qcoeff=2d0*coeff_G/3d0
              disc=1d0+4d0*qcoeff*(source+coeff_G*tbar_ij)
              a_tgt=source
              if(abs(qcoeff)>1d-30) then
                 if(disc>0d0) then
                    a_tgt=(-1d0+sqrt(disc))/(2d0*qcoeff)
                 else
                    a_tgt=-1d0/(2d0*qcoeff)
                 end if
              end if
              if(galileon_tracker .and. abs(qcoeff)>1d-30) then
                 ! With epsilon=1/3 and centered mixed derivatives,
                 ! the trace-free quadratic target is independent of
                 ! the centre cell.  Relax the analytic physical root;
                 ! a negative discriminant is projected to its double
                 ! root exactly as in the ECOSMOG prescription.
                 residual=lapl-a_tgt
                 jacobian=-6d0*dx2_inv
              else
                 residual = lapl + coeff_G * vain_term - source
                 ! Exact local Jacobian of the compatible Vainshtein stencil.
                 jacobian = -6d0*dx2_inv + coeff_G * &
                      & (-8d0*lapl*dx2_inv-2d0*dmix_du)
              end if

              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 ! Damped Newton; floor the clamp scale at 1% of the
                 ! delta=1 source so cells with delta~0 still relax
                 if(galileon_tracker) then
                    sclamp=0.5d0*dx2*max(abs(a_tgt),abs(source), &
                         & 1d-2*omega_m*aexp/abs(beta_G))
                 else
                    sclamp=0.5d0*dx2*max(abs(source), &
                         & 1d-2*omega_m*aexp/abs(beta_G))
                 end if
                 if(abs(delta_u) > sclamp) delta_u = sign(sclamp, delta_u)
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
              end if

              res_max = max(res_max, abs(residual))
              src_max = max(src_max, abs(source))

           end do
        end do
     end do

  end do

end subroutine galileon_gauss_seidel

!#########################################################
!#########################################################
!  COUPLED DARK ENERGY — FORCE MODIFIER
!#########################################################
!#########################################################

!=========================================================
! apply_coupled_de_force: multiply f() by (1 + 2β²)
! G_eff = G*(1+2β²) for DM component (approximation: apply to all)
!=========================================================
subroutine apply_coupled_de_force(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  real(dp)::enhancement
  integer,dimension(1:nvector)::ind_cell

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  enhancement = 1d0 + 2d0 * beta_cde**2

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_cell) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+active(ilevel)%igrid(igrid+i-1)
        end do

        do i=1,ngrid
           do idim=1,ndim
              f(ind_cell(i),idim) = f(ind_cell(i),idim) * enhancement
           end do
        end do
     end do
  end do

end subroutine apply_coupled_de_force

#ifdef HYDRO_CUDA
!#########################################################
!#########################################################
! GPU scalar-solver support (scalar_cuda_kernels.cu).
! The Newton-GS sweeps run on the GPU; boundary Dirichlet
! data and halo cell lists are prepared here once per
! level solve (the parent level is frozen during a solve).
!#########################################################
!#########################################################

!=========================================================
! scalar_lookup_icell: same-level cell lookup returning the
! cell index (index-returning variant of scalar_lookup_cell)
!=========================================================
subroutine scalar_lookup_icell(ilevel, ix_in, iy_in, iz_in, icell_out, found)
  use amr_commons
  use morton_hash
  use morton_keys, only: mkey_t, morton_encode
  implicit none
  integer,intent(in)::ilevel
  integer(8),intent(in)::ix_in,iy_in,iz_in
  integer,intent(out)::icell_out
  logical,intent(out)::found
  integer(8)::ix,iy,iz,ncx,ncy,ncz,gx,gy,gz
  integer::igrid,ind
  type(mkey_t)::key

  found=.false.
  icell_out=0
  if(.not. allocated(mort_table)) return
  if(ilevel < 1 .or. ilevel > size(mort_table)) return

  ncx=int(nx,8)*2_8**ilevel
  ncy=int(ny,8)*2_8**ilevel
  ncz=int(nz,8)*2_8**ilevel
  ix=modulo(ix_in,ncx)
  iy=modulo(iy_in,ncy)
  iz=modulo(iz_in,ncz)
  key=morton_encode(ix/2_8,iy/2_8,iz/2_8)
  igrid=morton_hash_lookup(mort_table(ilevel),key)
  if(igrid <= 0) return

  ind=1+int(modulo(ix,2_8))+2*int(modulo(iy,2_8))+4*int(modulo(iz,2_8))
  icell_out=ncoarse+(ind-1)*ngridmax+igrid
  found=.true.
end subroutine scalar_lookup_icell

!=========================================================
! build_scalar_halo_indices: flat emission/reception cell
! lists for the GPU scalar halo (same enumeration as
! make_virtual_fine_dp packing; see build_mg_halo_indices)
!=========================================================
subroutine build_scalar_halo_indices(ilevel)
  use amr_commons
  use scalar_gpu_commons
  use scalar_cuda_interface
  use iso_c_binding
  implicit none
  integer,intent(in)::ilevel
  integer::icpu,i,j,idx,iskip

  sgpu_n_emit = 0
  sgpu_n_recv = 0
  do icpu = 1, ncpu
     sgpu_n_emit = sgpu_n_emit + emission(icpu,ilevel)%ngrid * twotondim
     sgpu_n_recv = sgpu_n_recv + reception(icpu,ilevel)%ngrid * twotondim
  end do

  if(allocated(sgpu_emit_cells)) deallocate(sgpu_emit_cells)
  if(allocated(sgpu_recv_cells)) deallocate(sgpu_recv_cells)
  if(allocated(sgpu_emit_buf))   deallocate(sgpu_emit_buf)
  if(allocated(sgpu_recv_buf))   deallocate(sgpu_recv_buf)
  allocate(sgpu_emit_cells(1:max(sgpu_n_emit,1)))
  allocate(sgpu_recv_cells(1:max(sgpu_n_recv,1)))
  allocate(sgpu_emit_buf(1:max(sgpu_n_emit,1)))
  allocate(sgpu_recv_buf(1:max(sgpu_n_recv,1)))

  idx = 0
  do icpu = 1, ncpu
     if(emission(icpu,ilevel)%ngrid > 0) then
        do j = 1, twotondim
           iskip = ncoarse + (j-1)*ngridmax
           do i = 1, emission(icpu,ilevel)%ngrid
              idx = idx + 1
              sgpu_emit_cells(idx) = emission(icpu,ilevel)%igrid(i) + iskip
           end do
        end do
     end if
  end do

  idx = 0
  do icpu = 1, ncpu
     if(reception(icpu,ilevel)%ngrid > 0) then
        do j = 1, twotondim
           iskip = ncoarse + (j-1)*ngridmax
           do i = 1, reception(icpu,ilevel)%ngrid
              idx = idx + 1
              sgpu_recv_cells(idx) = reception(icpu,ilevel)%igrid(i) + iskip
           end do
        end do
     end if
  end do

  call cuda_scal_halo_setup_c( &
       sgpu_emit_cells, int(sgpu_n_emit, c_int), &
       sgpu_recv_cells, int(sgpu_n_recv, c_int))

end subroutine build_scalar_halo_indices

!=========================================================
! make_virtual_scalar_gpu: scalar_gr halo exchange with the
! field resident on the GPU (gather emission cells to host,
! standard MPI exchange, scatter reception cells back)
!=========================================================
subroutine make_virtual_scalar_gpu(ilevel)
  use amr_commons
  use poisson_commons
  use scalar_gpu_commons
  use scalar_cuda_interface
  use iso_c_binding
  implicit none
  integer,intent(in)::ilevel
  integer::i

  if(sgpu_n_emit == 0 .and. sgpu_n_recv == 0) return

  if(sgpu_n_emit > 0) then
     call cuda_scal_halo_gather_c(sgpu_emit_buf, int(sgpu_n_emit, c_int))
     do i = 1, sgpu_n_emit
        scalar_gr(sgpu_emit_cells(i)) = sgpu_emit_buf(i)
     end do
  end if

  call make_virtual_fine_dp(scalar_gr(1), ilevel)

  if(sgpu_n_recv > 0) then
     do i = 1, sgpu_n_recv
        sgpu_recv_buf(i) = scalar_gr(sgpu_recv_cells(i))
     end do
     call cuda_scal_halo_scatter_c(sgpu_recv_buf, int(sgpu_n_recv, c_int))
  end if

end subroutine make_virtual_scalar_gpu

!=========================================================
! scalar_gpu_begin: decide (collectively) to run this level
! solve on the GPU, prepare grid tables, coarse-fine
! Dirichlet blocks and halo lists, and upload the field.
! MPI-collective — every rank must call at the same point.
!=========================================================
subroutine scalar_gpu_begin(ilevel, need18, ok)
  use amr_commons
  use poisson_commons
  use scalar_gpu_commons
  use scalar_cuda_interface
  use cuda_commons
  use morton_keys, only: mkey_t, grid_to_morton, morton_decode
  use iso_c_binding
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  logical,intent(in)::need18
  logical,intent(out)::ok

  integer::ncache,ia,igrid_amr,jf,je,ind,io,slot,info,flag,flag_all
  integer::noff,k,nbnd,base
  integer::bx,by,bz,ox,oy,oz,tx,ty,tz,dxs,dys,dzs,nnz,gtab,idx2
  integer::icell_live
  integer(8)::gx,gy,gz,cx,cy,cz,txa,tya,tza
  logical::found,has_missing
  real(dp)::val
  integer(c_long_long)::ncell_c
  type(mkey_t)::key

  ok=.false.
  scal_gpu_active=.false.

  flag=0
  if(gpu_scalar .and. cuda_pool_is_initialized_c()/=0) flag=1
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(flag,flag_all,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  flag=flag_all
#endif
  if(flag==0) return

  ncache=active(ilevel)%ngrid
  noff=6
  if(need18) noff=18

  if(ncache > 0) then
     call vain_prepare_uniform_cache(ilevel)

     if(allocated(sgpu_face)) then
        if(size(sgpu_face) < 6*ncache) then
           deallocate(sgpu_face,sgpu_edge,sgpu_bnd_slot)
        end if
     end if
     if(.not.allocated(sgpu_face)) then
        allocate(sgpu_face(1:6*max(ncache,1)))
        allocate(sgpu_edge(1:12*max(ncache,1)))
        allocate(sgpu_bnd_slot(1:max(ncache,1)))
     end if

     ! Pack per-grid neighbor tables and count boundary grids
     nbnd=0
     do ia=1,ncache
        has_missing=.false.
        do jf=1,6
           sgpu_face((ia-1)*6+jf)=vain_face_grid(ia,jf)
           if(vain_face_grid(ia,jf)<=0) has_missing=.true.
        end do
        do je=1,4
           sgpu_edge((ia-1)*12+je)  =vain_xy_grid(ia,je)
           sgpu_edge((ia-1)*12+4+je)=vain_xz_grid(ia,je)
           sgpu_edge((ia-1)*12+8+je)=vain_yz_grid(ia,je)
        end do
        if(need18) then
           do je=1,4
              if(vain_xy_grid(ia,je)<=0) has_missing=.true.
              if(vain_xz_grid(ia,je)<=0) has_missing=.true.
              if(vain_yz_grid(ia,je)<=0) has_missing=.true.
           end do
        end if
        if(has_missing) then
           sgpu_bnd_slot(ia)=nbnd
           nbnd=nbnd+1
        else
           sgpu_bnd_slot(ia)=-1
        end if
     end do
     sgpu_nbnd=nbnd
     sgpu_noff=noff

     if(allocated(sgpu_bnd_live)) deallocate(sgpu_bnd_live)
     if(allocated(sgpu_bnd_val))  deallocate(sgpu_bnd_val)
     allocate(sgpu_bnd_live(1:max(nbnd*8*noff,1)))
     allocate(sgpu_bnd_val (1:max(nbnd*8*noff,1)))
     sgpu_bnd_live=0
     sgpu_bnd_val=scal_gpu_sentinel

     ! Fill boundary closures: live same-level cell index where the cell
     ! exists but the grid tables cannot reach it, frozen parent-CIC
     ! Dirichlet value otherwise (sentinel = zero-gradient fallback).
!$omp parallel do private(ia,slot,igrid_amr,key,gx,gy,gz,ind,bx,by,bz, &
!$omp& cx,cy,cz,io,ox,oy,oz,tx,ty,tz,dxs,dys,dzs,nnz,gtab,idx2,base,k, &
!$omp& txa,tya,tza,icell_live,found,val) schedule(dynamic)
     do ia=1,ncache
        slot=sgpu_bnd_slot(ia)
        if(slot<0) cycle
        igrid_amr=active(ilevel)%igrid(ia)
        key=grid_to_morton(igrid_amr,ilevel)
        call morton_decode(key,gx,gy,gz)
        base=slot*8*noff
        do ind=1,twotondim
           bx=mod(ind-1,2); by=mod((ind-1)/2,2); bz=(ind-1)/4
           cx=2_8*gx+int(bx,8); cy=2_8*gy+int(by,8); cz=2_8*gz+int(bz,8)
           do io=1,noff
              ox=sgpu_off(1,io); oy=sgpu_off(2,io); oz=sgpu_off(3,io)
              tx=bx+ox; ty=by+oy; tz=bz+oz
              dxs=0; if(tx<0) dxs=-1; if(tx>1) dxs=1
              dys=0; if(ty<0) dys=-1; if(ty>1) dys=1
              dzs=0; if(tz<0) dzs=-1; if(tz>1) dzs=1
              nnz=abs(dxs)+abs(dys)+abs(dzs)
              if(nnz==0) then
                 gtab=igrid_amr
              else if(nnz==1) then
                 if(dxs/=0) then
                    idx2=merge(1,2,dxs<0)
                 else if(dys/=0) then
                    idx2=merge(3,4,dys<0)
                 else
                    idx2=merge(5,6,dzs<0)
                 end if
                 gtab=sgpu_face((ia-1)*6+idx2)
              else
                 if(dzs==0) then
                    idx2=((dxs+1)/2)*2+((dys+1)/2)+1
                 else if(dys==0) then
                    idx2=4+((dxs+1)/2)*2+((dzs+1)/2)+1
                 else
                    idx2=8+((dys+1)/2)*2+((dzs+1)/2)+1
                 end if
                 gtab=sgpu_edge((ia-1)*12+idx2)
              end if
              if(gtab>0) cycle

              k=base+(ind-1)*noff+io
              txa=cx+int(ox,8); tya=cy+int(oy,8); tza=cz+int(oz,8)
              call scalar_lookup_icell(ilevel,txa,tya,tza,icell_live,found)
              if(found) then
                 sgpu_bnd_live(k)=icell_live
              else
                 call scalar_sample_offset(igrid_amr,ind,ilevel,ox,oy,oz, &
                      & scal_gpu_sentinel,val)
                 sgpu_bnd_val(k)=val
              end if
           end do
        end do
     end do
!$omp end parallel do

     call build_scalar_halo_indices(ilevel)

     ncell_c=int(ncoarse,c_long_long) &
          & +int(twotondim,c_long_long)*int(ngridmax,c_long_long)
     call cuda_scal_upload_c(scalar_gr, rho, ncell_c, &
          & active(ilevel)%igrid, sgpu_face, sgpu_edge, sgpu_bnd_slot, &
          & sgpu_bnd_live, sgpu_bnd_val, &
          & int(ncache,c_int), int(nbnd,c_int), int(noff,c_int))
  else
     sgpu_n_emit=0
     sgpu_n_recv=0
  end if

  ! Collective agreement on upload success (ranks without grids
  ! participate trivially in the halo/ALLREDUCE steps)
  flag=0
  if(ncache==0 .or. cuda_scal_is_ready_c()/=0) flag=1
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(flag,flag_all,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  flag=flag_all
#endif
  if(flag==0) then
     call cuda_scal_release_c()
     return
  end if

  ok=.true.
  scal_gpu_active=.true.
  scal_gpu_level=ilevel

end subroutine scalar_gpu_begin

!=========================================================
! scalar_gpu_sweep_halo: one red+black Newton-GS sweep on
! the GPU followed by the scalar halo exchange (mirrors the
! CPU sweep + make_virtual_fine_dp pair in the solve loops)
!=========================================================
subroutine scalar_gpu_sweep_halo(ilevel, model, params, tracker, &
     & res_max, src_max)
  use amr_commons
  use scalar_cuda_interface
  use iso_c_binding
  implicit none
  integer,intent(in)::ilevel,model,tracker
  real(dp),dimension(1:12),intent(in)::params
  real(dp),intent(out)::res_max,src_max
  real(c_double)::res_c,src_c

  res_max=0d0
  src_max=0d0
  if(active(ilevel)%ngrid > 0) then
     call cuda_scal_sweep_c(int(model,c_int), params, &
          & int(ngridmax,c_int), int(ncoarse,c_int), &
          & int(tracker,c_int), res_c, src_c)
     res_max=res_c
     src_max=src_c
  end if

  call make_virtual_scalar_gpu(ilevel)

end subroutine scalar_gpu_sweep_halo

!=========================================================
! scalar_gpu_end: download the converged field and release
! the per-solve GPU arrays
!=========================================================
subroutine scalar_gpu_end(ilevel)
  use amr_commons
  use poisson_commons
  use scalar_gpu_commons
  use scalar_cuda_interface
  use iso_c_binding
  implicit none
  integer,intent(in)::ilevel
  integer(c_long_long)::ncell_c

  if(scal_gpu_active .and. active(ilevel)%ngrid > 0) then
     ncell_c=int(ncoarse,c_long_long) &
          & +int(twotondim,c_long_long)*int(ngridmax,c_long_long)
     call cuda_scal_download_c(scalar_gr, ncell_c)
  end if
  call cuda_scal_release_c()
  scal_gpu_active=.false.
  scal_gpu_level=-1

end subroutine scalar_gpu_end
#endif
