program ramses_density_slab
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer, parameter :: nbase_x=512, nbase_y=512
  integer, parameter :: nzoom_x=708, nzoom_y=1000
  real(dp), parameter :: default_lbox_mpc=128.0_dp
  real(dp), parameter :: xmin_zoom_mpc=55.0_dp
  real(dp), parameter :: xmax_zoom_mpc=72.7_dp
  real(dp), parameter :: ymin_zoom_mpc=52.0_dp
  real(dp), parameter :: ymax_zoom_mpc=77.0_dp
  real(dp), parameter :: zmin_slab_mpc=55.75_dp
  real(dp), parameter :: zmax_slab_mpc=72.75_dp

  integer :: ierr, ilevel, nchar_env, env_status
  integer(kind=8) :: nleaf_local, nleaf_global
  integer(kind=8), allocatable :: nleaf_level_local(:), nleaf_level_global(:)
  real(dp) :: lbox_mpc
  real(dp), allocatable :: base_local(:,:), base_global(:,:)
  real(dp), allocatable :: zoom_local(:,:), zoom_global(:,:)
  character(len=512) :: output_prefix, env_value

  call read_params

  lbox_mpc=default_lbox_mpc
  env_value=' '
  call get_environment_variable('RAMSES_LBOX_MPC_H',env_value, &
       length=nchar_env,status=env_status)
  if(env_status==0 .and. nchar_env>0)read(env_value(1:nchar_env),*)lbox_mpc

  output_prefix='density_slab'
  env_value=' '
  call get_environment_variable('RAMSES_SLAB_PREFIX',env_value, &
       length=nchar_env,status=env_status)
  if(env_status==0 .and. nchar_env>0)output_prefix=env_value(1:nchar_env)

  if(ndim/=3)then
     if(myid==1)write(*,*)'ramses_density_slab requires NDIM=3'
     call clean_stop
  endif
  if(.not.pic .or. .not.poisson)then
     if(myid==1)write(*,*)'ramses_density_slab requires pic and poisson'
     call clean_stop
  endif
  if(nrestart<=0)then
     if(myid==1)write(*,*)'ramses_density_slab requires a restart snapshot'
     call clean_stop
  endif

  call init_amr
  call init_time
  if(poisson)call init_poisson
  if(pic)call init_part
  if(pic)call init_tree

  allocate(base_local(nbase_x,nbase_y),zoom_local(nzoom_x,nzoom_y))
  allocate(base_global(nbase_x,nbase_y),zoom_global(nzoom_x,nzoom_y))
  allocate(nleaf_level_local(levelmin:nlevelmax))
  allocate(nleaf_level_global(levelmin:nlevelmax))
  base_local=0.0_dp
  zoom_local=0.0_dp
  base_global=0.0_dp
  zoom_global=0.0_dp
  nleaf_level_local=0_8
  nleaf_level_global=0_8

  if(myid==1)then
     write(*,'(A,I0,A,I0,A,I0)') &
          'Density reconstruction: levelmin=',levelmin, &
          ', levelmax=',nlevelmax,', MPI ranks=',ncpu
     write(*,'(A,2F10.3)')'LOS slab [Mpc/h]: ',zmin_slab_mpc,zmax_slab_mpc
     write(*,'(A,4F10.3)')'Zoom xy [Mpc/h]: ',xmin_zoom_mpc, &
          xmax_zoom_mpc,ymin_zoom_mpc,ymax_zoom_mpc
  endif

  do ilevel=levelmin,nlevelmax
     if(numbtot(1,ilevel)==0)cycle
     call make_tree_fine(ilevel)
     call rho_fine(ilevel,2)
     call accumulate_level(ilevel,lbox_mpc,base_local,zoom_local, &
          nleaf_level_local(ilevel))
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  enddo

#ifdef WITHOUTMPI
  base_global=base_local
  zoom_global=zoom_local
  nleaf_level_global=nleaf_level_local
#else
  call MPI_REDUCE(base_local,base_global,nbase_x*nbase_y, &
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(zoom_local,zoom_global,nzoom_x*nzoom_y, &
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(nleaf_level_local,nleaf_level_global, &
       nlevelmax-levelmin+1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif

  if(myid==1)then
     nleaf_global=sum(nleaf_level_global)
     call write_products(trim(output_prefix),lbox_mpc,base_global, &
          zoom_global,nleaf_level_global,nleaf_global)
  endif

#ifndef WITHOUTMPI
  call MPI_FINALIZE(ierr)
#endif

contains

  subroutine accumulate_level(ilevel,lbox,map_base,map_zoom,nleaf)
    use amr_commons
    use poisson_commons
    implicit none
    integer, intent(in) :: ilevel
    real(dp), intent(in) :: lbox
    real(dp), intent(inout) :: map_base(:,:),map_zoom(:,:)
    integer(kind=8), intent(out) :: nleaf
    integer :: i,ind,iskip,icell,ix,iy,iz,nx_loc
    real(dp) :: dx_code,dx_mpc,scale
    real(dp) :: xcen,ycen,zcen,xlo,xhi,ylo,yhi,zlo,zhi
    real(dp) :: xc(1:twotondim,1:ndim)
    real(dp) :: skip_loc(1:3)

    nx_loc=icoarse_max-icoarse_min+1
    scale=boxlen/dble(nx_loc)
    skip_loc=0.0_dp
    skip_loc(1)=dble(icoarse_min)
    skip_loc(2)=dble(jcoarse_min)
    skip_loc(3)=dble(kcoarse_min)
    dx_code=0.5_dp**ilevel
    dx_mpc=dx_code*scale*lbox/boxlen

    do ind=1,twotondim
       iz=(ind-1)/4
       iy=(ind-1-4*iz)/2
       ix=(ind-1-2*iy-4*iz)
       xc(ind,1)=(dble(ix)-0.5_dp)*dx_code
       xc(ind,2)=(dble(iy)-0.5_dp)*dx_code
       xc(ind,3)=(dble(iz)-0.5_dp)*dx_code
    enddo

    nleaf=0_8
    do ind=1,twotondim
       iskip=ncoarse+(ind-1)*ngridmax
       do i=1,active(ilevel)%ngrid
          icell=active(ilevel)%igrid(i)+iskip
          xcen=(xg(active(ilevel)%igrid(i),1)+xc(ind,1)-skip_loc(1)) &
               *scale*lbox/boxlen
          ycen=(xg(active(ilevel)%igrid(i),2)+xc(ind,2)-skip_loc(2)) &
               *scale*lbox/boxlen
          zcen=(xg(active(ilevel)%igrid(i),3)+xc(ind,3)-skip_loc(3)) &
               *scale*lbox/boxlen
          xlo=xcen-0.5_dp*dx_mpc
          xhi=xcen+0.5_dp*dx_mpc
          ylo=ycen-0.5_dp*dx_mpc
          yhi=ycen+0.5_dp*dx_mpc
          zlo=zcen-0.5_dp*dx_mpc
          zhi=zcen+0.5_dp*dx_mpc

          if(ilevel==levelmin)then
             call deposit_cell(map_base,0.0_dp,lbox,0.0_dp,lbox, &
                  zmin_slab_mpc,zmax_slab_mpc,xlo,xhi,ylo,yhi,zlo,zhi, &
                  rho(icell))
          endif

          if(son(icell)==0 .and. xhi>xmin_zoom_mpc .and. &
               xlo<xmax_zoom_mpc .and. yhi>ymin_zoom_mpc .and. &
               ylo<ymax_zoom_mpc .and. zhi>zmin_slab_mpc .and. &
               zlo<zmax_slab_mpc)then
             nleaf=nleaf+1_8
             call deposit_cell(map_zoom,xmin_zoom_mpc,xmax_zoom_mpc, &
                  ymin_zoom_mpc,ymax_zoom_mpc,zmin_slab_mpc,zmax_slab_mpc, &
                  xlo,xhi,ylo,yhi,zlo,zhi,rho(icell))
          endif
       enddo
    enddo
  end subroutine accumulate_level

  subroutine deposit_cell(map,xmin,xmax,ymin,ymax,zmin,zmax, &
       xlo,xhi,ylo,yhi,zlo,zhi,density)
    implicit none
    real(dp), intent(inout) :: map(:,:)
    real(dp), intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(dp), intent(in) :: xlo,xhi,ylo,yhi,zlo,zhi,density
    integer :: i,j,i0,i1,j0,j1,nxm,nym
    real(dp) :: dxp,dyp,dzo,ox,oy,px0,px1,py0,py1

    dzo=max(0.0_dp,min(zhi,zmax)-max(zlo,zmin))
    if(dzo<=0.0_dp)return
    if(xhi<=xmin .or. xlo>=xmax .or. yhi<=ymin .or. ylo>=ymax)return

    nxm=size(map,1)
    nym=size(map,2)
    dxp=(xmax-xmin)/dble(nxm)
    dyp=(ymax-ymin)/dble(nym)
    i0=max(1,min(nxm,int(floor((max(xlo,xmin)-xmin)/dxp))+1))
    i1=max(1,min(nxm,int(ceiling((min(xhi,xmax)-xmin)/dxp))))
    j0=max(1,min(nym,int(floor((max(ylo,ymin)-ymin)/dyp))+1))
    j1=max(1,min(nym,int(ceiling((min(yhi,ymax)-ymin)/dyp))))

    do j=j0,j1
       py0=ymin+dble(j-1)*dyp
       py1=py0+dyp
       oy=max(0.0_dp,min(yhi,py1)-max(ylo,py0))
       if(oy<=0.0_dp)cycle
       do i=i0,i1
          px0=xmin+dble(i-1)*dxp
          px1=px0+dxp
          ox=max(0.0_dp,min(xhi,px1)-max(xlo,px0))
          if(ox<=0.0_dp)cycle
          map(i,j)=map(i,j)+density*dzo*ox*oy/(dxp*dyp)
       enddo
    enddo
  end subroutine deposit_cell

  subroutine write_products(prefix,lbox,map_base,map_zoom,nleaf_level,nleaf)
    implicit none
    character(len=*), intent(in) :: prefix
    real(dp), intent(in) :: lbox
    real(dp), intent(in) :: map_base(:,:),map_zoom(:,:)
    integer(kind=8), intent(in) :: nleaf_level(:),nleaf
    integer :: i,iu
    real(dp) :: mean_column,base_integral,zoom_integral

    mean_column=rho_tot*(zmax_slab_mpc-zmin_slab_mpc)
    base_integral=sum(map_base)*(lbox/dble(nbase_x)) &
         *(lbox/dble(nbase_y))
    zoom_integral=sum(map_zoom) &
         *((xmax_zoom_mpc-xmin_zoom_mpc)/dble(nzoom_x)) &
         *((ymax_zoom_mpc-ymin_zoom_mpc)/dble(nzoom_y))

    open(newunit=iu,file=trim(prefix)//'_base.bin',access='stream', &
         form='unformatted',status='replace')
    write(iu)map_base
    close(iu)
    open(newunit=iu,file=trim(prefix)//'_leaf.bin',access='stream', &
         form='unformatted',status='replace')
    write(iu)map_zoom
    close(iu)

    open(newunit=iu,file=trim(prefix)//'_meta.txt',form='formatted', &
         status='replace')
    write(iu,'(A,ES24.16)')'aexp=',aexp
    write(iu,'(A,ES24.16)')'redshift=',1.0_dp/aexp-1.0_dp
    write(iu,'(A,I0)')'levelmin=',levelmin
    write(iu,'(A,I0)')'levelmax=',nlevelmax
    write(iu,'(A,ES24.16)')'lbox_mpc_h=',lbox
    write(iu,'(A,2(ES24.16,1X))')'slab_z_mpc_h=', &
         zmin_slab_mpc,zmax_slab_mpc
    write(iu,'(A,4(ES24.16,1X))')'zoom_xy_mpc_h=',xmin_zoom_mpc, &
         xmax_zoom_mpc,ymin_zoom_mpc,ymax_zoom_mpc
    write(iu,'(A,2(I0,1X))')'base_shape=',nbase_x,nbase_y
    write(iu,'(A,2(I0,1X))')'leaf_shape=',nzoom_x,nzoom_y
    write(iu,'(A,ES24.16)')'rho_mean_code=',rho_tot
    write(iu,'(A,ES24.16)')'mean_column_code=',mean_column
    write(iu,'(A,ES24.16)')'base_integral=',base_integral
    write(iu,'(A,ES24.16)')'leaf_integral=',zoom_integral
    write(iu,'(A,I0)')'leaf_cells_intersecting=',nleaf
    do i=levelmin,nlevelmax
       write(iu,'(A,I0,A,I0)')'leaf_cells_level_',i,'=', &
            nleaf_level(i-levelmin+1)
    enddo
    close(iu)
  end subroutine write_products

end program ramses_density_slab
