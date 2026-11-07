module poisson_commons
  use amr_commons
  use poisson_parameters
  use qumond_composite

  ! Experimental fixed-hierarchy history path; the shared solver state remains
  ! authoritative unless explicitly enabled. This does not repair AMR source loss.
  logical :: qh_checked=.false.,qh_enabled=.false.,qh_publish=.false.,qh_evolve=.false.
  integer :: qh_kick_unit=0
  integer :: qh_stage=0
  type(composite_mesh_t) :: qh_composite_mesh
  real(dp),allocatable :: qh_current(:,:),qh_previous(:,:),qh_rho_newton(:)
  real(dp),allocatable :: qh_time(:,:),qh_old_time(:,:)
  integer,allocatable :: qh_epoch(:,:)
  logical,allocatable :: qh_valid(:,:)
  integer,allocatable :: qh_cell_level(:),qh_grid_level(:),qh_son(:)
  real(dp),allocatable :: qh_xg(:,:)

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:)  ::rho_star          ! Star density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  ! Uniform-grid BIMOND fields. Retain the twin force until both particle
  ! half-kicks have consumed the previous gravity epoch.
  real(dp),allocatable :: bimond_rho_twin(:),bimond_rho_ordinary(:)
  real(dp),allocatable :: bimond_rho_phantom(:),bimond_phi_difference(:)
  real(dp),allocatable :: bimond_f_difference(:,:),bimond_f_twin(:,:)
  logical :: bimond_fields_ready=.false.
  real(dp) :: bimond_acceleration_max=0d0

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Refined-MG exchange policy.  The communication graph is fixed during one
  ! solve, so auto selection is performed once per level and reused by every
  ! V-cycle.  Values are the EXCHANGE_* identifiers from dynamic_exchange.
  integer, allocatable, dimension(:) :: mg_forward_backend
  integer, allocatable, dimension(:) :: mg_reverse_backend
  real(dp), allocatable, dimension(:) :: mg_exchange_send_dp
  real(dp), allocatable, dimension(:) :: mg_exchange_recv_dp
  integer, allocatable, dimension(:) :: mg_exchange_send_int
  integer, allocatable, dimension(:) :: mg_exchange_recv_int

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! A restart restores phi from the gravity dump.  Consume it once per level
  ! as the first Poisson initial guess; fresh runs retain the normal predictor.
  logical, allocatable, dimension(:) :: phi_restart_available

  ! A checkpoint may advertise its restored potential as a safe warm-start
  ! guess only when every active level has completed a Poisson solve since its
  ! most recent topology change.  Legacy checkpoints have no such marker.
  logical, allocatable, dimension(:) :: phi_checkpoint_level_valid

  ! Pre-computed neighbor grids for fine-level MG solver
  integer, allocatable, dimension(:,:) :: nbor_grid_fine

  ! Pre-computed neighbor grids/cpus for coarse-level MG solver
  type mg_nbor_cache_t
     integer, allocatable, dimension(:,:) :: grid  ! (0:2*ndim, 1:ngrid)
     integer, allocatable, dimension(:,:) :: cpu   ! (0:2*ndim, 1:ngrid)
  end type mg_nbor_cache_t
  type(mg_nbor_cache_t), allocatable, dimension(:) :: nbor_mg_cache

  ! GPU MG halo exchange arrays (emission/reception cell indices and buffers)
  integer, allocatable :: mg_halo_emit_cells(:)
  integer, allocatable :: mg_halo_recv_cells(:)
  real(dp), allocatable :: mg_halo_emit_buf(:)
  real(dp), allocatable :: mg_halo_recv_buf(:)
  integer :: mg_halo_n_emit = 0
  integer :: mg_halo_n_recv = 0

  ! GPU MG restrict/interp flat arrays
  integer, allocatable :: mg_ri_flat_offset(:)
  integer :: mg_ri_total_coarse = 0
  real(dp), allocatable :: mg_ri_coarse_rhs(:)
  real(dp), allocatable :: mg_ri_coarse_phi(:)

  ! Modified gravity scalar field (f_R or phi_nDGP)
  real(dp),allocatable,dimension(:)  ::scalar_gr       ! Scalar field
  real(dp),allocatable,dimension(:)  ::scalar_gr_old   ! Previous step (initial guess)

  ! Uniform Vainshtein-stencil topology.  nDGP/Galileon nonlinear
  ! corrections reuse an unchanged active-grid layout many times within
  ! one scalar solve, so build these Morton neighbours once per solve.
  integer,allocatable,dimension(:,:) ::vain_face_grid  ! (active grid, 6 faces)
  integer,allocatable,dimension(:,:) ::vain_xy_grid    ! (active grid, 4 xy edges)
  integer,allocatable,dimension(:,:) ::vain_xz_grid    ! (active grid, 4 xz edges)
  integer,allocatable,dimension(:,:) ::vain_yz_grid    ! (active grid, 4 yz edges)
  integer ::vain_cache_level=-1
  integer ::vain_cache_ngrid=-1
  logical ::vain_cache_ready=.false.

  ! Fuzzy Dark Matter (FDM) wavefunction
  real(dp),allocatable,dimension(:)  ::psi_re          ! Re(psi)
  real(dp),allocatable,dimension(:)  ::psi_im          ! Im(psi)

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

contains

  subroutine qumond_history_configure()
    character(len=8)::env
    integer::status
    if(qh_checked)return
    env=''
    call get_environment_variable('RAR_QUMOND_STAGE_HISTORY',env,status=status)
    qh_enabled=status==0.and.trim(env)=='1'
    env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_PUBLISH',env,status=status)
    qh_publish=status==0.and.trim(env)=='1'
    env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_EVOLVE',env,status=status)
    qh_evolve=status==0.and.trim(env)=='1'
    if((qh_publish.and..not.qh_enabled).or.(qh_evolve.and..not.qh_publish))then
       write(*,*) 'QUMOND composite publication/evolution requires stage history and publication'
       call clean_stop
    endif
    qh_checked=.true.
  end subroutine qumond_history_configure

  subroutine qumond_record_kick(phase,cell,lev,source,receiver,halfdt,force,oldp,olde,newp,newe)
    integer,intent(in)::phase,cell,lev
    real(dp),intent(in)::source,receiver,halfdt,force(3),oldp(3),olde,newp(3),newe
    if(qh_kick_unit==0)return
!$omp critical(rar_gravity_kick_record)
    write(qh_kick_unit,'(3(I0,1X),14(ES24.16,1X))') phase,cell,lev,source,receiver,halfdt, &
         force,oldp,olde,newp,newe
!$omp end critical(rar_gravity_kick_record)
  end subroutine qumond_record_kick

  subroutine qumond_history_enter(stage,ilevel)
#include "amr_index.h"
    integer,intent(in)::stage,ilevel
    character(len=8)::env
    integer::status,lev,g,grid,child,cell
    integer,allocatable::levels(:)
    call qumond_history_configure()
    if(.not.qh_enabled)return
    if(stage<1.or.stage>2)call clean_stop
#ifdef HYDRO_CUDA
    write(*,*) 'QUMOND stage history: GPU build unsupported'
    call clean_stop
#endif
    if(ncpu/=1.or.cosmo.or.nrestart/=0.or.nremap/=0.or.(.not.static.and..not.qh_evolve).or. &
         (.not.use_mond.and..not.qh_publish).or.mond_type/=1.or.bimond_twin_matter.or. &
         nlevelmax/=levelmin+1.or.any(nsubcycle(levelmin:nlevelmax)/=1).or..not.use_fftw)then
       write(*,*) 'QUMOND stage history requires serial fixed two-level static CPU FFTW, no restart/remap/subcycling'
       call clean_stop
    endif
    if(qh_evolve)then
       if(static.or.remap_thresh>0d0.or.lb_force_remap.or.use_sgs.or.pressure_fix.or. &
            eeos_poly_coeff>0d0.or.simple_boundary)then
          write(*,*) 'QUMOND evolving gas requires nonstatic fixed topology, no remap/SGS/pressure fix/eEOS/boundary'
          call clean_stop
       endif
#ifdef SOLVERmhd
       write(*,*) 'QUMOND evolving gas does not support MHD'
       call clean_stop
#endif
    endif
    if(qh_publish)then
       if(.not.hydro.or.pic.or.star.or.sink.or.cooling.or.neq_chem.or.rt.or.T2_star>0d0.or. &
            use_fR.or.use_nDGP.or.use_symmetron.or.use_dilaton.or.use_galileon.or.use_adm.or.use_fdm.or. &
            mond_mu_type/=1)then
          write(*,*) 'QUMOND composite publication supports fixed gas/analytic source only; no evolution or extra gravity'
          call clean_stop
       endif
    endif
    allocate(levels(size(xg,1)));levels=0
    do lev=levelmin,nlevelmax
       do g=1,active(lev)%ngrid
          levels(active(lev)%igrid(g))=lev
       enddo
    enddo
    if(.not.allocated(qh_current))then
       allocate(qh_current(size(phi),2),qh_previous(size(phi),2),qh_rho_newton(size(phi)))
       allocate(qh_time(2,nlevelmax),qh_old_time(2,nlevelmax),qh_valid(2,nlevelmax),qh_epoch(2,nlevelmax))
       allocate(qh_cell_level(size(phi)),qh_grid_level(size(xg,1)),qh_son(size(son)),qh_xg(size(xg,1),ndim))
       qh_rho_newton=0d0
       qh_current=0d0;qh_previous=0d0;qh_time=0d0;qh_old_time=0d0;qh_valid=.false.
       qh_epoch=-1
       qh_grid_level=levels;qh_son=son;qh_xg=0d0;qh_cell_level=0
       do lev=levelmin,nlevelmax
          do g=1,active(lev)%ngrid
             grid=active(lev)%igrid(g)
             qh_xg(grid,:)=xg(grid,:)
             do child=1,twotondim
                cell=ICELL_OF(grid,child)
                qh_cell_level(cell)=lev
             enddo
          enddo
       enddo
    else
       if(size(phi)/=size(qh_current,1).or.size(xg,1)/=size(qh_grid_level))then
          write(*,*) 'QUMOND stage history: allocation changed; reconstruction not implemented'
          call clean_stop
       endif
       if(any(levels/=qh_grid_level))then
          write(*,*) 'QUMOND stage history: hierarchy changed; reconstruction not implemented'
          call clean_stop
       endif
       do grid=1,size(levels)
          if(levels(grid)==0)cycle
          if(any(xg(grid,:)/=qh_xg(grid,:)))then
             write(*,*) 'QUMOND stage history: grid relocation unsupported'
             call clean_stop
          endif
       enddo
       do cell=1,size(phi)
          if(qh_cell_level(cell)==0)cycle
          if(son(cell)/=qh_son(cell))then
             write(*,*) 'QUMOND stage history: refinement changed; reconstruction not implemented'
             call clean_stop
          endif
       enddo
    endif
    deallocate(levels)
    qh_stage=stage
  end subroutine qumond_history_enter

  subroutine qumond_history_capture(stage,ilevel)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer,intent(in)::stage,ilevel
    integer::cell
    logical::first
    real(dp)::gauge,amplitude
    if(.not.qh_enabled)return
    if(qh_stage/=stage)call clean_stop
    first=.not.qh_valid(stage,ilevel)
    if(.not.first)then
       if(t<qh_time(stage,ilevel))then
          write(*,*) 'QUMOND stage history: time reversal'
          call clean_stop
       endif
       if(t>qh_time(stage,ilevel))then
          where(qh_cell_level==ilevel)qh_previous(:,stage)=qh_current(:,stage)
          qh_old_time(stage,ilevel)=qh_time(stage,ilevel)
       endif
    endif
    do cell=1,size(phi)
       if(qh_cell_level(cell)/=ilevel)cycle
       if(.not.ieee_is_finite(phi(cell)))then
          write(*,*) 'QUMOND stage history: nonfinite solved potential'
          call clean_stop
       endif
       qh_current(cell,stage)=phi(cell)
       if(stage==1)qh_rho_newton(cell)=rho(cell)
       if(first)qh_previous(cell,stage)=phi(cell)
    enddo
    if(ilevel==levelmin)then
       gauge=sum(pack(qh_current(:,stage),qh_cell_level==ilevel))/ &
            real(count(qh_cell_level==ilevel),dp)
       amplitude=maxval(abs(pack(qh_current(:,stage),qh_cell_level==ilevel)))
       if(abs(gauge)>1d-10*max(amplitude,tiny(1d0)))then
          write(*,*) 'QUMOND stage history: base periodic gauge is not canonical',gauge,amplitude
          call clean_stop
       endif
    endif
    qh_epoch(stage,ilevel)=nstep_coarse
    if(first)qh_old_time(stage,ilevel)=t
    qh_time(stage,ilevel)=t;qh_valid(stage,ilevel)=.true.
    write(*,*) 'RAR_QUMOND_HISTORY_V1',stage,ilevel,t,qh_old_time(stage,ilevel),first
    if(stage==1.and.ilevel==nlevelmax)call qumond_composite_operator_check()
  end subroutine qumond_history_capture

  subroutine qumond_composite_operator_check()
    character(len=8)::env
    integer::status,ierr,d,f,l,r,ncheck,ninterface
    real(dp),allocatable::p(:),op(:),ph(:),gradient(:,:)
    real(dp)::vector(3),target(3),constant_error,affine_error,band,flux_balance,phantom_balance
    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,a0,external(3),fourpi
    env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_OPERATOR',env,status=status)
    if(status/=0.or.trim(env)/='1')return
    if(qh_composite_mesh%n==0)then
       call composite_native_build(qh_composite_mesh,ierr)
       if(ierr/=0)then
          write(*,*) 'QUMOND composite operator: native topology failed',ierr
          call clean_stop
       endif
    endif
    associate(mesh=>qh_composite_mesh)
    allocate(p(mesh%n),op(mesh%n),ph(mesh%n),gradient(3,mesh%n))
    p=1d0;call composite_apply(mesh,p,op,ierr)
    if(ierr/=0)call clean_stop
    constant_error=maxval(abs(op));affine_error=0d0;ncheck=0;ninterface=0
    band=2d0*maxval(mesh%width)
    do d=1,3
       p=mesh%center(d,:);call composite_slopes(mesh,p,gradient,ierr)
       if(ierr/=0)call clean_stop
       target=0d0;target(d)=1d0
       do f=1,mesh%nface
          l=mesh%left(f);r=mesh%right(f)
          if(any(mesh%center(:,l)<=band).or.any(mesh%center(:,l)>=mesh%box-band).or. &
               any(mesh%center(:,r)<=band).or.any(mesh%center(:,r)>=mesh%box-band))cycle
          call composite_face_gradient(mesh,p,gradient,f,vector)
          affine_error=max(affine_error,maxval(abs(vector-target)));ncheck=ncheck+1
          if(mesh%width(l)/=mesh%width(r))ninterface=ninterface+1
       enddo
    enddo
    if(constant_error>1d-12.or.affine_error>1d-9.or.ncheck==0.or.ninterface==0)then
       write(*,*) 'QUMOND composite operator: constant/affine failure',constant_error,affine_error,ninterface
       call clean_stop
    endif
    p=qh_current(mesh%cell,1)
    call composite_apply(mesh,p,op,ierr)
    if(ierr/=0)call clean_stop
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    a0=a0_mond*scale_t**2/scale_l;external=g_ext_mond*scale_t**2/scale_l
    fourpi=4d0*acos(-1d0)
    call composite_phantom(mesh,p,a0,fourpi,external,ph,ierr)
    if(ierr/=0)call clean_stop
    flux_balance=abs(sum(op*mesh%volume))/max(sum(abs(op)*mesh%volume),tiny(1d0))
    phantom_balance=abs(sum(ph*mesh%volume))/max(sum(abs(ph)*mesh%volume),tiny(1d0))
    if(flux_balance>1d-12.or.phantom_balance>1d-12)then
       write(*,*) 'QUMOND composite operator: shared flux balance failure',flux_balance,phantom_balance
       call clean_stop
    endif
    write(*,*) 'RAR_QUMOND_COMPOSITE_V1',nstep_coarse,mesh%n,mesh%nface,ninterface, &
         constant_error,affine_error,flux_balance,phantom_balance
    if(nstep_coarse==0)call qumond_composite_solver_check()
    end associate
  end subroutine qumond_composite_operator_check

  subroutine qumond_composite_solver_check()
    character(len=8)::env
    integer::status,ierr,iterations
    real(dp),allocatable::exact(:),rhs(:),answer(:),potential_n(:),phantom(:),work(:)
    real(dp)::relative,error,twopi,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
    real(dp)::a0,fourpi,external(3),removed_mean
    env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_SOLVER',env,status=status)
    if(status/=0.or.trim(env)/='1')return
    if(mond_mu_type/=1)then
       write(*,*) 'QUMOND composite solver: only simple nu is implemented'
       call clean_stop
    endif
    associate(mesh=>qh_composite_mesh)
    allocate(exact(mesh%n),rhs(mesh%n),answer(mesh%n),potential_n(mesh%n),phantom(mesh%n),work(mesh%n))
    twopi=2d0*acos(-1d0)
    exact=sin(twopi*mesh%center(1,:)/mesh%box)+0.3d0*cos(2d0*twopi*mesh%center(2,:)/mesh%box) &
         +0.2d0*sin(twopi*mesh%center(3,:)/mesh%box)
    exact=exact-sum(exact*mesh%volume)/sum(mesh%volume)
    call composite_apply(mesh,exact,rhs,ierr)
    if(ierr/=0)call clean_stop
    answer=0d0
    call composite_solve(mesh,rhs,answer,1d-9,1200,40,iterations,relative,ierr)
    error=sqrt(sum(mesh%volume*(answer-exact)**2)/sum(mesh%volume*exact**2))
    write(*,*) 'RAR_QUMOND_SOLVER_V1',0,mesh%n,iterations,relative,error,ierr
    if(ierr/=0.or.error>1d-6)then
       write(*,*) 'QUMOND composite solver: manufactured recovery failed'
       call clean_stop
    endif
    answer=0d0
    call composite_solve(mesh,rhs,answer,1d-12,1,1,iterations,relative,ierr)
    if(ierr/=2.or.relative<=1d-12)then
       write(*,*) 'QUMOND composite solver: nonconvergence was not rejected'
       call clean_stop
    endif
    work=1d0;answer=0d0
    call composite_solve(mesh,work,answer,1d-9,10,5,iterations,relative,ierr)
    if(ierr/=4)then
       write(*,*) 'QUMOND composite solver: incompatible source was not rejected'
       call clean_stop
    endif
    write(*,*) 'RAR_QUMOND_SOLVER_REJECTION_V1',2,4
    fourpi=4d0*acos(-1d0)
    removed_mean=sum(qh_rho_newton(mesh%cell)*mesh%volume)/sum(mesh%volume)
    rhs=fourpi*(qh_rho_newton(mesh%cell)-removed_mean)
    potential_n=qh_current(mesh%cell,1)
    call composite_solve(mesh,rhs,potential_n,1d-9,1200,40,iterations,relative,ierr)
    write(*,*) 'RAR_QUMOND_SOLVER_V1',1,mesh%n,iterations,relative,removed_mean,ierr
    if(ierr/=0)then
       write(*,*) 'QUMOND composite solver: Newtonian stage failed'
       call clean_stop
    endif
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    a0=a0_mond*scale_t**2/scale_l;external=g_ext_mond*scale_t**2/scale_l
    call composite_phantom(mesh,potential_n,a0,fourpi,external,phantom,ierr)
    if(ierr/=0)call clean_stop
    rhs=rhs+fourpi*phantom
    rhs=rhs-sum(rhs*mesh%volume)/sum(mesh%volume)
    answer=potential_n
    call composite_solve(mesh,rhs,answer,1d-9,1200,40,iterations,relative,ierr)
    write(*,*) 'RAR_QUMOND_SOLVER_V1',2,mesh%n,iterations,relative, &
         sum(phantom*mesh%volume),ierr
    if(ierr/=0)then
       write(*,*) 'QUMOND composite solver: QUMOND stage failed'
       call clean_stop
    endif
    end associate
  end subroutine qumond_composite_solver_check

  real(dp) function qumond_history_value(cell,level,tfrac) result(value)
    integer,intent(in)::cell,level
    real(dp),intent(in)::tfrac
    integer::stage
    if(.not.qh_enabled)then
       value=phi(cell)+(phi(cell)-phi_old(cell))*tfrac
       return
    endif
    stage=qh_stage
    if(stage<1.or.stage>2.or.level<levelmin.or.level>nlevelmax)call clean_stop
    if(cell<1.or.cell>size(phi))call clean_stop
    if(qh_cell_level(cell)/=level.or..not.qh_valid(stage,level))then
       write(*,*) 'QUMOND stage history: unavailable same-stage coarse value',stage,level,cell
       call clean_stop
    endif
    if(qh_epoch(stage,level)/=nstep_coarse)then
       write(*,*) 'QUMOND stage history: stale coarse epoch'
       call clean_stop
    endif
    value=qh_current(cell,stage)+(qh_current(cell,stage)-qh_previous(cell,stage))*tfrac
  end function qumond_history_value


  subroutine bimond_source_mean(ilevel)
#include "amr_index.h"
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer,intent(in)::ilevel
    integer::i,ind,cell,info
    real(dp)::local_sum(2),global_sum(2)
    local_sum=0d0
    do i=1,active(ilevel)%ngrid
       do ind=1,twotondim
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          local_sum(1)=local_sum(1)+rho(cell)
          local_sum(2)=local_sum(2)+1d0
       enddo
    enddo
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_sum,global_sum,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    global_sum=local_sum
#endif
    rho_tot=global_sum(1)/global_sum(2)
  end subroutine bimond_source_mean

  subroutine bimond_prepare_difference(ilevel)
#include "amr_index.h"
    integer,intent(in)::ilevel
    integer::i,ind,cell
    call make_virtual_reverse_dp(bimond_rho_twin(1),ilevel)
    call make_virtual_fine_dp(bimond_rho_twin(1),ilevel)
    call make_virtual_reverse_dp(bimond_rho_ordinary(1),ilevel)
    call make_virtual_fine_dp(bimond_rho_ordinary(1),ilevel)
    do i=1,active(ilevel)%ngrid
       do ind=1,twotondim
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          rho(cell)=bimond_rho_ordinary(cell)-bimond_rho_twin(cell)
       enddo
    enddo
    call make_virtual_fine_dp(rho(1),ilevel)
    call bimond_source_mean(ilevel)
  end subroutine bimond_prepare_difference

  subroutine bimond_prepare_ordinary(ilevel)
#include "amr_index.h"
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer,intent(in)::ilevel
    integer::i,ind,cell,info
    real(dp)::local_sum(2),global_sum(2)
    local_sum=0d0
    ! compute_mond_phantom_density has added rho_ph to the signed source.
    ! The Poisson backends restore their mean-subtracted rho before return.
    do i=1,active(ilevel)%ngrid
       do ind=1,twotondim
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          bimond_rho_phantom(cell)=rho(cell)-(bimond_rho_ordinary(cell)-bimond_rho_twin(cell))
          local_sum(1)=local_sum(1)+bimond_rho_phantom(cell)
          local_sum(2)=local_sum(2)+abs(bimond_rho_phantom(cell))
          rho(cell)=bimond_rho_ordinary(cell)+bimond_rho_phantom(cell)
       enddo
    enddo
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_sum,global_sum,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    global_sum=local_sum
#endif
    if(myid==1)write(*,'(A,2ES24.15)') ' BIMOND phantom source sum, absolute sum: ',global_sum
    if(abs(global_sum(1))>1d-10*max(global_sum(2),tiny(1d0)))then
       if(myid==1)write(*,'(A)') 'ERROR: BIMOND periodic phantom source has a nonzero mean'
       call clean_stop
    endif
    call make_virtual_fine_dp(rho(1),ilevel)
    call bimond_source_mean(ilevel)
  end subroutine bimond_prepare_ordinary

  subroutine bimond_finish_fields(ilevel)
#include "amr_index.h"
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer,intent(in)::ilevel
    integer::i,ind,cell,idim,info
    real(dp)::acc_local,acc_global
    bimond_f_twin=f-bimond_f_difference
    do idim=1,ndim
       call make_virtual_fine_dp(bimond_f_twin(1,idim),ilevel)
    enddo
    acc_local=0d0
    do i=1,active(ilevel)%ngrid
       do ind=1,twotondim
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          acc_local=max(acc_local,sqrt(sum(f(cell,:)**2)),sqrt(sum(bimond_f_twin(cell,:)**2)))
       enddo
    enddo
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(acc_local,acc_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
    acc_global=acc_local
#endif
    bimond_acceleration_max=acc_global
    bimond_fields_ready=.true.
  end subroutine bimond_finish_fields

end module poisson_commons
