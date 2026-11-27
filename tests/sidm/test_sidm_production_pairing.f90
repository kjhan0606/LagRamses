! Adapted from SHA-256 d37ac2bc2e895984c44c4fa62d4e9b3c978fd8e64eb7083be15154b4509d3386.
! The historical expected-events diagnostic is absent from this main version;
! expected-rate assertions and all 13 production scatter calls are retained.
module pm_commons
  implicit none
  integer,parameter :: IRandNumSize=4,MAX_PART=128,MAX_GRID=16
  integer(kind=1),parameter :: PTYPE_DM=0,PTYPE_ISIDM_EX1=10,PTYPE_ISIDM_EX2=11
  integer(kind=1),parameter :: PTYPE_STAR=1,PTYPE_SINK=2
  real(kind(1.0d0)) :: xp(MAX_PART,3),vp(MAX_PART,3),mp(MAX_PART)
  integer(kind=8) :: idp(MAX_PART)
  integer :: headp(MAX_GRID),nextp(MAX_PART),numbp(MAX_GRID)
  integer(kind=1) :: ptypep(MAX_PART)
contains
  integer function ptype_to_isidm_state(pt)
    integer(kind=1),intent(in)::pt
    write(*,'(A,I0)') 'Unexpected inelastic state lookup for ptype ',pt
    ptype_to_isidm_state=-99
    error stop 2
  end function ptype_to_isidm_state

  function isidm_state_to_ptype(istate) result(pt)
    integer,intent(in)::istate
    integer(kind=1)::pt
    write(*,'(A,I0)') 'Unexpected inelastic state update for state ',istate
    pt=-99
    error stop 2
  end function isidm_state_to_ptype
end module pm_commons

module amr_parameters
  implicit none
  integer,parameter :: twotondim=8
  integer :: amr_block_size=64
end module amr_parameters

module amr_commons
  use amr_parameters,only:twotondim,amr_block_size
  use pm_commons,only:MAX_GRID
  implicit none
  private :: MAX_GRID
  integer,parameter :: dp=kind(1.0d0),MAX_CELL=512
  integer :: icoarse_min,icoarse_max,jcoarse_min,jcoarse_max
  integer :: kcoarse_min,kcoarse_max,ncoarse,myid
  integer :: next(MAX_GRID),son(MAX_CELL)
  real(dp) :: xg(MAX_GRID,3),boxlen,dtnew(32),aexp
  real(dp) :: unit_scale_l,unit_scale_t,unit_scale_d,unit_scale_v
  real(dp) :: unit_scale_nH,unit_scale_T2
  integer :: sidm_npart_min,sidm_nstates
  real(dp) :: sidm_vrel_max,sidm_cross_section,sidm_v0,sidm_power
  real(dp) :: sidm_a_transition,sidm_sigma_ratio,sidm_a_width
  real(dp) :: sidm_epsilon,sidm_fdiss,sidm_mchi,sidm_energy(0:9)
  logical :: sidm_inelastic,sidm_excited_init,sidm_probability_guard=.false.
  real(dp) :: sidm_courant=0.1d0
  character(len=32) :: sidm_type,sidm_a_type,sidm_angular
end module amr_commons

module random
  implicit none
  integer,parameter :: MAX_RANDOM=4096
  real(kind(1.0d0)) :: scripted_random(MAX_RANDOM)
  integer :: scripted_count=0,scripted_next=1
contains
  subroutine rng_set(count,values)
    integer,intent(in)::count
    real(kind(1.0d0)),intent(in)::values(:)
    if(count<0.or.count>MAX_RANDOM.or.count>size(values)) error stop 3
    scripted_random=0.0d0
    scripted_random(1:count)=values(1:count)
    scripted_count=count
    scripted_next=1
  end subroutine rng_set

  subroutine ranf(seed,value)
    integer,intent(inout)::seed(:)
    real(kind(1.0d0)),intent(out)::value
    if(scripted_next>scripted_count) then
       write(*,'(A,I0)') 'Scripted RNG exhausted at draw ',scripted_next
       error stop 4
    end if
    value=scripted_random(scripted_next)
    scripted_next=scripted_next+1
  end subroutine ranf

  subroutine rng_require_consumed()
    if(scripted_next/=scripted_count+1) then
       write(*,'(A,I0,A,I0)') 'Unused scripted draws: consumed ', &
            scripted_next-1,' of ',scripted_count
       error stop 5
    end if
  end subroutine rng_require_consumed
end module random

program test_sidm_production_pairing
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use pm_commons
  use amr_commons
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: mpi_test_rank,mpi_test_size,mpi_test_ierr
#endif
  integer :: sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  integer :: sidm_n_vrel_reject
  integer :: sidm_n_leaf_occupied,sidm_n_leaf_active,sidm_n_dm_leaf
  real(dp) :: sidm_dp(3),sidm_dEk,sidm_dp_max,sidm_dEk_max
  real(dp) :: sidm_Pmax_level,sidm_dEdiss,sidm_expected_rate
  integer :: sidm_diag_enabled
  common /sidm_diag/ sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down, &
       sidm_n_vrel_reject
  common /sidm_cons/ sidm_dp,sidm_dEk,sidm_dp_max,sidm_dEk_max
  common /sidm_pmax_c/ sidm_Pmax_level
  common /sidm_diss/ sidm_dEdiss
  common /sidm_estimator_counts/ sidm_n_leaf_occupied, &
       sidm_n_leaf_active,sidm_n_dm_leaf
  common /sidm_estimator_rate/ sidm_expected_rate
  common /sidm_diag_config/ sidm_diag_enabled

  integer,parameter :: MAX_FIXTURE_PART=128
  integer :: grid_of(MAX_FIXTURE_PART),expected_child(MAX_FIXTURE_PART)
  logical :: expected_eligible(MAX_FIXTURE_PART)
  integer :: ncreated,production_calls
  real(dp) :: dx_mesh,dx_local,mesh_scale,particle_mass_code
  real(dp),parameter :: reject_draw=0.999999d0,accept_draw=1.0d-7
  integer :: seeds(IRandNumSize,0:0)
  real(dp) :: velocity_fixture(4,3)
  character(len=32) :: requested_case

  velocity_fixture(1,:)=(/ 0.10d0, 0.02d0,-0.01d0 /)
  velocity_fixture(2,:)=(/-0.05d0, 0.20d0, 0.04d0 /)
  velocity_fixture(3,:)=(/ 0.03d0,-0.12d0, 0.08d0 /)
  velocity_fixture(4,:)=(/ 0.22d0, 0.05d0,-0.03d0 /)

  call get_command_argument(1,requested_case)
  if(len_trim(requested_case)>0) then
#ifndef WITHOUTMPI
     if(trim(requested_case)=='mpi_invalid') then
        call MPI_INIT(mpi_test_ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_test_rank,mpi_test_ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_test_size,mpi_test_ierr)
        if(mpi_test_size/=2) call MPI_ABORT(MPI_COMM_WORLD,2,mpi_test_ierr)
        call setup_single_cell(2)
        myid=mpi_test_rank+1
        if(mpi_test_rank==1) then
           sidm_cross_section=1.0d0
        else
           sidm_cross_section=0.0d0
        endif
        call rng_set(4,(/0.90d0,0.10d0,0.20d0,0.30d0/))
        write(*,'(A,I0)') 'MPI_GUARD_READY_RANK=',mpi_test_rank
        flush(6)
        call MPI_BARRIER(MPI_COMM_WORLD,mpi_test_ierr)
        call run_production_scatter(16,myid,2,1,seeds,0,1)
        call MPI_BARRIER(MPI_COMM_WORLD,mpi_test_ierr)
        print '(A)', 'UNEXPECTED: invalid probability did not abort MPI'
        call MPI_FINALIZE(mpi_test_ierr)
        stop
     endif
#endif
     if(trim(requested_case)=='zero'.or.trim(requested_case)=='one') then
        call test_probability_boundary(trim(requested_case))
        stop
     endif
     call test_invalid_probability(trim(requested_case))
     print '(A)', 'UNEXPECTED: invalid probability was accepted'
     stop
  endif

  production_calls=0
  call test_membership_and_eligibility()
  call test_pairing_and_conservation(3)
  call test_pairing_and_conservation(4)
  call test_probability_normalization(3)
  call test_probability_normalization(4)
  if(production_calls/=13) error stop 15
  write(*,'(A,I0)') 'PRODUCTION_CALLS=',production_calls
  write(*,'(A)') 'SIDM production leaf selection and pairing: PASS'

contains
  subroutine reset_case()
    integer :: ip
    xp=0.0d0; vp=0.0d0; mp=0.0d0
    idp=0; ptypep=PTYPE_DM
    headp=0; nextp=0; numbp=0
    next=0; son=0; xg=0.0d0
    dtnew=0.0d0
    ncreated=0; grid_of=0; expected_child=0; expected_eligible=.false.
    boxlen=2.5d0
    icoarse_min=3; icoarse_max=6
    jcoarse_min=1; jcoarse_max=4
    kcoarse_min=0; kcoarse_max=3
    ncoarse=0; myid=1; aexp=1.0d0
    amr_block_size=64
    xg(2,:)=(/3.5d0,1.5d0,0.5d0/)
    xg(5,:)=(/4.5d0,2.5d0,1.5d0/)
    dx_mesh=0.5d0**16
    mesh_scale=boxlen/dble(icoarse_max-icoarse_min+1)
    dx_local=dx_mesh*mesh_scale
    unit_scale_l=1.0d0/dx_local
    ! Keep physical velocity and time conversions dimensionally consistent.
    unit_scale_t=1.0d4; unit_scale_d=1.0d0
    unit_scale_v=unit_scale_l/unit_scale_t
    unit_scale_nH=1.0d0; unit_scale_T2=1.0d0
    particle_mass_code=0.1d0/unit_scale_l**3
    dtnew(16)=1.0d-2
    sidm_npart_min=2
    sidm_vrel_max=0.0d0
    sidm_cross_section=1.0d-6
    sidm_v0=1.0d0; sidm_power=0.0d0
    sidm_a_transition=0.5d0; sidm_sigma_ratio=1.0d0
    sidm_a_width=0.1d0; sidm_epsilon=0.01d0
    sidm_fdiss=0.0d0; sidm_mchi=1.0d0
    sidm_nstates=2; sidm_energy=0.0d0
    sidm_inelastic=.false.; sidm_excited_init=.false.
    sidm_type='constant'; sidm_a_type='none'; sidm_angular='isotropic'
    sidm_n_scatter=0; sidm_n_pairs=0; sidm_n_up=0; sidm_n_down=0
    sidm_n_vrel_reject=0
    sidm_n_leaf_occupied=0; sidm_n_leaf_active=0; sidm_n_dm_leaf=0
    sidm_dp=0.0d0; sidm_dEk=0.0d0; sidm_dp_max=0.0d0; sidm_dEk_max=0.0d0
    sidm_Pmax_level=0.0d0; sidm_dEdiss=0.0d0
    sidm_expected_rate=0.0d0
    sidm_diag_enabled=1
    seeds(:,0)=(/11,23,37,41/)
    do ip=1,MAX_FIXTURE_PART
       grid_of(ip)=0; expected_child(ip)=0; expected_eligible(ip)=.false.
    end do
  end subroutine reset_case

  subroutine run_production_scatter(ilevel,icpu,kgrid,subnump,thread_seeds, &
       mythread_loc,nthreads_in)
    integer,intent(in)::ilevel,icpu,kgrid,subnump,mythread_loc,nthreads_in
    integer,intent(inout)::thread_seeds(IRandNumSize,0:nthreads_in-1)
    production_calls=production_calls+1
    call sub_sidm_scatter(ilevel,icpu,kgrid,subnump,thread_seeds, &
         mythread_loc,nthreads_in,dtnew(ilevel),dtnew(ilevel))
  end subroutine run_production_scatter

  subroutine add_particle(grid,child,ptype,particle_id,mesh_position,eligible)
    integer,intent(in)::grid,child,particle_id
    integer(kind=1),intent(in)::ptype
    real(dp),intent(in)::mesh_position(3)
    logical,intent(in)::eligible
    integer :: ip
    if(ncreated>=MAX_FIXTURE_PART) error stop 6
    ncreated=ncreated+1; ip=ncreated
    grid_of(ip)=grid; expected_child(ip)=child
    expected_eligible(ip)=eligible
    xp(ip,:)=(mesh_position-(/dble(icoarse_min),dble(jcoarse_min), &
         dble(kcoarse_min)/))*mesh_scale
    idp(ip)=particle_id; ptypep(ip)=ptype
    mp(ip)=particle_mass_code
    vp(ip,:)=(/1.0d-2*dble(ip),2.0d-3*dble(ip), &
         3.0d-3*dble(mod(ip,7))/)
  end subroutine add_particle

  subroutine add_center_particle(grid,child,ptype,particle_id,eligible)
    integer,intent(in)::grid,child,particle_id
    integer(kind=1),intent(in)::ptype
    logical,intent(in)::eligible
    integer :: bits(3),axis
    real(dp) :: position(3)
    bits(1)=mod(child-1,2)
    bits(2)=mod((child-1)/2,2)
    bits(3)=mod((child-1)/4,2)
    position=xg(grid,:)
    do axis=1,3
       position(axis)=position(axis)+dble(2*bits(axis)-1)*0.5d0*dx_mesh
    end do
    call add_particle(grid,child,ptype,particle_id,position,eligible)
  end subroutine add_center_particle

  subroutine link_particles(scramble)
    logical,intent(in)::scramble
    integer :: grid,ip,nlocal,k,left,right,local(MAX_FIXTURE_PART)
    integer :: order(MAX_FIXTURE_PART),tmp
    headp=0; nextp=0; numbp=0
    do grid=1,MAX_GRID
       nlocal=0
       do ip=1,ncreated
          if(grid_of(ip)==grid) then
             nlocal=nlocal+1; local(nlocal)=ip
          end if
       end do
       if(nlocal==0) cycle
       if(scramble) then
          left=1; right=nlocal
          do k=1,nlocal
             if(mod(k,2)==1) then
                order(k)=local(right); right=right-1
             else
                order(k)=local(left); left=left+1
             end if
          end do
       else
          order(1:nlocal)=local(1:nlocal)
       end if
       headp(grid)=order(1); numbp(grid)=nlocal
       do k=1,nlocal-1
          nextp(order(k))=order(k+1)
       end do
       nextp(order(nlocal))=0
    end do
  end subroutine link_particles

  integer function oracle_count(grid,child)
    integer,intent(in)::grid,child
    integer :: ip
    oracle_count=0
    do ip=1,ncreated
       if(grid_of(ip)==grid.and.expected_child(ip)==child.and.expected_eligible(ip)) &
            oracle_count=oracle_count+1
    end do
  end function oracle_count

  subroutine append_membership_sequence()
    integer :: grids(2),ig,grid,child,n,draw_count,seq_count
    real(dp) :: sequence(MAX_RANDOM)
    grids=(/2,5/); seq_count=0
    do ig=1,2
       grid=grids(ig)
       do child=1,8
          if(son(cell_index(grid,child))/=0) cycle
          n=oracle_count(grid,child)
          if(n<sidm_npart_min) cycle
          draw_count=n-1+n/2
          if(grid==2.and.child==2) then
             if(n/=2) error stop 7
             if(seq_count+4>MAX_RANDOM) error stop 8
             sequence(seq_count+1:seq_count+4)=(/0.75d0,accept_draw,0.5d0,0.0d0/)
             seq_count=seq_count+4
          else
             if(seq_count+draw_count>MAX_RANDOM) error stop 9
             sequence(seq_count+1:seq_count+draw_count)=reject_draw
             seq_count=seq_count+draw_count
          end if
       end do
    end do
    call rng_set(seq_count,sequence)
  end subroutine append_membership_sequence

  integer function cell_index(grid,child)
    integer,intent(in)::grid,child
    cell_index=ncoarse+((grid-1)/amr_block_size) &
         *(twotondim*amr_block_size)+(child-1)*amr_block_size &
         +mod(grid-1,amr_block_size)+1
  end function cell_index

  subroutine expected_membership_counters(noccupied,nactive,ndm,npairs)
    integer,intent(out)::noccupied,nactive,ndm,npairs
    integer :: grids(2),ig,grid,child,n
    grids=(/2,5/); noccupied=0; nactive=0; ndm=0; npairs=0
    do ig=1,2
       grid=grids(ig)
       do child=1,8
          if(son(cell_index(grid,child))/=0) cycle
          n=oracle_count(grid,child)
          if(n>0) then
             noccupied=noccupied+1; ndm=ndm+n
          end if
          if(n>=sidm_npart_min) then
             nactive=nactive+1; npairs=npairs+n/2
          end if
       end do
    end do
  end subroutine expected_membership_counters

  real(dp) function independent_probability(i,j,n,npairs)
    integer,intent(in)::i,j,n,npairs
    real(dp) :: m1,m2,vrel
    m1=mp(i)*unit_scale_d*unit_scale_l**3
    m2=mp(j)*unit_scale_d*unit_scale_l**3
    vrel=sqrt(sum(((vp(i,:)-vp(j,:))*unit_scale_v)**2))
    independent_probability=sidm_cross_section*0.5d0*(m1+m2)*vrel &
         *(dtnew(16)*unit_scale_t)/(dx_local*unit_scale_l)**3 &
         *dble(n*(n-1))/dble(2*npairs)
  end function independent_probability

  subroutine expected_membership_rates(events,pmax)
    real(dp),intent(out)::events,pmax
    integer :: grids(2),ig,grid,child,n,npairs,ip,particle_list(MAX_FIXTURE_PART),k
    real(dp) :: probability
    events=0.0d0; pmax=0.0d0; grids=(/2,5/)
    do ig=1,2
       grid=grids(ig)
       do child=1,8
          if(son(cell_index(grid,child))/=0) cycle
          n=oracle_count(grid,child)
          if(n<sidm_npart_min) cycle
          npairs=n/2; k=0; ip=headp(grid)
          do while(ip/=0)
             if(expected_child(ip)==child.and.expected_eligible(ip)) then
                k=k+1; particle_list(k)=ip
             end if
             ip=nextp(ip)
          end do
          if(k/=n) error stop 8
          do ip=1,npairs
             probability=independent_probability(particle_list(2*ip-1), &
                  particle_list(2*ip),n,npairs)
             events=events+probability; pmax=max(pmax,probability)
          end do
       end do
    end do
  end subroutine expected_membership_rates

  subroutine test_membership_and_eligibility()
    integer :: child,repeat,particle_id,axis,state,other_bits(3),bits(3),bit
    integer :: expected_child_id,noccupied,nactive,ndm,npairs,ip
    real(dp) :: position(3),saved_vp(MAX_PART,3),expected_vp(MAX_PART,3)
    real(dp) :: saved_xp(MAX_PART,3)
    real(dp) :: saved_mp(MAX_PART),events,pmax
    integer(kind=8) :: saved_idp(MAX_PART)
    integer(kind=1) :: saved_ptypep(MAX_PART)
    call reset_case()
    call assert_integer(cell_index(2,1),2,'L16 block child index for grid 2')
    call assert_integer(cell_index(5,1),5,'noncontiguous L16 grid index')
    call assert_integer(cell_index(5,8),453,'L16 block child stride')
    next(2)=5; next(5)=0
    call add_center_particle(2,1,PTYPE_DM,1,.true.)
    call add_center_particle(2,1,PTYPE_STAR,2,.false.)
    call add_center_particle(2,2,PTYPE_DM,3,.true.)
    call add_center_particle(2,2,PTYPE_ISIDM_EX1,4,.true.)
    do repeat=1,2
       call add_center_particle(2,3,PTYPE_DM,10+repeat,.true.)
    end do
    do repeat=1,3
       call add_center_particle(2,4,PTYPE_DM,20+repeat,.true.)
    end do
    do repeat=1,4
       call add_center_particle(2,5,PTYPE_DM,30+repeat,.true.)
    end do
    do repeat=1,2
       call add_center_particle(2,6,PTYPE_DM,40+repeat,.true.)
    end do
    do repeat=1,2
       call add_center_particle(2,7,PTYPE_DM,50+repeat,.true.)
    end do
    call add_center_particle(2,7,PTYPE_STAR,60,.false.)
    call add_center_particle(2,7,PTYPE_SINK,61,.false.)
    call add_center_particle(2,7,PTYPE_DM,0,.false.)
    call add_center_particle(2,7,PTYPE_DM,-1,.false.)
    call add_center_particle(2,8,PTYPE_DM,70,.true.)
    call add_center_particle(2,8,PTYPE_ISIDM_EX2,71,.true.)
    son(cell_index(2,6))=1

    do axis=1,3
       other_bits=0
       select case(axis)
       case(1); other_bits(2)=0; other_bits(3)=1
       case(2); other_bits(1)=1; other_bits(3)=0
       case(3); other_bits(1)=0; other_bits(2)=1
       end select
       do state=-1,1
          bits=other_bits
          if(state>=0) then
             bits(axis)=1
          else
             bits(axis)=0
          end if
          expected_child_id=1+bits(1)+2*bits(2)+4*bits(3)
          position=xg(5,:)
          do bit=1,3
             if(bit==axis) then
                if(state<0) position(bit)=xg(5,bit)-0.25d0*dx_mesh
                if(state>0) position(bit)=xg(5,bit)+0.25d0*dx_mesh
             else
                position(bit)=xg(5,bit)+dble(2*other_bits(bit)-1) &
                     *0.25d0*dx_mesh
             end if
          end do
          do repeat=1,2
             particle_id=100+10*axis+3*state+repeat
             call add_particle(5,expected_child_id,PTYPE_DM,particle_id,position,.true.)
          end do
       end do
    end do
    call add_center_particle(5,8,PTYPE_STAR,500,.false.)
    call add_center_particle(5,8,PTYPE_SINK,501,.false.)
    call link_particles(.true.)
    call append_membership_sequence()
    saved_vp=vp; saved_xp=xp; saved_mp=mp
    saved_idp=idp; saved_ptypep=ptypep
    expected_vp=saved_vp
    call reference_isotropic_equal_mass(saved_vp(3,:),saved_vp(4,:), &
         expected_vp(3,:),expected_vp(4,:))
    call run_production_scatter(16,1,2,2,seeds,0,1)
    call rng_require_consumed()
    call expected_membership_counters(noccupied,nactive,ndm,npairs)
    call assert_integer(sidm_n_leaf_occupied,noccupied,'occupied leaf count')
    call assert_integer(sidm_n_leaf_active,nactive,'active leaf count')
    call assert_integer(sidm_n_dm_leaf,ndm,'eligible DM count')
    call assert_integer(sidm_n_pairs,npairs,'sampled pair count')
    call assert_integer(sidm_n_scatter,1,'accepted DM/EX1 scatter count')
    call assert_integer(sidm_n_vrel_reject,0,'pathological velocity rejection count')
    call assert_true(all(ieee_is_finite(vp)).and. &
         all(ieee_is_finite(expected_vp)).and. &
         maxval(abs(vp-expected_vp))<2.0d-12, &
         'explicit DM/EX1 partner scatters; rejected and excluded particles stay unchanged')
    call assert_true(all(xp==saved_xp),'particle positions remain unchanged')
    call assert_true(all(mp==saved_mp),'particle masses remain unchanged')
    call assert_true(all(idp==saved_idp),'particle IDs remain unchanged')
    call assert_true(all(ptypep==saved_ptypep),'particle types remain unchanged')
    call expected_membership_rates(events,pmax)
    call assert_close(sidm_expected_rate,events/(dtnew(16)*unit_scale_t), &
         'expected rate')
    call assert_close(sidm_Pmax_level,pmax,'maximum pair probability')
    call assert_true(pmax<0.05d0,'all membership probabilities remain below 0.05')
  end subroutine test_membership_and_eligibility

  subroutine setup_single_cell(n)
    integer,intent(in)::n
    integer :: ip
    call reset_case()
    next(2)=0
    do ip=1,n
       call add_center_particle(2,1,PTYPE_DM,ip,.true.)
       vp(ip,:)=velocity_fixture(ip,:)
    end do
    call link_particles(.false.)
  end subroutine setup_single_cell

  subroutine test_invalid_probability(which)
    character(len=*),intent(in) :: which
    call setup_single_cell(2)
    select case(which)
    case('above_one')
       sidm_cross_section=1.0d0
    case('negative')
       sidm_cross_section=-1.0d-6
    case('nonfinite')
       sidm_cross_section=ieee_value(0.0d0,ieee_quiet_nan)
    case default
       error stop 'unknown probability fixture'
    end select
    call rng_set(4,(/0.90d0,0.10d0,0.20d0,0.30d0/))
    call run_production_scatter(16,1,2,1,seeds,0,1)
  end subroutine test_invalid_probability

  subroutine test_probability_boundary(which)
    character(len=*),intent(in) :: which
    call setup_single_cell(2)
    unit_scale_l=1.0d0; unit_scale_t=1.0d0
    unit_scale_d=1.0d0; unit_scale_v=1.0d0
    dtnew(16)=1.0d0
    mp(1:2)=1.0d0
    vp(1,:)=(/0.0d0,0.0d0,0.0d0/)
    vp(2,:)=(/1.0d0,0.0d0,0.0d0/)
    if(which=='zero') then
       sidm_cross_section=0.0d0
       call rng_set(2,(/0.90d0,0.10d0/))
    else
       ! In this two-particle setup, probability = sigma / cell volume.
       sidm_cross_section=dx_local**3
       call rng_set(4,(/0.90d0,0.10d0,0.20d0,0.30d0/))
    endif
    call run_production_scatter(16,1,2,1,seeds,0,1)
    call rng_require_consumed()
    if(which=='zero') then
       call assert_close(sidm_Pmax_level,0.0d0,'zero probability boundary')
       call assert_integer(sidm_n_scatter,0,'zero probability scattering')
    else
       call assert_close(sidm_Pmax_level,1.0d0,'unit probability boundary')
       call assert_integer(sidm_n_scatter,1,'unit probability scattering')
    endif
    print '(A,A)', 'Valid probability boundary passed: ',which
  end subroutine test_probability_boundary

  subroutine matching_case(n,case_id,order,shuffle,nshuffle)
    integer,intent(in)::n,case_id
    integer,intent(out)::order(4),nshuffle
    real(dp),intent(out)::shuffle(3)
    order=0; shuffle=0.0d0
    if(n==3) then
       nshuffle=2
       select case(case_id)
       case(1); order(1:3)=(/1,2,3/); shuffle(1:2)=(/0.90d0,0.75d0/)
       case(2); order(1:3)=(/1,3,2/); shuffle(1:2)=(/0.50d0,0.75d0/)
       case(3); order(1:3)=(/3,2,1/); shuffle(1:2)=(/0.10d0,0.75d0/)
       case default; error stop 9
       end select
    else if(n==4) then
       nshuffle=3
       select case(case_id)
       case(1); order=(/1,2,3,4/); shuffle=(/0.90d0,0.90d0,0.75d0/)
       case(2); order=(/1,3,2,4/); shuffle=(/0.90d0,0.50d0,0.75d0/)
       case(3); order=(/1,4,3,2/); shuffle=(/0.375d0,0.90d0,0.75d0/)
       case default; error stop 10
       end select
    else
       error stop 11
    end if
  end subroutine matching_case

  subroutine set_matching_random(n,nshuffle,shuffle,accept_pairs)
    integer,intent(in)::n,nshuffle
    real(dp),intent(in)::shuffle(3)
    logical,intent(in)::accept_pairs
    real(dp) :: sequence(MAX_RANDOM)
    integer :: count,k,pairs
    count=0
    do k=1,nshuffle
       count=count+1; sequence(count)=shuffle(k)
    end do
    pairs=n/2
    do k=1,pairs
       count=count+1
       if(accept_pairs) then
          sequence(count)=accept_draw
          count=count+1; sequence(count)=0.5d0
          count=count+1; sequence(count)=0.0d0
       else
          sequence(count)=reject_draw
       end if
    end do
    call rng_set(count,sequence)
  end subroutine set_matching_random

  subroutine selected_pair_expectation(n,order,events,pmax)
    integer,intent(in)::n,order(4)
    real(dp),intent(out)::events,pmax
    integer :: pair_index,i,j
    real(dp) :: probability
    events=0.0d0; pmax=0.0d0
    do pair_index=1,n/2
       i=order(2*pair_index-1); j=order(2*pair_index)
       probability=independent_probability(i,j,n,n/2)
       events=events+probability; pmax=max(pmax,probability)
    end do
  end subroutine selected_pair_expectation

  subroutine test_pairing_and_conservation(n)
    integer,intent(in)::n
    integer :: case_id,order(4),nshuffle,ip,pair_index,i,j
    real(dp) :: shuffle(3),saved_vp(MAX_PART,3),expected_vp(MAX_PART,3)
    real(dp) :: saved_xp(MAX_PART,3),saved_mp(MAX_PART)
    integer(kind=8) :: saved_idp(MAX_PART)
    integer(kind=1) :: saved_ptypep(MAX_PART)
    real(dp) :: events,pmax,dp_before(3),dp_after(3),energy_before,energy_after
    do case_id=1,3
       call setup_single_cell(n)
       call matching_case(n,case_id,order,shuffle,nshuffle)
       call set_matching_random(n,nshuffle,shuffle,.true.)
       saved_vp=vp; expected_vp=vp; saved_xp=xp; saved_mp=mp
       saved_idp=idp; saved_ptypep=ptypep
       call selected_pair_expectation(n,order,events,pmax)
       call assert_true(pmax<0.05d0,'accepted matching probability below 0.05')
       call run_production_scatter(16,1,2,1,seeds,0,1)
       call rng_require_consumed()
       do pair_index=1,n/2
          i=order(2*pair_index-1); j=order(2*pair_index)
          call reference_isotropic_equal_mass(saved_vp(i,:),saved_vp(j,:), &
               expected_vp(i,:),expected_vp(j,:))
       end do
       call assert_true(all(ieee_is_finite(vp)).and. &
            all(ieee_is_finite(expected_vp)).and. &
            maxval(abs(vp-expected_vp))<2.0d-12, &
            'production partner identities and elastic velocities')
       call assert_integer(sidm_n_pairs,n/2,'production sampled-pair count')
       call assert_integer(sidm_n_scatter,n/2,'accepted scatter count')
       call assert_integer(sidm_n_leaf_occupied,1,'single-cell occupied count')
       call assert_integer(sidm_n_leaf_active,1,'single-cell active count')
       call assert_integer(sidm_n_dm_leaf,n,'single-cell DM count')
       call assert_close(sidm_expected_rate,events/(dtnew(16)*unit_scale_t), &
            'selected pair rate')
       call assert_close(sidm_Pmax_level,pmax,'selected pair Pmax')
       call assert_true(all(xp==saved_xp),'elastic scatter leaves positions unchanged')
       call assert_true(all(mp==saved_mp),'elastic scatter leaves masses unchanged')
       call assert_true(all(idp==saved_idp),'elastic scatter leaves IDs unchanged')
       call assert_true(all(ptypep==saved_ptypep), &
            'elastic scatter leaves particle types unchanged')
       dp_before=0.0d0; dp_after=0.0d0
       energy_before=0.0d0; energy_after=0.0d0
       do ip=1,n
          dp_before=dp_before+saved_vp(ip,:)*unit_scale_v*mp(ip) &
               *unit_scale_d*unit_scale_l**3
          dp_after=dp_after+vp(ip,:)*unit_scale_v*mp(ip) &
               *unit_scale_d*unit_scale_l**3
          energy_before=energy_before+0.5d0*mp(ip)*unit_scale_d*unit_scale_l**3 &
               *sum((saved_vp(ip,:)*unit_scale_v)**2)
          energy_after=energy_after+0.5d0*mp(ip)*unit_scale_d*unit_scale_l**3 &
               *sum((vp(ip,:)*unit_scale_v)**2)
       end do
       call assert_true(all(ieee_is_finite(dp_before)).and. &
            all(ieee_is_finite(dp_after)).and. &
            ieee_is_finite(energy_before).and.ieee_is_finite(energy_after).and. &
            maxval(abs(dp_after-dp_before))<5.0d-12, &
            'independent elastic momentum conservation')
       call assert_close(energy_after,energy_before,'independent elastic energy conservation')
    end do
  end subroutine test_pairing_and_conservation

  subroutine reference_isotropic_equal_mass(v1,v2,out1,out2)
    real(dp),intent(in)::v1(3),v2(3)
    real(dp),intent(out)::out1(3),out2(3)
    real(dp) :: vcm(3),vrel(3),magnitude
    vcm=0.5d0*(v1+v2); vrel=v1-v2; magnitude=sqrt(sum(vrel**2))
    out1=vcm+0.5d0*magnitude*(/1.0d0,0.0d0,0.0d0/)
    out2=vcm-0.5d0*magnitude*(/1.0d0,0.0d0,0.0d0/)
  end subroutine reference_isotropic_equal_mass

  subroutine test_probability_normalization(n)
    integer,intent(in)::n
    integer :: case_id,order(4),nshuffle
    real(dp) :: shuffle(3),events,pmax,matched_average,all_pairs_sum
    real(dp) :: vrel,mass,dt,volume
    integer :: i,j
    matched_average=0.0d0
    do case_id=1,3
       call setup_single_cell(n)
       call matching_case(n,case_id,order,shuffle,nshuffle)
       call set_matching_random(n,nshuffle,shuffle,.false.)
       call selected_pair_expectation(n,order,events,pmax)
       call assert_true(pmax<0.05d0,'normalization probabilities below 0.05')
       call run_production_scatter(16,1,2,1,seeds,0,1)
       call rng_require_consumed()
       call assert_close(sidm_expected_rate,events/(dtnew(16)*unit_scale_t), &
            'normalization selected-pair rate')
       call assert_close(sidm_Pmax_level,pmax,'normalization selected-pair Pmax')
       call assert_integer(sidm_n_pairs,n/2,'normalization sampled-pair count')
       call assert_integer(sidm_n_scatter,0,'normalization case rejected all scatters')
       matched_average=matched_average+sidm_expected_rate*dtnew(16)*unit_scale_t/3.0d0
    end do
    mass=mp(1)*unit_scale_d*unit_scale_l**3
    dt=dtnew(16)*unit_scale_t
    volume=(dx_local*unit_scale_l)**3
    all_pairs_sum=0.0d0
    do i=1,n-1
       do j=i+1,n
          vrel=sqrt(sum(((velocity_fixture(i,:)-velocity_fixture(j,:)) &
               *unit_scale_v)**2))
          all_pairs_sum=all_pairs_sum+sidm_cross_section*mass*dt*vrel/volume
       end do
    end do
    call assert_close(matched_average,all_pairs_sum, &
         'average production matching rate equals independent all-pairs target')
  end subroutine test_probability_normalization

  subroutine assert_integer(actual,expected,label)
    integer,intent(in)::actual,expected
    character(len=*),intent(in)::label
    if(actual/=expected) then
       write(*,'(A,A,A,I0,A,I0)') 'FAIL: ',trim(label), &
            ' got ',actual,' expected ',expected
       error stop 12
    end if
  end subroutine assert_integer

  subroutine assert_true(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if(.not.condition) then
       write(*,'(A,A)') 'FAIL: ',trim(label)
       error stop 13
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,label)
    real(dp),intent(in)::actual,expected
    character(len=*),intent(in)::label
    real(dp) :: tolerance
    if(.not.ieee_is_finite(actual).or..not.ieee_is_finite(expected)) then
       write(*,'(A,A)') 'FAIL: non-finite value in ',trim(label)
       error stop 16
    end if
    tolerance=5.0d-12*max(1.0d0,abs(expected))
    if(abs(actual-expected)>tolerance) then
       write(*,'(A,A,A,ES24.16,A,ES24.16,A,ES12.4)') 'FAIL: ',trim(label), &
            ' got ',actual,' expected ',expected,' tolerance ',tolerance
       error stop 14
    end if
  end subroutine assert_close
end program test_sidm_production_pairing

subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons,only:dp,unit_scale_l,unit_scale_t,unit_scale_d, &
       unit_scale_v,unit_scale_nH,unit_scale_T2
  implicit none
  real(dp),intent(out)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  scale_l=unit_scale_l; scale_t=unit_scale_t; scale_d=unit_scale_d
  scale_v=unit_scale_v; scale_nH=unit_scale_nH; scale_T2=unit_scale_T2
end subroutine units

! The test isolates the real scatter routine from the wrapper.  Guarded
! subcycling is exercised by the integrated tests, not this legacy fixture.
subroutine sidm_probability_guard_abort(message)
  implicit none
  character(len=*),intent(in)::message
  write(*,*) trim(message)
  error stop 1
end subroutine sidm_probability_guard_abort
