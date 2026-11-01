program test_sidm_particle_geometry
  use, intrinsic :: iso_fortran_env, only: real64
  use sidm_particle_geometry, only: sidm_particle_subcell_index
  implicit none
  include 'mpif.h'

  integer,parameter::ndim=3,nx_loc=64
  real(real64),parameter::pi=3.1415926535897932384626433832795_real64
  real(real64)::boxlen,scale,skip_loc(ndim),xg(ndim),xp(ndim)
  real(real64)::xcoord(ndim),eps,split_delta
  real(real64)::trajectory_start(ndim),trajectory_velocity(ndim),trajectory_time
  integer::mask,axis,idim,case_id,got,expected,expected_mask,test_count
  integer::myid,nproc,info

  call MPI_INIT(info)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,info)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,info)

  boxlen=8.0_real64*pi
  scale=boxlen/real(nx_loc,real64)
  skip_loc=[4.0_real64,11.0_real64,23.0_real64]
  xg=[8.5_real64,18.5_real64,31.5_real64]
  eps=1.0e-12_real64
  test_count=0

  ! All eight octants in a non-unit box with nonzero global coarse offsets.
  do mask=0,7
     do idim=1,ndim
        if(btest(mask,idim-1)) then
           xcoord(idim)=xg(idim)+0.125_real64
        else
           xcoord(idim)=xg(idim)-0.125_real64
        end if
        xp(idim)=(xcoord(idim)-skip_loc(idim))*scale
     end do
     got=sidm_particle_subcell_index(xp,xg,scale,skip_loc,ndim)
     call require_rank_consistent(got,'octant classification')
     expected=oracle_index(xp,xg,scale,skip_loc,ndim)
     call require(got==mask+1,'octant classification')
     call require(got==expected,'independent integer oracle')
     test_count=test_count+1
  end do

  ! Exercise just below, exactly on, and just above every split plane.
  do axis=1,ndim
     xcoord=xg-0.125_real64
     do idim=1,ndim
        if(idim==axis) cycle
        xcoord(idim)=xg(idim)+0.125_real64
     end do
     do case_id=1,3
        select case(case_id)
        case(1)
           split_delta=-eps
        case(2)
           split_delta=0.0_real64
        case(3)
           split_delta=eps
        end select
        xcoord(axis)=xg(axis)+split_delta
        xp=(xcoord-skip_loc)*scale
        got=sidm_particle_subcell_index(xp,xg,scale,skip_loc,ndim)
        call require_rank_consistent(got,'split-plane boundary classification')
        expected=oracle_index(xp,xg,scale,skip_loc,ndim)
        call require(got==expected,'split-plane boundary classification')
        if(case_id==2) call require(btest(got-1,axis-1), &
             'center belongs to upper half')
        test_count=test_count+1
     end do
  end do

  ! Follow deterministic particle paths through each split plane in both
  ! directions.  At the exact crossing, the >= convention assigns the
  ! particle to the upper half; a reverse path then moves to the lower half.
  do axis=1,ndim
     trajectory_start=xg+0.25_real64
     trajectory_start(axis)=xg(axis)-0.25_real64
     trajectory_velocity=0.0_real64
     trajectory_velocity(axis)=0.5_real64
     do case_id=0,4
        trajectory_time=0.25_real64*real(case_id,real64)
        xcoord=trajectory_start+trajectory_time*trajectory_velocity
        xp=(xcoord-skip_loc)*scale
        got=sidm_particle_subcell_index(xp,xg,scale,skip_loc,ndim)
        call require_rank_consistent(got,'forward split-plane trajectory')
        expected_mask=0
        do idim=1,ndim
           if(xcoord(idim)>=xg(idim)) expected_mask=expected_mask+2**(idim-1)
        end do
        expected=expected_mask+1
        call require(got==expected,'forward split-plane trajectory')
        test_count=test_count+1
     end do

     trajectory_start=xg+0.25_real64
     trajectory_velocity=0.0_real64
     trajectory_velocity(axis)=-0.5_real64
     do case_id=0,4
        trajectory_time=0.25_real64*real(case_id,real64)
        xcoord=trajectory_start+trajectory_time*trajectory_velocity
        xp=(xcoord-skip_loc)*scale
        got=sidm_particle_subcell_index(xp,xg,scale,skip_loc,ndim)
        call require_rank_consistent(got,'reverse split-plane trajectory')
        expected_mask=0
        do idim=1,ndim
           if(xcoord(idim)>=xg(idim)) expected_mask=expected_mask+2**(idim-1)
        end do
        expected=expected_mask+1
        call require(got==expected,'reverse split-plane trajectory')
        test_count=test_count+1
     end do
  end do

  ! A diagonal trajectory crosses all three planes at the same sample time.
  trajectory_start=xg-0.25_real64
  trajectory_velocity=0.5_real64
  do case_id=0,4
     trajectory_time=0.25_real64*real(case_id,real64)
     xcoord=trajectory_start+trajectory_time*trajectory_velocity
     xp=(xcoord-skip_loc)*scale
     got=sidm_particle_subcell_index(xp,xg,scale,skip_loc,ndim)
     call require_rank_consistent(got,'diagonal octant trajectory')
     expected_mask=0
     do idim=1,ndim
        if(xcoord(idim)>=xg(idim)) expected_mask=expected_mask+2**(idim-1)
     end do
     expected=expected_mask+1
     call require(got==expected,'diagonal octant trajectory')
     test_count=test_count+1
  end do

  ! Periodic upper-edge records are wrapped to the lower edge before the
  ! production classifier, matching particle_tree's global-origin mapping.
  xg=skip_loc+0.5_real64
  xp=(real(nx_loc,real64)-skip_loc)*scale
  do axis=1,ndim
     xcoord=xp/scale+skip_loc
     if(xcoord(axis)>=real(nx_loc,real64)) then
        xp(axis)=xp(axis)-(real(nx_loc,real64)-skip_loc(axis))*scale
     end if
  end do
  got=sidm_particle_subcell_index(xp,xg,scale,skip_loc,ndim)
  call require_rank_consistent(got,'periodic edge classification')
  call require(got==1,'periodic upper edge maps to lower edge owner')
  expected=oracle_index(xp,xg,scale,skip_loc,ndim)
  call require(got==expected,'periodic edge integer oracle')
  test_count=test_count+1

  if(myid==0) write(*,'(A,I0,A,I0)') &
       'PASS sidm_particle_geometry cases=',test_count,' MPI ranks=',nproc
  call MPI_FINALIZE(info)

contains

  integer function oracle_index(position,grid_center,cell_scale,origin,nd) &
       result(cell_index)
    real(real64),intent(in)::position(:),grid_center(:),cell_scale,origin(:)
    integer,intent(in)::nd
    integer::d,child_bit
    real(real64)::global_coordinate

    cell_index=1
    do d=1,nd
       global_coordinate=position(d)/cell_scale+origin(d)
       child_bit=floor(2.0_real64*(global_coordinate-grid_center(d))+ &
            1.0_real64)
       if(child_bit==1) cell_index=cell_index+2**(d-1)
    end do
  end function oracle_index

  subroutine require(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) then
       write(*,'(A)') 'FAIL sidm_particle_geometry: '//trim(message)
       error stop 1
    end if
  end subroutine require

  subroutine require_rank_consistent(value,message)
    integer,intent(in)::value
    character(len=*),intent(in)::message
    integer::rank_min,rank_max,mpi_error
    call MPI_ALLREDUCE(value,rank_min,1,MPI_INTEGER,MPI_MIN, &
         MPI_COMM_WORLD,mpi_error)
    call MPI_ALLREDUCE(value,rank_max,1,MPI_INTEGER,MPI_MAX, &
         MPI_COMM_WORLD,mpi_error)
    call require(rank_min==rank_max,'rank-independent '//trim(message))
  end subroutine require_rank_consistent

end program test_sidm_particle_geometry
