! Physical sink bodies own gravity and motion; legacy clouds sample donor gas.
module qumond_composite_sinks
  use amr_commons
  use pm_commons
  use hydro_commons, only: uold
  use qumond_composite, only: composite_mesh_t
  use qumond_composite_clouds, only: native_cloud_map
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: native_sinks_source,native_sinks_force,native_sinks_cached
  public :: native_sinks_begin,native_sinks_drift,native_sinks_finish
  public :: native_sinks_reconcile,native_sinks_export,native_sinks_inventory,native_sinks_sync_clouds
  integer,allocatable :: identities(:),ptr(:),leaf(:)
  real(dp),allocatable :: masses(:),source_x(:,:),weights(:),accel(:,:)
  real(dp),allocatable :: start_x(:,:),start_v(:,:),start_a(:,:),half_v(:,:)
  integer :: nb=0
  real(dp) :: cloud_width=0d0,domain_width=0d0,displacement_limit=0d0
  logical :: initialized=.false.,step_open=.false.,drifted=.false.
contains
  subroutine fail(reason)
    character(len=*),intent(in) :: reason
    write(*,*) 'Native composite sinks failed: ',trim(reason)
    call clean_stop
  end subroutine

  subroutine initialize(mesh)
    type(composite_mesh_t),intent(in) :: mesh
    integer :: i,j,status,ios
    character(len=128) :: env
    if(nsink<0.or.nsink>nsinkmax.or..not.allocated(msink).or..not.allocated(idsink)) &
         call fail('invalid physical sink table')
    nb=nsink
    allocate(identities(nb),masses(nb),source_x(nb,3),accel(nb,3))
    allocate(start_x(nb,3),start_v(nb,3),start_a(nb,3),half_v(nb,3))
    identities=idsink(:nb);masses=msink(:nb)
    if(any(identities<=0))call fail('nonpositive physical sink identity')
    do i=1,nb
       do j=1,i-1
          if(identities(i)==identities(j))call fail('duplicate physical sink identity')
       enddo
    enddo
    if(.not.initialized)then
       cloud_width=mesh%box*0.5d0**levelmin
       env='';call get_environment_variable('RAR_NATIVE_SINK_CLOUD_WIDTH',env,status=status)
       if(status==-1)call fail('truncated physical cloud width')
       if(status==0.and.len_trim(env)>0)then
          read(env,*,iostat=ios)cloud_width
          if(ios/=0)call fail('invalid physical cloud width')
       endif
       if(.not.ieee_is_finite(cloud_width))call fail('nonfinite physical cloud width')
       if(cloud_width<mesh%box*0.5d0**levelmin.or.cloud_width>mesh%box) &
            call fail('physical sink cloud must span at least a base cell and at most the box')
       domain_width=mesh%box;displacement_limit=0.5d0*minval(mesh%width)
    endif
    source_x=-huge(1d0);accel=0d0;initialized=.true.
    call validate()
    write(*,*) 'RAR_NATIVE_SINK_MODEL_V1',nb,cloud_width
  end subroutine

  subroutine validate()
    if(.not.initialized.or.nsink/=nb)call fail('physical sink count changed outside source transaction')
    if(any(idsink(:nb)/=identities).or.any(msink(:nb)/=masses)) &
         call fail('physical sink identity or mass changed outside source transaction')
    if(any(masses<=0d0).or..not.all(ieee_is_finite(masses)).or. &
         .not.all(ieee_is_finite(xsink(:nb,:))).or..not.all(ieee_is_finite(vsink(:nb,:)))) &
         call fail('invalid physical sink state')
  end subroutine

  logical function native_sinks_cached() result(ok)
    ok=.false.
    if(.not.initialized)return
    call validate()
    ok=all(xsink(:nb,:)==source_x)
  end function

  subroutine native_sinks_source(mesh,density)
    type(composite_mesh_t),intent(in) :: mesh
    real(dp),intent(inout) :: density(:)
    real(dp),allocatable :: sink_rho(:)
    real(dp) :: max_error,integral,expected
    integer :: i,k,j
    if(.not.initialized)call initialize(mesh)
    call validate()
    if(mesh%box/=domain_width)call fail('physical sink domain changed')
    source_x=xsink(:nb,:)
    call native_cloud_map(mesh,source_x,cloud_width,ptr,leaf,weights,max_error)
    allocate(sink_rho(mesh%n));sink_rho=0d0
    do i=1,nb
       do k=ptr(i),ptr(i+1)-1
          j=leaf(k);sink_rho(j)=sink_rho(j)+masses(i)*weights(k)/mesh%volume(j)
       enddo
    enddo
    integral=sum(sink_rho*mesh%volume);expected=sum(masses)
    if(abs(integral-expected)>1d-12*max(expected,tiny(1d0)))call fail('physical sink deposit mass mismatch')
    density=density+sink_rho
    write(*,*) 'RAR_NATIVE_SINK_SOURCE_V1',nstep_coarse,t,nb,size(weights),expected,integral,max_error
  end subroutine

  subroutine native_sinks_force(mesh,force)
    type(composite_mesh_t),intent(in) :: mesh
    real(dp),intent(in) :: force(:,:)
    integer :: i,k
    if(.not.native_sinks_cached())call fail('stale physical sink source weights')
    accel=0d0
    do i=1,nb
       do k=ptr(i),ptr(i+1)-1
          accel(i,:)=accel(i,:)+weights(k)*force(mesh%cell(leaf(k)),:)
       enddo
    enddo
    if(.not.all(ieee_is_finite(accel)))call fail('nonfinite physical sink force')
  end subroutine

  subroutine native_sinks_begin(dt)
    real(dp),intent(inout) :: dt
    real(dp) :: speed,acc,bound
    integer :: i,d
    call validate()
    if(step_open)call fail('overlapping sink steps')
    start_x=xsink(:nb,:);start_v=vsink(:nb,:);start_a=accel
    do i=1,nb
       do d=1,3
          speed=abs(start_v(i,d));acc=abs(start_a(i,d))
          if(speed==0d0.and.acc==0d0)cycle
          bound=2d0*displacement_limit/(speed+sqrt(speed**2+2d0*acc*displacement_limit))
          dt=min(dt,bound)
       enddo
    enddo
    step_open=.true.;drifted=.false.
  end subroutine

  subroutine native_sinks_drift(dt)
    real(dp),intent(in) :: dt
    integer :: g,grid,p,j,i
    call validate()
    if(.not.step_open.or.drifted)call fail('invalid physical sink drift phase')
    if(any(xsink(:nb,:)/=start_x).or.any(vsink(:nb,:)/=start_v)) &
         call fail('legacy consumer changed physical sink during common step')
    half_v=start_v+0.5d0*dt*start_a
    if(nb>0)then
       if(maxval(abs(dt*half_v))>displacement_limit*(1d0+1d-12))call fail('physical sink displacement bound')
    endif
    do g=1,active(levelmin)%ngrid
       grid=active(levelmin)%igrid(g);p=headp(grid)
       do j=1,numbp(grid)
          if(ptypep(p)==PTYPE_SINK)then
             i=int(-idp(p));if(i<1.or.i>nb)call fail('invalid drifting cloud owner')
             xp(p,:)=xp(p,:)+dt*half_v(i,:);vp(p,:)=half_v(i,:)
          endif
          p=nextp(p)
       enddo
    enddo
    ! The ordinary native particle drift relinks all lists and wraps proxies.
    xsink(:nb,:)=modulo(start_x+dt*half_v,domain_width);vsink(:nb,:)=half_v
    drifted=.true.
  end subroutine

  subroutine native_sinks_finish(dt)
    real(dp),intent(in) :: dt
    integer :: i,unit,ios,status
    character(len=8) :: env
    character(len=80) :: filename
    call validate()
    if(.not.step_open.or..not.drifted)call fail('invalid physical sink endpoint kick')
    vsink(:nb,:)=half_v+0.5d0*dt*accel
    call native_sinks_sync_clouds(.true.)
    env='';call get_environment_variable('RAR_QUMOND_GRAVITY_KICKS',env,status=status)
    if(status==0.and.trim(env)=='1')then
    write(filename,'("sink_step_",I5.5,".tsv")')nstep_coarse
    open(newunit=unit,file=trim(filename),status='new',action='write',iostat=ios)
    if(ios/=0)call fail('sink step export exists or cannot be created')
    write(unit,'(A)') '# slot physical_id species mass dt x0[3] v0[3] a0[3] vhalf[3] x1[3] a1[3] v1[3]'
    do i=1,nb
       write(unit,'(3(I0,1X),23(ES24.16,1X))')i,identities(i),PTYPE_SINK,masses(i),dt, &
            start_x(i,:),start_v(i,:),start_a(i,:),half_v(i,:),xsink(i,:),accel(i,:),vsink(i,:)
    enddo
    close(unit)
    endif
    step_open=.false.;drifted=.false.
  end subroutine

  subroutine native_sinks_sync_clouds(canonical)
    logical,intent(in) :: canonical
    real(dp),allocatable :: cloud_mass(:)
    integer,allocatable :: counts(:)
    integer :: lev,g,grid,p,j,i,last_level
    allocate(cloud_mass(nsink),counts(nsink));cloud_mass=0d0;counts=0
    last_level=nlevelmax;if(canonical)last_level=levelmin
    do lev=levelmin,last_level
       do g=1,active(lev)%ngrid
          grid=active(lev)%igrid(g);p=headp(grid)
          do j=1,numbp(grid)
             if(ptypep(p)==PTYPE_SINK)then
                i=int(-idp(p));if(i<1.or.i>nsink)call fail('invalid sampling cloud owner')
                cloud_mass(i)=cloud_mass(i)+mp(p);counts(i)=counts(i)+1
                vp(p,:)=vsink(i,:)
             endif
             p=nextp(p)
          enddo
       enddo
    enddo
    if(any(counts==0))call fail('physical sink has no sampling clouds')
    if(any(abs(cloud_mass-msink(:nsink))>1d-12*msink(:nsink)))call fail('sampling cloud mass differs from physical body')
  end subroutine

  subroutine native_sinks_reconcile(mesh)
    type(composite_mesh_t),intent(in) :: mesh
    if(.not.initialized.or.step_open)call fail('sink source reconciliation inside open dynamics step')
    deallocate(identities,masses,source_x,accel,start_x,start_v,start_a,half_v)
    call initialize(mesh)
    call native_sinks_sync_clouds(.true.)
  end subroutine

  subroutine native_sinks_inventory(mesh,label)
    type(composite_mesh_t),intent(in) :: mesh
    character(len=*),intent(in) :: label
    real(dp) :: gas(5),body(5)
    integer :: d
    gas(1)=sum(uold(mesh%cell,1)*mesh%volume);body(1)=sum(msink(:nsink))
    do d=1,3
       gas(d+1)=sum(uold(mesh%cell,d+1)*mesh%volume)
       body(d+1)=sum(msink(:nsink)*vsink(:nsink,d))
    enddo
    gas(5)=sum(uold(mesh%cell,5)*mesh%volume)
    body(5)=0.5d0*sum(msink(:nsink)*sum(vsink(:nsink,:)**2,dim=2))
    write(*,*) 'RAR_NATIVE_SINK_INVENTORY_V1',trim(label),' ',nstep_coarse,t,nsink,gas,body
  end subroutine

  subroutine native_sinks_export(epoch)
    integer,intent(in) :: epoch
    integer :: unit,ios,i
    character(len=80) :: filename
    write(filename,'("sinks_epoch_",I5.5,".tsv")')epoch
    open(newunit=unit,file=trim(filename),status='new',action='write',iostat=ios)
    if(ios/=0)call fail('sink source export exists or cannot be created')
    write(unit,'(A,1X,ES24.16,1X,A)') '# cloud_width',cloud_width,'physical_sink'
    write(unit,'(A,1X,I0,1X,ES24.16)') '# epoch',nstep_coarse,t
    write(unit,'(A)') '# slot physical_id species mass x[3] v[3] acceleration[3]'
    do i=1,nb
       write(unit,'(3(I0,1X),10(ES24.16,1X))')i,identities(i),PTYPE_SINK,masses(i), &
            xsink(i,:),vsink(i,:),accel(i,:)
    enddo
    close(unit)
  end subroutine
end module qumond_composite_sinks
