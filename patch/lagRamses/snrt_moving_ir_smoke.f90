! Native interaction test with a prescribed conservative material callback.
! This exercises the real IR iteration, not a PAH physical-model calibration
! or a live RAMSES/AMR/production test. Compile with snrt_moving_scatter and
! snrt_dust_ir in an isolated directory; no Makefile integration is required.
program snrt_moving_ir_smoke
  use snrt_dust_ir
  implicit none
  integer,parameter::nb=2,ng=2,nd=6,nc=2
  real(dust_dp),parameter::c=2.99792458d10,ev=1.602176634d-12
  type(dust_ir_table)::table
  type(dust_ir_diagnostics)::diag
  real(dust_dp)::direction(3,nd),w(nd),e(ng,nd,nc),eold(ng,nd,nc),mat(nc),mat0(nc),temp(nc)
  real(dust_dp)::photons(ng,nc),pop(2,nc),oldpop(2,nc),rho(nb,nc),p(3,nb,nc),p0(3,nb,nc)
  real(dust_dp)::alpha(nb,ng,nc),sca(nb,ng,nc),work(nb,nc),heat_last(nb,ng,nc),events_last(nb,ng,nc)
  real(dust_dp)::emit_rate(nb,ng,nc),density(nc),primary(nc),capacity(nc),group_ev(ng)
  real(dust_dp)::dt,dx,chat,delta,dp(3),radp(3),total,ki,oldtemp(nc),oldphotons(ng,nc)
  integer::neighbors(6,nc),i,b,d,k,ierr,assertions,callbacks
  logical::blocked(6,nc),reject_material
  assertions=0;direction=0;w=1d0/nd;group_ev=[.01d0,.02d0]
  do k=1,3
     direction(k,2*k-1)=1;direction(k,2*k)=-1
  enddo
  call snrt_dust_ir_initialize(table,group_ev,[.005d0,.01d0],[.1d0,.2d0],[5d0,30d0,100d0],5d0,ierr)
  call require(ierr==0,'table initialization')

  call reset();e=0;e(:,1,:)=.006d0;p0=p;call save()
  call advance()
  call require(ierr==0,'beam absorption')
  call require(all(p(1,:,:)>0).and.sum(work)>0,'finite absorption recoil work')
  call check_closure('beam absorption')
  call require(sum(mat-mat0)<sum(eold-e)/nd,'material receives absorption minus work')
  call require(any(abs(heat_last-events_last*spread(spread(group_ev*ev,1,nb),3,nc))>1d-17), &
       'captured events are distinct from deposited heat')

  call reset();e=0;alpha=0
  do i=1,nc
     p(1,1,i)=rho(1,i)*c*.003d0;p(2,2,i)=rho(2,i)*c*(-.002d0)
  enddo
  emit_rate(1,1,:)=.001d0;emit_rate(2,2,:)=.002d0
  call save();call advance()
  call require(ierr==0,'moving phase emission')
  call require(sum(work)<0.and.sum(e)>0,'emission includes recoil and lab energy boost')
  call require(p(1,1,1)<p0(1,1,1).and.p(2,2,1)>p0(2,2,1),'distinct emitting-phase recoil')
  call check_closure('emission')

  call reset();e(:,1,:)=.004d0
  emit_rate(1,1,:)=.0001d0;emit_rate(2,2,:)=.0002d0
  do i=1,nc
     p(:,1,i)=rho(1,i)*c*[.002d0,-.001d0,0d0]
     p(:,2,i)=rho(2,i)*c*[-.001d0,.002d0,.001d0]
  enddo
  sca=.4d0;call save();call advance()
  call require(ierr==0,'absorption emission and scatter')
  call check_closure('combined interactions')
  call require(callbacks>1,'actual nonlinear material iteration exercised')

  call reset();alpha=0;sca=0;blocked=.false.
  neighbors(:,1)=2;neighbors(:,2)=1
  e=0;e(:,1,1)=.01d0;call save();call advance()
  call require(ierr==0,'spatial transport without interactions')
  call require(any(e/=eold),'transport actually changes radiation')
  call require(all(p==p0).and.all(work==0),'spatial transport is not a force')
  call check_closure('transport only')

  call reset();call save();reject_material=.true.
  oldtemp=temp;oldphotons=photons;oldpop=pop
  call advance()
  call require(ierr/=0,'material callback rejection')
  call require(all(e==eold).and.all(p==p0).and.all(mat==mat0).and.all(work==-123), &
       'radiation material momentum work rollback')
  call require(all(temp==oldtemp).and.all(photons==oldphotons).and.all(pop==oldpop),'population diagnostics rollback')
  call require(diag%mechanical_erg==-456,'diagnostic transaction rollback')

  call reset();call save();alpha(2,1,2)=-1;call advance()
  call require(ierr/=0.and.all(e==eold).and.all(p==p0).and.all(work==-123),'bad phase opacity rollback')
  write(*,'(A,I0)')'MOVING_IR_SMOKE_PASS assertions=',assertions
contains
  subroutine reset()
    density=1;primary=0;capacity=1;mat=.1d0;temp=30;photons=0;pop=1
    e=.001d0;p=0;rho=1d-20;rho(2,:)=2d-20
    alpha(1,:,:)=.2d0;alpha(2,:,:)=.4d0;sca=0;emit_rate=0
    dt=.2d0;dx=1;chat=1;blocked=.true.
    neighbors=0
    work=-123;diag=dust_ir_diagnostics();diag%mechanical_erg=-456
    reject_material=.false.;callbacks=0
  end subroutine
  subroutine save()
    eold=e;p0=p;mat0=mat
  end subroutine
  subroutine advance()
    call snrt_dust_ir_advance(table,direction,w,neighbors,dx,dt,chat,density,primary,e,temp,photons,diag,ierr, &
         1d-11,160,mat,capacity,blocked_face=blocked,population=pop,cell_absorption=sum(alpha,dim=1), &
         thin_reabsorption=.true., &
         phase_density=rho,phase_momentum=p,phase_absorption=alpha,phase_scattering=sca,phase_work=work, &
         moving_material_dispatch=material)
  end subroutine
  subroutine material(heat,events,step_dt,old,next,phase_rate,next_energy,next_t,status)
    real(dust_dp),intent(in)::heat(:,:,:),events(:,:,:),step_dt,old(:,:)
    real(dust_dp),intent(out)::next(:,:),phase_rate(:,:,:),next_energy(:),next_t(:)
    integer,intent(out)::status
    callbacks=callbacks+1;status=dust_ok
    if(reject_material)then
       status=dust_err_state;return
    endif
    next=old;phase_rate=emit_rate;next_t=30
    next_energy=mat0+sum(sum(heat-step_dt*phase_rate,dim=1),dim=1)
    heat_last=heat;events_last=events
  end subroutine
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       write(*,'(A,A,A,I0,A,I0)')'FAIL: ',label,' ierr=',ierr,' callbacks=',callbacks
       error stop 1
    endif
    assertions=assertions+1
  end subroutine
  subroutine check_closure(label)
    character(len=*),intent(in)::label
    total=sum(eold)/nd+sum(mat0);ki=0;radp=0;dp=0
    do i=1,nc
       do b=1,nb
          delta=dot_product(.5d0*(p(:,b,i)+p0(:,b,i))/rho(b,i),p(:,b,i)-p0(:,b,i))
          ki=ki+delta;dp=dp+p(:,b,i)-p0(:,b,i)
          call require(abs(delta-work(b,i))<1d-12*total,label//' per-phase work')
       enddo
       do d=1,nd
          radp=radp+sum(e(:,d,i)-eold(:,d,i))*w(d)*direction(:,d)/c
       enddo
    enddo
    delta=sum(e-eold)/nd+sum(mat-mat0)+ki
    call require(abs(delta)<2d-11*total,label//' energy closure')
    call require(maxval(abs(dp+radp))*c<2d-11*total,label//' momentum closure')
    call require(abs(diag%mechanical_erg-ki)<1d-12*total,label//' signed diagnostic work')
    call require(all(e>=0).and.all(mat>=0),label//' positive state')
  end subroutine
end program
