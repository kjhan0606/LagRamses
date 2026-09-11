program stellar_radioactive_runtime_test
  use amr_commons
  use hydro_commons
  use stellar_enrichment_config, only: configured_radioactive_model
  use stellar_radioactive_decay, only: radioactive_model
  use stellar_radioactive_transport, only: radioactive_constraints
  use stellar_radioactive_runtime
  implicit none
  real(dp)::saved(32,nvar),x(14),states(14,8),g(16),expected(2)
  real(dp)::qin(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar)
  real(dp)::qm(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim),qp(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim)
  real(dp)::u1(nvector,0:twondim,nvar),u2(nvector,twotondim,nvar)
  integer::idx(14),j,i,k,d,l
  configured_radioactive_model=radioactive_model
  active(1)%ngrid=1;active(1)%igrid(1)=1;son(2)=9;dtnew=1
  uold(:,1)=1;uold(:,5)=9;uold(:,8)=.3d0
  uold(:,9)=.6d0;uold(:,10)=.1d0;uold(:,15)=.04d0;uold(:,19)=.06d0
  uold(:,20)=.03d0;uold(:,21)=.02d0;unew=-7;saved=uold
  call radioactive_advance_level(1)
  if(fixture_stopped.or.fixture_virtual_calls/=4)stop 1
  if(any(uold(2,:)/=saved(2,:)).or.any(uold(9:,:)/=saved(9:,:)))stop 2
  if(any(unew/=-7).or.any(uold(:,1:14)/=saved(:,1:14)).or.any(uold(:,16:18)/=saved(:,16:18)))stop 3
  if(abs(uold(1,20)-.015d0)>1d-16.or.abs(uold(1,15)-.055d0)>1d-16)stop 4
  if(abs(uold(1,19)-uold(1,21)-.04d0)>1d-16)stop 5
  call radioactive_advance_level(1)
  if(fixture_stopped.or.abs(uold(1,20)-.0075d0)>1d-16)stop 6
  ! The callback stages ALL owned leaves before the first state write.
  uold(7,21)=1;saved=uold
  call radioactive_advance_level(1)
  if(.not.fixture_stopped.or.any(uold/=saved))stop 7
  fixture_stopped=.false.;configured_radioactive_model='none'
  call radioactive_advance_level(1)
  if(fixture_stopped.or.any(uold/=saved))stop 8
  configured_radioactive_model=radioactive_model
  print *, 'RADIOACTIVE_OWNED_LEAF_OLD_GAS_ONCE_COVERED_SKIP_NO_UNEW_OVERWRITE_ATOMIC_PASS'

  idx=[8,(j,j=9,19),20,21]
  x=0;x(1)=.3d0;x(2)=.6d0;x(3)=.1d0;x(8)=.04d0;x(12)=.06d0;x(13)=.03d0;x(14)=.02d0
  qin=0;qm=0;qp=0
  do k=ku1,ku2
     do j=ju1,ju2
        do i=iu1,iu2
           do l=1,nvector
              qin(l,i,j,k,idx)=x
              do d=1,ndim
                 qm(l,i,j,k,idx,d)=x;qp(l,i,j,k,idx,d)=x
                 qm(l,i,j,k,19,d)=.01d0;qp(l,i,j,k,19,d)=.11d0
                 qm(l,i,j,k,20,d)=.25d0;qp(l,i,j,k,20,d)=-.19d0
              enddo
           enddo
        enddo
     enddo
  enddo
  call radioactive_limit_faces(qin,qm,qp,nvector)
  if(fixture_stopped)stop 9
  ! Include halo-adjacent i/j/k=0: these are real Riemann donor faces.
  do k=0,3
     do j=0,3
        do i=0,3
           do d=1,ndim
              do l=1,nvector
                 g=radioactive_constraints(qm(l,i,j,k,idx,d));if(minval(g)<-1d-16)stop 10
                 g=radioactive_constraints(qp(l,i,j,k,idx,d));if(minval(g)<-1d-16)stop 11
              enddo
           enddo
        enddo
     enddo
  enddo
  if(qp(1,-1,-1,-1,20,1)/=-.19d0)stop 12
  u1=0;u2=0
  do l=1,nvector
     u1(l,0,idx)=x
     do j=1,8
        u2(l,j,idx)=x
        u2(l,j,19)=.06d0+(-1d0)**j*.05d0
        u2(l,j,20)=.03d0+(-1d0)**j*.22d0
     enddo
  enddo
  call radioactive_limit_children(u1,u2,nvector)
  if(fixture_stopped)stop 13
  do l=1,nvector
     if(maxval(abs(sum(u2(l,:,idx),dim=1)/8-x))>1d-16)stop 14
     do j=1,8
        g=radioactive_constraints(u2(l,j,idx));if(minval(g)<-1d-16)stop 15
     enddo
  enddo
  print *, 'RADIOACTIVE_ACTUAL_FACE_HALO_AND_CHILD_HOOKS_CONSERVATIVE_HOST_LIMITING_PASS'
end program
