program agn_feedback_deposition_smoke
  use amr_parameters, only: dp
  use agn_feedback_deposition
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  real(dp) :: row(5), before(5), deferred, dm, momentum(3), energy
  real(dp) :: volumes(3), offsets(3), total_mass, total_momentum(3), total_energy
  real(dp) :: weights(2), lobe_sum(2), axis(3)
  integer :: ierr, i

  row=[1d0,0d0,0d0,0d0,5d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],3d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==0 .and. row(5)==8d0 .and. deferred==0d0,'thermal below cap')
  row=[1d0,0d0,0d0,0d0,5d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],8d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==0 .and. row(5)==10d0 .and. deferred==9d0,'cap crossing and nonunit volume')
  call check(abs(row(5)*3d0+deferred-39d0)<1d-12,'gas plus deferred equals old plus input')
  row=[1d0,0d0,0d0,0d0,12d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],4d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==0 .and. row(5)==12d0 .and. deferred==12d0,'already hot gas is not cooled')

  ! A donor with rho=2, v=3, internal energy density 4 loses drho=.5.
  ! average_AGN preserves velocity/internal energy, so its post-removal row
  ! is [1.5,4.5,0,0,10.75]. Return at the captured DONOR velocity, not v_BH.
  row=[1.5d0,4.5d0,0d0,0d0,10.75d0]
  dm=0.5d0
  momentum=[1.5d0,0d0,0d0]
  energy=4d0+0.5d0*dm*3d0**2
  call agn_deposit_cell(row,dm,momentum,energy,2d0,2d0,1d0,100d0,deferred,ierr)
  call check(ierr==0 .and. maxval(abs(row-[2d0,6d0,0d0,0d0,17d0]))<1d-12 .and. &
       deferred==0d0,'moving donor fallback restores mass momentum and kinetic energy')

  ! Synthetic cylinder with unequal cell volumes, including its midplane.
  ! Kernel shape=1; each lobe has volume-weighted norm=2.
  volumes=[1d0,2d0,1d0]; offsets=[-1d0,0d0,1d0]
  lobe_sum=0d0
  do i=1,3
     call agn_jet_geometry([offsets(i),0d0,0d0],[1d0,0d0,0d0],1d0,weights,axis)
     lobe_sum=lobe_sum+weights*volumes(i)
  enddo
  total_mass=0d0; total_momentum=0d0; total_energy=0d0
  do i=1,3
     call agn_jet_geometry([offsets(i),0d0,0d0],[1d0,0d0,0d0],1d0,weights,axis)
     call agn_jet_delta(2d0,weights,lobe_sum,[0d0,0d0,0d0],axis, &
          sqrt(8d0),dm,momentum,energy)
     if(i==2)call check(all(momentum==0d0) .and. abs(energy-2d0)<1d-12, &
          'midplane retains opposed-stream energy with zero net momentum')
     total_mass=total_mass+dm*volumes(i)
     total_momentum=total_momentum+momentum*volumes(i)
     row=[1d0,0d0,0d0,0d0,1d0]
     call agn_deposit_cell(row,dm,momentum,energy,volumes(i),2d0,1d0,100d0,deferred,ierr)
     call check(ierr==0,'cylinder cell deposition')
     total_energy=total_energy+(row(5)-1d0)*volumes(i)+deferred
  end do
  call check(abs(total_mass-2d0)<1d-12 .and. maxval(abs(total_momentum))<1d-12 .and. &
       abs(total_energy-8d0)<1d-12,'cylinder sums loaded mass and supplied kinetic energy')

  row=[1d0,0d0,0d0,0d0,5d0]; before=row
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],ieee_value(0d0,ieee_quiet_nan), &
       3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==agn_deposit_invalid_source .and. all(row==before) .and. deferred==0d0, &
       'invalid source is non-mutating')
  row=[1d0,3d0,0d0,0d0,1d0]; before=row
  call agn_deposit_cell(row,0.5d0,[1.5d0,0d0,0d0],8d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==agn_deposit_invalid_receiver .and. all(row==before) .and. deferred==0d0, &
       'negative incoming internal energy rejects; do not silently lose returned mass')
  call check(agn_eddington_ratio(0d0,0d0)==0d0,'idle Eddington ratio is finite zero')
  call entrainment_event([0d0,0d0,0d0],.false.,100d0)
  call entrainment_event([17d0,-11d0,5d0],.false.,1d0)
  call entrainment_event([0d0,0d0,0d0],.true.,1d0)
  call layouts_and_payload()
  write(*,'(A)') 'AGN_NATIVE_CELL_COUPLING_SMOKE_OK'
contains
  subroutine entrainment_event(boost,missing_lobe,cap)
    real(dp),intent(in)::boost(3),cap
    logical,intent(in)::missing_lobe
    real(dp)::gas(9,4),initial(9,4),vol(4),pos(3,3),lw(2),norms(2),ax(3)
    real(dp)::mass,vel(3),fractions(3),donor_vel(3),receiver_vel(3),drho,dpvec(3),de
    real(dp)::saved,dsave,change(9),ejet,thermal,jet_speed,donor_specific
    real(dp)::packed(7),unpacked_mass,unpacked_vel(3),unpacked_frac(3),zero_pack(7)
    integer::fields(3),k,status
    fields=[6,8,9]; vol=[2d0,1d0,2d0,3d0]
    donor_vel=[3d0,1d0,-0.5d0]+boost
    receiver_vel=[-0.2d0,0.3d0,0.1d0]+boost
    gas(:,1)=[8d0,8d0*donor_vel,12d0+4d0*sum(donor_vel**2),1.6d0,7d0,4.8d0,0.08d0]
    do k=2,4
       gas(:,k)=[2d0,2d0*receiver_vel,2d0+sum(receiver_vel**2),0.1d0,3d0,1.4d0,0.002d0]
    enddo
    initial=gas
    pos(:,1)=[-0.5d0,0.2d0,0d0]
    pos(:,2)=[0d0,0.1d0,0d0]
    pos(:,3)=[0.5d0,0.6d0,0d0]
    if(missing_lobe)pos(1,:)=0.5d0
    norms=0d0
    do k=1,3
       call agn_jet_geometry(pos(:,k),[1d0,0d0,0d0],1d0,lw,ax)
       norms=norms+lw*vol(k+1)
    enddo
    call agn_withdraw_cell(gas(:,1),fields,1,8d0,vol(1),mass,vel,fractions,status)
    call check(status==0 .and. abs(mass-4d0)<1d-12,'withdrawal uses actual donor and 25 percent cap')
    call check(maxval(abs(gas(fields,1)/gas(1,1)-initial(fields,1)/initial(1,1)))<1d-12, &
         'donor fractions remain unchanged on mass removal')
    donor_specific=(gas(5,1)-0.5d0*sum(gas(2:4,1)**2)/gas(1,1))/gas(1,1)
    call check(abs(donor_specific-2d0)<1d-10,'cold-loading donor specific heat rises by 4/3')
    call agn_pack_load(mass,vel,fractions,packed)
    call agn_pack_load(0d0,[0d0,0d0,0d0],[0d0,0d0,0d0],zero_pack)
    call agn_unpack_load(packed+zero_pack,unpacked_mass,unpacked_vel,unpacked_frac)
    call check(unpacked_mass==mass .and. all(unpacked_vel==vel) .and. all(unpacked_frac==fractions), &
         'synthetic remote payload retains donor velocity and composition through packed sum')
    ejet=20d0; saved=0d0
    if(any(norms<=0d0))then
       drho=mass/vol(1); dpvec=drho*vel; de=ejet/vol(1)+0.5d0*drho*sum(vel**2)
       call agn_deposit_material(gas(:,1),fields,1,fractions,drho,dpvec,de,vol(1), &
            2d0,1d0,cap,dsave,status)
       call check(status==0,'missing lobe returns whole event at donor')
       saved=dsave
       call check(all(gas(:,2:4)==initial(:,2:4)),'missing lobe deposits no partial distributed jet')
       call check(maxval(abs(gas(fields,1)-initial(fields,1)))<1d-12,'fallback restores all transported species')
    else
       jet_speed=sqrt(2d0*0.6d0*ejet/mass)
       thermal=0.4d0*ejet/sum(vol(2:4))
       do k=1,3
          call agn_jet_geometry(pos(:,k),[1d0,0d0,0d0],1d0,lw,ax)
          call agn_jet_delta(unpacked_mass,lw,norms,unpacked_vel,ax,jet_speed,drho,dpvec,de)
          call agn_deposit_material(gas(:,k+1),fields,1,unpacked_frac,drho,dpvec,de+thermal, &
               vol(k+1),2d0,1d0,cap,dsave,status)
          call check(status==0,'asymmetric mixed kinetic and thermal receiver')
          saved=saved+dsave
       enddo
       call check(gas(9,2)/gas(1,2)>initial(9,2)/initial(1,2), &
            'receiver abundance changes toward donor; not merely global species invariance')
    endif
    change=0d0
    do k=1,4
       change=change+(gas(:,k)-initial(:,k))*vol(k)
    enddo
    call check(abs(change(1))<1d-11 .and. maxval(abs(change(2:4)))<1d-10, &
         'whole event conserves mass and vector momentum under boost')
    call check(abs(change(5)+saved-ejet)<1d-9,'whole event gas plus deferred energy equals EAGN under boost')
    if(cap<=1d0)call check(saved>0d0,'low-cap event actually exercises deferred energy')
    call check(maxval(abs(change(fields)))<1d-11 .and. all(gas(7,:)==initial(7,:)), &
         'whole event species conserved and delayed reservoir untouched')
  end subroutine entrainment_event

  subroutine layouts_and_payload()
    integer::fields(3),phase_fields(12),empty(0),status
    real(dp)::gas(9),snapshot(9),mass,vel(3),frac(3),dsave
    call agn_scalar_map(9,6,8,2,[7,0,0,0],fields,status)
    call check(status==0 .and. all(fields==[6,8,9]),'legacy chemical layout excludes delayed reservoir')
    call agn_scalar_map(18,6,8,11,[7,0,0,0],phase_fields,status)
    call check(status==0 .and. phase_fields(12)==18,'channel-resolved layout includes all eleven stored elements')
    call agn_scalar_map(5,0,6,0,[0,0,0,0],empty,status)
    call check(status==0,'empty chemistry layout is valid')
    call agn_scalar_map(9,6,7,2,[7,0,0,0],fields,status)
    call check(status/=0,'chemical reservoir overlap rejected')
    call agn_scalar_map(8,6,8,2,[7,0,0,0],fields,status)
    call check(status/=0,'chemical extent rejected')
    fields=[6,8,9]
    gas=[8d0,24d0,0d0,0d0,48d0,1.6d0,7d0,4.8d0,0.08d0]
    call agn_withdraw_cell(gas,fields,1,4d0,2d0,mass,vel,frac,status)
    call check(status==0 .and. gas(1)==6d0,'first shared donor withdrawal')
    call agn_withdraw_cell(gas,fields,1,4d0,2d0,mass,vel,frac,status)
    call check(status==0 .and. mass==3d0 .and. gas(1)==4.5d0 .and. &
         maxval(abs(frac-[0.2d0,0.6d0,0.01d0]))<1d-12,'second shared donor uses current state and cap')
    snapshot=gas; frac(2)=ieee_value(0d0,ieee_quiet_nan)
    call agn_deposit_material(gas,fields,1,frac,1d0,[0d0,0d0,0d0],3d0,2d0, &
         2d0,1d0,100d0,dsave,status)
    call check(status==agn_deposit_invalid_source .and. all(gas==snapshot), &
         'invalid composition commits neither hydro nor species')
    gas(8)=-1d0; snapshot=gas
    call agn_withdraw_cell(gas,fields,1,1d0,2d0,mass,vel,frac,status)
    call check(status==agn_deposit_invalid_receiver .and. all(gas==snapshot), &
         'invalid donor composition is rejected without withdrawal')
    call check(agn_contains_donor([0.5d0,0d0,0d0],1d0) .and. &
         .not.agn_contains_donor([-0.5d0,0d0,0d0],1d0),'half-open donor ownership removes cell-face tie')
  end subroutine layouts_and_payload

  subroutine check(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       write(*,'(A)') 'FAIL: '//label
       error stop 1
    endif
    write(*,'(A)') 'PASS: '//label
  end subroutine check
end program agn_feedback_deposition_smoke
