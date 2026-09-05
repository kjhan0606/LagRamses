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
  call accepted_accretion()
  call overlapping_events(.false.)
  call overlapping_events(.true.)
  write(*,'(A)') 'AGN_NATIVE_CELL_COUPLING_SMOKE_OK'
contains
  subroutine accepted_accretion()
    real(dp)::gas(9),initial(9),gross,net,erg,erg2,merged(2),reordered(2),gas1(4),gas2(5)
    integer::status,field(1)
    gas=[8d0,16d0,0d0,0d0,24d0,1.6d0,7d0,4.8d0,0.08d0]; initial=gas
    call agn_accretion_receipt(8d0,8d0,2d0,100d0,0.75d0,0.1d0,5d33,gross,net,erg,status)
    call check(status==0.and.gross==4d0.and.abs(net-3.6d0)<1d-14,'floor clips gross before retained mass and radiation')
    call agn_accrete_scalars(gas,[6,8,9],1,gross,2d0,status)
    call check(status==0.and.maxval(abs(gas([6,8,9])-0.75d0*initial([6,8,9])))<1d-14, &
         'actual accreted fraction removes total metal and all declared elements')
    call check(all(gas(1:5)==initial(1:5)).and.gas(7)==initial(7),'scalar helper leaves hydro and reserved convention to caller')
    ! Next accepted event uses its own efficiency and unit conversion.
    call agn_accretion_receipt(6d0,6d0,2d0,1d0,0.75d0,0.2d0,5d33,gross,net,erg2,status)
    call check(status==0.and.gross==1d0.and.net==0.8d0,'second event contemporaneous efficiency')
    call check(abs((erg+erg2)/(0.6d0*5d33*(2.99792458d10)**2)-1d0)<1d-14, &
         'accepted energies add without latest-efficiency reconstruction')
    call agn_merge_pending([erg,0d0,erg2],[2,1,2],merged,status)
    call check(status==0.and.merged(1)==0d0.and.merged(2)==erg+erg2,'new sink zero and merger sum accepted energy')
    call agn_merge_pending(merged,[2,1],reordered,status)
    call check(status==0.and.all(reordered==merged(2:1:-1)),'reorder preserves pending ownership')
    call agn_merge_pending([erg,erg2],[1,3],merged,status)
    call check(status/=0,'invalid merger mapping rejected before sink commit')
    call agn_accretion_receipt(6d0,8d0,2d0,100d0,0.75d0,0.2d0,5d33,gross,net,erg2,status)
    call check(status==0.and.gross==0d0.and.net==0d0.and.erg2==0d0,'exhausted floor supplies no photons')
    call agn_accretion_receipt(8d0,8d0,2d0,1d0,0.75d0,0d0,5d33,gross,net,erg2,status)
    call check(status==0.and.net==gross.and.erg2==0d0,'zero efficiency is valid zero radiative energy')
    call agn_accretion_receipt(8d0,8d0,2d0,1d0,0.75d0,1d0,5d33,gross,net,erg2,status)
    call check(status/=0.and.gross==0d0,'invalid efficiency rejected before donor mutation')
    gas=initial
    call agn_accrete_scalars(gas,[6,6],1,4d0,2d0,status)
    call check(status/=0.and.all(gas==initial),'invalid accretion scalar map non-mutating')
    gas=initial; gas(8)=-1d-14; initial=gas
    call agn_accrete_scalars(gas,[6,8,9],1,4d0,2d0,status)
    call check(status==2.and.all(gas==initial),'finite negative constituent requests whole-event skip without clamping')
    gas(8)=1d0; gas(6)=9d0; initial=gas
    call agn_accrete_scalars(gas,[6,8,9],1,4d0,2d0,status)
    call check(status==2.and.all(gas==initial),'metal above gas density requests event skip without mutation')
    call agn_accretion_receipt(8d0,0d0,2d0,100d0,0.75d0,0.1d0,5d33,gross,net,erg2,status)
    call check(status==0.and.gross==4d0.and.abs(net-3.6d0)<1d-14, &
         'zero initial reference uses current donor floor instead of emptying cell or aborting')
    call agn_scalar_map(4,4,0,0,[0],field,status,3)
    call check(status==0.and.field(1)==4,'1D metal index follows actual hydro extent')
    gas1=[8d0,16d0,24d0,1.6d0]
    call agn_accrete_scalars(gas1,field,1,4d0,2d0,status,3)
    call check(status==0.and.abs(gas1(4)-1.2d0)<1d-14,'1D scalar withdrawal leaves hydro slots intact')
    call agn_scalar_map(5,5,0,0,[0],field,status,4)
    call check(status==0.and.field(1)==5,'2D metal index follows actual hydro extent')
    gas2=[8d0,16d0,0d0,24d0,1.6d0]
    call agn_accrete_scalars(gas2,field,1,4d0,2d0,status,4)
    call check(status==0.and.abs(gas2(5)-1.2d0)<1d-14,'2D scalar withdrawal leaves hydro slots intact')
  end subroutine accepted_accretion

  subroutine overlapping_events(reverse_order)
    logical,intent(in)::reverse_order
    real(dp)::gas(9,3),initial(9,3),vol(3),rho_ref,pre_mass,norm(2),lw(2),ax(3),pos(3)
    real(dp)::loaded(2),vel(3,2),frac(3,2),drho,dpvec(3),de,dsave,saved,requested(4),budget(4)
    real(dp)::change(9),thermal_density,speed
    integer::k,j,event,jet,status
    vol=[2d0,1d0,3d0]; pos=[-0.5d0,0d0,0.5d0]
    do k=1,3
       gas(:,k)=[2d0,0d0,0d0,0d0,2d0,0.2d0,7d0,1d0,0.02d0]
    enddo
    gas(:,1)=[8d0,0d0,0d0,0d0,8d0,1.6d0,7d0,4.8d0,0.08d0]
    initial=gas; pre_mass=sum(gas(1,:)*vol)
    call agn_withdraw_cell(gas(:,1),[6,8,9],1,2d0,vol(1),loaded(1),vel(:,1),frac(:,1),status)
    call check(status==0,'overlap first donor loading')
    call agn_withdraw_cell(gas(:,1),[6,8,9],1,1d0,vol(1),loaded(2),vel(:,2),frac(:,2),status)
    call check(status==0,'overlap shared donor loading')
    call check(abs(sum(gas(1,:)*vol)-pre_mass)>1d0,'old thermal normalization actually differs after loading')
    norm=0d0
    do k=1,3
       call agn_jet_geometry([pos(k),0d0,0d0],[1d0,0d0,0d0],1d0,lw,ax)
       norm=norm+lw*vol(k)
    enddo
    ! Thermal, jet, saved thermal replay, second shared-donor jet.
    budget=[7d0,4d0,5d0,3d0]; saved=0d0; requested=0d0
    do k=1,3
       rho_ref=gas(1,k)
       if(k==1)rho_ref=rho_ref+sum(loaded)/vol(k)
       call check(abs(rho_ref-initial(1,k))<1d-14,'pre-loading density restored once before cell event loop')
       do j=1,4
          event=j
          if(reverse_order)event=5-j
          drho=0d0; dpvec=0d0
          if(event==1.or.event==3)then
             de=budget(event)*rho_ref/pre_mass
             requested(event)=requested(event)+de*vol(k)
             call agn_deposit_cell(gas(1:5,k),drho,dpvec,de,vol(k),2d0,1d0,1d0,dsave,status)
          else
             jet=event/2
             call agn_jet_geometry([pos(k),0d0,0d0],[1d0,0d0,0d0],1d0,lw,ax)
             speed=sqrt(budget(event)/loaded(jet))
             call agn_jet_delta(loaded(jet),lw,norm,vel(:,jet),ax,speed,drho,dpvec,de)
             thermal_density=0.5d0*budget(event)/sum(vol)
             de=de+thermal_density
             call agn_deposit_material(gas(:,k),[6,8,9],1,frac(:,jet),drho,dpvec,de, &
                  vol(k),2d0,1d0,1d0,dsave,status)
          endif
          call check(status==0,'overlapping native deposition accepted')
          saved=saved+dsave
       enddo
    enddo
    change=matmul(gas-initial,vol)
    call check(maxval(abs(requested([1,3])-budget([1,3])))<1d-12,'thermal and replay requested energies normalize separately')
    call check(abs(change(5)+saved-sum(budget))<1d-11.and.saved>0d0,'overlap gas plus deferred energy closes with active cap')
    call check(maxval(abs(change([1,2,3,4,6,8,9])))<1d-11,'overlap preserves loaded mass momentum and species')
  end subroutine overlapping_events

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
