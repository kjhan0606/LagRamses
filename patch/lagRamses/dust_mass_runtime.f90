module dust_mass_runtime
  use amr_commons
  use hydro_commons
  use dust_mass_physics
  use snrt_dust_contract
  use snrt_state, only: snrt_state_get_slot,snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii
  use snrt_thermochemistry, only: snrt_mean_molecular_weight
#include "amr_index.h"
  implicit none
contains
  subroutine dust_injection_specific_energy(value,ierr)
    real(dp),intent(out)::value ! erg / gram of dust
    integer,intent(out)::ierr
    integer::k,n
    real(dp)::w,t
    value=0;ierr=1;n=snrt_dust_contract_number_temperature;t=dust_injection_temperature
    if(.not.snrt_dust_contract_loaded.or.snrt_dust_contract_version/=4.or.n<2)return
    if(t<snrt_dust_contract_ir_background_k.or.t<snrt_dust_contract_temperature_k(1).or. &
       t>snrt_dust_contract_temperature_k(n).or.snrt_dust_contract_mass_per_h_g<=0)return
    k=1
    do while(k<n-1)
       if(t<=snrt_dust_contract_temperature_k(k+1))exit
       k=k+1
    enddo
    w=log(t/snrt_dust_contract_temperature_k(k))/ &
         log(snrt_dust_contract_temperature_k(k+1)/snrt_dust_contract_temperature_k(k))
    value=((1-w)*snrt_dust_contract_internal_energy_per_h_erg(k)+ &
         w*snrt_dust_contract_internal_energy_per_h_erg(k+1))/snrt_dust_contract_mass_per_h_g
    if(.not.ieee_is_finite(value).or.value<=0)return
    ierr=0
  end subroutine

  subroutine dust_mass_advance_level(ilevel)
    integer,intent(in)::ilevel
    integer::i,ind,cell,n,k,status,bad,all_bad,info,slot
    integer,allocatable::cells(:)
    real(dp),allocatable::stage(:,:)
    real(dp)::sl,st,sd,sv,snh,st2,rho,metal,dust,eg,ed,temp,a,b,next,q,new_ed,dt,mu
    include 'mpif.h'
    if(.not.dust_mass_enabled)return
    call units(sl,st,sd,sv,snh,st2)
    dt=dtnew(ilevel)*st
    n=active(ilevel)%ngrid*twotondim
    allocate(cells(n),stage(n,3));k=0;bad=0
    ! Prepare all local leaf updates before any write. One existing MPI
    ! collective rejects invalid states before the level is committed.
    do ind=1,twotondim
       do i=1,active(ilevel)%ngrid
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          if(son(cell)/=0)cycle
          k=k+1;cells(k)=cell
          rho=uold(cell,1);metal=uold(cell,imetal);dust=uold(cell,idust);ed=uold(cell,idust_energy)
          if(.not.all(ieee_is_finite(uold(cell,1:nvar))).or.rho<=0)then
             bad=1;cycle
          endif
          eg=uold(cell,ndim+2)-.5d0*sum(uold(cell,2:ndim+1)**2)/rho
#if NENER>0
          eg=eg-sum(uold(cell,inener:inener+NENER-1))
#endif
          mu=snrt_mean_molecular_weight(0d0,0d0,0d0)
          slot=snrt_state_get_slot(cell)
          if(slot>0)mu=snrt_mean_molecular_weight(snrt_hydrogen_ii(slot),snrt_helium_ii(slot),snrt_helium_iii(slot))
          temp=eg/rho*(gamma-1)*st2*mu
          call dust_mass_rates(rho*sd,metal*sd,temp,a,b,status)
          if(status/=0)then
             bad=1;cycle
          endif
          call dust_mass_step(metal,dust,dt,a,b,next,status)
          if(status/=0)then
             bad=1;cycle
          endif
          call dust_mass_exchange(dust,next,ed,eg,new_ed,q,status)
          if(status/=0)then
             bad=1;cycle
          endif
          stage(k,:)=[next,new_ed,uold(cell,ndim+2)-q]
       enddo
    enddo
    call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(all_bad/=0.or.info/=0)then
       if(myid==1)write(*,*)'ERROR: dust mass evolution rejected invalid mass/energy/rate state'
       call MPI_ABORT(MPI_COMM_WORLD,11,info)
    endif
    do i=1,k
       cell=cells(i)
       uold(cell,idust)=stage(i,1);uold(cell,idust_energy)=stage(i,2);uold(cell,ndim+2)=stage(i,3)
    enddo
    deallocate(cells,stage)
    ! Restrict updated leaves into covered coarse cells, including non-RT runs.
    call upload_fine(ilevel)
    call make_virtual_fine_dp(uold(1,idust),ilevel)
    call make_virtual_fine_dp(uold(1,idust_energy),ilevel)
    call make_virtual_fine_dp(uold(1,ndim+2),ilevel)
  end subroutine
end module
