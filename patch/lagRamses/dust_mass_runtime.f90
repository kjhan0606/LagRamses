module dust_mass_runtime
  use amr_commons
  use hydro_commons
  use dust_mass_physics
  use dust_composition_material
  use dust_composition_optics, only: d03_optics_binding
  use snrt_dust_ir, only: snrt_dust_material_temperature
  use stellar_native_units, only: solar_mass_cgs
  use snrt_dust_contract
  use snrt_state, only: snrt_state_get_slot,snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii
  use snrt_thermochemistry, only: snrt_mean_molecular_weight
  use snrt_atomic_cooling, only: atomic_temperature
#include "amr_index.h"
  implicit none
contains
  subroutine dust_injection_specific_energy(value,ierr,grains)
    real(dp),intent(out)::value ! erg / gram of dust
    integer,intent(out)::ierr
    real(dp),optional,intent(in)::grains(2)
    integer::k,n
    real(dp)::w,t
    real(dp)::composition(2),specific(1)
    value=0;ierr=1;n=snrt_dust_contract_number_temperature;t=dust_injection_temperature
    if(.not.snrt_dust_contract_loaded.or.snrt_dust_contract_version/=4.or.n<2)return
    if(dust_optics_enabled())then
       if(.not.snrt_dust_contract_scattering_enabled)return
       if(.not.d03_optics_binding( &
            snrt_dust_contract_group_edges_ev(1:snrt_dust_contract_number_groups+1), &
            snrt_dust_contract_absorption_mean_energy_ev(1:snrt_dust_contract_number_groups), &
            snrt_dust_contract_ir_energy_ev(1:snrt_dust_contract_number_ir), &
            dust_size_radius_cm,dust_size_density))return
    endif
    if(t<snrt_dust_contract_ir_background_k.or.t<snrt_dust_contract_temperature_k(1).or. &
       t>snrt_dust_contract_temperature_k(n).or.snrt_dust_contract_mass_per_h_g<=0)return
    if(dust_material_composition_enabled())then
       if(any(snrt_dust_contract_temperature_k(1:n)<dl01_t(1)).or. &
            any(snrt_dust_contract_temperature_k(1:n)>dl01_t(dl01_n)))return
       composition=1;if(present(grains))composition=grains
       call dust_composition_curve([t],composition,1d0,specific,ierr)
       if(ierr==0)value=specific(1)
       return
    endif
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
    real(dp)::grains(2),next_grains(2),bins(4),next_bins(4),bin_totals(2)
    real(dp)::shocked_bins(4),volume,cell_mass,sn_energy
    real(dp)::gas_x(11),ions(3)
    real(dp)::material_curve(snrt_dust_contract_max_temperature),td,unew_specific(1)
    integer::nt
    include 'mpif.h'
    if(.not.dust_mass_enabled)return
    call units(sl,st,sd,sv,snh,st2)
    nt=snrt_dust_contract_number_temperature
    dt=dtnew(ilevel)*st;grains=0;next_grains=0
    volume=(boxlen/2d0**ilevel)**ndim
    n=active(ilevel)%ngrid*twotondim
    allocate(cells(n),stage(n,9));k=0;bad=0
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
          if(dust_atomic_cooling_enabled())then
             call dust_gas_elements(uold(cell,ichem:ichem+10),uold(cell,idust_species:idust_species+1),gas_x,status)
             if(status/=0)then
                bad=1;cycle
             endif
             ions=0
             if(slot>0)ions=[snrt_hydrogen_ii(slot),snrt_helium_ii(slot),snrt_helium_iii(slot)]
             temp=atomic_temperature(rho*sd,gas_x/rho,ions,gamma,eg*sd*sv**2)
          endif
          if(dust_composition_enabled())then
             grains=uold(cell,idust_species:idust_species+1)
             if(abs(sum(grains)-dust)>128*epsilon(1d0)*max(dust,sum(grains),tiny(1d0)))then
                if(bad==0)write(*,*)'Dust aggregate/species mismatch: ',dust,grains,sum(grains)-dust
                bad=1;cycle
             endif
             if(dust_two_size_enabled())then
                bins=uold(cell,idust_bins:idust_bins+3)
                bin_totals=[sum(bins(1:2)),sum(bins(3:4))]
                if(any(abs(bin_totals-grains)>128*epsilon(1d0)*max(bin_totals,grains,tiny(1d0))))then
                   if(bad==0)write(*,*)'Dust species/size mismatch: ',grains,bins
                   bad=1;cycle
                endif
                shocked_bins=bins
                if(dust_sn_shocks)then
                   cell_mass=rho*volume*sd*sl**3/solar_mass_cgs
                   sn_energy=uold(cell,idust_shock)*volume*sd*sl**3*sv**2
                   call dust_shock_step(cell_mass,sn_energy,bins, &
                        dust_injection_bins(uold(cell,idust_fresh:idust_fresh+1)),shocked_bins,status)
                   if(status/=0)then
                      if(bad==0)write(*,*)'Dust SN/fresh transaction rejected: ',cell,sn_energy,bins
                      bad=1;cycle
                   endif
                endif
                call dust_size_step(rho,metal,uold(cell,ichem:ichem+10),shocked_bins,temp,sd,dt,next_bins,status)
                next_grains=[sum(next_bins(1:2)),sum(next_bins(3:4))]
             else
                call dust_species_step(rho,metal,uold(cell,ichem:ichem+10),grains,temp,sd,dt,next_grains,status)
             endif
             next=sum(next_grains)
          else
             call dust_mass_rates(rho*sd,metal*sd,temp,a,b,status)
             if(status/=0)then
                bad=1;cycle
             endif
             call dust_mass_step(metal,dust,dt,a,b,next,status)
          endif
          if(status/=0)then
             if(bad==0)write(*,*)'Dust rate/mass rejected: ',cell,rho,metal,dust,temp,grains,next_grains
             bad=1;cycle
          endif
          if(dust_material_composition_enabled().and.dust>0)then
             call dust_composition_curve(snrt_dust_contract_temperature_k(1:nt),grains,1d0, &
                  material_curve(1:nt),status)
             if(status==0)call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:nt), &
                  material_curve(1:nt),ed*sv**2/dust,td,status)
             if(status==0)call dust_composition_curve([td],next_grains,1d0,unew_specific,status)
             if(status==0)then
                new_ed=next*unew_specific(1)/sv**2;q=new_ed-ed
                if(.not.ieee_is_finite(new_ed).or.q>eg)status=1
             endif
          else
             call dust_mass_exchange(dust,next,ed,eg,new_ed,q,status)
          endif
          if(status/=0)then
             if(bad==0)write(*,*)'Dust material exchange rejected: ',dust,next,ed,eg
             bad=1;cycle
          endif
          stage(k,1:3)=[next,new_ed,uold(cell,ndim+2)-q]
          if(dust_composition_enabled())stage(k,4:5)=next_grains
          if(dust_two_size_enabled())stage(k,6:9)=next_bins
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
       if(dust_composition_enabled())uold(cell,idust_species:idust_species+1)=stage(i,4:5)
       if(dust_two_size_enabled())uold(cell,idust_bins:idust_bins+3)=stage(i,6:9)
       if(dust_sn_shocks)uold(cell,idust_shock:idust_fresh+1)=0d0
    enddo
    deallocate(cells,stage)
    ! Restrict updated leaves into covered coarse cells, including non-RT runs.
    call upload_fine(ilevel)
    call make_virtual_fine_dp(uold(1,idust),ilevel)
    call make_virtual_fine_dp(uold(1,idust_energy),ilevel)
    call make_virtual_fine_dp(uold(1,ndim+2),ilevel)
    if(dust_composition_enabled())then
       call make_virtual_fine_dp(uold(1,idust_species),ilevel)
       call make_virtual_fine_dp(uold(1,idust_species+1),ilevel)
    endif
    if(dust_two_size_enabled())then
       do i=idust_bins,idust_bins+3
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
    endif
    if(dust_sn_shocks)then
       do i=idust_shock,idust_fresh+1
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
    endif
  end subroutine
end module
