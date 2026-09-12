module dust_mass_runtime
#ifdef DUST_DYNAMICS
  use dust_phase_state, only: dust_phase_read,dust_phase_kinetic,dust_phase_erode
#endif
  use amr_commons
  use hydro_commons
  use dust_mass_physics
  use dust_composition_material
  use dust_iron_compare, only: iron_compare_curve,iron_compare_temperature,iron_mixture_enthalpy
  use dust_iron_optics, only: fe_radius_cm,fe_density
  use dust_sublimation_physics, only: dust_sublimation_step
  use dust_composition_optics, only: d03_optics_binding
  use snrt_dust_ir, only: snrt_dust_material_temperature
  use stellar_native_units, only: solar_mass_cgs
  use snrt_dust_contract
  use snrt_state, only: snrt_state_get_slot,snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii
  use snrt_thermochemistry, only: snrt_mean_molecular_weight
  use snrt_atomic_cooling, only: atomic_temperature
#ifdef SNRT_CHIMES
  use snrt_chimes_runtime
  use snrt_chimes, only: chimes_ns
#endif
#include "amr_index.h"
  implicit none
contains
  subroutine dust_expansion_level(ilevel,expansion_ratio)
    ! Pressureless grain enthalpy has no Hubble PdV term. Its physical
    ! specific energy stays fixed, whereas scale_v**2 varies as a**(-2).
    ! Update only the in-flight conserved reservoir, at the same global
    ! clock boundary used by update_cosmomag, before the hydro flux update.
    ! Mass carriers and gas total energy must NOT receive this factor.
    integer,intent(in)::ilevel
    real(dp),intent(in)::expansion_ratio
    real(dp)::factor
    integer::ind,i,cell,icpu
    if(.not.dust_mass_enabled)return
    if(.not.ieee_is_finite(expansion_ratio).or.expansion_ratio<=0)then
       write(*,*)'ERROR: invalid dust expansion ratio'
       call clean_stop
       return
    endif
    factor=expansion_ratio**2
    if(.not.ieee_is_finite(factor).or.factor<=0)then
       write(*,*)'ERROR: invalid dust expansion factor'
       call clean_stop
       return
    endif
    do ind=1,twotondim
       do i=1,active(ilevel)%ngrid
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          unew(cell,idust_energy)=unew(cell,idust_energy)*factor
       enddo
       do icpu=1,ncpu
          do i=1,reception(icpu,ilevel)%ngrid
             cell=ICELL_OF(reception(icpu,ilevel)%igrid(i),ind)
             unew(cell,idust_energy)=unew(cell,idust_energy)*factor
          enddo
       enddo
    enddo
  end subroutine

  subroutine dust_injection_specific_energy(value,ierr,grains,metallic_iron)
    real(dp),intent(out)::value ! erg / gram of dust
    integer,intent(out)::ierr
    real(dp),optional,intent(in)::grains(2)
    real(dp),optional,intent(in)::metallic_iron
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
       if(.not.dust_material_domain(snrt_dust_contract_temperature_k(1:n)))return
       composition=1;if(present(grains))composition=grains
       if(dust_iron_enabled().and.present(metallic_iron))then
          call iron_compare_curve([t],composition,metallic_iron,1d0,specific,ierr)
       else
          call dust_composition_curve([t],composition,1d0,specific,ierr)
       endif
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
    real(dp)::iron_mass,enthalpy_hi,solid_pah(2),iron_bins(2),next_iron(2)
    real(dp)::gas_x(11),ions(3)
    real(dp)::material_curve(snrt_dust_contract_max_temperature),td,unew_specific(2),material_weight
    real(dp)::evaporated_bins(4),evaporated_code(4),evaporated_ed,vapor_heat,phase_energy,olivine_latent
    integer::nt,material_interval,grain_bin
#ifdef DUST_DYNAMICS
    real(dp)::phase_mass(ndust_phase),phase_pd(3,ndust_phase),phase_rg,phase_pg(3),phase_ke
    real(dp)::phase_p(3,0:6),fe_phase_p(3,0:2),mixing_heat,phase_row(nvar)
    real(dp),allocatable::momentum_stage(:,:,:)
    integer::phase_b,phase_j,mass_phases
#endif
#ifdef SNRT_CHIMES
    real(dp)::chemical(chimes_ns),elements(11),grain_elements(11)
    real(dp),allocatable::chemical_stage(:,:)
#endif
    include 'mpif.h'
    if(.not.dust_mass_enabled)return
#ifdef SNRT_CHIMES
    if(dust_chimes_enabled())call chimes_prepare_level(ilevel)
#endif
    call units(sl,st,sd,sv,snh,st2)
    olivine_latent=0
    if(dust_silicate_sublimation_enabled())then
       call dust_olivine_phase_reference(olivine_latent,status)
       if(status/=0)call MPI_ABORT(MPI_COMM_WORLD,11,info)
    endif
    nt=snrt_dust_contract_number_temperature
    dt=dtnew(ilevel)*st;grains=0;next_grains=0
    volume=(boxlen/2d0**ilevel)**ndim
    n=active(ilevel)%ngrid*twotondim
    allocate(cells(n),stage(n,11));k=0;bad=0
#ifdef DUST_DYNAMICS
    mass_phases=4
    if(dust_iron_enabled())mass_phases=6
    if(dust_relative_motion)allocate(momentum_stage(3,mass_phases,n))
#endif
#ifdef SNRT_CHIMES
    if(dust_chimes_enabled())allocate(chemical_stage(chimes_ns,n))
#endif
    ! Prepare all local leaf updates before any write. One existing MPI
    ! collective rejects invalid states before the level is committed.
    do ind=1,twotondim
       do i=1,active(ilevel)%ngrid
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          if(son(cell)/=0)cycle
          k=k+1;cells(k)=cell
          rho=uold(cell,1);metal=uold(cell,imetal);dust=uold(cell,idust);ed=uold(cell,idust_energy)
          iron_mass=0;iron_bins=0;next_iron=0
          solid_pah=0
          if(dust_pah_enabled())solid_pah=dust_pah_inventory(uold(cell,idust_pah:idust_pah+dust_pah_nstate()-1))
          if(dust_iron_enabled())then
             iron_bins=uold(cell,idust_iron:idust_iron+1);iron_mass=sum(iron_bins);next_iron=iron_bins
          endif
          if(.not.all(ieee_is_finite(uold(cell,1:nvar))).or.rho<=0)then
             bad=1;cycle
          endif
          eg=uold(cell,ndim+2)-.5d0*sum(uold(cell,2:ndim+1)**2)/rho-magnetic_energy(uold(cell,:))
#ifdef DUST_DYNAMICS
          if(dust_relative_motion)then
             call dust_phase_read(uold(cell,:),phase_mass,phase_pd,phase_rg,phase_pg,phase_ke,status)
             if(status/=0)then
                if(bad==0)write(*,*)'Dust initial phase state rejected: ',cell,status
                bad=1;cycle
             endif
             eg=uold(cell,ndim+2)-phase_ke
             phase_p=0;phase_p(:,0)=phase_pg;phase_p(:,1:mass_phases)=phase_pd(:,1:mass_phases);mixing_heat=0
          endif
#endif
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
#ifdef SNRT_CHIMES
             grain_elements=uold(cell,ichem:ichem+10)
             if(dust_chimes_enabled())then
                call chimes_grain_budget(cell,grain_elements,status)
                if(status/=0)then
                   if(bad==0)write(*,*)'Dust atomic grain budget rejected: ',cell,status
                   bad=1;cycle
                endif
                temp=eg*sd*sv**2/chimes_live_capacity(uold(cell,ichimes:ichimes+chimes_ns-1),sd)
             endif
#endif
             grains=uold(cell,idust_species:idust_species+1)
             if(abs(sum(grains)+iron_mass-dust)>128*epsilon(1d0)*max(dust,sum(grains)+iron_mass,tiny(1d0)))then
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
#ifdef DUST_DYNAMICS
                if(dust_relative_motion)then
                   call dust_phase_erode(bins,shocked_bins,phase_p(:,0),phase_p(:,1:4),status)
                   if(status/=0)then
                      if(bad==0)write(*,*)'Dust shock phase transfer rejected: ',cell,status,bins,shocked_bins
                      bad=1;cycle
                   endif
                   if(dust_iron_enabled().or.dust_pah_enabled())then
                      call dust_size_step_reserved_iron(rho,metal,grain_elements,shocked_bins,iron_mass, &
                           temp,sd,dt,next_bins,status,pah_hc=solid_pah,phase_momentum=phase_p(:,0:4),mixing_heat=mixing_heat)
                   else
                      call dust_size_step(rho,metal,grain_elements,shocked_bins,temp,sd,dt,next_bins,status, &
                           phase_momentum=phase_p(:,0:4),mixing_heat=mixing_heat)
                   endif
                else
#endif
#ifdef SNRT_CHIMES
                if(dust_iron_enabled().or.dust_pah_enabled())then
                   call dust_size_step_reserved_iron(rho,metal,grain_elements,shocked_bins,iron_mass, &
                        temp,sd,dt,next_bins,status,pah_hc=solid_pah)
                else
                   call dust_size_step(rho,metal,grain_elements,shocked_bins,temp,sd,dt,next_bins,status)
                endif
#else
                call dust_size_step(rho,metal,uold(cell,ichem:ichem+10),shocked_bins,temp,sd,dt,next_bins,status)
#endif
#ifdef DUST_DYNAMICS
                endif
#endif
                next_grains=[sum(next_bins(1:2)),sum(next_bins(3:4))]
             else
                call dust_species_step(rho,metal,uold(cell,ichem:ichem+10),grains,temp,sd,dt,next_grains,status)
             endif
             if(status==0.and.dust_iron_enabled().and.dust_fe_kinetics)then
#ifdef SNRT_CHIMES
#ifdef DUST_DYNAMICS
                if(dust_relative_motion)then
                   fe_phase_p(:,0)=phase_p(:,0);fe_phase_p(:,1:2)=phase_p(:,5:6)
                   ! Fe has no molecular carrier in CHIMES. Use total nuclei
                   ! here: sputtering n_H must include hydrogen in molecules,
                   ! unlike the free-atom budget needed for C/S accretion.
                   call dust_fe_kinetics_step(rho,uold(cell,ichem:ichem+10),next_grains,iron_bins,temp,sd,dt, &
                        fe_radius_cm,fe_density(1),next_iron,status,solid_pah,fe_phase_p,mixing_heat)
                   if(status==0)then
                      phase_p(:,0)=fe_phase_p(:,0);phase_p(:,5:6)=fe_phase_p(:,1:2)
                   endif
                else
#endif
                   call dust_fe_kinetics_step(rho,uold(cell,ichem:ichem+10),next_grains,iron_bins,temp,sd,dt, &
                        fe_radius_cm,fe_density(1),next_iron,status,pah_hc=solid_pah)
#ifdef DUST_DYNAMICS
                endif
#endif
#else
                status=1 ! Fe requires the explicit CHIMES depletion path.
#endif
             endif
             next=sum(next_grains)+sum(next_iron)
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
#ifdef DUST_DYNAMICS
          if(dust_relative_motion)then
             phase_row=uold(cell,:);phase_row(idust_bins:idust_bins+3)=next_bins
             if(dust_iron_enabled())phase_row(idust_iron:idust_iron+1)=next_iron
             do phase_b=1,mass_phases
                phase_j=idust_momentum+3*(phase_b-1)
                phase_row(phase_j:phase_j+2)=phase_p(:,phase_b)
             enddo
             ! The rates conserve total E: all lost component KE is already
             ! thermal. Do not add mixing_heat to the total-energy field.
             eg=eg+phase_ke-dust_phase_kinetic(phase_row)
          endif
#endif
          if((dust_iron_enabled().or.dust_pah_enabled().or.dust_relative_motion.or.cosmo).and.dust>0)then
             ! Fe, relative IR and cosmological CMB use analytic mixture enthalpy
             ! (including phase plateaus), unlike the C/silicate U(log T) path.
             call iron_compare_temperature(grains,iron_mass,ed*sv**2,td,status)
             if(status==0)call iron_mixture_enthalpy(td,[next_grains,sum(next_iron)],new_ed,enthalpy_hi,status)
             if(status==0)then
                new_ed=new_ed/sv**2
                if(all(next_grains==grains).and.all(next_iron==iron_bins))new_ed=ed
                q=new_ed-ed
             endif
          else if(dust_material_composition_enabled().and.dust>0)then
             call dust_composition_curve(snrt_dust_contract_temperature_k(1:nt),grains,1d0, &
                  material_curve(1:nt),status)
             if(status==0)call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:nt), &
                  material_curve(1:nt),ed*sv**2/dust,td,status)
             if(status==0)then
                material_interval=1
                do while(material_interval<nt-1)
                   if(td<=snrt_dust_contract_temperature_k(material_interval+1))exit
                   material_interval=material_interval+1
                enddo
                call dust_composition_curve( &
                     snrt_dust_contract_temperature_k(material_interval:material_interval+1), &
                     next_grains,1d0,unew_specific,status)
                material_weight=log(td/snrt_dust_contract_temperature_k(material_interval))/ &
                     log(snrt_dust_contract_temperature_k(material_interval+1)/ &
                     snrt_dust_contract_temperature_k(material_interval))
                material_weight=max(0d0,min(1d0,material_weight))
             endif
             if(status==0)then
                ! Use the same U(log T) interpolation as the inverse above.
                ! Re-evaluating the analytic curve at td changes energy even
                ! when no mass process occurred between table knots.
                new_ed=next*((1-material_weight)*unew_specific(1)+material_weight*unew_specific(2))/sv**2
                if(all(next_grains==grains))new_ed=ed
                q=new_ed-ed
                if(.not.ieee_is_finite(new_ed))status=1
             endif
          else
             call dust_mass_exchange(dust,next,ed,eg,new_ed,q,status)
          endif
          if(status==0.and.dust_sublimation_enabled())then
             ! Ephase=L*(C_total-C_solid), derived from existing conserved
             ! carriers. Apply the SAME reference to growth, sputtering and
             ! SN destruction: releasing latent heat only for one operator
             ! would break the closed-cell phase-energy budget.
             q=q-dust_carbon_latent*(next_grains(1)-grains(1))/sv**2
             q=q-olivine_latent*(next_grains(2)-grains(2))/sv**2
          endif
          if(status==0)then
             if(.not.ieee_is_finite(q).or.q>eg)status=1
          endif
          if(status==0.and.dust_sublimation_enabled().and..not.dust_sublimation_rt_enabled())then
             call dust_sublimation_step(snrt_dust_contract_temperature_k(1:nt), &
                  snrt_dust_contract_ir_background_k,next_bins*sd,new_ed*sd*sv**2,dt, &
                  evaporated_bins,evaporated_ed,vapor_heat,phase_energy,status)
             if(status==0)then
                ! Transfer the represented survival fraction, not a cgs/code
                ! round trip of the absolute mass. Exact zero erosion must
                ! remain an identity; a one-ulp increase is not solid growth.
                evaporated_code=0
                do grain_bin=1,4
                   if(next_bins(grain_bin)<=0)cycle
                   if(next_bins(grain_bin)*sd<=0)then
                      status=1;exit
                   endif
                   if(evaporated_bins(grain_bin)==next_bins(grain_bin)*sd)then
                      evaporated_code(grain_bin)=next_bins(grain_bin)
                   else
                      evaporated_code(grain_bin)=next_bins(grain_bin)* &
                           (evaporated_bins(grain_bin)/(next_bins(grain_bin)*sd))
                   endif
                enddo
                if(status/=0)then
                   bad=1;cycle
                endif
#ifdef DUST_DYNAMICS
                if(dust_relative_motion)then
                   call dust_phase_erode(next_bins,evaporated_code,phase_p(:,0),phase_p(:,1:4),status)
                   if(status/=0)then
                      if(bad==0)write(*,*)'Dust sublimation phase transfer rejected: ',cell,status,next_bins,evaporated_code
                      if(bad==0)write(*,*)'Dust sublimation mass delta/nonfinite momentum: ', &
                           evaporated_code-next_bins,count(.not.ieee_is_finite(phase_p(:,0:4)))
                      bad=1;cycle
                   endif
                endif
#endif
                next_bins=evaporated_code
                next_grains=[sum(next_bins(1:2)),sum(next_bins(3:4))];next=sum(next_grains)
                new_ed=evaporated_ed/(sd*sv**2)
                q=q-vapor_heat/(sd*sv**2)
             endif
          endif
          if(status/=0)then
             if(bad==0)write(*,*)'Dust material exchange rejected: ',dust,next,ed,eg
             bad=1;cycle
          endif
          stage(k,1:3)=[next,new_ed,uold(cell,ndim+2)-q]
#ifdef SNRT_CHIMES
          if(dust_chimes_enabled())then
             call chimes_cell_state(cell,next_grains,chemical_stage(:,k),elements,status,metallic_iron=sum(next_iron))
             if(status/=0)then
                if(bad==0)write(*,*)'Dust final chemical reconciliation rejected: ',cell,status
                bad=1
             endif
          endif
#endif
          if(dust_composition_enabled())stage(k,4:5)=next_grains
          if(dust_two_size_enabled())stage(k,6:9)=next_bins
          if(dust_iron_enabled())stage(k,10:11)=next_iron
#ifdef DUST_DYNAMICS
          if(dust_relative_motion)then
             momentum_stage(:,:,k)=phase_p(:,1:mass_phases)
             phase_row=uold(cell,:);phase_row(idust_bins:idust_bins+3)=next_bins
             if(dust_iron_enabled())phase_row(idust_iron:idust_iron+1)=next_iron
             phase_row(ndim+2)=stage(k,3)
             do phase_b=1,mass_phases
                phase_j=idust_momentum+3*(phase_b-1)
                phase_row(phase_j:phase_j+2)=phase_p(:,phase_b)
             enddo
             eg=phase_row(ndim+2)-dust_phase_kinetic(phase_row)
             if(nener>0)eg=eg-sum(phase_row(inener:inener+nener-1))
             if(.not.ieee_is_finite(eg).or.eg<0)then
                if(bad==0)write(*,*)'Dust final phase energy rejected: ',cell,eg
                bad=1
             endif
          endif
#endif
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
       if(dust_iron_enabled())uold(cell,idust_iron:idust_iron+1)=stage(i,10:11)
#ifdef DUST_DYNAMICS
       if(dust_relative_motion)then
          do phase_b=1,mass_phases
             phase_j=idust_momentum+3*(phase_b-1)
             uold(cell,phase_j:phase_j+2)=momentum_stage(:,phase_b,i)
          enddo
       endif
#endif
       if(dust_sn_shocks)uold(cell,idust_shock:idust_fresh+1)=0d0
#ifdef SNRT_CHIMES
       if(dust_chimes_enabled())uold(cell,ichimes:ichimes+chimes_ns-1)=chemical_stage(:,i)
#endif
    enddo
    deallocate(cells,stage)
    ! Restrict updated leaves into covered coarse cells, including non-RT runs.
    call upload_fine(ilevel)
#ifdef DUST_DYNAMICS
    if(dust_relative_motion)then
       do i=idust_momentum,idust_momentum+3*ndust_phase-1
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
    endif
#endif
#ifdef SNRT_CHIMES
    if(dust_chimes_enabled())then
       do i=ichimes,ichimes+chimes_ns-1
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
    endif
#endif
    call make_virtual_fine_dp(uold(1,idust),ilevel)
    if(dust_iron_enabled())then
       call make_virtual_fine_dp(uold(1,idust_iron),ilevel)
       call make_virtual_fine_dp(uold(1,idust_iron+1),ilevel)
    endif
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
