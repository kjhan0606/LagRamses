! One radiation field shared by bulk C/silicate/Fe and non-equilibrium PAHs.
! Absolute IR convention; no subtraction of an untracked thermal bath.
module dust_pah_mixed
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dust_pah_live_model
  use dust_pah_radiation, only: pah_absorbed_step,pah_charged_absorbed_step
  use dust_pah_hydrogen, only: pah_hydrogen_charged_step,pah_h2_binding
  use dust_iron_radiation, only: iron_radiative_cell
  use dust_iron_optics, only: fe_six_opacity_basis
  use dust_composition_optics, only: d03_opacity_basis,d03_ng,d03_nir
  use dust_mass_physics, only: dust_iron_enabled,dust_fe_max_temperature,dust_pah_nbin,dust_pah_molecule_g, &
       dust_pah_charged,dust_pah_nstate,dust_pah_hydrogenated,dust_pah_charge_size,dust_pah_h2_enabled
  use snrt_dust_contract
  use snrt_dust_ir
  implicit none
  private
  public::pah_mixed_advance
contains
  subroutine pah_mixed_advance(table,direction,weight,neighbors,dx,dt,chat,bins,primary, &
       radiation,population,bulk_energy,temperature,photons,diag,ierr,ghosts,remote,blocked, &
       gas_energy,gas_capacity,conductance,transfer,primary_pah_heat,primary_pah_captures, &
       phase_density,phase_momentum,phase_absorption,phase_scattering,phase_work,gas_electrons,electron_capacity, &
       primary_population,gas_atomic_h,gas_molecular_h2)
    type(dust_ir_table),intent(in)::table
    real(dust_dp),intent(in)::direction(:,:),weight(:),dx,dt,chat,bins(:,:),primary(:,:)
    integer,intent(in)::neighbors(:,:)
    real(dust_dp),intent(inout)::radiation(:,:,:),population(:,:),bulk_energy(:),temperature(:),photons(:,:)
    type(dust_ir_diagnostics),intent(inout)::diag
    integer,intent(out)::ierr
    real(dust_dp),intent(in)::ghosts(:,:,:)
    integer,intent(in)::remote(:,:)
    logical,intent(in)::blocked(:,:)
    real(dust_dp),intent(inout)::gas_energy(:)
    real(dust_dp),intent(in)::gas_capacity(:),conductance(:)
    real(dust_dp),intent(out)::transfer(:)
    real(dust_dp),optional,intent(in)::primary_pah_heat(:,:),primary_pah_captures(:,:)
    real(dust_dp),optional,intent(in)::phase_density(:,:),phase_absorption(:,:,:),phase_scattering(:,:,:)
    real(dust_dp),optional,intent(inout)::phase_momentum(:,:,:),phase_work(:,:)
    real(dust_dp),optional,intent(inout)::gas_electrons(:)
    real(dust_dp),optional,intent(in)::electron_capacity
    real(dust_dp),optional,intent(in)::primary_population(:,:)
    real(dust_dp),optional,intent(inout)::gas_atomic_h(:)
    real(dust_dp),optional,intent(inout)::gas_molecular_h2(:)
    real(dust_dp)::pa(d03_ng,6),ps(d03_ng,6),pg(d03_ng,6),ia(d03_nir,6),isc(d03_nir,6),ig(d03_nir,6)
    real(dust_dp),allocatable::band(:,:,:),sigma(:,:),pf(:,:),irf(:,:),total(:),unit_density(:),capacity(:)
    real(dust_dp),allocatable::eb(:),qb(:),number(:),work_t(:),old_pop(:,:)
    real(dust_dp),allocatable::pfc(:,:,:),irfc(:,:,:),peheat(:),next_ne(:)
    real(dust_dp),allocatable::next_h(:)
    real(dust_dp),allocatable::next_h2(:)
    real(dust_dp)::charge_number(2),alpha(2)
    real(dust_dp)::x,occ,photon,denom,reference,phase(4)
    integer::ng,nt,nc,i,j,k,g,status,nbulk
    ierr=dust_err_config;transfer=0
    call pah_live_prepare(status)
    if(status/=0)return
    ng=snrt_dust_contract_number_ir;nt=snrt_dust_contract_number_temperature;nc=size(bulk_energy)
    if(ng/=d03_nir.or.any(shape(bins)/=[6,nc]))return
    if(any(shape(primary)/=[d03_ng,nc]).or.any(shape(population)/=[dust_pah_nstate(),nc]))return
    if(dust_pah_charged())then
       if(.not.present(gas_electrons).or..not.present(electron_capacity))return
       if(present(phase_density).or.present(primary_pah_heat))return
       if(size(gas_electrons)/=nc)return
       if(any(.not.ieee_is_finite(gas_electrons)).or.any(gas_electrons<0))return
       if(.not.ieee_is_finite(electron_capacity).or.electron_capacity<=0)return
    endif
    if(present(primary_population))then
       if(any(shape(primary_population)/=shape(population)))return
       if(any(.not.ieee_is_finite(primary_population)).or.any(primary_population<0))return
    endif
    if(dust_pah_hydrogenated())then
       if(.not.present(gas_atomic_h))return
       if(size(gas_atomic_h)/=nc)return
       if(any(.not.ieee_is_finite(gas_atomic_h)).or.any(gas_atomic_h<0))return
    endif
    if(dust_pah_h2_enabled())then
       if(.not.present(gas_molecular_h2))return
       if(size(gas_molecular_h2)/=nc)return
       if(any(.not.ieee_is_finite(gas_molecular_h2)).or.any(gas_molecular_h2<0))return
    else if(present(gas_molecular_h2))then
       return
    endif
    nbulk=merge(6,4,dust_iron_enabled())
    if(present(phase_density).neqv.present(phase_momentum))return
    if(present(phase_density).neqv.present(phase_absorption))return
    if(present(phase_density).neqv.present(phase_work))return
    if(present(phase_scattering).and..not.present(phase_density))return
    if(present(phase_density))then
       if(any(shape(phase_density)/=[nbulk+1,nc]))return
       if(any(shape(phase_absorption)/=[nbulk+1,ng,nc]))return
       if(any(primary>0).and..not.present(primary_pah_heat))return
    endif
    if(present(primary_pah_heat).neqv.present(primary_pah_captures))return
    if(present(primary_pah_heat))then
       if(any(shape(primary_pah_heat)/=shape(primary)).or.any(shape(primary_pah_captures)/=shape(primary)))return
       if(any(.not.ieee_is_finite(primary_pah_heat)).or.any(.not.ieee_is_finite(primary_pah_captures)))return
       if(any(primary_pah_heat<0).or.any(primary_pah_captures<0).or.any(primary_pah_heat>primary))return
    endif
    if(size(gas_energy)/=nc.or.size(gas_capacity)/=nc.or.size(conductance)/=nc.or.size(transfer)/=nc)return
    reference=snrt_dust_contract_mass_per_h_g
    pa=0;ps=0;pg=0;ia=0;isc=0;ig=0
    if(dust_iron_enabled())then
       call fe_six_opacity_basis(reference,pa,ps,pg,ia,isc,ig,status)
    else
       call d03_opacity_basis(reference,pa(:,1:4),ps(:,1:4),pg(:,1:4),ia(:,1:4),isc(:,1:4),ig(:,1:4),status)
    endif
    if(status/=0)return
    allocate(band(ng,nt,6),sigma(ng,nc),pf(d03_ng,nc),irf(ng,nc),total(nc),unit_density(nc), &
         capacity(nc),eb(nc),qb(nc),work_t(nc))
    number=sum(population,dim=1);old_pop=population
    allocate(pfc(d03_ng,2,nc),irfc(ng,2,nc),peheat(nc),next_ne(nc))
    pfc=0;irfc=0;peheat=0;next_ne=0
    if(dust_pah_charged())next_ne=gas_electrons
    if(dust_pah_hydrogenated())next_h=gas_atomic_h
    if(dust_pah_h2_enabled())next_h2=gas_molecular_h2
    unit_density=1;capacity=1;eb=bulk_energy;qb=0;work_t=temperature
    do i=1,nc
       total(i)=bulk_energy(i)+gas_energy(i)+dot_product(pah_level,population(:,i))
       charge_number=[sum(population(1:dust_pah_charge_size(),i)),0d0]
       if(dust_pah_charged())charge_number(2)=sum(population(dust_pah_charge_size()+1:,i))
       do g=1,ng
          alpha=charge_number*[pah_ir_sigma(g),pah_ir_ion_sigma(g)]
          sigma(g,i)=dot_product(ia(g,:),bins(:,i))/reference+sum(alpha)
          irf(g,i)=0
          if(sigma(g,i)>0)irfc(g,:,i)=alpha/sigma(g,i)
          irf(g,i)=sum(irfc(g,:,i))
          if(number(i)>0.and.pah_ir_supported(g)==0)then
             if(any(radiation(g,:,i)>0).or.any(ghosts(g,:,:)>0))return
          endif
       enddo
       do g=1,d03_ng
          ! Primary transport used the full-step initial charge mixture.
          ! IR subcycling must not repartition those already captured photons
          ! with a different opacity. IR itself uses each substep's mixture.
          if(present(primary_population))then
             charge_number=[sum(primary_population(1:dust_pah_charge_size(),i)),0d0]
             if(dust_pah_charged())charge_number(2)=sum(primary_population(dust_pah_charge_size()+1:,i))
          endif
          alpha=charge_number*[pah_primary_sigma(g),pah_primary_ion_sigma(g)]
          denom=dot_product(pa(g,:),bins(:,i))/reference+sum(alpha)
          pf(g,i)=0
          if(denom>0)pfc(g,:,i)=alpha/denom
          pf(g,i)=sum(pfc(g,:,i))
          if(number(i)>0.and.snrt_dust_contract_absorption_mean_energy_ev(g)>pah_max_primary_ev.and.primary(g,i)>0)return
       enddo
    enddo
    ! Same Planck quadrature as the existing bulk receiver; units per ref H.
    do j=1,nt
       do g=1,ng
          photon=snrt_dust_contract_ir_energy_ev(g)*1.602176634d-12
          x=photon/(1.380649d-16*snrt_dust_contract_temperature_k(j))
          if(x<1d-3)then
             occ=1/x-.5d0+x/12-x**3/720
          else
             occ=exp(-x)/(1-exp(-x))
          endif
          band(g,j,:)=8*acos(-1d0)*photon**3/(6.62607015d-27**3*2.99792458d10**2)*ia(g,:)* &
               snrt_dust_contract_ir_weight_ev(g)*1.602176634d-12*occ
       enddo
    enddo
    if(present(phase_density))then
    call snrt_dust_ir_advance(table,direction,weight,neighbors,dx,dt,chat,unit_density,sum(primary,dim=1)/dt, &
         radiation,work_t,photons,diag,ierr,1d-9,256,total,capacity,ghosts,remote,blocked, &
         thin_reabsorption=.true.,population=population,cell_absorption=sigma,phase_density=phase_density, &
         phase_momentum=phase_momentum,phase_absorption=phase_absorption,phase_scattering=phase_scattering, &
         phase_work=phase_work,moving_material_dispatch=moving_material)
    else
    call snrt_dust_ir_advance(table,direction,weight,neighbors,dx,dt,chat,unit_density,sum(primary,dim=1)/dt, &
         radiation,work_t,photons,diag,ierr,1d-9,256,total,capacity,ghosts,remote,blocked, &
         thin_reabsorption=.true.,population=population,population_dispatch=material,cell_absorption=sigma)
    endif
    if(ierr/=0)return
    bulk_energy=eb;gas_energy=gas_energy+peheat-qb;transfer=qb-peheat;temperature=work_t
    if(dust_pah_charged())gas_electrons=next_ne
    if(dust_pah_hydrogenated())gas_atomic_h=next_h
    if(dust_pah_h2_enabled())gas_molecular_h2=next_h2
  contains
    subroutine moving_material(ir_heat,ir_captured,step_dt,old,next,phase_rate,next_energy,next_t,status)
      real(dust_dp),intent(in)::ir_heat(:,:,:),ir_captured(:,:,:),step_dt,old(:,:)
      real(dust_dp),intent(out)::next(:,:),phase_rate(:,:,:),next_energy(:),next_t(:)
      integer,intent(out)::status
      real(dust_dp)::pah_rate(ng),bulk_rate(ng),pu,pt,over,heating,phase(4),bulk_alpha(nbulk),den
      real(dust_dp)::pah_heat(d03_ng),pah_count(d03_ng+ng)
      integer::cell,g
      next=old;phase_rate=0;next_energy=0;next_t=0;status=dust_err_state
      do cell=1,nc
         pu=0;pt=0;pah_rate=0;pah_heat=0;pah_count=0
         if(present(primary_pah_heat))then
            pah_heat=primary_pah_heat(:,cell)
            pah_count(1:d03_ng)=primary_pah_captures(:,cell)
         endif
         pah_count(d03_ng+1:)=ir_captured(nbulk+1,:,cell)
         call pah_absorbed_step(pah_live_model, &
              [snrt_dust_contract_absorption_mean_energy_ev(1:d03_ng),snrt_dust_contract_ir_energy_ev(1:ng)], &
              [pah_heat,ir_heat(nbulk+1,:,cell)],old(:,cell),step_dt, &
              next(:,cell),pah_rate,pu,pt,over,status,captured_photons=pah_count)
         if(status/=0)return
         heating=(sum(primary(:,cell)-pah_heat)+sum(ir_heat(1:nbulk,:,cell)))/step_dt
         ! Both material reservoirs start from their immutable substep state.
         call iron_radiative_cell(snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_temperature_k(1), &
              bins(:,cell),bulk_energy(cell),heating,step_dt,reference,band,gas_energy(cell),gas_capacity(cell), &
              conductance(cell),eb(cell),next_t(cell),phase,bulk_rate,qb(cell),status,absolute_emission=.true., &
              photon_ev=snrt_dust_contract_ir_energy_ev(1:ng))
         if(status/=0)return
         if(dust_iron_enabled().and.next_t(cell)>dust_fe_max_temperature)then
            status=dust_err_range;return
         endif
         do g=1,ng
            bulk_alpha=ia(g,1:nbulk)*bins(1:nbulk,cell)/reference;den=sum(bulk_alpha)
            if(den>0)then
               phase_rate(1:nbulk,g,cell)=bulk_rate(g)*(bulk_alpha/den)
            else if(bulk_rate(g)>0)then
               status=dust_err_state;return
            endif
         enddo
         phase_rate(nbulk+1,:,cell)=pah_rate
         next_energy(cell)=eb(cell)+pu+gas_energy(cell)-qb(cell)
      enddo
      status=dust_ok
    end subroutine moving_material

    subroutine material(ir_absorbed,step_dt,old,next,rate,next_energy,next_t,status)
      real(dust_dp),intent(in)::ir_absorbed(:,:),step_dt,old(:,:)
      real(dust_dp),intent(out)::next(:,:),rate(:,:),next_energy(:),next_t(:)
      integer,intent(out)::status
      real(dust_dp)::pah_rate(ng),bulk_rate(ng),pu,pt,over,heating,phase(4)
      real(dust_dp)::pah_heat(d03_ng),pah_count(d03_ng+ng)
      real(dust_dp)::charge_old(dust_pah_nbin,2),charge_next(dust_pah_nbin,2),captures(d03_ng+ng,2)
      real(dust_dp)::h_old(dust_pah_nbin,0:13,2),h_next(dust_pah_nbin,0:13,2)
      real(dust_dp)::energies(d03_ng+ng),ip_change,cv,eg
      integer::q
      integer::cell
      next=old;rate=0;next_energy=0;next_t=0;status=1
      do cell=1,nc
         pu=0;pt=0;pah_rate=0
         pah_heat=primary(:,cell)*pf(:,cell)
         cv=gas_capacity(cell);eg=gas_energy(cell)
         if(dust_pah_charged())then
            energies=[snrt_dust_contract_absorption_mean_energy_ev(1:d03_ng), &
                 snrt_dust_contract_ir_energy_ev(1:ng)]
            do q=1,2
               captures(:,q)=[primary(:,cell)*pfc(:,q,cell),ir_absorbed(:,cell)*irfc(:,q,cell)]/ &
                    (energies*1.602176634d-12)
            enddo
            peheat(cell)=0;ip_change=0
            if(dust_pah_hydrogenated())then
               h_old=reshape(old(:,cell),[dust_pah_nbin,14,2]);h_next=h_old
               next_ne(cell)=gas_electrons(cell);next_h(cell)=gas_atomic_h(cell)
               if(dust_pah_h2_enabled())then
                  next_h2(cell)=gas_molecular_h2(cell)
                  call pah_hydrogen_charged_step(pah_h_models,energies,captures,h_old,step_dt,eg/cv, &
                       h_next,next_ne(cell),next_h(cell),peheat(cell),pah_rate,status,gas_h2=next_h2(cell))
                  cv=cv+electron_capacity*(next_h2(cell)-gas_molecular_h2(cell))
               else
                  call pah_hydrogen_charged_step(pah_h_models,energies,captures,h_old,step_dt,eg/cv, &
                       h_next,next_ne(cell),next_h(cell),peheat(cell),pah_rate,status)
               endif
               if(status/=0)return
               next(:,cell)=reshape(h_next,[dust_pah_nstate()])
               cv=cv+electron_capacity*(next_h(cell)-gas_atomic_h(cell))
            else
            charge_old=reshape(old(:,cell),[dust_pah_nbin,2]);charge_next=charge_old
            call pah_charged_absorbed_step(pah_charge_models,energies,captures,charge_old,gas_electrons(cell), &
                 eg/cv,step_dt,charge_next,next_ne(cell),peheat(cell),ip_change,pah_rate,status)
            if(status/=0)return
            next(:,cell)=reshape(charge_next,[dust_pah_nstate()])
            endif
            pu=dot_product(pah_level,next(:,cell))
            eg=eg+peheat(cell);cv=cv+electron_capacity*(next_ne(cell)-gas_electrons(cell))
            if(eg<=0.or.cv<=0)then
               status=dust_err_state;return
            endif
         else if(present(primary_pah_heat))then
            pah_heat=primary_pah_heat(:,cell)
            pah_count(1:d03_ng)=primary_pah_captures(:,cell)
            pah_count(d03_ng+1:)=ir_absorbed(:,cell)*irf(:,cell)/ &
                 (snrt_dust_contract_ir_energy_ev(1:ng)*1.602176634d-12)
            call pah_absorbed_step(pah_live_model, &
                 [snrt_dust_contract_absorption_mean_energy_ev(1:d03_ng),snrt_dust_contract_ir_energy_ev(1:ng)], &
                 [pah_heat,ir_absorbed(:,cell)*irf(:,cell)],old(:,cell),step_dt, &
                 next(:,cell),pah_rate,pu,pt,over,status,captured_photons=pah_count)
         else
         call pah_absorbed_step(pah_live_model, &
              [snrt_dust_contract_absorption_mean_energy_ev(1:d03_ng),snrt_dust_contract_ir_energy_ev(1:ng)], &
              [primary(:,cell)*pf(:,cell),ir_absorbed(:,cell)*irf(:,cell)],old(:,cell),step_dt, &
              next(:,cell),pah_rate,pu,pt,over,status)
         endif
         if(status/=0)return
         heating=(sum(primary(:,cell)-pah_heat)+sum(ir_absorbed(:,cell)*(1-irf(:,cell))))/step_dt
         call iron_radiative_cell(snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_temperature_k(1), &
              bins(:,cell),bulk_energy(cell),heating,step_dt,reference,band,eg,cv, &
              conductance(cell),eb(cell),next_t(cell),phase,bulk_rate,qb(cell),status,absolute_emission=.true., &
              photon_ev=snrt_dust_contract_ir_energy_ev(1:ng))
         if(status/=0)then
            write(*,'(A,2I6,*(ES23.15,1X))')' PAH mixed bulk rejected cell/status, bins,U,H,dt,Eg,Cv,K: ', &
                 cell,status,bins(:,cell),bulk_energy(cell),heating,step_dt,gas_energy(cell),gas_capacity(cell),conductance(cell)
            return
         endif
         if(dust_iron_enabled().and.next_t(cell)>dust_fe_max_temperature)then
            status=dust_err_range;return
         endif
         rate(:,cell)=bulk_rate+pah_rate
         next_energy(cell)=eb(cell)+pu+eg-qb(cell)
         ! Reference chemical energy to the INITIAL H2 in this substep.
         ! Dissociating captured H2 costs D0; this is not gas thermal heat
         ! and not an additional advected variable or a negative energy bath.
         if(dust_pah_h2_enabled())next_energy(cell)=next_energy(cell)+ &
              (gas_molecular_h2(cell)-next_h2(cell))*pah_h2_binding
      enddo
      status=0
    end subroutine
  end subroutine
end module
