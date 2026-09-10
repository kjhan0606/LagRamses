! Real neutral-PAH + six bulk components through the moving IR callback.
! Requires existing SNRT_DUST_CONTRACT and SNRT_PAH_NEUTRAL_TABLE inputs.
! Native bounded-model test only; no RAMSES evolution or production claim.
program snrt_moving_pah_smoke
  use dust_pah_mixed
  use dust_pah_live_model
  use dust_iron_material,only:iron_mixture_enthalpy
  use dust_iron_optics,only:fe_nir,fe_six_opacity_basis
  use dust_mass_physics,only:dust_mass_enabled,dust_iron_model,dust_pah_nbin,dust_pah_molecule_g
  use snrt_dust_contract
  use snrt_dust_ir
  implicit none
  integer,parameter::nd=6,nb=7,nc=1
  real(dust_dp),parameter::c=2.99792458d10,ev=1.602176634d-12
  type(dust_ir_table)::tab
  type(dust_ir_diagnostics)::diag
  real(dust_dp)::rays(3,nd),weight(nd),field(fe_nir,nd,nc),old_field(fe_nir,nd,nc),ph(fe_nir,nc)
  real(dust_dp)::pop(dust_pah_nbin,nc),old_pop(dust_pah_nbin,nc),bins(6,nc),primary(9,nc)
  real(dust_dp)::ed(nc),old_ed(nc),td(nc),gas(nc),cv(nc),conductance(nc),transfer(nc)
  real(dust_dp)::rho(nb,nc),p(3,nb,nc),p0(3,nb,nc),work(nb,nc),alpha(nb,fe_nir,nc),sca(nb,fe_nir,nc)
  real(dust_dp)::pa(9,6),ps(9,6),pg(9,6),ia(fe_nir,6),isc(fe_nir,6),ig(fe_nir,6)
  real(dust_dp)::pah_heat(9,nc),captures(9,nc),ghost(fe_nir,nd,0),hi,reference,den,delta,kin,scale,dp(3)
  integer::neighbors(6,nc),remote(6,nc),status,b,g,d
  logical::blocked(6,nc)
  character(len=2048)::path
  call get_environment_variable('SNRT_DUST_CONTRACT',path)
  call snrt_dust_contract_load(trim(path),status)
  call require(status==0,'contract')
  dust_mass_enabled=.true.;dust_iron_model='fe_electric_compare_v1'
  call pah_live_prepare(status);call require(status==0,'PAH preparation')
  reference=snrt_dust_contract_mass_per_h_g
  call fe_six_opacity_basis(reference,pa,ps,pg,ia,isc,ig,status)
  call require(status==0,'physical phase optics')
  call snrt_dust_ir_initialize(tab,snrt_dust_contract_ir_energy_ev(1:fe_nir), &
       snrt_dust_contract_ir_weight_ev(1:fe_nir),snrt_dust_contract_ir_absorption_per_h_cm2(1:fe_nir), &
       snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature),10d0,status)
  call require(status==0,'IR table')
  bins=1d-18;pop=0;pop(1,:)=1d2;old_pop=pop
  call iron_mixture_enthalpy(20d0,[2d-18,2d-18,2d-18],ed(1),hi,status)
  call require(status==0,'material energy');old_ed=ed;td=20
  rho(1:6,:)=bins;rho(7,:)=sum(pop,dim=1)*dust_pah_molecule_g
  alpha=0;sca=0;p=0
  do b=1,6
     alpha(b,:,1)=ia(:,b)*bins(b,1)/reference
     sca(b,:,1)=(isc(:,b)-ig(:,b))*bins(b,1)/reference
     p(:,b,1)=rho(b,1)*c*[.0002d0*b,-.0001d0,0d0]
  enddo
  alpha(7,:,1)=pah_ir_sigma*sum(pop(:,1));p(2,7,1)=rho(7,1)*c*.002d0;p0=p
  primary=0;primary(2,1)=1d-12;pah_heat=0;captures=0
  den=dot_product(pa(2,:),bins(:,1))/reference+sum(pop(:,1))*pah_primary_sigma(2)
  pah_heat(2,1)=primary(2,1)*sum(pop(:,1))*pah_primary_sigma(2)/den
  captures(2,1)=pah_heat(2,1)*1.001d0/(snrt_dust_contract_absorption_mean_energy_ev(2)*ev)
  rays=0;weight=1d0/nd
  do d=1,3
     rays(d,2*d-1)=1;rays(d,2*d)=-1
  enddo
  field=0;ph=0
  g=minloc(abs(snrt_dust_contract_ir_energy_ev(1:fe_nir)-.03d0),dim=1)
  ! Resolve the impulse against finite-precision absolute phase momenta.
  ! The very weak 2e-11 fixture had dp below their ulp, so it could not
  ! measure mechanical recoil even though the material solve converged.
  field(g,1,1)=2d-6;old_field=field
  gas=1d-12;cv=1d-14;conductance=0;transfer=0
  neighbors=0;remote=0;blocked=.true.;work=-123
  call pah_mixed_advance(tab,rays,weight,neighbors,1d12,1d0,3d8,bins,primary,field,pop,ed,td,ph,diag,status, &
       ghost,remote,blocked,gas,cv,conductance,transfer,primary_pah_heat=pah_heat,primary_pah_captures=captures, &
       phase_density=rho,phase_momentum=p,phase_absorption=alpha,phase_scattering=sca,phase_work=work)
  call require(status==0,'moving mixed material iteration')
  call require(sum(pop(2:,:))>0.and.sum(ph)>0,'actual PAH excitation and emission')
  call require(abs(sum(pop)-sum(old_pop))<1d-11*sum(old_pop),'PAH molecule conservation')
  kin=0;dp=0
  do b=1,nb
     kin=kin+dot_product(.5d0*(p(:,b,1)+p0(:,b,1))/rho(b,1),p(:,b,1)-p0(:,b,1))
     dp=dp+p(:,b,1)-p0(:,b,1)
  enddo
  do d=1,nd
     dp=dp+sum(field(:,d,1)-old_field(:,d,1))*weight(d)*rays(:,d)/c
  enddo
  delta=sum(ed-old_ed)+sum(field-old_field)/nd+dot_product(pah_level,pop(:,1)-old_pop(:,1))+kin-sum(primary)
  scale=sum(primary)+sum(old_field)/nd
  write(*,'(A,4ES16.7)')'mixed diagnostic energy/momentum/dK/reportedWork: ',delta,maxval(abs(dp))*c,kin,sum(work)
  call require(abs(delta)<5d-9*scale,'radiation material PAH kinetic closure')
  call require(maxval(abs(dp))*c<5d-9*scale,'combined phase momentum closure')
  call require(abs(kin-sum(work))<5d-9*scale,'reported mechanical work')
  call require(any(work/=0),'nonzero physical IR work')
  write(*,'(A,3ES15.6)')'MOVING_PAH_SMOKE_PASS energy/momentum/work: ',abs(delta)/scale,maxval(abs(dp))*c/scale,sum(work)
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       write(*,'(A,A,A,I0)')'FAIL: ',label,' status=',status
       error stop 1
    endif
  end subroutine
end program
