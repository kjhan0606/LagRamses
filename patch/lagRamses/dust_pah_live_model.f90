! Shared definition for neutral, fixed-H charged and H-state PAH comparisons.
! Source, transport, chemistry and restart use the same state-dependent mass,
! excitation/binding energy and charge-resolved optical projections. H loss
! is optional; carbon-skeleton destruction is not represented.
module dust_pah_live_model
  use dust_mass_physics, only: dust_pah_nbin,dust_pah_nstate,dust_pah_charged, &
       dust_pah_molecule_g,dust_pah_hc,dust_injection_temperature,dust_pah_hydrogenated, &
       dust_pah_charge_size,dust_pah_state_mass,dust_pah_h2_enabled
  use dust_pah_radiation
  use dust_pah_hydrogen
  use dust_stochastic_physics, only: dust_pah_modes,dust_vibrational_curve
  use snrt_dust_contract
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type(pah_radiative_model),public,save::pah_live_model
  type(pah_radiative_model),public,save::pah_charge_models(2)
  type(pah_hydrogen_model),public,save::pah_h_models(2)
  real(real64),allocatable,public,save::pah_level(:),pah_injection(:)
  real(real64),public,save::pah_injection_specific_u,pah_max_primary_ev=4d0
  real(real64),allocatable,public,save::pah_ir_ion_sigma(:),pah_primary_ion_sigma(:)
  real(real64),allocatable,public,save::pah_ir_sigma(:),pah_primary_sigma(:),pah_ir_supported(:)
  logical,save::ready=.false.
  logical,save::charged=.false.
  logical,save::hydrogenated=.false.
  logical,save::molecular_capture=.false.
  public::pah_live_prepare,pah_live_identity,pah_mass_excitation,pah_primary_alpha
contains
  function pah_primary_alpha(population) result(alpha)
    real(real64),intent(in)::population(:)
    real(real64)::alpha(size(pah_primary_sigma))
    alpha=sum(population(1:dust_pah_charge_size()))*pah_primary_sigma
    if(charged)alpha=alpha+sum(population(dust_pah_charge_size()+1:))*pah_primary_ion_sigma
  end function
  real(real64) function pah_mass_excitation(mass) result(energy)
    ! Input g/cm3 per state; output erg/cm3. Same for g -> erg.
    real(real64),intent(in)::mass(:)
    integer::k
    if(.not.hydrogenated)then
       energy=dot_product(mass,pah_level)/dust_pah_molecule_g
    else
       energy=0
       do k=1,size(mass)
          energy=energy+mass(k)*pah_level(k)/dust_pah_state_mass(k)
       enddo
    endif
  end function

  subroutine pah_live_prepare(ierr)
    integer,intent(out)::ierr
    real(real64),allocatable::e(:),w(:),cs(:),ei(:),wi(:),ci(:)
    real(real64)::x,f,modes(102),u(1),cv(1)
    real(real64)::base(dust_pah_nbin),injection(dust_pah_nbin),bond(0:13),attach(0:13)
    integer::i,k,n,nt,status,q,nh,offset
    character(len=2048)::path
    ierr=0
    if(ready)then
       if(charged.neqv.dust_pah_charged())ierr=1
       if(hydrogenated.neqv.dust_pah_hydrogenated())ierr=1
       if(molecular_capture.neqv.dust_pah_h2_enabled())ierr=1
       return
    endif
    ierr=1
    charged=dust_pah_charged();hydrogenated=dust_pah_hydrogenated()
    molecular_capture=dust_pah_h2_enabled()
    pah_max_primary_ev=merge(13.6d0,4d0,charged)
    if(.not.snrt_dust_contract_loaded.or.snrt_dust_contract_version/=4)return
    call get_environment_variable('SNRT_PAH_NEUTRAL_TABLE',path,status=status)
    if(status/=0.or.len_trim(path)==0)return
    call pah_neutral_optics(trim(path),24,e,w,cs,status)
    if(status/=0)return
    if(charged)then
       call get_environment_variable('SNRT_PAH_ION_TABLE',path,status=status)
       if(status/=0.or.len_trim(path)==0)return
       call pah_charge_optics(trim(path),24,1,ei,wi,ci,status)
       if(status/=0)return
       if(size(ei)/=size(e))return
       if(any(ei/=e).or.any(wi/=w))return
    endif
    n=snrt_dust_contract_number_ir;nt=snrt_dust_contract_number_temperature
    if(n<2.or.nt<2)return
    if(allocated(pah_ir_sigma))deallocate(pah_ir_sigma,pah_ir_supported,pah_primary_sigma, &
         pah_ir_ion_sigma,pah_primary_ion_sigma,pah_level,pah_injection)
    allocate(pah_ir_sigma(n),pah_ir_supported(n),pah_primary_sigma(snrt_dust_contract_number_groups))
    allocate(pah_ir_ion_sigma(n),pah_primary_ion_sigma(size(pah_primary_sigma)))
    allocate(pah_level(dust_pah_nstate()),pah_injection(dust_pah_nstate()))
    pah_ir_ion_sigma=0;pah_primary_ion_sigma=0
    pah_ir_supported=1
    do i=1,n
       x=snrt_dust_contract_ir_energy_ev(i)
       if(x<e(1))then
          ! Explicit Rayleigh/Drude asymptote (LD01 eq12), anchored at the
          ! 1000 micron table boundary. This is a named model continuation,
          ! NOT claimed to be measured table data at centimetre wavelengths.
          pah_ir_sigma(i)=cs(1)*(x/e(1))**2
          if(charged)pah_ir_ion_sigma(i)=ci(1)*(x/e(1))**2
       else if(x>e(size(e)))then
          ! No assumed X-ray continuation. The live caller rejects occupied
          ! unsupported bands before any state commit; zero is a mask only.
          pah_ir_sigma(i)=0;pah_ir_supported(i)=0
       else
          pah_ir_sigma(i)=sample(x,cs)
          if(charged)pah_ir_ion_sigma(i)=sample(x,ci)
       endif
       if(charged.and.x>pah_max_primary_ev)then
          pah_ir_sigma(i)=0;pah_ir_ion_sigma(i)=0;pah_ir_supported(i)=0
       endif
    enddo
    pah_primary_sigma=0
    do i=1,size(pah_primary_sigma)
       x=snrt_dust_contract_absorption_mean_energy_ev(i)
       if(x>=e(1).and.x<=pah_max_primary_ev)then
          pah_primary_sigma(i)=sample(x,cs)
          if(charged)pah_primary_ion_sigma(i)=sample(x,ci)
       endif
    enddo
    pah_level(1)=0
    do i=2,dust_pah_nbin
       pah_level(i)=1.602176634d-12*.012d0*(64d0/.012d0)**(real(i-2,real64)/(dust_pah_nbin-2))
    enddo
    call dust_pah_modes(24,12,modes,status)
    if(status/=0)return
    call dust_vibrational_curve(modes,[dust_injection_temperature],u,cv,status)
    if(status/=0.or.u(1)>pah_level(dust_pah_nbin))return
    k=1
    do while(k<dust_pah_nbin-1)
       if(u(1)<=pah_level(k+1))exit
       k=k+1
    enddo
    f=(u(1)-pah_level(k))/(pah_level(k+1)-pah_level(k))
    pah_injection=0;pah_injection(k)=1-f;pah_injection(k+1)=f
    pah_injection_specific_u=u(1)/dust_pah_molecule_g
    call pah_radiative_prepare(pah_live_model,24,12,pah_level(1:dust_pah_nbin), &
         snrt_dust_contract_ir_energy_ev(1:n),snrt_dust_contract_ir_weight_ev(1:n),pah_ir_sigma,status)
    if(status/=0)return
    if(charged)then
       pah_charge_models(1)=pah_live_model
       call pah_radiative_prepare(pah_charge_models(2),24,12,pah_level(1:dust_pah_nbin), &
            snrt_dust_contract_ir_energy_ev(1:n),snrt_dust_contract_ir_weight_ev(1:n),pah_ir_ion_sigma,status)
       if(status/=0)return
       ! Ionization energy is carried by the cation mass, not by gas heat or
       ! a duplicate passive energy field. Injection remains neutral only.
       if(.not.hydrogenated)pah_level(dust_pah_nbin+1:)=pah_level(1:dust_pah_nbin)+7.02d0*1.602176634d-12
    endif
    if(hydrogenated)then
       base=pah_level(1:dust_pah_nbin);injection=pah_injection(1:dust_pah_nbin);pah_injection=0
       pah_injection(12*dust_pah_nbin+1:13*dust_pah_nbin)=injection
       pah_injection_specific_u=u(1)/dust_pah_state_mass(12*dust_pah_nbin+1)
       do q=1,2
          call pah_hydrogen_parameters(q-1,bond,attach,status)
          if(status/=0)return
          call pah_hydrogen_prepare(pah_h_models(q),pah_charge_models(q),q-1,.0005d0*1.602176634d-12,status, &
               h2_capture=molecular_capture)
          if(status/=0)return
          do nh=0,13
             offset=((q-1)*14+nh)*dust_pah_nbin
             pah_level(offset+1:offset+dust_pah_nbin)=base+bond(nh)+(q-1)*7.02d0*1.602176634d-12
          enddo
       enddo
    endif
    ready=.true.;ierr=0
  contains
    real(real64) function sample(x,cross) result(sigma)
      real(real64),intent(in)::x
      real(real64),intent(in)::cross(:)
      integer::j
      real(real64)::a
      j=1
      do while(j<size(e)-1)
         if(x<=e(j+1))exit
         j=j+1
      enddo
      a=log(x/e(j))/log(e(j+1)/e(j))
      sigma=exp((1-a)*log(cross(j))+a*log(cross(j+1)))
    end function
  end subroutine

  subroutine pah_live_identity(values,ierr)
    real(real64),allocatable,intent(out)::values(:)
    integer,intent(out)::ierr
    real(real64),allocatable::h1(:),h2(:)
    call pah_live_prepare(ierr)
    if(ierr/=0)return
    ! Bind actual projected coefficients, energy/temperature axes, quadrature,
    ! representation and mass conventions; never a path-only restart match.
    ! Version 2 binds the zero-T DL01/Fe bulk and cold Planck continuation
    ! used by the shared absolute-radiation callback (also when Fe is off).
    values=[2d0,24d0,12d0,dust_pah_molecule_g,dust_pah_hc,2d0,4d0,64d0, &
         real(dust_pah_nbin,real64),pah_level,pah_injection,pah_injection_specific_u, &
         pah_ir_sigma,pah_ir_supported,pah_primary_sigma, &
         snrt_dust_contract_ir_energy_ev(1:size(pah_ir_sigma)), &
         snrt_dust_contract_ir_weight_ev(1:size(pah_ir_sigma))]
    if(charged)values=[3d0,values(2:),2d0,7.02d0,13.6d0,1d-5,300d0,10d0,1d4, &
         pah_ir_ion_sigma,pah_primary_ion_sigma]
    if(hydrogenated)then
       call pah_hydrogen_identity(pah_h_models(1),h1,ierr)
       if(ierr/=0)return
       call pah_hydrogen_identity(pah_h_models(2),h2,ierr)
       if(ierr/=0)return
       values=[4d0,values(2:),14d0,1.66d-24,h1,h2]
       if(molecular_capture)values(1)=5d0
    endif
  end subroutine
end module
