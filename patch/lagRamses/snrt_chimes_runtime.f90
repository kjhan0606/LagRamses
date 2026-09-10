! Native live state adapter. Species are transported as m_H*n_i, in code
! density units; their sum is never added to rho, metals or the gas energy.
module snrt_chimes_runtime
  use iso_c_binding, only: c_null_char,c_ptr,c_null_ptr
  use amr_commons
  use hydro_commons
  use snrt_chimes
  use snrt_spectral_contract, only: snrt_chimes_band_enabled,snrt_chimes_bank_sha256
  use snrt_thermochemistry, only: snrt_secondary_tables_loaded,snrt_secondary_tables_load_from_environment
  use snrt_atomic_cooling, only: atomic_mh
  use dust_phase_state, only: dust_dynamics_enabled,dust_phase_kinetic
  use dust_mass_physics, only: dust_chimes_enabled,dust_gas_elements,olivine_fraction,dust_iron_enabled, &
       dust_pah_enabled,dust_pah_nstate,dust_pah_hc,dust_pah_charged,dust_pah_solid_charge,dust_pah_inventory
  use dust_mass_physics, only: dust_condensation,dust_growth,dust_sputtering,dust_sn_shocks, &
       dust_sublimation,dust_coagulation,dust_shattering
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public::chimes_live_initialize,chimes_live_identity,chimes_prepare_level
  public::chimes_initialize_hierarchy
  public::chimes_consistent_carriers
  public::chimes_cell_state,chimes_grain_budget,chimes_live_capacity,chimes_live_stage
  public::chimes_live_band_stage
  real(dp),parameter::mass_number(11)=[1d0,4d0,12d0,14d0,16d0,20d0,24d0,28d0,32d0,40d0,56d0]
  logical,save::ready=.false.
  type(c_ptr),save::band_handle=c_null_ptr
contains
  real(dp) function cell_solid_charge(cell) result(q)
    integer,intent(in)::cell
    q=0
    if(dust_pah_charged())q=dust_pah_solid_charge(uold(cell,idust_pah:idust_pah+dust_pah_nstate()-1),atomic_mh)
  end function
  subroutine chimes_consistent_carriers(values)
    ! Consistent multi-species advection: one common normalization for all
    ! independent chemical and dust mass carriers; derive total element
    ! carriers from exactly those fluxes/states. Valid for signed face flux.
    real(dp),intent(inout)::values(nvar)
    real(dp)::elements(11),total,factor,hc(2)
    integer::status
    if(.not.dust_chimes_enabled())return
    status=chimes_nuclear_sums(values(ichimes:ichimes+chimes_ns-1),elements)
    if(status/=0)then
       call clean_stop;return
    endif
    total=sum(elements*mass_number)+sum(values(idust_bins:idust_bins+3))
    if(dust_iron_enabled())total=total+sum(values(idust_iron:idust_iron+1))
    if(dust_pah_enabled())total=total+sum(values(idust_pah:idust_pah+dust_pah_nstate()-1))
    if(total==0)then
       if(values(1)/=0)call clean_stop
       return
    endif
    factor=values(1)/total
    if(factor<0.or..not.ieee_is_finite(factor))then
       call clean_stop;return
    endif
    values(ichimes:ichimes+chimes_ns-1)=factor*values(ichimes:ichimes+chimes_ns-1)
    values(idust:idust+10)=factor*values(idust:idust+10)
    values(idust_species)=sum(values(idust_bins:idust_bins+1))
    values(idust_species+1)=sum(values(idust_bins+2:idust_bins+3))
    values(idust)=sum(values(idust_species:idust_species+1))
    elements=factor*elements*mass_number+olivine_fraction*values(idust_species+1)
    elements(3)=elements(3)+values(idust_species)
    if(dust_iron_enabled())then
       values(idust_iron:idust_iron+1)=factor*values(idust_iron:idust_iron+1)
       values(idust)=values(idust)+sum(values(idust_iron:idust_iron+1))
       elements(11)=elements(11)+sum(values(idust_iron:idust_iron+1))
    endif
    values(ichem:ichem+10)=elements
    if(dust_pah_enabled())then
       values(idust_pah:idust_pah+dust_pah_nstate()-1)=factor*values(idust_pah:idust_pah+dust_pah_nstate()-1)
       hc=dust_pah_inventory(values(idust_pah:idust_pah+dust_pah_nstate()-1))
       elements(1)=elements(1)+hc(1)
       elements(3)=elements(3)+hc(2)
       values(ichem:ichem+10)=elements
    endif
    values(imetal)=sum(elements(3:11))
  end subroutine

  subroutine chimes_live_initialize(ierr)
    integer,intent(out)::ierr
    character(len=500)::main,dir,paths(9)
    character(len=2)::number
    integer::i,n,status,nreaction,nshell
    character(len=64)::bank_hash
    real(dp)::bank_identity(32)
    ierr=0
    if(ready.or..not.dust_chimes_enabled())return
    ! Fresh hierarchy initialization and restart identity reads can precede
    ! the first RT step. Load the same pinned FS contract before the C
    ! receiver is initialized; never lazy-load from an OpenMP cell callback.
    if(.not.snrt_secondary_tables_loaded)then
       call snrt_secondary_tables_load_from_environment(ierr)
       if(ierr/=0)return
    endif
    ierr=1
    call get_environment_variable('SNRT_CHIMES_MAIN_DATA',main,length=n,status=status)
    if(status/=0.or.n<1.or.n>=len(main))return
    call get_environment_variable('SNRT_CHIMES_GROUP_DIR',dir,length=n,status=status)
    if(status/=0.or.n<1.or.n+15>=len(dir))return
    do i=1,9
       write(number,'(I2.2)')i
       paths(i)=trim(dir)//'/group_'//number//'.hdf5'//c_null_char
    enddo
    ierr=chimes_initialize(trim(main)//c_null_char,9,paths,500)
    if(ierr==0.and.dust_pah_charged())then
       if(chimes_charge_supported()/=1)ierr=1
    endif
    if(ierr==0.and.snrt_chimes_band_enabled())then
       ierr=1
       ! This receiver is atomic and dust-free, not a molecular/dust model.
       if(dust_iron_enabled().or.dust_pah_enabled().or.dust_dynamics_enabled())return
       if(any(dust_condensation/=0).or.dust_growth.or.dust_sputtering.or.dust_sn_shocks)return
       if(dust_coagulation.or.dust_shattering.or.trim(dust_sublimation)/='none')return
       call get_environment_variable('SNRT_CHIMES_BAND_TABLE',main,length=n,status=status)
       if(status/=0.or.n<1.or.n>=len(main))return
       ierr=chimes_band_load(trim(main)//c_null_char,band_handle,nreaction,nshell,bank_identity)
       if(ierr/=0)return
       do i=1,32
          write(bank_hash(2*i-1:2*i),'(Z2.2)')nint(bank_identity(i))
       enddo
       ! Hex from Z editing is uppercase; normalize before exact binding.
       do i=1,64
          n=iachar(bank_hash(i:i))
          if(n>=iachar('A').and.n<=iachar('F'))bank_hash(i:i)=achar(n+32)
       enddo
       if(bank_hash/=snrt_chimes_bank_sha256.or.nreaction/=311.or.nshell/=682)then
          call chimes_band_free(band_handle);band_handle=c_null_ptr;ierr=1;return
       endif
    endif
    if(ierr==0)ready=.true.
  end subroutine

  subroutine chimes_live_identity(values,ierr)
    real(dp),intent(out)::values(324)
    integer,intent(out)::ierr
    values=0
    call chimes_live_initialize(ierr)
    if(ierr/=0)return
    ierr=chimes_identity(values(1:320))
    ! Version fixes local shielding, first-order transport/chemistry split,
    ! translational EOS, zero additional CR ionization, optically thin
    ! unresolved fluorescent/cascade cooling (not injected into dust IR).
    ! Version 4 adds up to two tighter-tolerance retries of nucleus/charge
    ! failures, always from the original cell input and with unchanged gates.
    values(321:324)=[4d0,real(ichimes,dp),atomic_mh,chimes_boltzmann()]
    if(snrt_chimes_band_enabled())values(321)=5d0
  end subroutine

  subroutine chimes_cell_state(cell,grains,state,elements,ierr,previous_state,metallic_iron,pah_hc)
    integer,intent(in)::cell
    real(dp),intent(in)::grains(2)
    real(dp),intent(out)::state(chimes_ns),elements(11)
    integer,intent(out)::ierr
    real(dp)::gas(11),solid_fe,solid_pah(2)
    real(dp),optional,intent(in)::previous_state(chimes_ns)
    real(dp),optional,intent(in)::metallic_iron
    real(dp),optional,intent(in)::pah_hc(2)
    state=0;elements=0;ierr=1
    if(.not.ready.or.ichimes<1.or.ichimes+chimes_ns-1>nvar)return
    solid_fe=0
    if(dust_iron_enabled())then
       if(any(.not.ieee_is_finite(uold(cell,idust_iron:idust_iron+1))))return
       if(any(uold(cell,idust_iron:idust_iron+1)<0))return
       solid_fe=sum(uold(cell,idust_iron:idust_iron+1))
    endif
    if(present(metallic_iron))solid_fe=metallic_iron
    solid_pah=0
    if(dust_pah_enabled())then
       if(any(.not.ieee_is_finite(uold(cell,idust_pah:idust_pah+dust_pah_nstate()-1))))return
       if(any(uold(cell,idust_pah:idust_pah+dust_pah_nstate()-1)<0))return
       solid_pah=dust_pah_inventory(uold(cell,idust_pah:idust_pah+dust_pah_nstate()-1))
    endif
    if(present(pah_hc))solid_pah=pah_hc
    call dust_gas_elements(uold(cell,ichem:ichem+10),grains,gas,ierr,solid_fe,solid_pah)
    if(ierr/=0)return
    elements=gas/mass_number
    if(elements(1)<=0)then
       ierr=1;return
    endif
    if(all(uold(cell,ichimes:ichimes+chimes_ns-1)==0d0).and.nrestart>0)then
       ierr=1;return
    endif
    ! Feedback adds neutral atoms; dust destruction returns neutral atoms.
    ! For an evolved cell, preserve the transported ion/molecular state.
    if(present(previous_state))then
       ierr=chimes_reconcile_charged(elements,previous_state,cell_solid_charge(cell),state)
    else
       ierr=chimes_reconcile_charged(elements,uold(cell,ichimes:ichimes+chimes_ns-1),cell_solid_charge(cell),state)
    endif
  end subroutine

  subroutine chimes_grain_budget(cell,elements,ierr)
    integer,intent(in)::cell
    real(dp),intent(out)::elements(11)
    integer,intent(out)::ierr
    real(dp)::locked(11)
    elements=uold(cell,ichem:ichem+10)
    ierr=chimes_locked(uold(cell,ichimes:ichimes+chimes_ns-1),locked)
    if(ierr/=0)return
    elements=elements-locked*mass_number
    if(any(elements<0))ierr=1
  end subroutine

  real(dp) function chimes_live_capacity(state,density_unit) result(cv)
    real(dp),intent(in)::state(chimes_ns),density_unit
    cv=1.5d0*chimes_boltzmann()*sum(state)*density_unit/atomic_mh
  end function

  subroutine chimes_prepare_level(ilevel,synchronize)
    integer,intent(in)::ilevel
    logical,optional,intent(in)::synchronize
    integer::ind,i,cell,k,n,bad,all_bad,status,info
    integer,allocatable::cells(:)
    real(dp),allocatable::stage(:,:)
    real(dp)::elements(11)
    include 'mpif.h'
    if(.not.dust_chimes_enabled())return
    call chimes_live_initialize(status)
    bad=merge(0,1,status==0);k=0
    n=active(ilevel)%ngrid*twotondim
    allocate(cells(n),stage(chimes_ns,n))
    if(bad==0)then
       do ind=1,twotondim
          do i=1,active(ilevel)%ngrid
             cell=ICELL_OF(active(ilevel)%igrid(i),ind)
             if(son(cell)/=0)cycle
             k=k+1;cells(k)=cell
             call chimes_cell_state(cell,uold(cell,idust_species:idust_species+1),stage(:,k),elements,status)
             if(status/=0)bad=1
          enddo
       enddo
    endif
    call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(all_bad/=0.or.info/=0)then
       if(myid==1)write(*,*)'ERROR: CHIMES element/transport/source reconciliation failed'
       call MPI_ABORT(MPI_COMM_WORLD,31,info)
    endif
    do i=1,k
       uold(cells(i),ichimes:ichimes+chimes_ns-1)=stage(:,i)
    enddo
    deallocate(cells,stage)
    if(present(synchronize))then
       if(synchronize)then
          do i=ichimes,ichimes+chimes_ns-1
             call make_virtual_fine_dp(uold(1,i),ilevel)
          enddo
       endif
    endif
  end subroutine

  subroutine chimes_initialize_hierarchy()
    integer::ilevel,i,ind,cell,status,bad,all_bad,info
    real(dp)::elements(11),measured(11),state(chimes_ns),q
    include 'mpif.h'
    if(.not.dust_chimes_enabled())return
    if(nrestart>0)then
       ! A checkpoint already contains the physical state and coarse/fine
       ! update phase. Validate it without an extra reconciliation or hydro
       ! restriction: either can change the restart trajectory.
       call chimes_live_initialize(status)
       bad=merge(0,1,status==0)
       if(bad==0)then
          do ilevel=levelmin,nlevelmax
             do ind=1,twotondim
                do i=1,active(ilevel)%ngrid
                   cell=ICELL_OF(active(ilevel)%igrid(i),ind)
                   if(son(cell)/=0)cycle
                   call chimes_cell_state(cell,uold(cell,idust_species:idust_species+1),state,elements,status)
                   if(status/=0)bad=1
                   status=chimes_budget(uold(cell,ichimes:ichimes+chimes_ns-1),measured,q)
                   if(status/=0.or.any(abs(measured-elements)>1d-8*max(elements,1d-30)).or. &
                        abs(q+cell_solid_charge(cell))>1d-10*max(elements(1),1d-30))bad=1
                enddo
             enddo
          enddo
       endif
       call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
       if(all_bad/=0.or.info/=0)then
          if(myid==1)write(*,*)'ERROR: restart chemical state violates its element/charge ledger'
          call MPI_ABORT(MPI_COMM_WORLD,31,info)
       endif
       return
    endif
    do ilevel=nlevelmax,levelmin,-1
       call chimes_prepare_level(ilevel)
       call upload_fine(ilevel)
       do i=ichimes,ichimes+chimes_ns-1
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
    enddo
  end subroutine

  subroutine chimes_live_band_stage(cell,sd,sv,dt,length,td,chat,nd,number,radiation_energy, &
       state,energy,next_number,next_energy,ledger,ierr)
    integer,intent(in)::cell,nd
    real(dp),intent(in)::sd,sv,dt,length,td,chat,number(nd,9),radiation_energy(nd,9)
    real(dp),intent(out)::state(chimes_ns),energy,next_number(nd,9),next_energy(nd,9),ledger(9)
    integer,intent(out)::ierr
    real(dp)::old(chimes_ns),elements(11),abundance(chimes_ns),new(chimes_ns),controls(9)
    real(dp)::nonthermal,nh,t,cv
    state=0;energy=uold(cell,ndim+2);next_number=number;next_energy=radiation_energy;ledger=0;ierr=1
    if(.not.ready.or..not.snrt_chimes_band_enabled())return
    if(any(uold(cell,idust:idust+10)/=0))return
    call chimes_cell_state(cell,[0d0,0d0],old,elements,ierr)
    if(ierr/=0)return
    state=old
    nonthermal=.5d0*sum(uold(cell,2:ndim+1)**2)/uold(cell,1)+magnetic_energy(uold(cell,:))
#if NENER>0
    nonthermal=nonthermal+sum(uold(cell,inener:inener+NENER-1))
#endif
    cv=chimes_live_capacity(old,sd)
    ierr=1
    if(cv<=0)return
    t=(energy-nonthermal)*sd*sv**2/cv
    nh=elements(1)*sd/atomic_mh;abundance=old/elements(1)
    controls=[nh,t,td,dt,length,0d0,1d0,0d0,chat]
    ierr=chimes_cell_band_hot_atomic(band_handle,nd,controls,elements/elements(1),abundance,0d0, &
         number,radiation_energy,t,new,next_number,next_energy,ledger)
    if(ierr/=0)return
    state=new*elements(1)
    energy=nonthermal+chimes_live_capacity(state,sd)*t/(sd*sv**2)
    if(.not.ieee_is_finite(energy))ierr=1
  end subroutine

  subroutine chimes_live_stage(cell,sd,sv,dt,length,td,area,chat,photons,state,energy,next_photons,ierr,staged_row)
    integer,intent(in)::cell
    real(dp),intent(in)::sd,sv,dt,length,td,area,chat,photons(9)
    real(dp),intent(out)::state(chimes_ns),energy,next_photons(9)
    integer,intent(out)::ierr
    ! Optional hydro state AFTER staged mechanical/thermal updates. Read E,
    ! component kinetic energy and CR/nonthermal carriers from this SAME row
    ! so chemistry preserves prior radiation work. Chemical abundances and
    ! grain inventories still come from cell; this is not a composition stage.
    real(dp),optional,intent(in)::staged_row(:)
    real(dp)::old(chimes_ns),elements(11),controls(9),nh,cv,t,e,nonthermal,surface
    real(dp)::abundance(chimes_ns),next_abundance(chimes_ns),reactive_dust
    real(dp)::row(nvar)
    logical,save::failure_reported=.false.
    row=uold(cell,1:nvar);energy=row(ndim+2);next_photons=photons;state=0;ierr=1
    if(present(staged_row))then
       if(size(staged_row)/=nvar)return
       if(any(.not.ieee_is_finite(staged_row)).or.staged_row(1)<=0)return
       row=staged_row;energy=row(ndim+2)
    endif
    call chimes_cell_state(cell,uold(cell,idust_species:idust_species+1),old,elements,ierr)
    state=old
    if(ierr/=0)return
    if(dust_dynamics_enabled())then
       nonthermal=dust_phase_kinetic(row)
    else
       nonthermal=.5d0*sum(row(2:ndim+1)**2)/row(1)
    endif
#if NENER>0
    nonthermal=nonthermal+sum(row(inener:inener+NENER-1))
#endif
    nonthermal=nonthermal+magnetic_energy(uold(cell,:))
    e=(energy-nonthermal)*sd*sv**2
    cv=chimes_live_capacity(old,sd)
    if(cv<=0.or.e<=0)then
       ierr=1;return
    endif
    nh=elements(1)*sd/atomic_mh
    abundance=old/elements(1);t=e/cv
    ! Richings+2014a section 2.4, eq. 2.12 adopts 1e-21 cm^2/H
    ! for Cazaux--Tielens formation at solar dust abundance. Replace that
    ! area with the actual two-size geometric area, not the optical Q.
    ! D/D_MW = (dust/H mass)/0.01; the product D*surface = area/1e-21
    ! is independent of this bookkeeping MW mass normalization. Applying
    ! the same area scaling to grain recombination is a comparison closure.
    surface=1d0
    reactive_dust=uold(cell,idust)
    if(dust_iron_enabled().or.dust_pah_charged())reactive_dust=sum(uold(cell,idust_species:idust_species+1))
    ! The Fe comparison does NOT assign C/silicate catalytic formation or
    ! grain recombination coefficients to metallic surfaces.
    if(reactive_dust>0)surface=area/(1d-21*reactive_dust/(elements(1)*.01d0))
    controls=[nh,t,td,dt,length,reactive_dust/(elements(1)*.01d0),surface,0d0,chat]
    ierr=chimes_cell_charged(controls,elements/elements(1),abundance,cell_solid_charge(cell)/elements(1), &
         photons,t,next_abundance,next_photons)
    if(ierr/=0)then
!$omp critical(chimes_failure_report)
       if(.not.failure_reported)then
          write(*,*)'CHIMES rejected cell/status: ',cell,ierr
          write(*,*)'CHIMES controls [nH,Tgas,Tdust,dt,dx,D,area,CR,chat]: ',controls
          failure_reported=.true.
       endif
!$omp end critical(chimes_failure_report)
       return
    endif
    state=next_abundance*elements(1)
    energy=nonthermal+chimes_live_capacity(state,sd)*t/(sd*sv**2)
    if(.not.ieee_is_finite(energy))ierr=1
  end subroutine
end module
