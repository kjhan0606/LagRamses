! Native live state adapter. Species are transported as m_H*n_i, in code
! density units; their sum is never added to rho, metals or the gas energy.
module snrt_chimes_runtime
  use iso_c_binding, only: c_null_char,c_ptr,c_null_ptr
  use amr_commons
  use hydro_commons
  use snrt_chimes
  use dust_iron_photons, only: fe_uv_step
  use snrt_spectral_contract, only: snrt_chimes_band_enabled,snrt_chimes_bank_sha256, &
       snrt_chimes_cold_enabled,snrt_chimes_molecular_sha256,snrt_chimes_transition_enabled
  use snrt_thermochemistry, only: snrt_secondary_tables_loaded,snrt_secondary_tables_load_from_environment
  use snrt_atomic_cooling, only: atomic_mh
  use dust_composition_optics, only: d03_band_abs,d03_band_transport,d03_radius_cm,d03_solid_density
  use snrt_moving_scatter, only: snrt_moving_scatter_cell,snrt_scatter_c
  use dust_phase_state, only: dust_dynamics_enabled,dust_phase_kinetic
  use dust_mass_physics, only: dust_chimes_enabled,dust_gas_elements,olivine_fraction,dust_iron_enabled, &
       dust_pah_enabled,dust_pah_nstate,dust_pah_hc,dust_pah_charged,dust_pah_solid_charge,dust_pah_inventory
  use dust_mass_physics, only: dust_condensation,dust_growth,dust_sputtering,dust_sn_shocks, &
       dust_sublimation,dust_coagulation,dust_shattering,dust_two_size_enabled
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public::chimes_live_initialize,chimes_live_identity,chimes_prepare_level
  public::chimes_initialize_hierarchy
  public::chimes_consistent_carriers
  public::chimes_cell_state,chimes_grain_budget,chimes_live_capacity,chimes_live_stage
  public::chimes_live_band_stage
  public::chimes_live_molecular_stage
  public::chimes_live_cold_stage
  public::chimes_live_grain_scatter
  public::chimes_live_fe_uv_stage
  real(dp),parameter::mass_number(11)=[1d0,4d0,12d0,14d0,16d0,20d0,24d0,28d0,32d0,40d0,56d0]
  logical,save::ready=.false.
  type(c_ptr),save::band_handle=c_null_ptr
  type(c_ptr),save::molecular_handle=c_null_ptr
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
       ! Old spectral modes keep fixed masses. Kind7 also admits the existing
       ! mass/size operator, before opacity is reconstructed.
       if(dust_iron_enabled().or.dust_pah_enabled())return
       if(.not.snrt_chimes_transition_enabled())then
          if(dust_dynamics_enabled())return
          if(any(dust_condensation/=0).or.dust_sn_shocks.or.trim(dust_sublimation)/='none')return
          if(dust_growth.or.dust_sputtering.or.dust_coagulation.or.dust_shattering)return
       endif
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
       if(snrt_chimes_cold_enabled())then
          if(.not.dust_two_size_enabled())then
             ierr=1;return
          endif
          call get_environment_variable('SNRT_CHIMES_MOLECULAR_TABLE',main,length=n,status=status)
          if(status/=0.or.n<1.or.n>=len(main))then
             ierr=1;return
          endif
          ierr=chimes_molecular_load(trim(main)//c_null_char,molecular_handle,bank_identity)
          if(ierr/=0)return
          do i=1,32
             write(bank_hash(2*i-1:2*i),'(Z2.2)')nint(bank_identity(i))
          enddo
          do i=1,64
             n=iachar(bank_hash(i:i))
             if(n>=iachar('A').and.n<=iachar('F'))bank_hash(i:i)=achar(n+32)
          enddo
          if(bank_hash/=snrt_chimes_molecular_sha256)then
             call chimes_molecular_free(molecular_handle);molecular_handle=c_null_ptr;ierr=1;return
          endif
       endif
    endif
    if(ierr==0.and.snrt_chimes_transition_enabled())then
       if(chimes_transition_supported()/=1)ierr=1
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
    if(snrt_chimes_cold_enabled())values(321)=6d0
    if(snrt_chimes_transition_enabled())values(321)=7d0
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
    ! H is not a C/olivine growth donor. Size exchange uses total hydrogen
    ! nuclei, including H2, not the molecular-depleted accretion reservoir.
    ! Apply to every CHIMES mass path, not only the spectral comparison.
    locked(1)=0
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
    ! Serial level boundary: update the proper CMB before private cell solves.
    ! This is derived from the restored AMR clock, not independent restart state.
    if(status==0)status=chimes_set_expansion(merge(aexp,1d0,cosmo))
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

  subroutine chimes_live_cold_stage(cell,sd,sv,dt,length,td,chat,nd,number,radiation_energy, &
       state,energy,next_number,next_energy,ledger,grain_number,grain_energy,ierr,event_info, &
       staged_row,directions,phase_energy,phase_moment)
    integer,intent(in)::cell,nd
    real(dp),intent(in)::sd,sv,dt,length,td,chat,number(nd,9),radiation_energy(nd,9)
    real(dp),intent(inout)::state(chimes_ns),energy,next_number(nd,9),next_energy(nd,9),ledger(11)
    real(dp),intent(inout)::grain_number(9),grain_energy(9)
    integer,intent(out)::ierr
    real(dp),optional,intent(inout)::event_info(2)
    real(dp),optional,intent(in)::staged_row(nvar),directions(3,nd)
    real(dp),optional,intent(inout)::phase_energy(9,4),phase_moment(3,4)
    ierr=1
    if(.not.snrt_chimes_cold_enabled())return
    call chimes_live_molecular_stage(cell,band_handle,molecular_handle,sd,sv,dt,length,td,chat,nd, &
         number,radiation_energy,state,energy,next_number,next_energy,ledger,grain_number,grain_energy,ierr,event_info, &
         staged_row,directions,phase_energy,phase_moment)
  end subroutine

  subroutine chimes_live_molecular_stage(cell,atomic_bank,molecular_bank,sd,sv,dt,length,td,chat,nd, &
       number,radiation_energy,state,energy,next_number,next_energy,ledger,grain_number,grain_energy,ierr,event_info, &
       staged_row,directions,phase_energy,phase_moment)
    ! Caller-owned handles must be loaded/bound before threaded cell work.
    ! This adapter does not select a model or publish uold/radiation. The
    ! driver must disable transport absorption and stage IR before commit.
    integer,intent(in)::cell,nd
    type(c_ptr),intent(in)::atomic_bank,molecular_bank
    real(dp),intent(in)::sd,sv,dt,length,td,chat,number(nd,9),radiation_energy(nd,9)
    real(dp),intent(inout)::state(chimes_ns),energy,next_number(nd,9),next_energy(nd,9),ledger(11)
    real(dp),intent(inout)::grain_number(9),grain_energy(9)
    integer,intent(out)::ierr
    real(dp),optional,intent(inout)::event_info(2)
    real(dp),optional,intent(in)::staged_row(nvar),directions(3,nd)
    real(dp),optional,intent(inout)::phase_energy(9,4),phase_moment(3,4)
    real(dp)::alpha_phase(128,9,4),phase_e(9,4),phase_p(3,4)
    real(dp)::old(chimes_ns),elements(11),new(chimes_ns),controls(9),bins(4),grains(2),alpha(128,9)
    real(dp)::pn(nd,9),pe(nd,9),budget(11),gn(9),ge(9),row(nvar),staged(chimes_ns)
    real(dp)::nh,cv,t,nonthermal,thermal,area,ratio,surface,proposed_energy,events(2)
    integer::s,k,j
    ierr=1
    if(.not.ready.or..not.allocated(uold))return
    if(cell<lbound(uold,1).or.cell>ubound(uold,1).or.nd<1.or.nd>720)return
    if(any(.not.ieee_is_finite([sd,sv,dt,length,td,chat])).or.sd<=0.or.sv<=0)return
    if(abs(gamma-5d0/3d0)>1d-12)return
    if(.not.dust_two_size_enabled())return
    if(dust_iron_enabled().or.dust_pah_enabled())return
    if(dust_dynamics_enabled().and..not.snrt_chimes_transition_enabled())return
    if(present(staged_row).neqv.dust_dynamics_enabled())return
    if(present(directions).neqv.dust_dynamics_enabled())return
    if(present(phase_energy).neqv.dust_dynamics_enabled())return
    if(present(phase_moment).neqv.dust_dynamics_enabled())return
    if(idust<1.or.idust_bins<1.or.idust_bins+3>nvar.or.idust_species<1.or.idust_species+1>nvar)return
    row=uold(cell,1:nvar)
    if(present(staged_row))row=staged_row
    if(any(.not.ieee_is_finite(row)).or.row(1)<=0)return
    bins=row(idust_bins:idust_bins+3)
    if(any(bins<0))return
    grains=[sum(bins(1:2)),sum(bins(3:4))]
    if(any(abs(grains-row(idust_species:idust_species+1))>1d-8*max(grains,1d-30)))return
    if(abs(sum(grains)-row(idust))>1d-8*max(sum(grains),1d-30))return
    call chimes_cell_state(cell,grains,old,elements,ierr)
    if(ierr/=0)return
    ierr=1
    nonthermal=.5d0*sum(row(2:ndim+1)**2)/row(1)+magnetic_energy(row)
    if(dust_dynamics_enabled())nonthermal=dust_phase_kinetic(row)+magnetic_energy(row)
#if NENER>0
    nonthermal=nonthermal+sum(row(inener:inener+NENER-1))
#endif
    thermal=(row(ndim+2)-nonthermal)*sd*sv**2
    cv=chimes_live_capacity(old,sd)
    if(cv<=0.or.thermal<=0)return
    nh=elements(1)*sd/atomic_mh;t=thermal/cv
    alpha=0;area=0
    do s=1,2
       do k=1,2
          j=2*(s-1)+k
          alpha_phase(:,:,j)=bins(j)*sd*d03_band_abs(:,:,j)
          alpha=alpha+alpha_phase(:,:,j)
          area=area+bins(j)*sd*.75d0/(d03_solid_density(s)*d03_radius_cm(k))
       enddo
    enddo
    ! Same C/silicate geometric-area normalization as the existing grey
    ! live receiver. Optical Q is not the catalytic cross section.
    area=area/nh;ratio=sum(grains)/(elements(1)*.01d0);surface=1
    if(ratio>0)surface=area/(1d-21*ratio)
    controls=[nh,t,td,dt,length,ratio,surface,0d0,chat]
    if(dust_dynamics_enabled())then
    ierr=chimes_cell_band_cold_molecular(atomic_bank,molecular_bank,nd,controls,elements/elements(1), &
         old/elements(1),alpha,number,radiation_energy,t,new,pn,pe,budget,gn,ge, &
         transition=snrt_chimes_transition_enabled(),event_info=events,phase_alpha=alpha_phase, &
         directions=directions,phase_energy=phase_e,phase_moment=phase_p)
    else
    ierr=chimes_cell_band_cold_molecular(atomic_bank,molecular_bank,nd,controls,elements/elements(1), &
         old/elements(1),alpha,number,radiation_energy,t,new,pn,pe,budget,gn,ge, &
         transition=snrt_chimes_transition_enabled(),event_info=events)
    endif
    if(ierr/=0)return
    staged=new*elements(1)
    proposed_energy=nonthermal+chimes_live_capacity(staged,sd)*t/(sd*sv**2)
    ierr=1
    if(.not.ieee_is_finite(proposed_energy).or.any(.not.ieee_is_finite(staged)))return
    state=staged;energy=proposed_energy;next_number=pn;next_energy=pe;ledger=budget
    grain_number=gn;grain_energy=ge;ierr=0
    if(present(event_info))event_info=events
    if(present(phase_energy))phase_energy=phase_e
    if(present(phase_moment))phase_moment=phase_p
  end subroutine

  subroutine chimes_live_grain_scatter(nd,number,energy,momentum,mass,directions,weights,dt,chat,work,ierr)
    integer,intent(in)::nd
    real(dp),intent(inout)::number(nd,9),energy(nd,9),momentum(3,4),work(4)
    real(dp),intent(in)::mass(4),directions(3,nd),weights(nd),dt,chat
    integer,intent(out)::ierr
    real(dp),allocatable::nn(:,:),ee(:,:)
    real(dp)::tau(4,9),p(3,4),w(4),gn(nd,9),ge(nd,9),group_n
    integer::b,g,k,j
    real(dp),parameter::ev_erg=1.602176634d-12
    ierr=1
    if(.not.ready.or..not.snrt_chimes_transition_enabled())return
    if(nd<1.or.nd>720.or.any(.not.ieee_is_finite([mass,dt,chat])).or.any(mass<0).or.dt<0.or.chat<=0.or.chat>1)return
    p=momentum;w=0;tau=0
    ! Explicit group-grey moving closure: photon-weighted D03 transport Q
    ! from the current per-ray node reconstruction, frozen for this call.
    ! The existing nine-group moving solver conserves N and kinetic+lab E;
    ! this is not a new node-resolved or cross-group Doppler transport model.
    ! Exact zero incident N AND E: no node reconstruction is needed.
    ! Still run the receiver below to validate phase state/directions and
    ! preserve its transactional work/momentum contract. This is local
    ! scattering only, never permission to skip transport or dark chemistry.
    if(any(number/=0).or.any(energy/=0))then
    allocate(nn(nd,128*9),ee(nd,128*9))
    ierr=chimes_band_nodes(band_handle,nd,number,energy,nn,ee)
    if(ierr/=0)return
    do g=1,9
       group_n=sum(number(:,g))
       if(group_n<=0)cycle
       do k=1,128
          j=k+128*(g-1)
          do b=1,4
             tau(b,g)=tau(b,g)+sum(nn(:,j))/group_n*mass(b)*d03_band_transport(k,g,b)*snrt_scatter_c*chat*dt
          enddo
       enddo
    enddo
    endif
    gn=number;ge=energy*ev_erg
    call snrt_moving_scatter_cell(gn,ge,p,mass,tau,directions,weights,snrt_scatter_c,w,ierr)
    if(ierr/=0)return
    number=gn;energy=ge/ev_erg;momentum=p;work=w
  end subroutine chimes_live_grain_scatter

  subroutine chimes_live_fe_uv_stage(cell,sd,sv,dt,td,chat,photons,state,energy,solid,ledger,ierr)
    ! Operate on the ALREADY staged CHIMES state; never restore pre-chemistry
    ! abundances here. All arguments remain unchanged if any check fails.
    ! Here chat is cm/s, NOT the dimensionless light-speed fraction passed
    ! to chimes_live_stage's external CHIMES controls.
    integer,intent(in)::cell
    real(dp),intent(in)::sd,sv,dt,td,chat
    real(dp),intent(inout)::photons(9),state(chimes_ns),energy,solid,ledger(8)
    integer,intent(out)::ierr
    real(dp)::elements(11),after(11),q0,q1,nh,nonthermal,e,u,n(9),s(chimes_ns),receipt(8),bins(6)
    ierr=1
    if(.not.ready.or.min(sd,sv)<=0)return
    ierr=chimes_budget(state,elements,q0)
    if(ierr/=0)return
    ierr=1
    if(elements(1)<=0)return
    nh=elements(1)*sd/atomic_mh
    nonthermal=.5d0*sum(uold(cell,2:ndim+1)**2)/uold(cell,1)+magnetic_energy(uold(cell,:))
#if NENER>0
    nonthermal=nonthermal+sum(uold(cell,inener:inener+NENER-1))
#endif
    e=(energy-nonthermal)*sd*sv**2;u=solid;n=photons;s=state/elements(1);receipt=0
    bins=[uold(cell,idust_bins:idust_bins+3),uold(cell,idust_iron:idust_iron+1)]*sd
    call fe_uv_step(dt,chat,nh,td,bins,n,s,e,u,fe_secondary,receipt,ierr)
    if(ierr/=0)return
    s=s*elements(1)
    ierr=chimes_budget(s,after,q1)
    if(ierr/=0)return
    ierr=9
    if(any(abs(after-elements)>2d-12*max(elements,tiny(1d0))))return
    if(abs(q1-q0)>2d-12*max(sum(state),tiny(1d0)))return
    state=s;photons=n;energy=nonthermal+e/(sd*sv**2);solid=u;ledger=receipt;ierr=0
  contains
    subroutine fe_secondary(electron_ev,abundance,f,status)
      real(dp),intent(in)::electron_ev,abundance(157)
      real(dp),intent(out)::f(5)
      integer,intent(out)::status
      f=0
      status=chimes_secondary_partition(electron_ev,abundance,f)
    end subroutine
  end subroutine chimes_live_fe_uv_stage

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
