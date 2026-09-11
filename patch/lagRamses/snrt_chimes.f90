! ISO C binding to the native CHIMES receiver. Integration builds opt in;
! the no-CHIMES production baseline never needs this external dependency.
module snrt_chimes
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_thermochemistry, only: snrt_secondary_fractions,snrt_secondary_tables_loaded, &
       snrt_secondary_energy_slice,snrt_secondary_raw_grid,snrt_secondary_nenergy
  implicit none
  private
  integer,parameter,public::chimes_ns=157,chimes_ng=9
  public::chimes_initialize,chimes_neutral,chimes_budget,chimes_cell
  public::chimes_reconcile,chimes_locked,chimes_identity
  public::chimes_boltzmann
  public::chimes_nuclear_sums
  public::chimes_group_binding
  public::chimes_secondary_partition
  public::chimes_reconcile_charged,chimes_cell_charged
  public::chimes_charge_supported
  public::chimes_band_load,chimes_band_free,chimes_band_reactions,chimes_band_moments
  public::chimes_band_photo_step
  public::chimes_band_photo_dust_step
  public::chimes_molecular_load,chimes_molecular_free,chimes_band_photo_molecular_step
  public::chimes_molecular_factors,chimes_cell_band_cold_molecular
  public::chimes_molecular_temperature_max
  public::chimes_band_photo_molecular_groups
  public::chimes_band_nodes,chimes_band_photo_molecular_phases
  public::chimes_cell_band_hot_atomic
  public::chimes_atomize,chimes_transition_supported,chimes_cell_transition_dark
  public::chimes_round_subnormal_survivors
  interface
     integer(c_int) function chimes_band_nodes(handle,nd,number,energy,nn,ne) bind(C,name='snrt_chimes_band_nodes')
       import
       type(c_ptr),value::handle
       integer(c_int),value::nd
       real(c_double),intent(in)::number(nd,9),energy(nd,9)
       real(c_double),intent(inout)::nn(nd,128,9),ne(nd,128,9)
     end function
     integer(c_int) function chimes_band_photo_molecular_phases(handle,mol,nd,nh,dt,chat,alpha,shield, &
          pumping,old,number,energy,new,next_number,next_energy,ledger,grain_number,grain_energy, &
          phase_alpha,directions,phase_energy,phase_moment) bind(C,name='snrt_chimes_band_photo_molecular_phases')
       import
       type(c_ptr),value::handle,mol
       integer(c_int),value::nd
       real(c_double),value::nh,dt,chat,pumping
       real(c_double),intent(in)::alpha(128,9),shield(2),old(chimes_ns),number(nd,9),energy(nd,9)
       real(c_double),intent(in)::phase_alpha(128,9,4),directions(3,nd)
       real(c_double),intent(inout)::new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(10)
       real(c_double),intent(inout)::grain_number(9),grain_energy(9),phase_energy(9,4),phase_moment(3,4)
     end function
     integer(c_int) function chimes_transition_supported() bind(C,name='snrt_chimes_transition_supported')
       import
     end function
     integer(c_int) function chimes_atomize(t,old,new,tnext,cost) bind(C,name='snrt_chimes_atomize')
       import
       real(c_double),value::t
       real(c_double),intent(in)::old(chimes_ns)
       real(c_double),intent(inout)::new(chimes_ns),tnext,cost
     end function
     integer(c_int) function chimes_cell_transition_dark(controls,elements,old,atomic,t,new,elapsed) &
          bind(C,name='snrt_chimes_cell_transition_dark')
       import
       real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns)
       integer(c_int),value::atomic
       real(c_double),intent(out)::t,new(chimes_ns),elapsed
     end function
     real(c_double) function chimes_molecular_temperature_max() bind(C,name='snrt_chimes_molecular_temperature_max')
       import
     end function
     integer(c_int) function chimes_band_photo_molecular_groups(handle,mol,nd,nh,dt,chat,alpha,shield, &
          pumping,old,number,energy,new,next_number,next_energy,ledger,grain_number,grain_energy) &
          bind(C,name='snrt_chimes_band_photo_molecular_groups')
       import
       type(c_ptr),value::handle,mol
       integer(c_int),value::nd
       real(c_double),value::nh,dt,chat,pumping
       real(c_double),intent(in)::alpha(128,9),shield(2),old(chimes_ns),number(nd,9),energy(nd,9)
       real(c_double),intent(inout)::new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(10)
       real(c_double),intent(inout)::grain_number(9),grain_energy(9)
     end function
     integer(c_int) function chimes_cell_cold_dark(controls,elements,old,t,new) bind(C,name='snrt_chimes_cell_cold_dark')
       import
       real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns)
       real(c_double),intent(out)::t,new(chimes_ns)
     end function
     integer(c_int) function chimes_molecular_factors(t,nh,length,state,factors) bind(C,name='snrt_chimes_molecular_factors')
       import
       real(c_double),value::t,nh,length
       real(c_double),intent(in)::state(chimes_ns)
       real(c_double),intent(inout)::factors(3)
     end function
     integer(c_int) function chimes_molecular_load(path,handle,identity) bind(C,name='snrt_chimes_molecular_load')
       import
       character(c_char),intent(in)::path(*)
       type(c_ptr),intent(inout)::handle
       real(c_double),intent(inout)::identity(32)
     end function
     subroutine chimes_molecular_free(handle) bind(C,name='snrt_chimes_molecular_free')
       import
       type(c_ptr),value::handle
     end subroutine
     integer(c_int) function chimes_band_photo_molecular_step(handle,mol,nd,nh,dt,chat,alpha,shield, &
          pumping,old,number,energy,new,next_number,next_energy,ledger) bind(C,name='snrt_chimes_band_photo_molecular_step')
       import
       type(c_ptr),value::handle,mol
       integer(c_int),value::nd
       real(c_double),value::nh,dt,chat,pumping
       real(c_double),intent(in)::alpha(128,9),shield(2),old(chimes_ns),number(nd,9),energy(nd,9)
       real(c_double),intent(inout)::new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(10)
     end function
     integer(c_int) function chimes_band_photo_dust_step(handle,nd,nh,dt,chat,alpha,old,number,energy,new, &
          next_number,next_energy,ledger) bind(C,name='snrt_chimes_band_photo_dust_step')
       import
       type(c_ptr),value::handle
       integer(c_int),value::nd
       real(c_double),value::nh,dt,chat
       real(c_double),intent(in)::alpha(128,9),old(chimes_ns),number(nd,9),energy(nd,9)
       real(c_double),intent(inout)::new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(10)
     end function
     integer(c_int) function chimes_band_photo_step(handle,nd,nh,dt,chat,old,number,energy,new, &
          next_number,next_energy,ledger) bind(C,name='snrt_chimes_band_photo_step')
       import
       type(c_ptr),value::handle
       integer(c_int),value::nd
       real(c_double),value::nh,dt,chat
       real(c_double),intent(in)::old(chimes_ns),number(nd,9),energy(nd,9)
       real(c_double),intent(inout)::new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(8)
     end function
     ! Atomic spectral building block only; this does not enable the live
     ! CHIMES+N/E receiver. Load/free outside threaded cell loops.
     integer(c_int) function chimes_band_load(path,handle,nr,ns,identity) bind(C,name='snrt_chimes_band_load')
       import
       character(c_char),intent(in)::path(*)
       type(c_ptr),intent(inout)::handle
       integer(c_int),intent(inout)::nr,ns
       real(c_double),intent(inout)::identity(32)
     end function
     subroutine chimes_band_free(handle) bind(C,name='snrt_chimes_band_free')
       import
       type(c_ptr),value::handle
     end subroutine
     integer(c_int) function chimes_band_reactions(handle,nr,mapping) bind(C,name='snrt_chimes_band_reactions')
       import
       type(c_ptr),value::handle
       integer(c_int),value::nr
       integer(c_int),intent(inout)::mapping(5,nr)
     end function
     integer(c_int) function chimes_band_moments(handle,nd,nr,number,energy,moments) &
          bind(C,name='snrt_chimes_band_moments')
       import
       type(c_ptr),value::handle
       integer(c_int),value::nd,nr
       real(c_double),intent(in)::number(nd,9),energy(nd,9)
       ! N*sigma, N*sigma*E, N*sigma*(E-binding); not accepted photons.
       real(c_double),intent(inout)::moments(nr,3,9)
     end function
     integer(c_int) function chimes_charge_supported() bind(C,name='snrt_chimes_charge_supported')
       import
     end function
     ! solid_charge is signed elementary charges in the SAME normalization
     ! as old (per H for cell; number density also allowed for reconcile).
     ! Frozen solid charge here is not a PAH charging evolution algorithm.
     integer(c_int) function chimes_reconcile_charged(elements,old,solid_charge,new) &
          bind(C,name='snrt_chimes_reconcile_charged')
       import
       real(c_double),intent(in)::elements(11),old(chimes_ns)
       real(c_double),value::solid_charge
       real(c_double),intent(out)::new(chimes_ns)
     end function
     integer(c_int) function chimes_cell_charged(controls,elements,old,solid_charge,photons,t,new,next_photons) &
          bind(C,name='snrt_chimes_cell_charged')
       import
       real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns),photons(chimes_ng)
       real(c_double),value::solid_charge
       real(c_double),intent(out)::t,new(chimes_ns),next_photons(chimes_ng)
     end function
     integer(c_int) function chimes_secondary_partition(energy,state,fractions) &
          bind(C,name='snrt_chimes_secondary_partition')
       import
       real(c_double),value::energy
       real(c_double),intent(in)::state(chimes_ns)
       real(c_double),intent(out)::fractions(5)
     end function
     integer(c_int) function chimes_group_binding(n,edges,means) bind(C,name='snrt_chimes_group_binding')
       import
       integer(c_int),value::n
       real(c_double),intent(in)::edges(*),means(*)
     end function
     integer(c_int) function chimes_nuclear_sums(state,elements) bind(C,name='snrt_chimes_nuclear_sums')
       import
       real(c_double),intent(in)::state(chimes_ns)
       real(c_double),intent(out)::elements(11)
     end function
     real(c_double) function chimes_boltzmann() bind(C,name='snrt_chimes_boltzmann')
       import
     end function
     integer(c_int) function chimes_identity(values) bind(C,name='snrt_chimes_identity')
       import
       real(c_double),intent(out)::values(320)
     end function
     integer(c_int) function chimes_reconcile(elements,old,new) bind(C,name='snrt_chimes_reconcile')
       import
       real(c_double),intent(in)::elements(11),old(chimes_ns)
       real(c_double),intent(out)::new(chimes_ns)
     end function
     integer(c_int) function chimes_locked(abundance,elements) bind(C,name='snrt_chimes_locked')
       import
       real(c_double),intent(in)::abundance(chimes_ns)
       real(c_double),intent(out)::elements(11)
     end function
     integer(c_int) function chimes_initialize(path,n,paths,stride) bind(C,name='snrt_chimes_initialize')
       import
       character(c_char),intent(in)::path(*),paths(*)
       integer(c_int),value::n,stride
     end function
     integer(c_int) function chimes_neutral(elements,abundance) bind(C,name='snrt_chimes_neutral')
       import
       real(c_double),intent(in)::elements(11)
       real(c_double),intent(out)::abundance(chimes_ns)
     end function
     integer(c_int) function chimes_budget(abundance,elements,charge) bind(C,name='snrt_chimes_budget')
       import
       real(c_double),intent(in)::abundance(chimes_ns)
       real(c_double),intent(out)::elements(11),charge
     end function
     integer(c_int) function chimes_cell(controls,elements,old,photons,t,new,next_photons) bind(C,name='snrt_chimes_cell')
       import
       real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns),photons(chimes_ng)
       real(c_double),intent(out)::t,new(chimes_ns),next_photons(chimes_ng)
     end function
  end interface
contains
  integer function chimes_round_subnormal_survivors(scale,incoming_energy,number,energy) result(status)
    ! FP32 storage of photon N cannot retain subnormal survivors with an
    ! independent FP64 energy. Round ONLY this storage tail as a paired N/E
    ! packet, preserving its mean energy. Bound the absolute energy change
    ! to FP64 roundoff of this cell's incoming radiation, not a physical floor.
    real(c_double),intent(in)::scale,incoming_energy
    real(c_double),intent(inout)::number(:,:),energy(:,:)
    real(c_double)::n(size(number,1),size(number,2)),e(size(number,1),size(number,2)),q,change
    integer::i,g
    status=1
    if(any(shape(number)/=shape(energy)).or..not.ieee_is_finite(scale).or.scale<=0)return
    if(.not.ieee_is_finite(incoming_energy).or.incoming_energy<0)return
    if(any(.not.ieee_is_finite(number)).or.any(.not.ieee_is_finite(energy)))return
    if(any(number<0).or.any(energy<0).or.any((number==0).neqv.(energy==0)))return
    n=number;e=energy;change=0
    do g=1,size(n,2)
       do i=1,size(n,1)
          if(n(i,g)<=0.or.n(i,g)/scale>=real(tiny(0.0_c_float),c_double))cycle
          q=real(real(n(i,g)/scale,c_float),c_double)*scale
          e(i,g)=e(i,g)*(q/n(i,g));n(i,g)=q
          change=change+abs(e(i,g)-energy(i,g))
       enddo
    enddo
    if(change>64*epsilon(1d0)*incoming_energy)return
    number=n;energy=e;status=0
  end function

  integer function chimes_cell_band_cold_molecular(handle,mol,nd,controls,elements,old,alpha,number,energy, &
       temperature,new,next_number,next_energy,ledger,grain_number,grain_energy,transition,event_info, &
       phase_alpha,directions,phase_energy,phase_moment) result(status)
    ! Cold split, optionally extended by the explicit energy-aware transition.
    ! Molecular carriers convert only through the budgeted atomization helper;
    ! this is NOT permission to extrapolate molecular rates or raise Tmol_K.
    ! Shielding/pump factors are frozen at entry (first-order coefficient split).
    type(c_ptr),intent(in)::handle,mol
    integer,intent(in)::nd
    real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns),alpha(128,9),number(nd,9),energy(nd,9)
    real(c_double),intent(inout)::temperature,new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(11)
    real(c_double),optional,intent(inout)::grain_number(9),grain_energy(9)
    logical,optional,intent(in)::transition
    real(c_double),optional,intent(inout)::event_info(2) ! event count, cost eV/cm3
    real(c_double),optional,intent(in)::phase_alpha(128,9,4),directions(3,nd)
    real(c_double),optional,intent(inout)::phase_energy(9,4),phase_moment(3,4)
    real(c_double)::phase_e(9,4),phase_p(3,4)
    real(c_double)::photo(chimes_ns),chem(chimes_ns),pn(nd,9),pe(nd,9),budget(11),ctl(9),factors(3)
    real(c_double)::measured(11),charge,temp,gn(9),ge(9),post_photo,tmax
    real(c_double)::incoming(chimes_ns),projected(chimes_ns),tin,elapsed,root_time,cost,entry_cost,events(2)
    logical::general,atomic_remainder
    real(c_double),parameter::ev_erg=1.602176634d-12
    status=2
    if(present(phase_alpha).neqv.present(directions))return
    if(present(phase_alpha).neqv.present(phase_energy))return
    if(present(phase_alpha).neqv.present(phase_moment))return
    tmax=chimes_molecular_temperature_max()
    general=.false.;if(present(transition))general=transition
    if(general.and.chimes_transition_supported()/=1)return
    if(nd<1.or.nd>720.or..not.c_associated(handle).or..not.c_associated(mol))return
    if(any(.not.ieee_is_finite(controls)).or.any(controls<0).or.controls(1)<=0)return
    if(controls(2)<10.or.controls(2)>merge(1d9,tmax,general).or.controls(3)<1.or.controls(3)>1d4)return
    if(controls(9)<=0.or.controls(9)>1.or.any(.not.ieee_is_finite(elements)).or.any(elements<0))return
    if(2.99792458d10*controls(9)*controls(4)>controls(5)*(1+1d-12))return
    status=chimes_budget(old,measured,charge)
    if(status/=0)return
    status=2
    if(any(abs(measured-elements)>1d-8*max(elements,1d-20)).or.abs(charge)>1d-10)return
    incoming=old;tin=controls(2);entry_cost=0;events=0;atomic_remainder=.false.
    if(general.and.tin>=tmax.and.controls(4)>0)then
       status=chimes_atomize(tin,incoming,projected,temp,cost)
       if(status/=0)return
       events=[merge(1d0,0d0,any(incoming(138:157)>0)),cost*controls(1)]
       incoming=projected;tin=temp;entry_cost=cost*controls(1);atomic_remainder=.true.
    endif
    if(atomic_remainder.or.(general.and.tin>=tmax))then
       factors=[1d0,1d0,0d0];status=0
    else
       status=chimes_molecular_factors(tin,controls(1),controls(5),incoming,factors)
    endif
    if(status/=0)return
    if(present(phase_alpha))then
       status=chimes_band_photo_molecular_phases(handle,mol,nd,controls(1),controls(4),2.99792458d10*controls(9), &
            alpha,factors(1:2),factors(3),incoming,number,energy,photo,pn,pe,budget(1:10),gn,ge, &
            phase_alpha,directions,phase_e,phase_p)
    else
       status=chimes_band_photo_molecular_groups(handle,mol,nd,controls(1),controls(4),2.99792458d10*controls(9), &
            alpha,factors(1:2),factors(3),incoming,number,energy,photo,pn,pe,budget(1:10),gn,ge)
    endif
    if(status/=0)return
    if(controls(4)==0)then
       temperature=controls(2);new=old;next_number=number;next_energy=energy;ledger=0
       if(present(grain_number))grain_number=0
       if(present(grain_energy))grain_energy=0
       if(present(event_info))event_info=0
       if(present(phase_energy))phase_energy=0
       if(present(phase_moment))phase_moment=0
       return
    endif
    status=8
    if(sum(photo)<=0)return
    ctl=controls
    ctl(2)=(tin*sum(incoming)+budget(1)*ev_erg/(1.5d0*chimes_boltzmann()*controls(1)))/sum(photo)
    if(.not.ieee_is_finite(ctl(2)).or.ctl(2)<10.or.ctl(2)>merge(1d9,tmax,general))return
    post_photo=1.5d0*chimes_boltzmann()*controls(1)*sum(photo)*ctl(2)
    if(general)then
       if(ctl(2)>=tmax)then
          status=chimes_atomize(ctl(2),photo,projected,temp,cost)
          if(status/=0)return
          events=events+[merge(1d0,0d0,any(photo(138:157)>0)),cost*controls(1)]
          photo=projected;ctl(2)=temp;atomic_remainder=.true.
       endif
       status=chimes_cell_transition_dark(ctl,elements,photo,merge(1,0,atomic_remainder),temp,chem,elapsed)
       if(status==51)then
          ! A CVODE root, not a failed Newton trial or bisection on status50.
          status=8
          if(abs(temp-tmax)>1d-7*tmax.or.elapsed<0.or.elapsed>ctl(4))return
          root_time=elapsed
          status=chimes_atomize(temp,chem,projected,tin,cost)
          if(status/=0)return
          events=events+[merge(1d0,0d0,any(chem(138:157)>0)),cost*controls(1)]
          ctl(2)=tin;ctl(4)=max(0d0,ctl(4)-root_time)
          status=chimes_cell_transition_dark(ctl,elements,projected,1,temp,chem,elapsed)
       endif
    else
       status=chimes_cell_cold_dark(ctl,elements,photo,temp,chem)
    endif
    if(status/=0)return
    status=8
    if(temp<10.or.temp>merge(1d9,tmax,general))return
    budget(11)=(1.5d0*chimes_boltzmann()*controls(1)*sum(chem)*temp-post_photo)/ev_erg-entry_cost
    if(.not.ieee_is_finite(budget(11)))return
    temperature=temp;new=chem;next_number=pn;next_energy=pe;ledger=budget;status=0
    if(present(grain_number))grain_number=gn
    if(present(grain_energy))grain_energy=ge
    if(present(event_info))event_info=events
    if(present(phase_energy))phase_energy=phase_e
    if(present(phase_moment))phase_moment=phase_p
  end function

  integer function chimes_cell_band_hot_atomic(handle,nd,controls,elements,old,solid_q,number,energy, &
       temperature,new,next_number,next_energy,ledger) result(status)
    ! First-order photo -> nonradiative chemistry split on the existing
    ! CHIMES hot atomic network. NOT a live cold/molecular/dust receiver.
    ! ledger(9) is signed dark-step gas thermal change, not escaped radiation.
    type(c_ptr),intent(in)::handle
    integer,intent(in)::nd
    real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns),solid_q,number(nd,9),energy(nd,9)
    real(c_double),intent(inout)::temperature,new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(9)
    real(c_double)::photo(chimes_ns),chem(chimes_ns),pn(nd,9),pe(nd,9),budget(9),ctl(9)
    real(c_double)::measured(11),charge,zero_photons(9),dark_photons(9),before,after,temp,post_photo_energy
    real(c_double),parameter::ev_erg=1.602176634d-12,tmol=1d5
    status=2
    if(nd<1.or.nd>720.or..not.c_associated(handle))return
    if(any(.not.ieee_is_finite(controls)).or.any(controls<0).or.controls(1)<=0)return
    if(.not.ieee_is_finite(solid_q).or.any(.not.ieee_is_finite(elements)))return
    if(any(elements<0).or.controls(3)<1.or.controls(3)>1d4.or.controls(9)<=0.or.controls(9)>1)return
    ! The upstream hot network removes all molecular states (not H-/C-/O-).
    ! Never silently discard those carriers between the two operators.
    if(controls(2)<=tmol.or.controls(2)>1d9.or.controls(6)/=0.or.solid_q/=0)return
    if(any(old(138:157)/=0))return
    status=chimes_budget(old,measured,charge)
    if(status/=0)return
    status=2
    if(any(abs(measured-elements)>1d-8*max(elements,1d-20)).or.abs(charge)>1d-10)return
    status=chimes_band_photo_step(handle,nd,controls(1),controls(4),2.99792458d10*controls(9), &
         old,number,energy,photo,pn,pe,budget(1:8))
    if(status/=0)return
    if(controls(4)==0)then
       temperature=controls(2);new=old;next_number=number;next_energy=energy;ledger=0
       return
    endif
    before=sum(old);after=sum(photo)
    status=8
    if(after<=0)return
    ctl=controls
    ctl(2)=(controls(2)*before+budget(1)*ev_erg/(1.5d0*chimes_boltzmann()*controls(1)))/after
    if(.not.ieee_is_finite(ctl(2)).or.ctl(2)<=tmol.or.ctl(2)>1d9)return
    if(any(photo(138:157)/=0))return
    post_photo_energy=1.5d0*chimes_boltzmann()*controls(1)*after*ctl(2)
    zero_photons=0
    status=chimes_cell_charged(ctl,elements,photo,solid_q,zero_photons,temp,chem,dark_photons)
    if(status/=0)return
    status=8
    if(temp<=tmol.or.any(chem(138:157)/=0).or.any(dark_photons/=0))return
    budget(9)=(1.5d0*chimes_boltzmann()*controls(1)*sum(chem)*temp-post_photo_energy)/ev_erg
    if(.not.ieee_is_finite(budget(9)))return
    temperature=temp;new=chem;next_number=pn;next_energy=pe;ledger=budget;status=0
  end function

  integer(c_int) function fs_grid(n,energies,xi,values) bind(C,name='snrt_chimes_fs_grid')
    integer(c_int),value::n
    real(c_double),intent(out)::energies(n),xi(14),values(6,14,n)
    integer::status
    fs_grid=1
    if(n/=snrt_secondary_nenergy)return
    call snrt_secondary_raw_grid(energies,xi,values,status)
    fs_grid=status
  end function
  integer(c_int) function fs_samples(energy,values,xi) bind(C,name='snrt_chimes_fs_samples')
    real(c_double),value::energy
    real(c_double),intent(out)::values(6,14),xi(14)
    integer::status
    call snrt_secondary_energy_slice(energy,values,xi,status)
    fs_samples=status
  end function
  integer(c_int) function fs_ready() bind(C,name='snrt_chimes_fs_ready')
    fs_ready=merge(1,0,snrt_secondary_tables_loaded)
  end function
  integer(c_int) function fs_fractions(energy,xi,fractions) bind(C,name='snrt_chimes_fs_fractions')
    real(c_double),value::energy,xi
    real(c_double),intent(out)::fractions(5)
    integer::status
    call snrt_secondary_fractions(energy,xi,fractions(1),fractions(2), &
         fractions(3),fractions(4),fractions(5),status)
    fs_fractions=status
  end function
end module
