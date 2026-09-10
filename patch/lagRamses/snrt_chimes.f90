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
  public::chimes_cell_band_hot_atomic
  interface
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
  integer function chimes_cell_band_cold_molecular(handle,mol,nd,controls,elements,old,alpha,number,energy, &
       temperature,new,next_number,next_energy,ledger,grain_number,grain_energy) result(status)
    ! Bounded cold-cell split; NOT permission to raise upstream Tmol_K.
    ! No species are deleted to force a temperature transition to succeed.
    ! Shielding/pump factors are frozen at entry (first-order coefficient split).
    type(c_ptr),intent(in)::handle,mol
    integer,intent(in)::nd
    real(c_double),intent(in)::controls(9),elements(11),old(chimes_ns),alpha(128,9),number(nd,9),energy(nd,9)
    real(c_double),intent(inout)::temperature,new(chimes_ns),next_number(nd,9),next_energy(nd,9),ledger(11)
    real(c_double),optional,intent(inout)::grain_number(9),grain_energy(9)
    real(c_double)::photo(chimes_ns),chem(chimes_ns),pn(nd,9),pe(nd,9),budget(11),ctl(9),factors(3)
    real(c_double)::measured(11),charge,temp,gn(9),ge(9),post_photo,tmax
    real(c_double),parameter::ev_erg=1.602176634d-12
    status=2
    tmax=chimes_molecular_temperature_max()
    if(nd<1.or.nd>720.or..not.c_associated(handle).or..not.c_associated(mol))return
    if(any(.not.ieee_is_finite(controls)).or.any(controls<0).or.controls(1)<=0)return
    if(controls(2)<10.or.controls(2)>tmax.or.controls(3)<1.or.controls(3)>1d4)return
    if(controls(9)<=0.or.controls(9)>1.or.any(.not.ieee_is_finite(elements)).or.any(elements<0))return
    if(2.99792458d10*controls(9)*controls(4)>controls(5)*(1+1d-12))return
    status=chimes_budget(old,measured,charge)
    if(status/=0)return
    status=2
    if(any(abs(measured-elements)>1d-8*max(elements,1d-20)).or.abs(charge)>1d-10)return
    status=chimes_molecular_factors(controls(2),controls(1),controls(5),old,factors)
    if(status/=0)return
    status=chimes_band_photo_molecular_groups(handle,mol,nd,controls(1),controls(4),2.99792458d10*controls(9), &
         alpha,factors(1:2),factors(3),old,number,energy,photo,pn,pe,budget(1:10),gn,ge)
    if(status/=0)return
    if(controls(4)==0)then
       temperature=controls(2);new=old;next_number=number;next_energy=energy;ledger=0
       if(present(grain_number))grain_number=0
       if(present(grain_energy))grain_energy=0
       return
    endif
    status=8
    if(sum(photo)<=0)return
    ctl=controls
    ctl(2)=(controls(2)*sum(old)+budget(1)*ev_erg/(1.5d0*chimes_boltzmann()*controls(1)))/sum(photo)
    if(.not.ieee_is_finite(ctl(2)).or.ctl(2)<10.or.ctl(2)>tmax)return
    post_photo=1.5d0*chimes_boltzmann()*controls(1)*sum(photo)*ctl(2)
    status=chimes_cell_cold_dark(ctl,elements,photo,temp,chem)
    if(status/=0)return
    status=8
    if(temp<10.or.temp>tmax)return
    budget(11)=(1.5d0*chimes_boltzmann()*controls(1)*sum(chem)*temp-post_photo)/ev_erg
    if(.not.ieee_is_finite(budget(11)))return
    temperature=temp;new=chem;next_number=pn;next_energy=pe;ledger=budget;status=0
    if(present(grain_number))grain_number=gn
    if(present(grain_energy))grain_energy=ge
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
