! Native callback fixture only: tiny owned-leaf topology, no MPI or RAMSES.
module amr_parameters
  implicit none
  integer,parameter::dp=kind(1d0),ndim=3,nvector=2,twondim=6,twotondim=8,amr_block_size=1
end module
module hydro_parameters
  use amr_parameters
  implicit none
  integer,parameter::nvar=21,iu1=-1,iu2=4,ju1=-1,ju2=4,ku1=-1,ku2=4
  integer::imetal=8,ichem=9,iradioactive=20
end module
module amr_commons
  use amr_parameters
  implicit none
  type::active_t
     integer::ngrid=0,igrid(2)=0
  end type
  type(active_t)::active(2)
  integer::ncoarse=0,myid=1,son(32)=0
  real(dp)::dtnew(2)=0,aexp=1
  integer::fixture_virtual_calls=0
  logical::fixture_stopped=.false.
end module
module hydro_commons
  use hydro_parameters
  implicit none
  real(dp)::uold(32,nvar)=0,unew(32,nvar)=0
end module
subroutine units(sl,st,sd,sv,snh,st2)
  use stellar_radioactive_decay, only: radioactive_half_life_s
  implicit none
  real(8),intent(out)::sl,st,sd,sv,snh,st2
  sl=1;st=radioactive_half_life_s(1);sd=1;sv=1;snh=1;st2=1
end subroutine
subroutine make_virtual_fine_dp(values,ilevel)
  use amr_commons, only: fixture_virtual_calls
  implicit none
  real(8)::values(32)
  integer::ilevel
  fixture_virtual_calls=fixture_virtual_calls+1
end subroutine
subroutine clean_stop
  use amr_commons, only: fixture_stopped
  implicit none
  fixture_stopped=.true.
end subroutine
