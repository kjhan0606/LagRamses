program stellar_source_smoke
  use amr_parameters, only: dp
  use snrt_spectral_contract
  use snrt_stellar_source
  use snrt_agn_locator, only: snrt_agn_find_local_leaf
  use amr_commons, only: levelmin,nlevelmax,ngridmax,ncoarse,amr_block_size,boxlen, &
       icoarse_min,icoarse_max,jcoarse_min,kcoarse_min,active,xg,son
#include "amr_index.h"
  implicit none
  integer::ierr,cell,level
  real(dp)::whole(9),left(9),right(9),expected
  call snrt_spectral_contract_load_from_environment(ierr)
  if(ierr/=0)stop 1
  call stellar_sed_load(ierr)
  if(ierr/=0.or..not.stellar_sed_enabled)stop 2
  call stellar_photon_interval(0d0,2d0,0.1d0,3d0,whole,ierr)
  if(ierr/=0)stop 3
  expected=4.5d40*31557600d6*(1.5d0+2d0-1d0/99999d0)
  if(abs(whole(5)/expected-1d0)>1d-14)stop 4
  call stellar_photon_interval(0d0,0.75d0,0.1d0,3d0,left,ierr)
  if(ierr/=0)stop 5
  call stellar_photon_interval(0.75d0,2d0,0.1d0,3d0,right,ierr)
  if(ierr/=0)stop 6
  if(maxval(abs(whole-left-right))/maxval(whole)>1d-14)stop 7
  call stellar_photon_interval(-1d0,0d0,0.1d0,3d0,left,ierr)
  if(ierr/=0.or.any(left/=0d0))stop 8
  call stellar_photon_interval(0d0,100001d0,0.1d0,3d0,left,ierr)
  if(ierr==0)stop 9
  call stellar_photon_interval(0d0,1d0,0.21d0,3d0,left,ierr)
  if(ierr==0)stop 10
  write(*,*)'PASS native stellar photon integration, step splitting and bounds'
  levelmin=2;nlevelmax=2;ngridmax=2;ncoarse=1;amr_block_size=1;boxlen=1d0
  icoarse_min=0;icoarse_max=0;jcoarse_min=0;kcoarse_min=0
  allocate(active(2),xg(2,3),son(17))
  active(2)%ngrid=2
  allocate(active(2)%igrid(2));active(2)%igrid=[1,2]
  xg(1,:)=[0.25d0,0.25d0,0.25d0];xg(2,:)=[0.75d0,0.25d0,0.25d0];son=0
  call snrt_agn_find_local_leaf([0.375d0,0.125d0,0.125d0],cell,level)
  if(cell/=ICELL_OF(1,2).or.level/=2)stop 11
  call snrt_agn_find_local_leaf([0.375d0,0.125d0,0.125d0],cell,level,2,2)
  if(cell/=0)stop 12 ! No negative-fraction truncation into the right grid.
  call snrt_agn_find_local_leaf([-0.125d0,0.125d0,0.125d0],cell,level)
  if(cell/=0)stop 13
  call snrt_agn_find_local_leaf([0.5d0,0.125d0,0.125d0],cell,level)
  if(cell/=ICELL_OF(2,1).or.level/=2)stop 14 ! Half-open shared face.
  write(*,*)'PASS native source cell layout, half-open ownership and constant-time stellar hint'
end program
