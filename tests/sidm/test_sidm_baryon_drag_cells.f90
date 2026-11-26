! Minimal single-process fixture for the actual sidm_baryon_drag routine.
! One gas-filled child is selected at a time; the other children are empty.
module amr_parameters
  implicit none
  integer,parameter :: twotondim=8,amr_block_size=64
end module amr_parameters

module pm_commons
  implicit none
  integer,parameter :: dp_pm=kind(1.0d0)
  integer(kind=1),parameter :: PTYPE_DM=0,PTYPE_STAR=1,PTYPE_SINK=2
  real(dp_pm) :: xp(8,3),vp(8,3),mp(8)
  integer(kind=8) :: idp(8)
  integer(kind=1) :: ptypep(8)
  integer :: numbp(1),headp(1),nextp(8)
end module pm_commons

module amr_commons
  use amr_parameters
  implicit none
  integer,parameter :: dp=kind(1.0d0)
  integer,parameter :: MPI_INTEGER=1,MPI_DOUBLE_PRECISION=2,MPI_SUM=1,MPI_MAX=2,MPI_COMM_WORLD=0
  integer :: icoarse_min,icoarse_max,jcoarse_min,kcoarse_min,ncoarse,myid
  integer :: numbtot(1,3),headl(1,3),numbl(1,3),next(1),son(512)
  real(dp) :: xg(1,3),boxlen,dtnew(3),sidm_baryon_sigma,sidm_baryon_power
  logical :: sidm_baryon,hydro
end module amr_commons

module hydro_commons
  use amr_commons,only:dp
  implicit none
  integer,parameter :: nvar=5
  real(dp) :: uold(512,5),gamma
end module hydro_commons

program test_sidm_baryon_drag_cells
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
  integer :: child,ip,axis,icell,sign_bit
  real(dp) :: dx,scale,volume,dm_mass,expected_v,expected_momentum,expected_energy
  real(dp),parameter :: tol=1.0d-13
  real(dp) :: skip(3)

  boxlen=2.5d0
  icoarse_min=3; icoarse_max=6; jcoarse_min=1; kcoarse_min=0
  ncoarse=0; myid=1
  skip=(/3.0d0,1.0d0,0.0d0/)
  xg(1,:)=(/3.5d0,1.5d0,0.5d0/)
  scale=boxlen/dble(icoarse_max-icoarse_min+1)
  dx=0.5d0**3
  volume=(dx*scale)**3
  dm_mass=volume*0.1d0
  numbtot=0; numbtot(1,3)=1
  headl=0; headl(1,3)=1
  numbl=0; numbl(1,3)=1
  next=0; son=0
  numbp=8; headp=1
  do ip=1,8
     nextp(ip)=ip+1
     idp(ip)=ip
     ptypep(ip)=PTYPE_DM
     do axis=1,3
        sign_bit=ibits(ip-1,axis-1,1)
        xp(ip,axis)=(xg(1,axis)+(2*sign_bit-1)*dx*0.25d0-skip(axis))*scale
     enddo
  enddo
  nextp(8)=0
  sidm_baryon=.true.; hydro=.true.
  sidm_baryon_sigma=0.01d0; sidm_baryon_power=0.0d0
  dtnew=0.0d0; dtnew(3)=0.1d0
  gamma=5.0d0/3.0d0
  expected_v=0.999d0
  expected_momentum=-dm_mass*(expected_v-1.0d0)/volume
  expected_energy=-0.5d0*dm_mass*(expected_v**2-1.0d0)/volume

  do child=1,8
     vp=0.0d0; vp(:,1)=1.0d0; mp=dm_mass
     uold=0.0d0
     icell=1+(child-1)*amr_block_size
     uold(icell,1)=1.0d0
     call sidm_baryon_drag(3)
     do ip=1,8
        if(ip==child) then
           if(abs(vp(ip,1)-expected_v)>tol) error stop 'wrong particle drag kick'
        else
           if(abs(vp(ip,1)-1.0d0)>tol) error stop 'drag leaked to another child'
        endif
     enddo
     if(abs(uold(icell,2)-expected_momentum)>tol) error stop 'gas momentum mismatch'
     if(abs(uold(icell,5)-expected_energy)>tol) error stop 'gas energy mismatch'
  enddo
  print '(A)', 'SIDM baryon drag child selection and exchange: PASS'
end program test_sidm_baryon_drag_cells

subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  implicit none
  real(kind(1.0d0)),intent(out) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  scale_l=1.0d0; scale_t=1.0d0; scale_d=1.0d0
  scale_v=1.0d0; scale_nH=1.0d0; scale_T2=1.0d0
end subroutine units
