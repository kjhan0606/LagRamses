module pm_commons
  implicit none
  integer,parameter,private::dp=kind(1.0d0)
  real(dp)::xp(8,3)
end module pm_commons

module amr_commons
  implicit none
  integer,parameter::dp=kind(1.0d0)
  real(dp)::xg(1,3),boxlen
  integer::icoarse_min,icoarse_max,jcoarse_min,jcoarse_max
  integer::kcoarse_min,kcoarse_max
end module amr_commons

program test_sidm_cell_classification
  use pm_commons
  use amr_commons
  implicit none
  integer,external::cell_index_from_part
  integer::got,ix,iy,iz,ip
  real(dp)::scale,skip(3),dx

  boxlen=2.5d0
  icoarse_min=3
  icoarse_max=6
  jcoarse_min=1
  jcoarse_max=4
  kcoarse_min=0
  kcoarse_max=3
  xg(1,:)=(/3.5d0,1.5d0,0.5d0/)
  scale=boxlen/dble(icoarse_max-icoarse_min+1)
  dx=0.5d0**16
  skip=(/dble(icoarse_min),dble(jcoarse_min),dble(kcoarse_min)/)

  ! Exercise all eight level-16 child regions with non-unit physical
  ! scaling and nonzero local mesh offsets.
  do iz=0,1
     do iy=0,1
        do ix=0,1
           ip=1+ix+2*iy+4*iz
           xp(ip,:)=(xg(1,:)+(/dble(ix)-0.5d0, &
                dble(iy)-0.5d0,dble(iz)-0.5d0/)*dx-skip)*scale
           got=cell_index_from_part(ip,1,16)
           if(got/=ip)error stop 'one of the eight child classifications is wrong'
        end do
     end do
  end do

  ! The RAMSES particle tree computes a nonnegative integer child offset, so
  ! an exact split-plane coordinate belongs to the upper child.
  xp(8,:)=(xg(1,:)-skip)*scale
  got=cell_index_from_part(8,1,16)
  if(got/=8)error stop 'midplane tie must follow the particle-tree convention'

  write(*,'(A)')'SIDM production child classification: PASS'
end program test_sidm_cell_classification
