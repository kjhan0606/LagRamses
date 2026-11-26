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
  integer::got,ix,iy,iz,ip,axis,expected,ymax,zmax,setup
  real(dp)::scale,skip(3),dx,xmesh(3)

  ymax=0
  zmax=0
#if NDIM>1
  ymax=1
#endif
#if NDIM==3
  zmax=1
#endif
  dx=0.5d0**16

  do setup=1,2
     ! Exercise both the default map and non-unit scale with local offsets.
     if(setup==1)then
        boxlen=4d0
        icoarse_min=0; icoarse_max=3
        jcoarse_min=0; jcoarse_max=3
        kcoarse_min=0; kcoarse_max=3
        xg(1,:)=(/0.5d0,0.5d0,0.5d0/)
     else
        boxlen=2.5d0
        icoarse_min=3; icoarse_max=6
        jcoarse_min=1; jcoarse_max=4
        kcoarse_min=0; kcoarse_max=3
        xg(1,:)=(/3.5d0,1.5d0,0.5d0/)
     endif
     scale=boxlen/dble(icoarse_max-icoarse_min+1)
     skip=(/dble(icoarse_min),dble(jcoarse_min),dble(kcoarse_min)/)

     ! Every child in the supported dimensionality, safely off the split.
     do iz=0,zmax
        do iy=0,ymax
           do ix=0,1
              ip=1+ix+2*iy+4*iz
              xmesh=xg(1,:)+(/dble(ix)-0.5d0, &
                   dble(iy)-0.5d0,dble(iz)-0.5d0/)*dx
              xp(ip,:)=(xmesh-skip)*scale
              got=cell_index_from_part(ip,1,16)
              if(got/=ip)error stop 'child classification is wrong'
           enddo
        enddo
     enddo

     ! A particle exactly on one split belongs to its upper child while
     ! the remaining active dimensions stay on their lower side.
     do axis=1,NDIM
        xmesh=xg(1,:)-0.5d0*dx
        xmesh(axis)=xg(1,axis)
        xp(1,:)=(xmesh-skip)*scale
        expected=1+2**(axis-1)
        got=cell_index_from_part(1,1,16)
        if(got/=expected)error stop 'single-axis split-plane tie is wrong'
     enddo

     xp(1,:)=(xg(1,:)-skip)*scale
     got=cell_index_from_part(1,1,16)
     if(got/=2**NDIM)error stop 'multiple split-plane tie is wrong'
  enddo

  write(*,'(A,I0,A)') 'SIDM production child classification NDIM=',NDIM,': PASS'
end program test_sidm_cell_classification
