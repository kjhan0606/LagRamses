!=======================================================================
! FDM light-cone density output for lagRamses.
!
! Mirrors output_cone_hydro (light_cone.hydro2.f90) but projects the
! Fuzzy-Dark-Matter density |psi|^2 = psi_re^2 + psi_im^2 onto the light
! cone instead of CDM N-body particles (which do not exist under FDM).
!
! Per-leaf-cell payload written to getc(1:nhvar):
!   getc(1) = cell size dx (coarse units)
!   getc(2) = rho_fdm = psi_re^2 + psi_im^2
!   getc(3) = phi (gravitational potential)
!
! Velocity is written as zero (gvel=0). The physical FDM (Madelung)
! velocity v = (hbar/m/a^2) Im(grad psi / psi) requires a neighbour
! gradient gather and can be added later if needed.
!
! The cone-geometry helpers (perform_my_selection, compute_replica_box,
! init_cosmo_cone, coord_distance, myint, ...) are shared external
! subroutines defined in light_cone.hydro2.f90.
!=======================================================================
subroutine output_cone_fdm(obs)
  use amr_commons
  use pm_commons
  use poisson_commons
#include "amr_index.h"
  implicit none

#ifndef WITHOUTMPI
#include "mpif.h"
#endif

  integer::info,dummy_io,info2
  integer,parameter::tag=1118

  character(len=5) :: istep_str
  character(len=150) :: conedir, conecmd, conefile, infofile

  character(LEN=150)::fileloc
  character(LEN=5)::nchar
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: Lobserver(3)
  integer::igrid,idim,icpu,ilevel
  integer::i,j,k,l,m,ncell

  integer,dimension(:),allocatable::ind_grid
  integer(kind=8),dimension(:),allocatable::ii8,ii8out,ii8_out
  real(kind=8),dimension(:,:),allocatable::gpos,gvel,getc
  real(kind=8),dimension(:,:),allocatable::gposout,gvelout,getcout
  real(kind=8),dimension(:,:),allocatable::gpos_out,gvel_out,getc_out
  real(kind=8),dimension(:),allocatable::gzout,gz_out
  integer::nleaf,ncache,ibound,ngout,istart,iglun,ind,nhvar,ivar
  integer::end_tag,print_mark,mncell,tngout,elongated_axis_cone,obs,nprint
  real(kind=8) :: cpi(8,3),dx,coord_distance,Omega0,OmegaL,OmegaR,coverH0
  real(kind=8) :: dist1,dist2,lboxz(3),minboxr(3),maxboxr(3)
  real(kind=8) :: minboxr_cone(3),maxboxr_cone(3)
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::nfine,dist10,dist20,daexp,z10,z20

  logical::opened
  opened=.false.
  end_tag=-1
  minboxr=1
  maxboxr=0
  tngout=0

  if(nstep_coarse.lt.2) return
  nfine=1.
  do ilevel=levelmin,nlevelmax
     if(numbtot(1,ilevel) .gt. 0) then
        nfine=nfine*nsubcycle(ilevel-1)
     endif
  enddo

  daexp=max((aexp-aexp_old2)/nfine,aexp-aexp_old_fine)

  z20=1./aexp_old2-1.
  z10=1./aexp-1.

  z2=1./(aexp_old2-daexp)-1.
  z1=1./(aexp+daexp)-1.

  if(z1<0.) z1=0.
  if(z2.gt.zmax_cone) return
  if(abs(z2-z1)<1d-6) return

  ! FDM payload: dx, rho_fdm, phi
  nhvar=3
  om0in=omega_m
  omLin=omega_l
  hubin=h0/100.
  Lbox=boxlen_ini/hubin
  if(obs==1) then
     elongated_axis_cone=elongated_axis_cone1
     do idim=1,3
        Lobserver(idim)=observer_cone1(idim)*Lbox
        minboxr_cone(idim)=minboxr_cone1(idim)
        maxboxr_cone(idim)=maxboxr_cone1(idim)
     enddo
  else
     elongated_axis_cone=elongated_axis_cone2
     do idim=1,3
        Lobserver(idim)=observer_cone2(idim)*Lbox
        minboxr_cone(idim)=minboxr_cone2(idim)
        maxboxr_cone(idim)=maxboxr_cone2(idim)
     enddo
  endif

  if(myid==1 .and. obs==1) write(*,*)'Computing and dumping FDM lightcone (|psi|^2 leaf cells)'
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)
  dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)
  dist10=coord_distance(z10,Omega0,OmegaL,OmegaR,coverH0)
  dist20=coord_distance(z20,Omega0,OmegaL,OmegaR,coverH0)
  if(myid==1 .and. obs==1)  write(*,*)'FDM lightcone redshifts',z10,z20
  if(myid==1 .and. obs==1)  write(*,*)'Distance (Mpc)',dist10,dist20

  iglun=3*ncpu+myid+10

  call title(nstep_coarse, istep_str)
  if(obs==1) then
     conedir = "light_cone/cone_" // trim(istep_str) // "/observer1/"
  else
     conedir = "light_cone/cone_" // trim(istep_str) // "/observer2/"
  endif
  conecmd = "mkdir -p " // trim(conedir)
  if(.not.withoutmkdir) then
     if (myid==1) call execute_command_line(conecmd,wait=.true.)
  endif

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  call title(myid,nchar)

  conefile = trim(conedir)//'cone_fdm_'//trim(istep_str)//'.out'
  fileloc=TRIM(conefile)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZECONE>0) then
     if (mod(myid-1,IOGROUPSIZECONE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  mncell=65536
  allocate(gpos(1:ndim,1:mncell))
  allocate(gvel(1:ndim,1:mncell))
  allocate(getc(1:nhvar,1:mncell))
  allocate(ii8(1:mncell))
  allocate(gpos_out(1:mncell,1:ndim))
  allocate(gvel_out(1:mncell,1:ndim))
  allocate(getc_out(1:mncell,1:nhvar))
  allocate(ii8_out(1:mncell))
  allocate(gz_out(1:mncell))

  cpi(1,1:3)=(/0,0,0/)
  cpi(2,1:3)=(/1,0,0/)
  cpi(3,1:3)=(/0,1,0/)
  cpi(4,1:3)=(/1,1,0/)
  cpi(5,1:3)=(/0,0,1/)
  cpi(6,1:3)=(/1,0,1/)
  cpi(7,1:3)=(/0,1,1/)
  cpi(8,1:3)=(/1,1,1/)

  do idim=1,ndim
     maxboxr(idim)=maxboxr_cone(idim)*Lbox
     minboxr(idim)=minboxr_cone(idim)*Lbox
     lboxz(idim)=maxboxr(idim)-minboxr(idim)
     Lobserver(idim)=Lobserver(idim)-minboxr(idim)
  enddo
  nprint=0
  j=0
  print_mark=0
  do ibound=1,nboundary+ncpu
     if(ibound == myid) then
        do ilevel=1,nlevelmax
           dx=0.5**ilevel
           nleaf=0
           if(ibound<=ncpu) then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           endif
           if(ncache>0) then
              allocate(ind_grid(1:ncache))
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              enddo

              do ind=1,twotondim
                 do i=1,ncache
                    if(son(ICELL_OF(ind_grid(i),ind))==0) then
                       nleaf=nleaf+1
                    endif
                 enddo
              enddo
              j=0
              if(nleaf .gt. 0) then
                 do ind=1,twotondim
                    do i=1,ncache
                       if(son(ICELL_OF(ind_grid(i),ind))==0) then
                          j=j+1
                          do idim=1,ndim
                             gpos(idim,j)=(xg(ind_grid(i),idim)+(cpi(ind,idim)-0.5)*dx)*Lbox-minboxr(idim)
                             gvel(idim,j)=0.0d0
                          enddo
                          getc(1,j)=dx
                          getc(2,j)=psi_re(ICELL_OF(ind_grid(i),ind))**2 &
                               &   + psi_im(ICELL_OF(ind_grid(i),ind))**2
                          getc(3,j)=phi(ICELL_OF(ind_grid(i),ind))
                       endif
                       if(j==mncell) then
                          print_mark=1
                       endif
                       if(ind==twotondim .and. i==ncache) then
                          print_mark=1
                       endif

                       if (print_mark==1) then
                          ! Count number of leaf cells within the redshift range
                          call perform_my_selection(.true.,z1,z2, &
                               &                        om0in,omLin,hubin,lboxz, &
                               &                        Lobserver,elongated_axis_cone, &
                               &                        mncell, nhvar, &
                               &                        ii8,gpos,gvel,getc,j, &
                               &                        ii8out,gposout,gvelout,getcout,gzout,ngout,.false.)

                          if(ngout > 0) then
                             allocate(gposout(1:ndim,1:ngout))
                             allocate(gvelout(1:ndim,1:ngout))
                             allocate(getcout(1:nhvar,1:ngout))
                             allocate(ii8out(1:ngout))
                             allocate(gzout(1:ngout))
                             tngout=tngout+ngout

                             ! Perform actual selection
                             call perform_my_selection(.false.,z1,z2, &
                                  &                        om0in,omLin,hubin,lboxz, &
                                  &                        Lobserver,elongated_axis_cone, &
                                  &                        mncell, nhvar, &
                                  &                        ii8,gpos,gvel,getc,j, &
                                  &                        ii8out,gposout,gvelout,getcout,gzout,ngout,.false.)
                             k=0
                             do while(k .lt. ngout)
                                if(nprint+ngout-k .le. mncell) then
                                   ncell=nprint+ngout-k
                                else
                                   ncell=mncell
                                endif

                                do m=nprint+1,ncell
                                   k=k+1
                                   do idim=1,ndim
                                      gpos_out(m,idim)=gposout(idim,k)/Lbox
                                      gvel_out(m,idim)=gvelout(idim,k)
                                   enddo
                                   gz_out(m)=gzout(k)
                                   do l=1,nhvar
                                      getc_out(m,l)=getcout(l,k)
                                   enddo
                                enddo
                                nprint=ncell
                                if(nprint .eq. mncell) print_mark=2

                                if(print_mark==2) then
                                   if(.not.opened) then
                                      open(iglun,file=TRIM(fileloc),form='unformatted')
                                      rewind(iglun)
                                      write(iglun)ncpu
                                      write(iglun)ndim
                                      write(iglun)nhvar
                                      write(iglun)boxlen_ini
                                      write(iglun)aexp_old2-daexp
                                      write(iglun)aexp+daexp
                                      write(iglun)dist2/Lbox
                                      write(iglun)dist1/Lbox
                                      write(iglun)aexp_old2
                                      write(iglun)aexp
                                      write(iglun)dist20/Lbox
                                      write(iglun)dist10/Lbox
                                      opened=.true.
                                   endif
                                   write(iglun)nprint
                                   do idim=1,ndim
                                      write(iglun)gpos_out(1:nprint,idim)
                                   end do
                                   do idim=1,ndim
                                      write(iglun)gvel_out(1:nprint,idim)
                                   end do
                                   write(iglun)gz_out(1:nprint)
                                   do l=1,nhvar
                                      write(iglun)getc_out(1:nprint,l)
                                   end do
                                   nprint=0
                                   print_mark=0
                                endif
                             enddo
                          endif
                          if(allocated(gposout)) deallocate(gposout)
                          if(allocated(gvelout)) deallocate(gvelout)
                          if(allocated(getcout)) deallocate(getcout)
                          if(allocated(ii8out)) deallocate(ii8out)
                          if(allocated(gzout)) deallocate(gzout)

                          if(j .eq. mncell) j=0
                          print_mark=0
                       endif
                    enddo
                 enddo
              endif
              deallocate(ind_grid)
           endif
        enddo
        if(nprint .gt. 0) then
           if(.not.opened) then
              open(iglun,file=TRIM(fileloc),form='unformatted')
              rewind(iglun)
              write(iglun)ncpu
              write(iglun)ndim
              write(iglun)nhvar
              write(iglun)boxlen_ini
              write(iglun)aexp_old2-daexp
              write(iglun)aexp+daexp
              write(iglun)dist2/Lbox
              write(iglun)dist1/Lbox
              write(iglun)aexp_old2
              write(iglun)aexp
              write(iglun)dist20/Lbox
              write(iglun)dist10/Lbox
              opened=.true.
           endif
           write(iglun)nprint
           do idim=1,ndim
              write(iglun)gpos_out(1:nprint,idim)
           end do
           do idim=1,ndim
              write(iglun)gvel_out(1:nprint,idim)
           end do
           write(iglun)gz_out(1:nprint)
           do k=1,nhvar
              write(iglun)getc_out(1:nprint,k)
           end do
           nprint=0
        endif
     endif
  enddo
  if(opened)write(iglun)end_tag
  if(opened)close(iglun)

  deallocate(gpos)
  deallocate(gvel)
  deallocate(getc)
  deallocate(ii8)
  deallocate(gpos_out)
  deallocate(gvel_out)
  deallocate(getc_out)
  deallocate(gz_out)
  deallocate(ii8_out)

  if (tngout>0) then
     open(iglun,file=TRIM(fileloc)//".txt",form='formatted')
     rewind(iglun)
     write(iglun,*) ncpu
     write(iglun,*) tngout
     close(iglun)
  endif

  if (myid == 1) then
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     infofile = trim(conedir)//'info_fdm_'//trim(istep_str)
     open(iglun,file=TRIM(infofile)//".txt",form='formatted')
     rewind(iglun)
     write(iglun,'("a-da          =",E23.15)') aexp_old2-daexp
     write(iglun,'("a             =",E23.15)') aexp+daexp
     write(iglun,'("Dist(a-da)    =",E23.15)') dist2/Lbox
     write(iglun,'("Dist(a)       =",E23.15)') dist1/Lbox
     write(iglun,'("a0-da0        =",E23.15)') aexp_old2
     write(iglun,'("a0            =",E23.15)') aexp
     write(iglun,'("Dist(a0-da0)  =",E23.15)') dist20/Lbox
     write(iglun,'("Dist(a0)      =",E23.15)') dist10/Lbox
     write(iglun,'("h0            =",E23.15)') hubin
     write(iglun,'("unit_l        =",E23.15)') scale_l
     write(iglun,'("unit_d        =",E23.15)') scale_d
     write(iglun,'("unit_t        =",E23.15)') scale_t
     write(iglun,'("unit_v        =",E23.15)') scale_v
     write(iglun,'("unit_nH       =",E23.15)') scale_nH
     write(iglun,'("unit_T2       =",E23.15)') scale_T2
     write(iglun,'("Lbox (cMpc/h) =",E23.15)') boxlen_ini
     close(iglun)
  endif

#ifndef WITHOUTMPI
  if(IOGROUPSIZECONE>0) then
     if(mod(myid,IOGROUPSIZECONE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

  if((opened.and.(tngout==0)).or.((.not.opened).and.(tngout>0))) then
     write(*,*)'Error in output_cone_fdm'
     write(*,*)'tngout=',tngout,'opened=',opened
     stop
  endif

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

end subroutine output_cone_fdm
