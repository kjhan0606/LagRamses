subroutine file_descriptor_hydro(filename)
  use amr_commons
  use hydro_commons
#ifdef DUST_LIVE
  use dust_mass_physics, only: dust_pah_enabled,dust_pah_nstate
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif

  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ivar,ilun

  if(verbose)write(*,*)'Entering file_descriptor_hydro'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
#ifdef SOLVERmhd
  ! HDF5 MHD stores raw conserved fields, not the legacy binary primitives.
  write(ilun,'("nvar        =",I11)')nvar_all
  write(ilun,'(A)')'variable #1: density'
  write(ilun,'(A)')'variable #2: momentum_x'
  write(ilun,'(A)')'variable #3: momentum_y'
  write(ilun,'(A)')'variable #4: momentum_z'
  write(ilun,'(A)')'variable #5: total_energy_including_magnetic'
  do ivar=6,8
     write(ilun,'(A,I0,A,I0)')'variable #',ivar,': magnetic_left_axis_',ivar-5
  enddo
  do ivar=9,nvar
     write(ilun,'(A,I0,A)')'variable #',ivar,': conserved_nonthermal_or_passive'
  enddo
  do ivar=nvar+1,nvar+3
     write(ilun,'(A,I0,A,I0)')'variable #',ivar,': magnetic_right_axis_',ivar-nvar
  enddo
  close(ilun)
  return
#endif

  ! Write run parameters
  write(ilun,'("nvar        =",I11)')nvar
  ivar=1
  write(ilun,'("variable #",I2,": density")')ivar
  ivar=2
  write(ilun,'("variable #",I2,": velocity_x")')ivar
  if(ndim>1)then
     ivar=3
     write(ilun,'("variable #",I2,": velocity_y")')ivar
  endif
  if(ndim>2)then
     ivar=4
     write(ilun,'("variable #",I2,": velocity_z")')ivar
  endif
#if NENER>0
  ! Non-thermal pressures
  do ivar=ndim+2,ndim+1+nener
     write(ilun,'("variable #",I2,": non_thermal_pressure_",I1)')ivar,ivar-ndim-1
  end do
#endif
  ivar=ndim+2+nener
  write(ilun,'("variable #",I2,": thermal_pressure")')ivar
#if NVAR>NDIM+2+NENER
  ! Passive scalars
  do ivar=ndim+3+nener,nvar
#ifdef DUST_LIVE
     if(dust_relative_motion.and.ivar>=idust_momentum.and.ivar<idust_momentum+3*ndust_phase)then
        write(ilun,'("variable #",I0,": dust_absolute_momentum_phase_",I0,"_axis_",I0)') &
             ivar,(ivar-idust_momentum)/3+1,mod(ivar-idust_momentum,3)+1
        cycle
     endif
     if(dust_pah_enabled().and.ivar>=idust_pah.and.ivar<idust_pah+dust_pah_nstate())then
        if(dust_pah_nstate()>1000)then
           write(ilun,'("variable #",I0,": pah_mass_state_",I4.4)')ivar,ivar-idust_pah
        else
           write(ilun,'("variable #",I0,": pah_mass_state_",I3.3)')ivar,ivar-idust_pah
        endif
        cycle
     endif
#endif
     ! CHIMES/PAH extend well beyond the old two-/one-digit display widths.
     write(ilun,'("variable #",I0,": passive_scalar_",I0)')ivar,ivar-ndim-2-nener
  end do
#endif
  
  close(ilun)

end subroutine file_descriptor_hydro

subroutine backup_hydro(filename)
#ifdef DUST_DYNAMICS
  use dust_phase_state, only: dust_phase_kinetic
#endif
  use amr_commons
  use hydro_commons
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif

  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,ilun,istart,ibound,irad
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1121
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZEOUT>0) then
     if (mod(myid-1,IOGROUPSIZEOUT)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif
  
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              do ivar=1,ndim+1
                 if(ivar==1)then
                    ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ICELL_OF(ind_grid(i),ind),1)
                    end do
                 else if(ivar>=2.and.ivar<=ndim+1)then
                    ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ICELL_OF(ind_grid(i),ind),ivar)/max(uold(ICELL_OF(ind_grid(i),ind),1),smallr)
                    end do
                 endif
                 write(ilun)xdp
              end do
#if NENER>0
              ! Write non-thermal pressures
              do ivar=ndim+3,ndim+2+nener
                 do i=1,ncache
                    xdp(i)=(gamma_rad(ivar-ndim-2)-1d0)*uold(ICELL_OF(ind_grid(i),ind),ivar)
                 end do
                 write(ilun)xdp
              end do
#endif
              ! Write thermal pressure
              do i=1,ncache
                 xdp(i)=uold(ICELL_OF(ind_grid(i),ind),ndim+2)
                 xdp(i)=xdp(i)-0.5d0*uold(ICELL_OF(ind_grid(i),ind),2)**2/max(uold(ICELL_OF(ind_grid(i),ind),1),smallr)
#if NDIM>1
                 xdp(i)=xdp(i)-0.5d0*uold(ICELL_OF(ind_grid(i),ind),3)**2/max(uold(ICELL_OF(ind_grid(i),ind),1),smallr)
#endif
#if NDIM>2
                 xdp(i)=xdp(i)-0.5d0*uold(ICELL_OF(ind_grid(i),ind),4)**2/max(uold(ICELL_OF(ind_grid(i),ind),1),smallr)
#endif
#if NENER>0
                 do irad=1,nener
                    xdp(i)=xdp(i)-uold(ICELL_OF(ind_grid(i),ind),ndim+2+irad)
                 end do
#endif
#ifdef DUST_DYNAMICS
                 if(dust_relative_motion)then
                    xdp(i)=uold(ICELL_OF(ind_grid(i),ind),ndim+2)-dust_phase_kinetic(uold(ICELL_OF(ind_grid(i),ind),:))
                    if(nener>0)xdp(i)=xdp(i)-sum(uold(ICELL_OF(ind_grid(i),ind),inener:inener+nener-1))
                 endif
#endif
                 xdp(i)=(gamma-1d0)*xdp(i)
              end do
              write(ilun)xdp
#if NVAR>NDIM+2+NENER
              ! Write passive scalars
              do ivar=ndim+3+nener,nvar
                 do i=1,ncache
                    xdp(i)=uold(ICELL_OF(ind_grid(i),ind),ivar)/max(uold(ICELL_OF(ind_grid(i),ind),1),smallr)
                    if(dust_relative_motion.and.ivar>=idust_momentum.and.ivar<idust_momentum+3*ndust_phase) &
                         xdp(i)=uold(ICELL_OF(ind_grid(i),ind),ivar)
                 end do
                 write(ilun)xdp
              end do
#endif
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
   
  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZEOUT>0) then
     if(mod(myid,IOGROUPSIZEOUT)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif
  
  
end subroutine backup_hydro
