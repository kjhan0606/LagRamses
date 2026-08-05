! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic gas dynamics using either
!              MUSCL-HANCOCK scheme or Collela's PLMDE scheme
!              with various slope limiters.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
! ----------------------------------------------------------------
subroutine unsplit(uin,gravin,flux,tmp,dx,dy,dz,dt,ngrid,uouter)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 
  real(dp),dimension(1:nvector,1:2,0:3,0:3,1:nvar,1:ndim),intent(in)::uouter
  real(dp),dimension(1:nvector,1:2,0:3,0:3,1:nvar,1:ndim)::qouter

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp 
#ifndef _OPENMP
  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2       ),save::cin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp
  
  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   ),save::tx
#else
  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::qin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2       )::cin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp
  
  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   )::tx
#endif

  ! Velocity divergence
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::divu_kjhan

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi


  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  ! Translate to primative variables, compute sound speeds  
  call ctoprim(uin,qin,cin,gravin,dt,ngrid)
  if(scheme=='weno5'.or.scheme=='weno5ppm'.or.scheme=='ppm')then
     call ctoprim_outer(uouter,qouter,ngrid)
  else
     qouter=zero
  end if

  ! Compute TVD slopes
  call uslope(qin,dq,dx,dt,ngrid)



  ! Compute 3D traced-states in all three directions
  if(scheme=='muscl'.or.scheme=='weno3'.or.scheme=='weno5'.or. &
       & scheme=='weno5ppm'.or.scheme=='ppm')then
#if NDIM==1
     call trace1d(qin,dq,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d(qin,dq,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d(qin,dq,qm,qp,dx,dy,dz,dt,ngrid,qouter)
#endif
  endif
  if(scheme=='plmde')then
#if NDIM==1
     call tracex  (qin,dq,cin,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call tracexy (qin,dq,cin,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call tracexyz(qin,dq,cin,qm,qp,dx,dy,dz,dt,ngrid)
#endif
  endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,fx,tx,ngrid)
  ! Save flux in output array
  do i=if1,if2
  do j=jlo,jhi
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dt/dx
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dt/dx
        end do
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,fx,tx,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jf1,jf2
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dt/dy
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dt/dy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,fx,tx,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf1,kf2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dt/dz
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dt/dz
        end do
     end do
  end do
  end do
  end do
#endif

  if(difmag>0.0)then
    call cmpdivu(qin,divu_kjhan,dx,dy,dz,ngrid)
    call consup(uin,flux,divu_kjhan,dt,ngrid)
  endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d(q,dq,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! Local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, ip, irad
  real(dp)::dtdx
  real(dp)::r, u, p, a
  real(dp)::drx, dux, dpx, dax
  real(dp)::sr0, su0, sp0, sa0
#if NENER>0
  real(dp),dimension(1:nener)::e, dex, se0
#endif
  
  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=3

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              p   =  q(l,i,j,k,ip)
#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
              end do
#endif
              ! TVD slopes in X direction
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dpx = dq(l,i,j,k,ip,1)
#if NENER>0
              do irad=1,nener
                 dex(irad) = dq(l,i,j,k,ip+irad,1)
              end do
#endif
              
              ! Source terms (including transverse derivatives)
              sr0 = -u*drx - (dux)*r
              sp0 = -u*dpx - (dux)*gamma*p
              su0 = -u*dux - (dpx)/r
#if NENER>0
              do irad=1,nener
                 su0 = su0 - (dex(irad))/r
                 se0(irad) = -u*dex(irad) &
                      & - (dux)*gamma_rad(irad)*e(irad)
              end do
#endif

              ! Right state
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
!              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Left state
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
!              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  ! Passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 sa0 = -u*dax             ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine trace1d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d(q,dq,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip, irad
  real(dp)::dtdx, dtdy
  real(dp)::r, u, v, p, a
  real(dp)::drx, dux, dvx, dpx, dax
  real(dp)::dry, duy, dvy, dpy, day
  real(dp)::sr0, su0, sv0, sp0, sa0
#if NENER>0
  real(dp),dimension(1:nener)::e, dex, dey, se0
#endif
  
  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              v   =  q(l,i,j,k,iv)
              p   =  q(l,i,j,k,ip)
#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
              end do
#endif

              ! TVD slopes in all directions
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dpx = dq(l,i,j,k,ip,1)
#if NENER>0
              do irad=1,nener
                 dex(irad) = dq(l,i,j,k,ip+irad,1)
              end do
#endif
              
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dpy = dq(l,i,j,k,ip,2)
#if NENER>0
              do irad=1,nener
                 dey(irad) = dq(l,i,j,k,ip+irad,2)
              end do
#endif
              
              ! source terms (with transverse derivatives)
              sr0 = -u*drx-v*dry - (dux+dvy)*r
              sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma*p
              su0 = -u*dux-v*duy - (dpx    )/r
              sv0 = -u*dvx-v*dvy - (dpy    )/r
#if NENER>0
              do irad=1,nener
                 su0 = su0 - (dex(irad))/r
                 sv0 = sv0 - (dey(irad))/r
                 se0(irad) = -u*dex(irad)-v*dey(irad) &
                      & - (dux+dvy)*gamma_rad(irad)*e(irad)
              end do
#endif

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
!              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
!              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
              qp(l,i,j,k,iu,2) = u - half*duy + su0*dtdy*half
              qp(l,i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
              qp(l,i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
!              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              if(qp(l,i,j,k,ir,2)<smallr)qp(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,2) = e(irad) - half*dey(irad) + se0(irad)*dtdy*half
              end do
#endif

              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
              qm(l,i,j,k,iu,2) = u + half*duy + su0*dtdy*half
              qm(l,i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
              qm(l,i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
!              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              if(qm(l,i,j,k,ir,2)<smallr)qm(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,2) = e(irad) + half*dey(irad) + se0(irad)*dtdy*half
              end do
#endif

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  ! passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sa0 = -u*dax-v*day       ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half   ! Top state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half   ! Bottom state
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d(q,dq,qm,qp,dx,dy,dz,dt,ngrid,qouter)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,1:2,0:3,0:3,1:nvar,1:ndim)::qouter

  ! declare local variables
  integer ::i, j, k, l, n, ivar_rec, idim_trace
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, iw, ip, irad
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::r, u, v, w, p, a, aleft, aright
  real(dp)::vmm,vm,v0,vp,vpp
  logical::troubled_x,troubled_y,troubled_z
  logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::troubled_x_cache
  logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::troubled_y_cache
  logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::troubled_z_cache
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0
  real(dp),dimension(1:nvar,1:ndim)::qedge_left,qedge_right
#if NENER>0
  real(dp),dimension(1:nener)::e, dex, dey, dez, se0
#endif
  
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              v   =  q(l,i,j,k,iv)
              w   =  q(l,i,j,k,iw)
              p   =  q(l,i,j,k,ip)
#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
              end do
#endif

              ! TVD slopes in all 3 directions
              drx = dq(l,i,j,k,ir,1)
              dpx = dq(l,i,j,k,ip,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dwx = dq(l,i,j,k,iw,1)
#if NENER>0
              do irad=1,nener
                 dex(irad) = dq(l,i,j,k,ip+irad,1)
              end do
#endif
              
              dry = dq(l,i,j,k,ir,2)
              dpy = dq(l,i,j,k,ip,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,j,k,iw,2)
#if NENER>0
              do irad=1,nener
                 dey(irad) = dq(l,i,j,k,ip+irad,2)
              end do
#endif
              
              drz = dq(l,i,j,k,ir,3)
              dpz = dq(l,i,j,k,ip,3)
              duz = dq(l,i,j,k,iu,3)
              dvz = dq(l,i,j,k,iv,3)
              dwz = dq(l,i,j,k,iw,3)
#if NENER>0
              do irad=1,nener
                 dez(irad) = dq(l,i,j,k,ip+irad,3)
              end do
#endif

              ! Reconstruct primitive states at each cell edge.  WENO5 and
              ! PPM use five cell averages along each coordinate line; only
              ! the i/j/k=0 and 3 trace cells touch the extra axial slabs.
              ! In WENO5-PPM, discontinuities in any hydro primitive select
              ! the same monotone PPM fallback for all reconstructed fields.
              troubled_x=.false.; troubled_y=.false.; troubled_z=.false.
              if(scheme=='weno5ppm')then
                 do ivar_rec=1,ndim+2+nener
                    if(i==0)then
                       vmm=qouter(l,1,j,k,ivar_rec,1)
                    else
                       vmm=q(l,i-2,j,k,ivar_rec)
                    end if
                    vm=q(l,i-1,j,k,ivar_rec); v0=q(l,i,j,k,ivar_rec)
                    vp=q(l,i+1,j,k,ivar_rec)
                    if(i==3)then
                       vpp=qouter(l,2,j,k,ivar_rec,1)
                    else
                       vpp=q(l,i+2,j,k,ivar_rec)
                    end if
                    call detect_five_trouble(vmm,vm,v0,vp,vpp,troubled_x)

                    if(j==0)then
                       vmm=qouter(l,1,i,k,ivar_rec,2)
                    else
                       vmm=q(l,i,j-2,k,ivar_rec)
                    end if
                    vm=q(l,i,j-1,k,ivar_rec); v0=q(l,i,j,k,ivar_rec)
                    vp=q(l,i,j+1,k,ivar_rec)
                    if(j==3)then
                       vpp=qouter(l,2,i,k,ivar_rec,2)
                    else
                       vpp=q(l,i,j+2,k,ivar_rec)
                    end if
                    call detect_five_trouble(vmm,vm,v0,vp,vpp,troubled_y)

                    if(k==0)then
                       vmm=qouter(l,1,i,j,ivar_rec,3)
                    else
                       vmm=q(l,i,j,k-2,ivar_rec)
                    end if
                    vm=q(l,i,j,k-1,ivar_rec); v0=q(l,i,j,k,ivar_rec)
                    vp=q(l,i,j,k+1,ivar_rec)
                    if(k==3)then
                       vpp=qouter(l,2,i,j,ivar_rec,3)
                    else
                       vpp=q(l,i,j,k+2,ivar_rec)
                    end if
                    call detect_five_trouble(vmm,vm,v0,vp,vpp,troubled_z)
                 end do
              end if
              troubled_x_cache(l,i,j,k)=troubled_x
              troubled_y_cache(l,i,j,k)=troubled_y
              troubled_z_cache(l,i,j,k)=troubled_z
              do ivar_rec=1,nvar
                 if(scheme=='weno3')then
                    call reconstruct_weno3(q(l,i-1,j,k,ivar_rec), &
                         & q(l,i,j,k,ivar_rec),q(l,i+1,j,k,ivar_rec), &
                         & qedge_left(ivar_rec,1),qedge_right(ivar_rec,1))
                    call reconstruct_weno3(q(l,i,j-1,k,ivar_rec), &
                         & q(l,i,j,k,ivar_rec),q(l,i,j+1,k,ivar_rec), &
                         & qedge_left(ivar_rec,2),qedge_right(ivar_rec,2))
                    call reconstruct_weno3(q(l,i,j,k-1,ivar_rec), &
                         & q(l,i,j,k,ivar_rec),q(l,i,j,k+1,ivar_rec), &
                         & qedge_left(ivar_rec,3),qedge_right(ivar_rec,3))
                 else if(scheme=='weno5'.or.scheme=='weno5ppm'.or.scheme=='ppm')then
                    if(i==0)then
                       vmm=qouter(l,1,j,k,ivar_rec,1)
                    else
                       vmm=q(l,i-2,j,k,ivar_rec)
                    end if
                    vm=q(l,i-1,j,k,ivar_rec); v0=q(l,i,j,k,ivar_rec)
                    vp=q(l,i+1,j,k,ivar_rec)
                    if(i==3)then
                       vpp=qouter(l,2,j,k,ivar_rec,1)
                    else
                       vpp=q(l,i+2,j,k,ivar_rec)
                    end if
                    call reconstruct_five(vmm,vm,v0,vp,vpp, &
                         & qedge_left(ivar_rec,1),qedge_right(ivar_rec,1),troubled_x)

                    if(j==0)then
                       vmm=qouter(l,1,i,k,ivar_rec,2)
                    else
                       vmm=q(l,i,j-2,k,ivar_rec)
                    end if
                    vm=q(l,i,j-1,k,ivar_rec); v0=q(l,i,j,k,ivar_rec)
                    vp=q(l,i,j+1,k,ivar_rec)
                    if(j==3)then
                       vpp=qouter(l,2,i,k,ivar_rec,2)
                    else
                       vpp=q(l,i,j+2,k,ivar_rec)
                    end if
                    call reconstruct_five(vmm,vm,v0,vp,vpp, &
                         & qedge_left(ivar_rec,2),qedge_right(ivar_rec,2),troubled_y)

                    if(k==0)then
                       vmm=qouter(l,1,i,j,ivar_rec,3)
                    else
                       vmm=q(l,i,j,k-2,ivar_rec)
                    end if
                    vm=q(l,i,j,k-1,ivar_rec); v0=q(l,i,j,k,ivar_rec)
                    vp=q(l,i,j,k+1,ivar_rec)
                    if(k==3)then
                       vpp=qouter(l,2,i,j,ivar_rec,3)
                    else
                       vpp=q(l,i,j,k+2,ivar_rec)
                    end if
                    call reconstruct_five(vmm,vm,v0,vp,vpp, &
                         & qedge_left(ivar_rec,3),qedge_right(ivar_rec,3),troubled_z)
                 else
                    qedge_left (ivar_rec,1)=q(l,i,j,k,ivar_rec)-half*dq(l,i,j,k,ivar_rec,1)
                    qedge_right(ivar_rec,1)=q(l,i,j,k,ivar_rec)+half*dq(l,i,j,k,ivar_rec,1)
                    qedge_left (ivar_rec,2)=q(l,i,j,k,ivar_rec)-half*dq(l,i,j,k,ivar_rec,2)
                    qedge_right(ivar_rec,2)=q(l,i,j,k,ivar_rec)+half*dq(l,i,j,k,ivar_rec,2)
                    qedge_left (ivar_rec,3)=q(l,i,j,k,ivar_rec)-half*dq(l,i,j,k,ivar_rec,3)
                    qedge_right(ivar_rec,3)=q(l,i,j,k,ivar_rec)+half*dq(l,i,j,k,ivar_rec,3)
                 end if
              end do

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
              sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
              su0 = -u*dux-v*duy-w*duz - (dpx        )/r
              sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
              sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r
#if NENER>0
              do irad=1,nener
                 su0 = su0 - (dex(irad))/r
                 sv0 = sv0 - (dey(irad))/r
                 sw0 = sw0 - (dez(irad))/r
                 se0(irad) = -u*dex(irad)-v*dey(irad)-w*dez(irad) & 
                      & - (dux+dvy+dwz)*gamma_rad(irad)*e(irad)
              end do
#endif

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = qedge_left(ir,1) + sr0*dtdx*half
              qp(l,i,j,k,ip,1) = qedge_left(ip,1) + sp0*dtdx*half
              qp(l,i,j,k,iu,1) = qedge_left(iu,1) + su0*dtdx*half
              qp(l,i,j,k,iv,1) = qedge_left(iv,1) + sv0*dtdx*half
              qp(l,i,j,k,iw,1) = qedge_left(iw,1) + sw0*dtdx*half
!              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,1) = qedge_left(ip+irad,1) + se0(irad)*dtdx*half
              end do
#endif

              ! Left state at left interface
              qm(l,i,j,k,ir,1) = qedge_right(ir,1) + sr0*dtdx*half
              qm(l,i,j,k,ip,1) = qedge_right(ip,1) + sp0*dtdx*half
              qm(l,i,j,k,iu,1) = qedge_right(iu,1) + su0*dtdx*half
              qm(l,i,j,k,iv,1) = qedge_right(iv,1) + sv0*dtdx*half
              qm(l,i,j,k,iw,1) = qedge_right(iw,1) + sw0*dtdx*half
!              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,1) = qedge_right(ip+irad,1) + se0(irad)*dtdx*half
              end do
#endif

              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = qedge_left(ir,2) + sr0*dtdy*half
              qp(l,i,j,k,ip,2) = qedge_left(ip,2) + sp0*dtdy*half
              qp(l,i,j,k,iu,2) = qedge_left(iu,2) + su0*dtdy*half
              qp(l,i,j,k,iv,2) = qedge_left(iv,2) + sv0*dtdy*half
              qp(l,i,j,k,iw,2) = qedge_left(iw,2) + sw0*dtdy*half
!              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              if(qp(l,i,j,k,ir,2)<smallr)qp(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,2) = qedge_left(ip+irad,2) + se0(irad)*dtdy*half
              end do
#endif

              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = qedge_right(ir,2) + sr0*dtdy*half
              qm(l,i,j,k,ip,2) = qedge_right(ip,2) + sp0*dtdy*half
              qm(l,i,j,k,iu,2) = qedge_right(iu,2) + su0*dtdy*half
              qm(l,i,j,k,iv,2) = qedge_right(iv,2) + sv0*dtdy*half
              qm(l,i,j,k,iw,2) = qedge_right(iw,2) + sw0*dtdy*half
!              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              if(qm(l,i,j,k,ir,2)<smallr)qm(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,2) = qedge_right(ip+irad,2) + se0(irad)*dtdy*half
              end do
#endif

              ! Back state at front interface
              qp(l,i,j,k,ir,3) = qedge_left(ir,3) + sr0*dtdz*half
              qp(l,i,j,k,ip,3) = qedge_left(ip,3) + sp0*dtdz*half
              qp(l,i,j,k,iu,3) = qedge_left(iu,3) + su0*dtdz*half
              qp(l,i,j,k,iv,3) = qedge_left(iv,3) + sv0*dtdz*half
              qp(l,i,j,k,iw,3) = qedge_left(iw,3) + sw0*dtdz*half
!              qp(l,i,j,k,ir,3) = max(smallr, qp(l,i,j,k,ir,3))
              if(qp(l,i,j,k,ir,3)<smallr)qp(l,i,j,k,ir,3)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,3) = qedge_left(ip+irad,3) + se0(irad)*dtdz*half
              end do
#endif

              ! Front state at back interface
              qm(l,i,j,k,ir,3) = qedge_right(ir,3) + sr0*dtdz*half
              qm(l,i,j,k,ip,3) = qedge_right(ip,3) + sp0*dtdz*half
              qm(l,i,j,k,iu,3) = qedge_right(iu,3) + su0*dtdz*half
              qm(l,i,j,k,iv,3) = qedge_right(iv,3) + sv0*dtdz*half
              qm(l,i,j,k,iw,3) = qedge_right(iw,3) + sw0*dtdz*half
!              qm(l,i,j,k,ir,3) = max(smallr, qm(l,i,j,k,ir,3))
              if(qm(l,i,j,k,ir,3)<smallr)qm(l,i,j,k,ir,3)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,3) = qedge_right(ip+irad,3) + se0(irad)*dtdz*half
              end do
#endif

              ! The stencil clamp below is applied before the Hancock
              ! predictor.  Guard the two thermodynamic variables again
              ! afterwards; this is a robustness check, not a proof of a
              ! discrete maximum principle for the full update.
              if(scheme=='weno5ppm')then
                 do idim_trace=1,ndim
                    if(qp(l,i,j,k,ir,idim_trace)<smallr) &
                         & qp(l,i,j,k,ir,idim_trace)=r
                    if(qm(l,i,j,k,ir,idim_trace)<smallr) &
                         & qm(l,i,j,k,ir,idim_trace)=r
                    if(qp(l,i,j,k,ip,idim_trace)< &
                         & max(qp(l,i,j,k,ir,idim_trace),smallr)*smallc**2/gamma) &
                         & qp(l,i,j,k,ip,idim_trace)=p
                    if(qm(l,i,j,k,ip,idim_trace)< &
                         & max(qm(l,i,j,k,ir,idim_trace),smallr)*smallc**2/gamma) &
                         & qm(l,i,j,k,ip,idim_trace)=p
                 end do
              end if

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  ! Passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 troubled_x=troubled_x_cache(l,i,j,k)
                 troubled_y=troubled_y_cache(l,i,j,k)
                 troubled_z=troubled_z_cache(l,i,j,k)
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 w   = q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sa0 = -u*dax-v*day-w*daz     ! Source terms
                 if(scheme=='weno3')then
                    call reconstruct_weno3(q(l,i-1,j,k,n),a,q(l,i+1,j,k,n),aleft,aright)
                 else if(scheme=='weno5'.or.scheme=='weno5ppm'.or.scheme=='ppm')then
                    if(i==0)then
                       vmm=qouter(l,1,j,k,n,1)
                    else
                       vmm=q(l,i-2,j,k,n)
                    end if
                    vm=q(l,i-1,j,k,n); vp=q(l,i+1,j,k,n)
                    if(i==3)then
                       vpp=qouter(l,2,j,k,n,1)
                    else
                       vpp=q(l,i+2,j,k,n)
                    end if
                    call reconstruct_five(vmm,vm,a,vp,vpp,aleft,aright,troubled_x)
                 else
                    aleft=a-half*dax
                    aright=a+half*dax
                 end if
                 qp(l,i,j,k,n,1) = aleft  + sa0*dtdx*half  ! Right state
                 qm(l,i,j,k,n,1) = aright + sa0*dtdx*half  ! Left state
                 if(scheme=='weno3')then
                    call reconstruct_weno3(q(l,i,j-1,k,n),a,q(l,i,j+1,k,n),aleft,aright)
                 else if(scheme=='weno5'.or.scheme=='weno5ppm'.or.scheme=='ppm')then
                    if(j==0)then
                       vmm=qouter(l,1,i,k,n,2)
                    else
                       vmm=q(l,i,j-2,k,n)
                    end if
                    vm=q(l,i,j-1,k,n); vp=q(l,i,j+1,k,n)
                    if(j==3)then
                       vpp=qouter(l,2,i,k,n,2)
                    else
                       vpp=q(l,i,j+2,k,n)
                    end if
                    call reconstruct_five(vmm,vm,a,vp,vpp,aleft,aright,troubled_y)
                 else
                    aleft=a-half*day
                    aright=a+half*day
                 end if
                 qp(l,i,j,k,n,2) = aleft  + sa0*dtdy*half  ! Bottom state
                 qm(l,i,j,k,n,2) = aright + sa0*dtdy*half  ! Upper state
                 if(scheme=='weno3')then
                    call reconstruct_weno3(q(l,i,j,k-1,n),a,q(l,i,j,k+1,n),aleft,aright)
                 else if(scheme=='weno5'.or.scheme=='weno5ppm'.or.scheme=='ppm')then
                    if(k==0)then
                       vmm=qouter(l,1,i,j,n,3)
                    else
                       vmm=q(l,i,j,k-2,n)
                    end if
                    vm=q(l,i,j,k-1,n); vp=q(l,i,j,k+1,n)
                    if(k==3)then
                       vpp=qouter(l,2,i,j,n,3)
                    else
                       vpp=q(l,i,j,k+2,n)
                    end if
                    call reconstruct_five(vmm,vm,a,vp,vpp,aleft,aright,troubled_z)
                 else
                    aleft=a-half*daz
                    aright=a+half*daz
                 end if
                 qp(l,i,j,k,n,3) = aleft  + sa0*dtdz*half  ! Front state
                 qm(l,i,j,k,n,3) = aright + sa0*dtdz*half  ! Back state
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reconstruct_weno3(qminus,qzero,qplus,qleft,qright)
  use amr_parameters, only:dp
  implicit none

  real(dp),intent(in)::qminus,qzero,qplus
  real(dp),intent(out)::qleft,qright
  real(dp)::beta_minus,beta_plus,eps
  real(dp)::alpha0,alpha1,weight0,weight1
  real(dp)::candidate0,candidate1,qmin,qmax

  beta_minus=(qzero-qminus)**2
  beta_plus =(qplus-qzero)**2
  eps=1.0d-12*max(1.0d0,qminus*qminus,qzero*qzero,qplus*qplus)

  ! State at the left edge of the central cell.
  alpha0=(2.0d0/3.0d0)/(eps+beta_minus)**2
  alpha1=(1.0d0/3.0d0)/(eps+beta_plus )**2
  weight0=alpha0/(alpha0+alpha1)
  weight1=alpha1/(alpha0+alpha1)
  candidate0=0.5d0*(qminus+qzero)
  candidate1=1.5d0*qzero-0.5d0*qplus
  qleft=weight0*candidate0+weight1*candidate1

  ! State at the right edge of the central cell.
  alpha0=(1.0d0/3.0d0)/(eps+beta_minus)**2
  alpha1=(2.0d0/3.0d0)/(eps+beta_plus )**2
  weight0=alpha0/(alpha0+alpha1)
  weight1=alpha1/(alpha0+alpha1)
  candidate0=1.5d0*qzero-0.5d0*qminus
  candidate1=0.5d0*(qzero+qplus)
  qright=weight0*candidate0+weight1*candidate1

  ! Preserve local bounds near strong contacts and shocks.
  qmin=min(qminus,qzero,qplus)
  qmax=max(qminus,qzero,qplus)
  qleft =min(qmax,max(qmin,qleft ))
  qright=min(qmax,max(qmin,qright))

end subroutine reconstruct_weno3
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reconstruct_five(qmm,qm,qzero,qp,qpp,qleft,qright,use_ppm)
  use amr_parameters, only:dp
  use hydro_parameters, only:scheme
  implicit none

  real(dp),intent(in)::qmm,qm,qzero,qp,qpp
  real(dp),intent(out)::qleft,qright
  logical,intent(in)::use_ppm
  real(dp)::qmin,qmax

  if(scheme=='ppm'.or.use_ppm)then
     call reconstruct_ppm(qmm,qm,qzero,qp,qpp,qleft,qright)
  else
     call reconstruct_weno5(qmm,qm,qzero,qp,qpp,qleft,qright)
  end if

  ! The hybrid is explicitly local-bound limited.  This does not make the
  ! complete Hancock update maximum-principle preserving, hence the scheme
  ! name describes its WENO5--PPM composition rather than claiming "BP".
  if(scheme=='weno5ppm')then
     qmin=min(qmm,qm,qzero,qp,qpp)
     qmax=max(qmm,qm,qzero,qp,qpp)
     qleft =min(qmax,max(qmin,qleft ))
     qright=min(qmax,max(qmin,qright))
  end if

end subroutine reconstruct_five
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine detect_five_trouble(qmm,qm,qzero,qp,qpp,troubled)
  use amr_parameters, only:dp
  implicit none

  real(dp),intent(in)::qmm,qm,qzero,qp,qpp
  logical,intent(inout)::troubled
  real(dp)::qmin,qmax,qrange,qscale,third_difference

  qmin=min(qmm,qm,qzero,qp,qpp)
  qmax=max(qmm,qm,qzero,qp,qpp)
  qrange=qmax-qmin
  qscale=max(1.0d-30,abs(qmm),abs(qm),abs(qzero),abs(qp),abs(qpp))
  third_difference=max(abs(qp-3.0d0*qzero+3.0d0*qm-qmm), &
       & abs(qpp-3.0d0*qp+3.0d0*qzero-qm))

  ! A third difference distinguishes an unresolved jump or kink from a
  ! resolved quadratic extremum (whose third difference vanishes).  Once
  ! any hydro primitive is troubled, retain that flag for the direction.
  troubled=troubled .or. (qrange>1.0d-10*qscale .and. &
       & third_difference>0.2d0*qrange)

end subroutine detect_five_trouble
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reconstruct_weno5(qmm,qm,qzero,qp,qpp,qleft,qright)
  use amr_parameters, only:dp
  implicit none

  real(dp),intent(in)::qmm,qm,qzero,qp,qpp
  real(dp),intent(out)::qleft,qright
  real(dp)::qmin,qmax,qscale

  qmin=min(qmm,qm,qzero,qp,qpp)
  qmax=max(qmm,qm,qzero,qp,qpp)
  qscale=max(1.0d0,abs(qmm),abs(qm),abs(qzero),abs(qp),abs(qpp))
  if(qmax-qmin<=1.0d-14*qscale)then
     qleft=qzero
     qright=qzero
     return
  end if

  ! WENO-Z nonlinear weights.  Reversing the stencil gives the state at
  ! the left edge of the same cell without a second set of coefficients.
  call weno5z_face(qmm,qm,qzero,qp,qpp,qright)
  call weno5z_face(qpp,qp,qzero,qm,qmm,qleft)

end subroutine reconstruct_weno5
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine weno5z_face(qmm,qm,qzero,qp,qpp,qface)
  use amr_parameters, only:dp
  implicit none

  real(dp),intent(in)::qmm,qm,qzero,qp,qpp
  real(dp),intent(out)::qface
  real(dp)::b0,b1,b2,tau5,eps
  real(dp)::a0,a1,a2,asum,w0,w1,w2
  real(dp)::p0,p1,p2,qscale

  b0=(13.0d0/12.0d0)*(qmm-2.0d0*qm+qzero)**2 &
       & +0.25d0*(qmm-4.0d0*qm+3.0d0*qzero)**2
  b1=(13.0d0/12.0d0)*(qm-2.0d0*qzero+qp)**2 &
       & +0.25d0*(qm-qp)**2
  b2=(13.0d0/12.0d0)*(qzero-2.0d0*qp+qpp)**2 &
       & +0.25d0*(3.0d0*qzero-4.0d0*qp+qpp)**2
  tau5=abs(b0-b2)
  qscale=max(1.0d0,qmm*qmm,qm*qm,qzero*qzero,qp*qp,qpp*qpp)
  eps=1.0d-14*qscale

  a0=0.1d0*(1.0d0+(tau5/(b0+eps))**2)
  a1=0.6d0*(1.0d0+(tau5/(b1+eps))**2)
  a2=0.3d0*(1.0d0+(tau5/(b2+eps))**2)
  asum=a0+a1+a2
  w0=a0/asum; w1=a1/asum; w2=a2/asum

  p0=(1.0d0/3.0d0)*qmm-(7.0d0/6.0d0)*qm+(11.0d0/6.0d0)*qzero
  p1=-(1.0d0/6.0d0)*qm+(5.0d0/6.0d0)*qzero+(1.0d0/3.0d0)*qp
  p2=(1.0d0/3.0d0)*qzero+(5.0d0/6.0d0)*qp-(1.0d0/6.0d0)*qpp
  qface=w0*p0+w1*p1+w2*p2

end subroutine weno5z_face
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reconstruct_ppm(qmm,qm,qzero,qp,qpp,qleft,qright)
  use amr_parameters, only:dp
  implicit none

  real(dp),intent(in)::qmm,qm,qzero,qp,qpp
  real(dp),intent(out)::qleft,qright
  real(dp)::dqcell,q6

  ! Fourth-order interface interpolation followed by the original PPM
  ! monotonicity constraints.  Contact steepening and shock flattening are
  ! deliberately left out of this first controlled comparison.
  qleft =(7.0d0/12.0d0)*(qm+qzero)-(1.0d0/12.0d0)*(qmm+qp)
  qright=(7.0d0/12.0d0)*(qzero+qp)-(1.0d0/12.0d0)*(qm+qpp)

  qleft =min(max(qm,qzero),max(min(qm,qzero),qleft))
  qright=min(max(qzero,qp),max(min(qzero,qp),qright))

  if((qright-qzero)*(qzero-qleft)<=0.0d0)then
     qleft=qzero
     qright=qzero
  else
     dqcell=qright-qleft
     q6=6.0d0*qzero-3.0d0*(qleft+qright)
     if(dqcell*q6>dqcell*dqcell)then
        qleft=3.0d0*qzero-2.0d0*qright
     else if(dqcell*q6<-(dqcell*dqcell))then
        qright=3.0d0*qzero-2.0d0*qleft
     end if
  end if

end subroutine reconstruct_ppm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &            flx,tmp,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  
  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp)::entho
#ifndef _OPENMP
  real(dp),dimension(1:nvector,1:nvar),save::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1),save::fgdnv
#else
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1)::fgdnv
#endif

  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Mass density
           do l = 1, ngrid
              qleft (l,1) = qm(l,i,j,k,1,xdim)
              qright(l,1) = qp(l,i,j,k,1,xdim)
           end do
           
           ! Normal velocity
           do l = 1, ngrid
              qleft (l,2) = qm(l,i,j,k,ln,xdim)
              qright(l,2) = qp(l,i,j,k,ln,xdim)
           end do
           
           ! Pressure
           do l = 1, ngrid
              qleft (l,3) = qm(l,i,j,k,ndim+2,xdim)
              qright(l,3) = qp(l,i,j,k,ndim+2,xdim)
           end do
           
           ! Tangential velocity 1
#if NDIM>1
           do l = 1, ngrid
              qleft (l,4) = qm(l,i,j,k,lt1,xdim)
              qright(l,4) = qp(l,i,j,k,lt1,xdim)
           end do
#endif
           ! Tangential velocity 2
#if NDIM>2
           do l = 1, ngrid
              qleft (l,5) = qm(l,i,j,k,lt2,xdim)
              qright(l,5) = qp(l,i,j,k,lt2,xdim)
           end do
#endif           
#if NVAR > NDIM + 2
           ! Other advected quantities
           do n = ndim+3, nvar
              do l = 1, ngrid
                 qleft (l,n) = qm(l,i,j,k,n,xdim)
                 qright(l,n) = qp(l,i,j,k,n,xdim)
              end do
           end do
#endif          
           ! Solve Riemann problem
           if(riemann.eq.'acoustic')then
              call riemann_acoustic(qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'exact')then
              call riemann_approx  (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'llf')then
              call riemann_llf     (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'hllc')then
              call riemann_hllc    (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'hll')then
              call riemann_hll     (qleft,qright,fgdnv,ngrid)
           else
              write(*,*)'unknown Riemann solver'
              stop
           end if
           
           ! Compute fluxes
           
           ! Mass density
           do l = 1, ngrid 
              flx(l,i,j,k,1) = fgdnv(l,1)
           end do
           
           ! Normal momentum
           do l = 1, ngrid
              flx(l,i,j,k,ln) = fgdnv(l,2)
           end do

           ! Transverse momentum 1
#if NDIM>1
           do l = 1, ngrid
              flx(l,i,j,k,lt1) = fgdnv(l,4)
           end do
#endif
           ! Transverse momentum 2
#if NDIM>2
           do l = 1, ngrid
              flx(l,i,j,k,lt2) = fgdnv(l,5)
           end do
#endif           
           ! Total energy
           do l = 1, ngrid
              flx(l,i,j,k,ndim+2) = fgdnv(l,3)
           end do

#if NVAR > NDIM + 2
           ! Other advected quantities
           do n = ndim+3, nvar
              do l = 1, ngrid
                 flx(l,i,j,k,n) = fgdnv(l,n)
              end do
           end do
#endif
           ! Normal velocity
           do l = 1, ngrid
              tmp(l,i,j,k,1) = half*(qleft(l,2)+qright(l,2))
           end do
           ! Internal energy flux
           do l = 1,ngrid
              tmp(l,i,j,k,2) = fgdnv(l,nvar+1)
           end do

        end do
     end do
  end do
  
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,c,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  integer ::i, j, k, l, n, idim, irad
  real(dp)::eint, smalle, smalle_poly, dtxhalf, oneoverrho
  real(dp)::eken, erad

  smalle = smallc**2/gamma/(gamma-one)
  dtxhalf = dt*half

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           do l = 1, ngrid

              ! Compute density
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)

              ! Compute velocities
              oneoverrho = one/q(l,i,j,k,1)
              q(l,i,j,k,2) = uin(l,i,j,k,2)*oneoverrho
#if NDIM>1
              q(l,i,j,k,3) = uin(l,i,j,k,3)*oneoverrho
#endif
#if NDIM>2
              q(l,i,j,k,4) = uin(l,i,j,k,4)*oneoverrho
#endif

              ! Compute specific kinetic energy
              eken = half*q(l,i,j,k,2)*q(l,i,j,k,2)
#if NDIM>1
              eken = eken + half*q(l,i,j,k,3)*q(l,i,j,k,3)
#endif
#if NDIM>2
              eken = eken + half*q(l,i,j,k,4)*q(l,i,j,k,4)
#endif
              ! Compute non-thermal pressure
              erad = zero
#if NENER>0
              do irad = 1,nener
                 q(l,i,j,k,ndim+2+irad) = (gamma_rad(irad)-one)*uin(l,i,j,k,ndim+2+irad)
                 erad = erad+uin(l,i,j,k,ndim+2+irad)*oneoverrho
              enddo
#endif
              ! Compute thermal pressure with polytropic floor (eEOS)
              smalle_poly = smalle
              if(eeos_poly_coeff > 0d0) then
                 smalle_poly = max(smalle, eeos_poly_coeff * q(l,i,j,k,1)**(eeos_poly_alpha-1d0))
              end if
              eint = MAX(uin(l,i,j,k,ndim+2)*oneoverrho-eken-erad, smalle_poly)
              ! Write back floored energy to conserved variable for consistency
              uin(l,i,j,k,ndim+2) = q(l,i,j,k,1)*(eint+eken+erad)
              q(l,i,j,k,ndim+2) = (gamma-one)*q(l,i,j,k,1)*eint

              ! Compute sound speed (c_eff includes SGS turbulent pressure)
              c(l,i,j,k)=gamma*q(l,i,j,k,ndim+2)
#if NENER>0
              do irad=1,nener
                 c(l,i,j,k)=c(l,i,j,k)+gamma_rad(irad)*q(l,i,j,k,ndim+2+irad)
              enddo
#endif
              if(use_sgs .and. sgs_hydro .and. isgs>0) then
                 ! BUG FIX: q(isgs) not yet set (passive scalars converted below)
                 ! Use uin(isgs) directly: uin(isgs) = rho*e_sgs, c is rho*c^2
                 c(l,i,j,k)=c(l,i,j,k)+(2d0/3d0)*max(uin(l,i,j,k,isgs),0d0)
              end if
              c(l,i,j,k)=sqrt(c(l,i,j,k)*oneoverrho)

              ! Gravity predictor step
              q(l,i,j,k,2) = q(l,i,j,k,2) + gravin(l,i,j,k,1)*dtxhalf
#if NDIM>1
              q(l,i,j,k,3) = q(l,i,j,k,3) + gravin(l,i,j,k,2)*dtxhalf
#endif
#if NDIM>2
              q(l,i,j,k,4) = q(l,i,j,k,4) + gravin(l,i,j,k,3)*dtxhalf
#endif

           end do
        end do
     end do
  end do



#if NVAR > NDIM + 2 + NENER
  ! Passive scalar
  do n = ndim+nener+3, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 oneoverrho = one/q(l,i,j,k,1)
                 q(l,i,j,k,n) = uin(l,i,j,k,n)*oneoverrho
              end do
           end do
        end do
     end do
  end do
#endif
 
end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim_outer(uin,q,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer,intent(in)::ngrid
  real(dp),dimension(1:nvector,1:2,0:3,0:3,1:nvar,1:ndim),intent(in)::uin
  real(dp),dimension(1:nvector,1:2,0:3,0:3,1:nvar,1:ndim),intent(out)::q

  integer::l,idim,iside,it1,it2,n,irad
  real(dp)::rho,oneoverrho,eken,erad,eint,smalle,smalle_poly

  smalle=smallc**2/gamma/(gamma-one)
  do idim=1,ndim
     do iside=1,2
        do it2=0,3
           do it1=0,3
              do l=1,ngrid
                 rho=max(uin(l,iside,it1,it2,1,idim),smallr)
                 q(l,iside,it1,it2,1,idim)=rho
                 oneoverrho=one/rho
                 q(l,iside,it1,it2,2,idim)= &
                      & uin(l,iside,it1,it2,2,idim)*oneoverrho
#if NDIM>1
                 q(l,iside,it1,it2,3,idim)= &
                      & uin(l,iside,it1,it2,3,idim)*oneoverrho
#endif
#if NDIM>2
                 q(l,iside,it1,it2,4,idim)= &
                      & uin(l,iside,it1,it2,4,idim)*oneoverrho
#endif
                 eken=half*q(l,iside,it1,it2,2,idim)**2
#if NDIM>1
                 eken=eken+half*q(l,iside,it1,it2,3,idim)**2
#endif
#if NDIM>2
                 eken=eken+half*q(l,iside,it1,it2,4,idim)**2
#endif
                 erad=zero
#if NENER>0
                 do irad=1,nener
                    q(l,iside,it1,it2,ndim+2+irad,idim)= &
                         & (gamma_rad(irad)-one)* &
                         & uin(l,iside,it1,it2,ndim+2+irad,idim)
                    erad=erad+uin(l,iside,it1,it2,ndim+2+irad,idim)*oneoverrho
                 end do
#endif
                 smalle_poly=smalle
                 if(eeos_poly_coeff>0.0d0)then
                    smalle_poly=max(smalle,eeos_poly_coeff*rho**(eeos_poly_alpha-1.0d0))
                 end if
                 eint=max(uin(l,iside,it1,it2,ndim+2,idim)*oneoverrho &
                      & -eken-erad,smalle_poly)
                 q(l,iside,it1,it2,ndim+2,idim)=(gamma-one)*rho*eint
              end do
           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  do idim=1,ndim
     do n=ndim+nener+3,nvar
        do iside=1,2
           do it2=0,3
              do it1=0,3
                 do l=1,ngrid
                    q(l,iside,it1,it2,n,idim)=uin(l,iside,it1,it2,n,idim) &
                         & /q(l,iside,it1,it2,1,idim)
                 end do
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine ctoprim_outer
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(q,dq,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi





  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

#if NDIM==1
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
                 do l = 1, ngrid
                    dlft = MIN(slope_type,2)*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = MIN(slope_type,2)*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/MIN(slope_type,2)
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else if(slope_type==4)then ! superbee
                 do l = 1, ngrid
                    dcen = q(l,i,j,k,2)*dt/dx
                    dlft = two/(one+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                    drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                    dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dsgn = sign(one, dlft)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                 end do
              else if(slope_type==5)then ! ultrabee
                 if(n==1)then
                    do l = 1, ngrid
                       dcen = q(l,i,j,k,2)*dt/dx
                       if(dcen>=0)then
                          dlft = two/(zero+dcen+1d-10)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                          drgt = two/(one -dcen      )*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       else
                          dlft = two/(one +dcen      )*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                          drgt = two/(zero-dcen+1d-10)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       endif
                       dsgn = sign(one, dlft)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else if(slope_type==6)then ! unstable
                 if(n==1)then
                    do l = 1, ngrid
                       dlft = (q(l,i,j,k,n)-q(l,i-1,j,k,n))
                       drgt = (q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       slop = 0.5*(dlft+drgt)
                       dlim = slop
                       dq(l,i,j,k,n,1) = dlim
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else if(slope_type==7)then ! van Leer
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else
                 write(*,*)'Unknown slope type'
                 stop
              end if
           end do
        end do
     end do     
  end do
#endif

#if NDIM==2              
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                    dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                    dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                    dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                    dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                    dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                    dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                    dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                    dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)
                    
                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    
                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))
                    
                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif
                    
                    dlim = slop
                    
                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy

                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==7)then ! van Leer
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0*dlft*drgt/(dlft+drgt))
                       end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2)=zero
                    else
                       dq(l,i,j,k,n,2)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
#endif

#if NDIM==3
  if(slope_type==1)then  ! minmod
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i  ,j,k,n) - q(l,i-1,j,k,n)
                    drgt = q(l,i+1,j,k,n) - q(l,i  ,j,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,1) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,1) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j  ,k,n) - q(l,i,j-1,k,n)
                    drgt = q(l,i,j+1,k,n) - q(l,i,j  ,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,2) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,2) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j,k  ,n) - q(l,i,j,k-1,n)
                    drgt = q(l,i,j,k+1,n) - q(l,i,j,k  ,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,3) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,3) = max(dlft,drgt)
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==2)then ! moncen
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                    dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                    dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                    dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                    dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                    dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                    dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                    dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                    dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                    dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                    dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                    dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                    dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                    dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                    dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                    dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                    dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                    dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                    dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                    dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                    dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                    dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                    dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                    dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                    dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                    dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                    dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)
                    
                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    
                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                    dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
                    
                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif
                    
                    dlim = slop
                    
                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy
                    dq(l,i,j,k,n,3) = dlim*dfz

                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==7)then ! van Leer
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2)=zero
                    else
                       dq(l,i,j,k,n,2)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = (q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3)=zero
                    else
                       dq(l,i,j,k,n,3)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = (q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif     
#endif
  
end subroutine uslope
