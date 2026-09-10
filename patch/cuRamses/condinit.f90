!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use amr_commons, only: boxlen
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar_all)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp)::phase,by,bz,rho,pres
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)
  u=0d0
#ifdef SOLVERmhd
  if(mhd_initial_condition/='uniform')then
     do i=1,nn
        select case(trim(mhd_initial_condition))
        case('alfven_x')
           ! Right-going circularly polarized wave: rho=Bx=1, v_perp=-B_perp.
           phase=2d0*acos(-1d0)*x(i,1)/boxlen
           by=.1d0*cos(phase);bz=.1d0*sin(phase)
           u(i,1)=1d0;u(i,3)=-by;u(i,4)=-bz
           u(i,6:8)=[1d0,by,bz]
           u(i,5)=.1d0/(gamma-1d0)+.5d0+by**2+bz**2
        case('brio_wu_x')
           rho=1d0;pres=1d0;by=1d0
           if(x(i,1)>=.5d0*boxlen)then
              rho=.125d0;pres=.1d0;by=-1d0
           endif
           u(i,1)=rho;u(i,6:8)=[.75d0,by,0d0]
           u(i,5)=pres/(gamma-1d0)+.5d0*(.75d0**2+by**2)
        end select
        ! Tangential fields depend only on x, normal field is constant:
        ! opposing faces are identical and discrete div(B) is exactly zero.
        u(i,nvar+1:nvar+3)=u(i,6:8)
     enddo
     return
  endif
#endif

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,nhydro+ivar)=q(1:nn,nhydro+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,nhydro+ivar)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=nhydro+1+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
