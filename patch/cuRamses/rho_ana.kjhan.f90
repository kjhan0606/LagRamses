!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::dmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,dd
  real(dp)::fourpi,star_norm,gas_norm,halo_norm,concentration
  real(dp)::radius_cyl,star_height,gas_height,halo_x

  if(gravity_type==-2.or.gravity_type==-3)then
     ! Fixed, cell-centred mass model for the controlled disc force pilot.
     ! gravity_params(1:6) = Mstar,Rstar,hstar,Mgas,Rgas,hgas.
     ! gravity_params(7:9) = Mhalo,rs,rt for gravity_type=-3.
     if(gravity_params(1)<0d0.or.gravity_params(4)<0d0.or. &
          gravity_params(2)<=0d0.or.gravity_params(3)<=0d0.or. &
          gravity_params(5)<=0d0.or.gravity_params(6)<=0d0)then
        write(*,'(A)') 'ERROR: invalid analytic disc mass or scale'
        call clean_stop
     end if
     if(gravity_type==-3)then
        if(gravity_params(7)<0d0.or.gravity_params(8)<=0d0.or. &
             gravity_params(9)<=gravity_params(8))then
           write(*,'(A)') 'ERROR: invalid truncated NFW halo mass or scale'
           call clean_stop
        end if
     end if
     fourpi=4d0*acos(-1d0)
     star_norm=gravity_params(1)/(fourpi*gravity_params(2)**2*gravity_params(3))
     gas_norm=gravity_params(4)/(fourpi*gravity_params(5)**2*gravity_params(6))
     halo_norm=0d0
     if(gravity_type==-3)then
        concentration=gravity_params(9)/gravity_params(8)
        halo_norm=gravity_params(7)/(fourpi*gravity_params(8)**3* &
             (log(1d0+concentration)-concentration/(1d0+concentration)))
     end if
!$omp parallel do private(i,rx,ry,rz,rr,dd,radius_cyl,star_height,gas_height,halo_x)
     do i=1,ncell
        rx=x(i,1)-0.5d0*boxlen
        ry=x(i,2)-0.5d0*boxlen
        rz=x(i,3)-0.5d0*boxlen
        radius_cyl=sqrt(rx**2+ry**2)
        star_height=1d0/cosh(rz/gravity_params(3))**2
        gas_height=1d0/cosh(rz/gravity_params(6))**2
        dd=star_norm*exp(-radius_cyl/gravity_params(2))*star_height &
             +gas_norm*exp(-radius_cyl/gravity_params(5))*gas_height
        if(gravity_type==-3)then
           rr=sqrt(radius_cyl**2+rz**2)
           if(rr<=gravity_params(9))then
              halo_x=max(rr,tiny(rr))/gravity_params(8)
              dd=dd+halo_norm/(halo_x*(1d0+halo_x)**2)
           end if
        end if
        d(i)=dd
     end do
     return
  end if

  emass=gravity_params(1) ! Softening length
  xmass=gravity_params(2) ! Point mass coordinates
  ymass=gravity_params(3)
  zmass=gravity_params(4)
  dmass=1.0/(emass*(1.0+emass)**2)

!$omp parallel do private(i,rx,ry,rz,rr,dd)
  do i=1,ncell
     rx=0.0d0; ry=0.0d0; rz=0.0d0
     rx=x(i,1)-xmass
#if NDIM>1
     ry=x(i,2)-ymass
#endif
#if NDIM>2
     rz=x(i,3)-zmass
#endif
     rr=sqrt(rx**2+ry**2+rz**2)
     dd=1.0/(rr*(1.0+rr)**2)
     d(i)=MIN(dd,dmass)
  end do

end subroutine rho_ana
