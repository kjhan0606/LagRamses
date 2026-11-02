!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
  real(dp)::source_cutoff,taper_inner,taper_width,taper_r,taper_z
  integer,parameter::max_profile=4096
  integer,save::nstar=0,ngas=0
  real(dp),save::rstar(max_profile),sstar(max_profile),rgas(max_profile),sgas(max_profile)
  logical,save::profiles_loaded=.false.
  real(dp)::star_sigma,gas_sigma,gas_taper,star_average,gas_average

  if(gravity_type==-5)then
     ! Fixed tabulated surface density, normalized to code units.
     ! Parameters: star/gas exponential heights, gas taper start/end.
     if(any(.not.ieee_is_finite(gravity_params(1:4))))then
        write(*,'(A)') 'ERROR: nonfinite tabulated profile parameters'
        call clean_stop
     end if
     if(gravity_params(1)<=0d0.or.gravity_params(2)<=0d0.or. &
          gravity_params(3)<0d0.or.gravity_params(4)<=gravity_params(3))then
        write(*,'(A)') 'ERROR: invalid tabulated profile parameters'
        call clean_stop
     end if
!$omp critical(ngc2403_profile_initialization)
     if(.not.profiles_loaded)then
        call read_profile('star_profile_code.dat',nstar,rstar,sstar)
        call read_profile('hi_profile_code.dat',ngas,rgas,sgas)
        profiles_loaded=.true.
     end if
!$omp end critical(ngc2403_profile_initialization)
!$omp parallel do private(i,rx,ry,rz,radius_cyl,star_sigma,gas_sigma,gas_taper,star_average,gas_average)
     do i=1,ncell
        rx=x(i,1)-0.5d0*boxlen
        ry=x(i,2)-0.5d0*boxlen
        rz=x(i,3)-0.5d0*boxlen
        radius_cyl=sqrt(rx*rx+ry*ry)
        star_sigma=interpolate_profile(radius_cyl,nstar,rstar,sstar)
        if(radius_cyl>rstar(nstar))star_sigma=0d0
        gas_sigma=interpolate_profile(radius_cyl,ngas,rgas,sgas)
        gas_taper=min(1d0,max(0d0,(gravity_params(4)-radius_cyl)/ &
             (gravity_params(4)-gravity_params(3))))
        star_average=exponential_cell_average(rz,dx,gravity_params(1))
        gas_average=exponential_cell_average(rz,dx,gravity_params(2))
        d(i)=star_sigma*star_average+gas_sigma*gas_taper*gas_average
     end do
!$omp end parallel do
     return
  end if

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
     if(gravity_params(10)<0d0)then
        write(*,'(A)') 'ERROR: negative analytic disc cutoff radius'
        call clean_stop
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
     source_cutoff=gravity_params(10)
     taper_inner=0.8d0*source_cutoff
     taper_width=0.2d0*source_cutoff
!$omp parallel do private(i,rx,ry,rz,rr,dd,radius_cyl,star_height,gas_height,halo_x,taper_r,taper_z)
     do i=1,ncell
        rx=x(i,1)-0.5d0*boxlen
        ry=x(i,2)-0.5d0*boxlen
        rz=x(i,3)-0.5d0*boxlen
        radius_cyl=sqrt(rx**2+ry**2)
        star_height=1d0/cosh(rz/gravity_params(3))**2
        gas_height=1d0/cosh(rz/gravity_params(6))**2
        dd=star_norm*exp(-radius_cyl/gravity_params(2))*star_height &
             +gas_norm*exp(-radius_cyl/gravity_params(5))*gas_height
        ! Optional compact support leaves the measured inner disc unchanged.
        ! A cosine taper from 0.8*cutoff to cutoff matches the same finite
        ! baryonic source on periodic boxes of different physical lengths.
        if(source_cutoff>0d0)then
           taper_r=0.5d0*(1d0+cos((fourpi/4d0)* &
                min(1d0,max(0d0,(radius_cyl-taper_inner)/taper_width))))
           taper_z=0.5d0*(1d0+cos((fourpi/4d0)* &
                min(1d0,max(0d0,(abs(rz)-taper_inner)/taper_width))))
           dd=dd*taper_r*taper_z
        end if
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

contains

  subroutine read_profile(filename,n,radius,surface)
    character(len=*),intent(in)::filename
    integer,intent(out)::n
    real(dp),intent(out)::radius(max_profile),surface(max_profile)
    integer::unit,ios,j
    open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0)then
       write(*,'(A,A)') 'ERROR: cannot open density profile: ',filename
       call clean_stop
    end if
    read(unit,*,iostat=ios)n
    if(ios/=0)then
       write(*,'(A,A)') 'ERROR: invalid density profile count: ',filename
       call clean_stop
    end if
    if(n<2.or.n>max_profile)then
       write(*,'(A,A)') 'ERROR: density profile count out of range: ',filename
       call clean_stop
    end if
    do j=1,n
       read(unit,*,iostat=ios)radius(j),surface(j)
       if(ios/=0)then
          write(*,'(A,A)') 'ERROR: invalid density profile row: ',filename
          call clean_stop
       end if
       if(.not.ieee_is_finite(radius(j)).or..not.ieee_is_finite(surface(j)))then
          write(*,'(A,A)') 'ERROR: nonfinite density profile row: ',filename
          call clean_stop
       end if
       if(radius(j)<0d0.or.surface(j)<0d0)then
          write(*,'(A,A)') 'ERROR: negative density profile value: ',filename
          call clean_stop
       end if
       if(j>1)then
          if(radius(j)<=radius(j-1))then
             write(*,'(A,A)') 'ERROR: unordered density profile radius: ',filename
             call clean_stop
          end if
       end if
    end do
    close(unit)
  end subroutine read_profile

  real(dp) function interpolate_profile(r,n,radius,surface) result(value)
    integer,intent(in)::n
    real(dp),intent(in)::r,radius(max_profile),surface(max_profile)
    integer::lo,hi,mid
    real(dp)::weight
    if(r<=radius(1))then
       value=surface(1)
    else if(r>=radius(n))then
       value=surface(n)
    else
       lo=1;hi=n
       do while(hi-lo>1)
          mid=(lo+hi)/2
          if(r<radius(mid))then
             hi=mid
          else
             lo=mid
          end if
       end do
       weight=(r-radius(lo))/(radius(hi)-radius(lo))
       value=surface(lo)*(1d0-weight)+surface(hi)*weight
    end if
  end function interpolate_profile

  real(dp) function exponential_cell_average(z,cell,h) result(value)
    real(dp),intent(in)::z,cell,h
    real(dp)::left,right,fl,fr
    left=z-cell/2d0;right=z+cell/2d0
    fl=sign(0.5d0,left)*(1d0-exp(-abs(left)/h))
    fr=sign(0.5d0,right)*(1d0-exp(-abs(right)/h))
    value=(fr-fl)/cell
  end function exponential_cell_average
end subroutine rho_ana
