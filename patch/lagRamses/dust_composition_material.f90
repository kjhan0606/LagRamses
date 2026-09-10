! Pure bulk graphite/olivine vibrational material comparison. A common dust
! temperature, not separate grain temperatures, stochastic heating or optics.
module dust_composition_material
  use dust_mass_physics
  implicit none
  include 'dust_dl01_composition_data.inc'
  real(real64),parameter :: dl01_hot_max=3000d0
  real(real64),parameter :: dust_olivine_molar_mass=24.305d0+55.845d0+28.085d0+4*15.999d0
contains
  subroutine dust_composition_cold(t,curve,ierr)
    ! Zero-point-subtracted DL01 bulk limit below the first 5 K knot.
    ! Integral x^2/(exp(x)-1) dx = 2*zeta(3); the 3D integral is pi^4/15.
    ! theta/T >= 100 here: the omitted exponential tail is negligible.
    ! Anchor each material to its existing knot (quadrature roundoff only).
    ! This is NOT a finite-molecule PAH heat capacity or a temperature floor.
    real(real64),intent(in)::t
    real(real64),intent(out)::curve(2)
    integer,intent(out)::ierr
    real(real64)::r,s2,s3
    curve=0;ierr=1
    if(.not.ieee_is_finite(t))return
    if(t<0.or.t>dl01_t(1))return
    r=t/dl01_t(1)
    s2=8d0*1.2020569031595942854d0/500d0**2
    s3=acos(-1d0)**4*dl01_t(1)/(5d0*1500d0**3)
    curve(1)=dl01_carbon(1)*r**3
    curve(2)=dl01_silicate(1)*r**3*(s2+s3*r)/(s2+s3)
    ierr=0
  end subroutine

  subroutine dust_olivine_phase_reference(latent,ierr)
    ! An explicit ideal endmember mixture, NOT activation energy from a
    ! kinetic fit. RH95 USGS Bulletin 2131 pp32--33: Hf(Fa)=-1478.2 and
    ! Hf(Fo)=-2173.0 kJ/mol. Atomic gas Hf at 298.15 K [Mg,Fe,Si,O]
    ! = [147.10,416.3,450.00,249.18] kJ/mol (NIST/NBS compilations).
    ! Calibrate constant phase energy to that atomization enthalpy at 298 K
    ! using this receiver's DL01 solid U and ideal monatomic gas H=5/2 NkT.
    ! Ignore excess solid-solution enthalpy, solid pV, electronic excitation
    ! and melting; those are declared approximations, not kinetic barriers.
    real(real64),intent(out)::latent
    integer,intent(out)::ierr
    real(real64)::solid_u(1)
    real(real64),parameter::r_gas=8.31446261815324d7,reference_t=298.15d0
    real(real64),parameter::atomization_kj=147.10d0+416.3d0+450d0+4*249.18d0+(1478.2d0+2173d0)/2
    call dust_composition_curve([reference_t],[0d0,1d0],1d0,solid_u,ierr)
    latent=0
    if(ierr/=0)return
    latent=(atomization_kj*1d10-2.5d0*7*r_gas*reference_t)/dust_olivine_molar_mass+solid_u(1)
    if(.not.ieee_is_finite(latent).or.latent<=0)ierr=1
  end subroutine

  subroutine dust_olivine_phase_identity(values,ierr)
    real(real64),intent(out)::values(15)
    integer,intent(out)::ierr
    real(real64)::latent
    call dust_olivine_phase_reference(latent,ierr)
    values=[1d0,20.06d0,16.29d0,22.23d0,84780d0,74420d0,90590d0, &
         1478.2d0,2173d0,147.10d0,416.3d0,450d0,249.18d0,dust_olivine_molar_mass,latent]
  end subroutine

  logical function dust_material_domain(nodes) result(ok)
    real(real64),intent(in)::nodes(:)
    ok=all(ieee_is_finite(nodes))
    ok=ok.and.all(nodes>=dl01_t(1)).and.all(nodes<=dl01_hot_max)
  end function

  real(real64) function dl01_warm_mode(t,theta,dimension) result(u)
    ! DL01 normalized bulk Debye mode: dimension*y**(dimension-1) dy.
    ! Sixteen-point Gauss-Legendre is used only at T>=300 K (theta/T<=8.35).
    ! No zero-point energy, latent heat, PAH modes or phase stability implied.
    real(real64),intent(in)::t,theta
    integer,intent(in)::dimension
    real(real64),parameter::x(8)=[.095012509837637440d0,.281603550779258913d0, &
         .458016777657227386d0,.617876244402643748d0,.755404408355003034d0, &
         .865631202387831744d0,.944575023073232576d0,.989400934991649933d0]
    real(real64),parameter::w(8)=[.189450610455068496d0,.182603415044923589d0, &
         .169156519395002538d0,.149595988816576733d0,.124628971255533872d0, &
         .095158511682492785d0,.062253523938647893d0,.027152459411754095d0]
    real(real64)::y,z
    integer::j,s
    u=0
    do j=1,8
       do s=-1,1,2
          y=(1+s*x(j))/2;z=theta*y/t
          u=u+w(j)*y**dimension/(exp(z)-1)
       enddo
    enddo
    u=u*.5d0*dimension*dust_kb*theta
  end function

  subroutine dust_material_identity(hot,v)
    ! Preserve the original cold checkpoint byte-for-byte. A hot contract
    ! binds the same sized identity sampled across its extended material domain.
    logical,intent(in)::hot
    real(real64),intent(out)::v(1+3*dl01_n)
    real(real64)::t(dl01_n),c(dl01_n),s(dl01_n)
    integer::j,status
    v=[1d0,dl01_t,dl01_carbon,dl01_silicate]
    if(.not.hot)return
    do j=1,dl01_n
       t(j)=dl01_t(1)*(dl01_hot_max/dl01_t(1))**(real(j-1,real64)/(dl01_n-1))
    enddo
    t(dl01_n)=dl01_hot_max
    call dust_composition_curve(t,[1d0,0d0],1d0,c,status)
    if(status/=0)error stop 'invalid hot carbon material identity'
    call dust_composition_curve(t,[0d0,1d0],1d0,s,status)
    if(status/=0)error stop 'invalid hot silicate material identity'
    v=[2d0,t,c,s]
  end subroutine

  subroutine dust_composition_curve(nodes,grains,normalization,curve,ierr)
    ! grains may be masses/densities/fractions; only their ratio is used.
    ! normalization=1 for erg/g, reference mass/H for erg/reference H.
    real(real64),intent(in)::nodes(:),grains(2),normalization
    real(real64),intent(out)::curve(size(nodes))
    integer,intent(out)::ierr
    real(real64)::f,w,uc,us
    integer::i,k
    curve=0;ierr=1
    if(.not.all(ieee_is_finite(nodes)).or..not.all(ieee_is_finite(grains)))return
    if(any(grains<0).or..not.ieee_is_finite(normalization).or.normalization<=0)return
    if(.not.dust_material_domain(nodes))return
    f=.5d0;if(sum(grains)>0)f=grains(1)/sum(grains)
    do i=1,size(nodes)
       if(nodes(i)>dl01_t(dl01_n))then
          uc=(dl01_warm_mode(nodes(i),863d0,2)+2*dl01_warm_mode(nodes(i),2504d0,2))/ &
               (12.011d0*1.66053906660d-24)
          us=(2*dl01_warm_mode(nodes(i),500d0,2)+dl01_warm_mode(nodes(i),1500d0,3))/ &
               (((24.305d0+55.845d0+28.085d0+4*15.999d0)/7)*1.66053906660d-24)
          curve(i)=normalization*(f*uc+(1-f)*us)
          cycle
       endif
       k=1
       do while(k<dl01_n-1)
          if(nodes(i)<=dl01_t(k+1))exit
          k=k+1
       enddo
       w=log(nodes(i)/dl01_t(k))/log(dl01_t(k+1)/dl01_t(k))
       uc=(1-w)*dl01_carbon(k)+w*dl01_carbon(k+1)
       us=(1-w)*dl01_silicate(k)+w*dl01_silicate(k+1)
       curve(i)=normalization*(f*uc+(1-f)*us)
    enddo
    if(any(.not.ieee_is_finite(curve)).or.any(curve<=0))return
    ierr=0
  end subroutine

  subroutine dust_composition_area(bins,normalization,area,ierr)
    ! Geometric pi*a^2 per reference dust mass: independent of optical Q.
    real(real64),intent(in)::bins(4),normalization
    real(real64),intent(out)::area
    integer,intent(out)::ierr
    integer::i,j
    area=0;ierr=1
    if(.not.dust_mass_parameters_ok().or.any(.not.ieee_is_finite(bins)))return
    if(any(bins<0).or..not.ieee_is_finite(normalization).or.normalization<=0)return
    if(sum(bins)>0)then
       do i=1,2
          do j=1,2
             area=area+.75d0*bins(2*i+j-2)/sum(bins)/(dust_size_density(i)*dust_size_radius_cm(j))
          enddo
       enddo
       area=area*normalization
    endif
    if(.not.ieee_is_finite(area))return
    ierr=0
  end subroutine
end module
