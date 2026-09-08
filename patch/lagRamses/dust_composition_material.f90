! Pure bulk graphite/olivine vibrational material comparison. A common dust
! temperature, not separate grain temperatures, stochastic heating or optics.
module dust_composition_material
  use dust_mass_physics
  implicit none
  include 'dust_dl01_composition_data.inc'
contains
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
    if(any(nodes<dl01_t(1)).or.any(nodes>dl01_t(dl01_n)))return
    f=.5d0;if(sum(grains)>0)f=grains(1)/sum(grains)
    do i=1,size(nodes)
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
