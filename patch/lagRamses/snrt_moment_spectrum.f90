! Bounded two-measure ANGULAR closure for a finite spectral band.
! A/B are endpoint weights, not a two-line spectral absorption model.
! The existing intra-band SED closure still receives angular N and E.
module snrt_moment_spectrum
  use snrt_moment_transport
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: mn_band_split,mn_band_split_angular,mn_band_join,mn_band_reconstruct
contains
  subroutine mn_band_split_angular(n,e,lo,hi,a,b,ierr)
    ! The material receiver independently sums positive N and E packets.
    ! At a band endpoint those two rounded sums may differ by a few ulps.
    ! Admit ONLY that arithmetic envelope, preserving N at the endpoint.
    ! The caller must include the resulting E difference in its numerical
    ! projection receipt (never in gas/dust heating). This is NOT a repair
    ! for an infeasible angular moment vector or a physical out-of-band SED.
    real(mn_dp),intent(in)::n(:),e(:),lo,hi
    real(mn_dp),intent(out)::a(:),b(:)
    integer,intent(out)::ierr
    real(mn_dp)::bound(size(n))
    ierr=mn_bad_input
    if(size(n)<1.or.size(n)/=size(e).or.size(a)/=size(n).or.size(b)/=size(n))return
    if(.not.all(ieee_is_finite([lo,hi])).or.lo<=0.or.hi<=lo)return
    if(any(.not.ieee_is_finite(n)).or.any(.not.ieee_is_finite(e)))return
    if(any(n<0).or.any(e<0).or.any(n==0.and.e/=0))return
    a=(e-lo*n)/(hi-lo);b=(hi*n-e)/(hi-lo)
    bound=8*epsilon(1d0)*max(e,hi*n)/(hi-lo)
    if(any(.not.ieee_is_finite(a)).or.any(.not.ieee_is_finite(b)).or.any(.not.ieee_is_finite(bound)))return
    if(any(a < -bound).or.any(b < -bound))return
    where(a<0)
       a=0;b=n
    endwhere
    where(b<0)
       b=0;a=n
    endwhere
    ierr=mn_ok
  end subroutine

  subroutine mn_band_split(n,e,lo,hi,a,b,ierr)
    real(mn_dp),intent(in)::n(:),e(:),lo,hi
    real(mn_dp),intent(out)::a(:),b(:)
    integer,intent(out)::ierr
    ierr=mn_bad_input
    if(size(n)<1.or.size(n)/=size(e).or.size(a)/=size(n).or.size(b)/=size(n))return
    if(.not.all(ieee_is_finite([lo,hi])).or.lo<=0.or.hi<=lo)return
    if(any(.not.ieee_is_finite(n)).or.any(.not.ieee_is_finite(e)))return
    a=(e-lo*n)/(hi-lo);b=(hi*n-e)/(hi-lo)
    if(a(1)<0.or.b(1)<0)return
    ierr=mn_ok
  end subroutine

  subroutine mn_band_join(a,b,lo,hi,n,e)
    real(mn_dp),intent(in)::a(:),b(:),lo,hi
    real(mn_dp),intent(out)::n(:),e(:)
    n=a+b;e=hi*a+lo*b
  end subroutine

  subroutine mn_band_reconstruct(basis,n,e,lo,hi,angular_n,angular_e,ierr)
    type(mn_basis),intent(in)::basis
    real(mn_dp),intent(in)::n(:),e(:),lo,hi
    real(mn_dp),intent(out)::angular_n(:),angular_e(:)
    integer,intent(out)::ierr
    real(mn_dp)::a(basis%nm),b(basis%nm),ia(basis%nq),ib(basis%nq),hint(basis%nm-1)
    call mn_band_split(n,e,lo,hi,a,b,ierr)
    if(ierr/=mn_ok)return
    hint=0
    call mn_reconstruct(basis,a,ia,ierr,dual_hint=hint)
    if(ierr/=mn_ok)return
    call mn_reconstruct(basis,b,ib,ierr,dual_hint=hint)
    if(ierr/=mn_ok)call mn_reconstruct(basis,b,ib,ierr)
    if(ierr/=mn_ok)return
    call mn_band_join(ia,ib,lo,hi,angular_n,angular_e)
  end subroutine
end module
