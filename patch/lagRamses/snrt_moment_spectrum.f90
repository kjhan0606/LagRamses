! Bounded two-measure ANGULAR closure for a finite spectral band.
! A/B are endpoint weights, not a two-line spectral absorption model.
! The existing intra-band SED closure still receives angular N and E.
module snrt_moment_spectrum
  use snrt_moment_transport
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: mn_band_split,mn_band_split_angular,mn_band_join,mn_band_reconstruct
  integer,parameter,public :: mn_band_err_nonfinite=11,mn_band_err_negative=12
  integer,parameter,public :: mn_band_err_zero_n_energy=13,mn_band_err_outside_edges=14
  ! Bound accumulation error for the receiver's 128-node dust spectrum plus
  ! the products and final N/E combination (four ULP/node plus scalar work).
  integer,parameter :: mn_band_roundoff_ulps=528
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
    real(mn_dp)::lower(size(n)),upper(size(n)),lower_gap(size(n)),upper_gap(size(n))
    real(mn_dp)::scale(size(n)),roundoff(size(n))
    ierr=mn_bad_input
    if(size(n)<1.or.size(n)/=size(e).or.size(a)/=size(n).or.size(b)/=size(n))return
    a=0._mn_dp;b=0._mn_dp
    if(.not.all(ieee_is_finite([lo,hi])).or.lo<=0.or.hi<=lo)return
    if(any(.not.ieee_is_finite(n)).or.any(.not.ieee_is_finite(e)))then
       ierr=mn_band_err_nonfinite;return
    endif
    if(any(n<0).or.any(e<0))then
       ierr=mn_band_err_negative;return
    endif
    if(any(n==0.and.e/=0))then
       ierr=mn_band_err_zero_n_energy;return
    endif
    lower=lo*n;upper=hi*n
    if(any(.not.ieee_is_finite(lower)).or.any(.not.ieee_is_finite(upper)))then
       ierr=mn_band_err_nonfinite;return
    endif
    lower_gap=e-lower;upper_gap=upper-e
    scale=max(abs(e),abs(lower),abs(upper))
    ! The material receiver sums up to 128 positive spectral-node contributions
    ! independently for N and E. The operation-count envelope covers the
    ! observed 9-ULP endpoint cancellation. Its energy correction is booked as
    ! a projection residual, never deposited as heat. ULP spacing remains
    ! nonzero when epsilon*scale would underflow for subnormal packets.
    roundoff=real(mn_band_roundoff_ulps,mn_dp)*(scale-nearest(scale,-1._mn_dp))
    if(any(.not.ieee_is_finite(roundoff)))then
       ierr=mn_band_err_nonfinite;return
    endif
    if(any(lower_gap < -roundoff).or.any(upper_gap < -roundoff))then
       ierr=mn_band_err_outside_edges;return
    endif
    a=lower_gap/(hi-lo)
    if(any(.not.ieee_is_finite(a)))then
       ierr=mn_band_err_nonfinite;a=0._mn_dp;b=0._mn_dp;return
    endif
    a=max(0._mn_dp,min(n,a))
    b=n-a
    where(lower_gap<0._mn_dp)
       a=0._mn_dp;b=n
    endwhere
    where(upper_gap<0._mn_dp)
       a=n;b=0._mn_dp
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
