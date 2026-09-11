! One convex reconstruction limiter for the hosts AND their isotope subsets.
! Ordering: total metal, H He C N O Ne Mg Si S Ca Fe, Al26 Fe60.
module stellar_radioactive_transport
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: radioactive_limit_states,radioactive_constraints
contains
  function radioactive_constraints(x) result(g)
    real(real64),intent(in)::x(14)
    real(real64)::g(16)
    g(:14)=x
    g(15)=x(12)-x(14)
    g(16)=x(1)-sum(x(4:12))-x(13)
  end function

  subroutine radioactive_limit_states(center,states,ierr)
    real(real64),intent(in)::center(14)
    real(real64),intent(inout)::states(:,:) ! 14 x faces or children
    integer,intent(out)::ierr
    real(real64)::g0(16),gf(16),theta,tol
    integer::f,j
    ierr=1
    if(size(states,1)/=14)return
    if(.not.all(ieee_is_finite(center)).or..not.all(ieee_is_finite(states)))return
    g0=radioactive_constraints(center)
    tol=64*epsilon(1d0)*maxval(abs(center))
    if(minval(center)<0.or.minval(g0)<-tol)return
    theta=1
    do f=1,size(states,2)
       gf=radioactive_constraints(states(:,f))
       do j=1,16
          if(gf(j)>=0)cycle
          if(g0(j)<=0)then
             theta=0
          else
             theta=min(theta,g0(j)/(g0(j)-gf(j)))
          endif
       enddo
    enddo
    if(theta<1)then
       ! Round toward the valid center, not an isotope mass floor or clip.
       theta=theta*(1-16*epsilon(1d0))
       do f=1,size(states,2)
          states(:,f)=center+theta*(states(:,f)-center)
       enddo
    endif
    ierr=0
  end subroutine
end module stellar_radioactive_transport
