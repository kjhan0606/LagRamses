! Phase 0 cumulative source increment.
!
! This module converts cumulative SSP quantities into the incremental source
! for one hydrodynamic interval.  It is intentionally independent of the AMR
! deposition code; the returned stellar_source_t is the only object that the
! deposition layer needs to consume.

module stellar_source_increment
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dust_mass_physics, only: dust_gas_elements,olivine_fraction,dust_pah_hc,dust_injection_bins
  use stellar_enrichment_config, only: stellar_dp
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, stellar_source_t, clear_cumulative, clear_source, &
       cumulative_difference
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_ssp_sources, only: integrate_ssp_channel, ssp_source_ok
  implicit none

  private
  integer, parameter, public :: source_increment_ok = 0
  integer, parameter, public :: source_increment_err_argument = 1
  integer, parameter, public :: source_increment_err_ssp = 2
  integer, parameter, public :: source_increment_err_nonfinite = 4
  integer, parameter, public :: source_increment_err_negative = 8

  public :: integrate_ssp_channel_increment
  public :: source_condensed_donors,source_condensed_momentum

contains

  subroutine integrate_ssp_channel_increment(table, population, channel_id, &
       previous_age_gyr, current_age_gyr, mass_min, mass_max, n_mass_bins, &
       source, ierr, current_cumulative)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    integer, intent(in) :: channel_id, n_mass_bins
    real(stellar_dp), intent(in) :: previous_age_gyr, current_age_gyr
    real(stellar_dp), intent(in) :: mass_min, mass_max
    type(stellar_source_t), intent(out) :: source
    integer, intent(out) :: ierr
    type(stellar_cumulative_t), intent(out), optional :: current_cumulative

    type(stellar_cumulative_t) :: previous, current
    integer :: ssp_ierr

    call clear_source(source)
    if (present(current_cumulative)) call clear_cumulative(current_cumulative)
    ierr = source_increment_ok

    if (.not. ieee_is_finite(previous_age_gyr) .or. &
         .not. ieee_is_finite(current_age_gyr) .or. &
         .not. ieee_is_finite(mass_min) .or. .not. ieee_is_finite(mass_max) .or. &
         previous_age_gyr < 0.0_stellar_dp .or. &
         current_age_gyr < previous_age_gyr .or. &
         mass_min <= 0.0_stellar_dp .or. mass_max <= mass_min .or. &
         n_mass_bins <= 0) then
       ierr = source_increment_err_argument
       return
    end if

    call integrate_ssp_channel(table, population, channel_id, previous_age_gyr, &
         mass_min, &
         mass_max, n_mass_bins, previous, ssp_ierr)
    if (ssp_ierr /= ssp_source_ok) then
       if(table%high_mass_ready)write(*,*) 'High-mass previous cumulative failed: ',channel_id,ssp_ierr
       ierr = source_increment_err_ssp
       return
    end if

    call integrate_ssp_channel(table, population, channel_id, &
         current_age_gyr, mass_min, mass_max, n_mass_bins, current, ssp_ierr)
    if (ssp_ierr /= ssp_source_ok) then
       if(table%high_mass_ready)write(*,*) 'High-mass current cumulative failed: ',channel_id,ssp_ierr
       ierr = source_increment_err_ssp
       return
    end if

    call cumulative_difference(current, previous, source)
    if (.not. source_values_finite(source)) then
       call clear_source(source)
       ierr = source_increment_err_nonfinite
       return
    end if
    if (.not. source_values_nonnegative(source)) then
       call clear_source(source)
       ierr = source_increment_err_negative
       return
    end if
    ! Condensation was already evaluated on source-node segments before SSP
    ! integration. Retain that result here, before channels are summed.
    source%channel_condensed_mass(channel_id,1:2)=source%dust_species
    if (present(current_cumulative)) current_cumulative = current
  end subroutine integrate_ssp_channel_increment

  subroutine source_condensed_donors(source,iron_mass,pah_mass,ierr)
    ! Preserve the existing aggregate condensation yields exactly. Only
    ! resolve their ephemeral mass provenance for component momentum.
    ! Fe uses residual Fe after olivine. The aggregate PAH minimum consumes
    ! H and residual C from the already mixed ejecta: sample each element's
    ! contributing channels proportionally, NOT each channel's total mass.
    type(stellar_source_t),intent(inout)::source
    real(stellar_dp),intent(in)::iron_mass,pah_mass
    integer,intent(out)::ierr
    real(stellar_dp)::work(size(source%channel_returned_mass),4)
    real(stellar_dp)::available(size(source%channel_returned_mass)),part(size(source%channel_returned_mass))
    real(stellar_dp)::target(4),tol,total
    integer::m,c
    ierr=source_increment_err_negative;work=source%channel_condensed_mass;work(:,3:4)=0
    target=[source%dust_species,iron_mass,pah_mass]
    if(any(.not.ieee_is_finite(target)).or.any(target<0))return
    if(any(.not.ieee_is_finite(work)).or.any(work<0))return
    do m=1,2
       total=sum(work(:,m));tol=256*epsilon(1d0)*max(total,target(m),tiny(1d0))
       if(abs(total-target(m))>tol)return
       call distribute(target(m),work(:,m),part,ierr)
       if(ierr/=0)return
       work(:,m)=part
    enddo
    available=source%channel_ejected_mass(:,11)-olivine_fraction(11)*work(:,2)
    available(4:)=0 ! No direct Ia/P(P)ISN Fe condensation prescription.
    tol=256*epsilon(1d0)*max(maxval(abs(source%channel_ejected_mass(:,11))),tiny(1d0))
    if(any(available < -tol))then
       ierr=source_increment_err_negative;return
    endif
    call distribute(iron_mass,max(available,0d0),part,ierr)
    if(ierr/=0)return
    work(:,3)=part
    do m=1,2
       if(m==1)then
          available=source%channel_ejected_mass(:,1)
       else
          available=source%channel_ejected_mass(:,3)-work(:,1)
       endif
       available(4:)=0 ! Same ordinary-channel donors used by aggregate PAH.
       tol=256*epsilon(1d0)*max(maxval(abs(source%channel_ejected_mass(:,2*m-1))),tiny(1d0))
       if(any(available < -tol))then
          ierr=source_increment_err_negative;return
       endif
       call distribute(pah_mass*dust_pah_hc(m),max(available,0d0),part,ierr)
       if(ierr/=0)return
       work(:,4)=work(:,4)+part
    enddo
    ierr=source_increment_err_negative
    do c=1,size(work,1)
       tol=512*epsilon(1d0)*max(source%channel_returned_mass(c),sum(work(c,:)),tiny(1d0))
       if(sum(work(c,:))>source%channel_returned_mass(c)+tol)return
    enddo
    source%channel_condensed_mass=work;ierr=source_increment_ok
  contains
    subroutine distribute(amount,reservoir,donated,status)
      real(stellar_dp),intent(in)::amount,reservoir(:)
      real(stellar_dp),intent(out)::donated(:)
      integer,intent(out)::status
      real(stellar_dp)::s,tolerance
      integer::largest
      status=source_increment_err_negative;donated=0;s=sum(reservoir)
      if(any(.not.ieee_is_finite(reservoir)).or.any(reservoir<0))return
      tolerance=512*epsilon(1d0)*max(amount,s,tiny(1d0))
      if(amount>s+tolerance)return
      if(amount>0)then
         if(s<=0)return
         donated=(reservoir/s)*amount
         largest=maxloc(reservoir,dim=1)
         donated(largest)=donated(largest)+amount-sum(donated)
      endif
      status=source_increment_ok
    end subroutine
  end subroutine source_condensed_donors

  subroutine source_condensed_momentum(source,bulk_velocity,scale_mass,scale_momentum,volume, &
       has_iron,has_pah,phase_mass,phase_momentum,ierr)
    type(stellar_source_t),intent(in)::source
    real(stellar_dp),intent(in)::bulk_velocity(3),scale_mass,scale_momentum,volume,phase_mass(:)
    logical,intent(in)::has_iron,has_pah
    real(stellar_dp),intent(inout)::phase_momentum(:,:)
    integer,intent(out)::ierr
    real(stellar_dp)::bins(size(phase_mass),size(source%channel_returned_mass))
    real(stellar_dp)::velocity(3),work(3,size(phase_mass)),s,tol,donor_mass
    integer::c,b,nb,k
    ierr=source_increment_err_argument;nb=4+merge(2,0,has_iron)+merge(1,0,has_pah)
    if(size(phase_mass)/=nb.or.any(shape(phase_momentum)/=[3,nb]))return
    if(min(scale_mass,scale_momentum,volume)<=0)return
    if(any(.not.ieee_is_finite([bulk_velocity,scale_mass,scale_momentum,volume])))return
    if(any(.not.ieee_is_finite(phase_mass)).or.any(phase_mass<0))return
    ierr=source_increment_err_negative;bins=0;work=0
    if(any(.not.ieee_is_finite(source%channel_condensed_mass)).or.any(source%channel_condensed_mass<0))return
    do c=1,size(bins,2)
       bins(1:4,c)=dust_injection_bins(source%channel_condensed_mass(c,1:2));k=4
       if(has_iron)then
          bins(6,c)=source%channel_condensed_mass(c,3);k=6
       endif
       if(has_pah)bins(k+1,c)=source%channel_condensed_mass(c,4)
    enddo
    bins=bins/scale_mass/volume
    do b=1,nb
       s=sum(bins(b,:));tol=1024*epsilon(1d0)*max(s,phase_mass(b),tiny(1d0))
       if(abs(s-phase_mass(b))>tol)return
       ! Only roundoff closure against the unchanged legacy field increment.
       if(s>0)bins(b,:)=bins(b,:)*(phase_mass(b)/s)
    enddo
    do c=1,size(bins,2)
       if(.not.ieee_is_finite(source%channel_returned_mass(c)))return
       if(any(.not.ieee_is_finite(source%channel_momentum(c,:))))return
       donor_mass=source%channel_returned_mass(c)/scale_mass
       if(donor_mass<=0)then
          if(any(bins(:,c)/=0).or.any(source%channel_momentum(c,:)/=0))return
          cycle
       endif
       velocity=bulk_velocity+(source%channel_momentum(c,:)/scale_momentum)/donor_mass
       do b=1,nb
          work(:,b)=work(:,b)+bins(b,c)*velocity
       enddo
    enddo
    if(any(.not.ieee_is_finite(work)))return
    phase_momentum=work;ierr=source_increment_ok
  end subroutine source_condensed_momentum

  logical function source_values_finite(source)
    type(stellar_source_t), intent(in) :: source
    integer :: i, j

    source_values_finite = all(ieee_is_finite(source%dust_species)).and. &
         ieee_is_finite(source%returned_mass) .and. &
         ieee_is_finite(source%energy)
    do i = 1, 3
       source_values_finite = source_values_finite .and. &
            ieee_is_finite(source%momentum(i))
    end do
    do i = 1, size(source%ejected_mass)
       source_values_finite = source_values_finite .and. &
            ieee_is_finite(source%ejected_mass(i)) .and. &
            ieee_is_finite(source%net_yield(i))
    end do
    do i = 1, size(source%channel_returned_mass)
       source_values_finite = source_values_finite .and. &
            ieee_is_finite(source%channel_returned_mass(i)) .and. &
            ieee_is_finite(source%channel_energy(i))
       do j = 1, 3
          source_values_finite = source_values_finite .and. &
               ieee_is_finite(source%channel_momentum(i,j))
       end do
       do j = 1, size(source%channel_ejected_mass, 2)
          source_values_finite = source_values_finite .and. &
               ieee_is_finite(source%channel_ejected_mass(i,j)) .and. &
               ieee_is_finite(source%channel_net_yield(i,j))
       end do
    end do
  end function source_values_finite

  logical function source_values_nonnegative(source)
    type(stellar_source_t), intent(in) :: source
    real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp
    real(stellar_dp) :: scale, tracked_mass, channel_tracked_mass
    real(stellar_dp)::gas(11)
    integer::dust_status
    integer :: channel, element

    source_values_nonnegative = source%returned_mass >= -tolerance .and. &
         source%energy >= -tolerance .and. &
         minval(source%ejected_mass) >= -tolerance .and. &
         minval(source%channel_returned_mass) >= -tolerance .and. &
         minval(source%channel_energy) >= -tolerance .and. &
         minval(source%channel_ejected_mass) >= -tolerance
    if (.not. source_values_nonnegative) return
    if(any(source%dust_species/=0d0))then
       call dust_gas_elements(source%ejected_mass,source%dust_species,gas,dust_status)
       if(dust_status/=0)then
          source_values_nonnegative=.false.;return
       endif
    endif

    tracked_mass = sum(source%ejected_mass)
    scale = max(1.0_stellar_dp, abs(source%returned_mass), abs(tracked_mass))
    if (tracked_mass > source%returned_mass + tolerance * scale) then
       source_values_nonnegative = .false.
       return
    end if
    if (abs(sum(source%channel_returned_mass) - source%returned_mass) > &
         tolerance * scale) then
       source_values_nonnegative = .false.
       return
    end if
    do element = 1, size(source%ejected_mass)
       scale = max(1.0_stellar_dp, abs(source%ejected_mass(element)), &
            abs(sum(source%channel_ejected_mass(:,element))))
       if (abs(sum(source%channel_ejected_mass(:,element)) - &
            source%ejected_mass(element)) > tolerance * scale) then
          source_values_nonnegative = .false.
          return
       end if
    end do
    do channel = 1, size(source%channel_returned_mass)
       channel_tracked_mass = sum(source%channel_ejected_mass(channel,:))
       scale = max(1.0_stellar_dp, &
            abs(source%channel_returned_mass(channel)), &
            abs(channel_tracked_mass))
       if (channel_tracked_mass > source%channel_returned_mass(channel) + &
            tolerance * scale) then
          source_values_nonnegative = .false.
          return
       end if
    end do
  end function source_values_nonnegative

end module stellar_source_increment
