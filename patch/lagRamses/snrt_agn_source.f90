module snrt_agn_source
  ! Photon budgets for the S_N AGN source path.
  !
  ! delta_inflow_mass_code is an increment of supplied inflow over the
  ! current radiation update, expressed in code mass units.  It must never be
  ! a sink seed mass, retained BH mass, or an unbounded Bondi supply.  Sink
  ! bookkeeping and AMR-cell deposition are intentionally kept outside this
  ! module.
  use, intrinsic :: iso_c_binding, only: c_float
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_parameters, only: dp
  implicit none

  real(dp), parameter :: snrt_c_cgs = 2.99792458d10
  real(dp), parameter :: snrt_ev_to_erg = 1.602176634d-12

contains

  subroutine snrt_agn_photon_budget(delta_inflow_mass_code, mass_unit_g, &
       delta_t_s, radiative_efficiency, ionizing_fraction, &
       mean_photon_energy_ev, luminosity_erg_s, emitted_photons)
    real(dp), intent(in) :: delta_inflow_mass_code, mass_unit_g, delta_t_s
    real(dp), intent(in) :: radiative_efficiency, ionizing_fraction
    real(dp), intent(in) :: mean_photon_energy_ev
    real(dp), intent(out) :: luminosity_erg_s, emitted_photons
    real(dp) :: radiated_energy_erg

    luminosity_erg_s = 0.0d0
    emitted_photons = 0.0d0

    if (.not. ieee_is_finite(delta_inflow_mass_code) .or. &
         .not. ieee_is_finite(mass_unit_g) .or. .not. ieee_is_finite(delta_t_s) .or. &
         .not. ieee_is_finite(radiative_efficiency) .or. &
         .not. ieee_is_finite(ionizing_fraction) .or. &
         .not. ieee_is_finite(mean_photon_energy_ev)) return
    if (delta_inflow_mass_code <= 0.0d0) return
    if (mass_unit_g <= 0.0d0 .or. delta_t_s <= 0.0d0) return
    ! The shared efficiency helper admits only [0,1); keep this public source
    ! boundary fail-closed as well so a future caller cannot bypass that
    ! convention by passing unity or a super-unity coefficient directly.
    if (radiative_efficiency <= 0.0d0 .or. radiative_efficiency >= 1.0d0 .or. &
         ionizing_fraction <= 0.0d0) return
    if (mean_photon_energy_ev <= 0.0d0) return

    radiated_energy_erg = radiative_efficiency * delta_inflow_mass_code * &
         mass_unit_g * snrt_c_cgs**2
    call snrt_agn_photon_budget_energy(radiated_energy_erg,delta_t_s,ionizing_fraction, &
         mean_photon_energy_ev,luminosity_erg_s,emitted_photons)
  end subroutine snrt_agn_photon_budget

  subroutine snrt_agn_photon_budget_energy(energy_erg,delta_t_s,fraction,mean_ev,luminosity,photons,ierr)
    ! Fractions are escaped fractions of the full accepted bolometric energy.
    real(dp), intent(in) :: energy_erg,delta_t_s,fraction,mean_ev
    real(dp), intent(out) :: luminosity,photons
    integer, optional, intent(out) :: ierr
    luminosity=0d0; photons=0d0
    if(present(ierr))ierr=1
    if(.not.all(ieee_is_finite([energy_erg,delta_t_s,fraction,mean_ev])))return
    if(energy_erg<0d0.or.delta_t_s<=0d0.or.fraction<0d0.or.fraction>1d0.or.mean_ev<=0d0)return
    luminosity=energy_erg/delta_t_s
    photons=(fraction*energy_erg)/(mean_ev*snrt_ev_to_erg)
    if(.not.all(ieee_is_finite([luminosity,photons])))then
       luminosity=0d0; photons=0d0
       return
    endif
    if(present(ierr))ierr=0
  end subroutine snrt_agn_photon_budget_energy

  pure subroutine snrt_agn_source_commit(pending_erg,source_ok)
    real(dp), intent(inout) :: pending_erg
    logical, intent(in) :: source_ok
    ! The caller invokes this only after the enclosing source -> RT/chemistry
    ! -> dust transaction has passed its global commit.  Until then the
    ! accepted event remains pending in memory; this routine clears the whole
    ! accepted event, rather than debiting a partial amount.
    if(source_ok)pending_erg=0d0
  end subroutine snrt_agn_source_commit

  subroutine snrt_agn_isotropic_packet(emitted_photons, angular_weights, &
       directional_photons)
    ! Split an angle-integrated photon count using the S_N quadrature.
    ! The caller retains emitted_photons as the conservation ledger; this
    ! routine only provides its directional representation for transport.
    real(dp), intent(in) :: emitted_photons
    real(dp), intent(in) :: angular_weights(:)
    real(dp), intent(out) :: directional_photons(size(angular_weights))
    real(dp) :: total_weight

    directional_photons = 0.0d0
    if (.not. ieee_is_finite(emitted_photons) .or. &
         any(.not. ieee_is_finite(angular_weights))) return
    if (emitted_photons <= 0.0d0) return

    total_weight = sum(max(angular_weights, 0.0d0))
    if (.not. ieee_is_finite(total_weight) .or. total_weight <= 0.0d0) return

    directional_photons = emitted_photons * max(angular_weights, 0.0d0) / &
         total_weight
  end subroutine snrt_agn_isotropic_packet

  subroutine snrt_agn_photons_to_density_code(emitted_photons, &
       cell_volume_code, length_unit_cm, n_h_unit_cm3, photon_density_code)
    ! Convert a cell-integrated physical photon count to the S_N state
    ! variable n_gamma / n_H,unit.  This keeps the c_float transport state
    ! in a usable dynamic range while leaving the integral photon ledger in
    ! double precision at the source boundary.
    real(dp), intent(in) :: emitted_photons, cell_volume_code
    real(dp), intent(in) :: length_unit_cm, n_h_unit_cm3
    real(dp), intent(out) :: photon_density_code

    photon_density_code = 0.0d0
    if (.not. ieee_is_finite(emitted_photons) .or. &
         .not. ieee_is_finite(cell_volume_code) .or. &
         .not. ieee_is_finite(length_unit_cm) .or. &
         .not. ieee_is_finite(n_h_unit_cm3)) return
    if (emitted_photons <= 0.0d0) return
    if (cell_volume_code <= 0.0d0 .or. length_unit_cm <= 0.0d0) return
    if (n_h_unit_cm3 <= 0.0d0) return

    photon_density_code = emitted_photons / (cell_volume_code * &
         length_unit_cm**3 * n_h_unit_cm3)
  end subroutine snrt_agn_photons_to_density_code

  subroutine snrt_agn_deposit_isotropic(state, slot, group, emitted_photons, &
       cell_volume_code, length_unit_cm, n_h_unit_cm3, angular_weights, &
       deposited_density_code, ierr)
    ! State layout is (direction, group, slot).  The source is supplied as
    ! an integrated physical photon count and is converted at this boundary
    ! to n_gamma / n_H,unit before its angular split.
    real(c_float), intent(inout) :: state(:,:,:)
    integer, intent(in) :: slot, group
    real(dp), intent(in) :: emitted_photons, cell_volume_code
    real(dp), intent(in) :: length_unit_cm, n_h_unit_cm3
    real(dp), intent(in) :: angular_weights(:)
    real(dp), intent(out) :: deposited_density_code
    integer, intent(out) :: ierr
    real(dp), allocatable :: directional_density(:)
    real(dp) :: state_limit

    ierr = 0
    deposited_density_code = 0.0d0
    if (slot < 1 .or. slot > size(state,3)) then
       ierr = 1
       return
    end if
    if (group < 1 .or. group > size(state,2)) then
       ierr = 2
       return
    end if
    if (size(state,1) /= size(angular_weights)) then
       ierr = 3
       return
    end if

    ! This public single-group entry point is retained for compatibility, but
    ! its production callers must not silently turn malformed inputs into a
    ! zero source.  Error 5 is a validation/range failure; no state mutation
    ! has occurred when it is returned.
    if (.not. ieee_is_finite(emitted_photons) .or. emitted_photons < 0.0d0 .or. &
         .not. ieee_is_finite(cell_volume_code) .or. cell_volume_code <= 0.0d0 .or. &
         .not. ieee_is_finite(length_unit_cm) .or. length_unit_cm <= 0.0d0 .or. &
         .not. ieee_is_finite(n_h_unit_cm3) .or. n_h_unit_cm3 <= 0.0d0 .or. &
         any(.not. ieee_is_finite(angular_weights)) .or. &
         any(angular_weights < 0.0d0) .or. &
         .not. all(ieee_is_finite(real(state(:,group,slot),dp))) ) then
       ierr = 5
       return
    end if
    if (.not. ieee_is_finite(sum(angular_weights)) .or. &
         sum(angular_weights) <= 0.0d0) then
       ierr = 5
       return
    end if

    call snrt_agn_photons_to_density_code(emitted_photons, cell_volume_code, &
         length_unit_cm, n_h_unit_cm3, deposited_density_code)
    if (.not. ieee_is_finite(deposited_density_code) .or. &
         deposited_density_code < 0.0d0) then
       ierr = 5
       return
    end if
    if (deposited_density_code == 0.0d0) return

    allocate(directional_density(size(angular_weights)))
    call snrt_agn_isotropic_packet(deposited_density_code, angular_weights, &
         directional_density)
    state_limit = real(huge(0.0_c_float), dp) * 0.5d0
    if (.not. all(ieee_is_finite(directional_density)) .or. &
         maxval(abs(real(state(:,group,slot),dp) + directional_density)) > &
         state_limit) then
       ierr = 4
       return
    end if
    state(:,group,slot) = state(:,group,slot) + real(directional_density,c_float)
  end subroutine snrt_agn_deposit_isotropic

  subroutine snrt_agn_deposit_transaction(state, slot, emitted_photons, &
       cell_volume_code, length_unit_cm, n_h_unit_cm3, angular_weights, &
       deposited_density_code, ierr, persistent_energy_shift, group_mean_energy_ev,quantize_source,source_mean_energy_ev)
    ! Prepare and commit all spectral groups for one source atomically.
    !
    ! ``state`` is (direction, group, slot).  Every validation and conversion
    ! is performed against temporary double-precision arrays first.  Thus a
    ! failure in a later group cannot leave an earlier group deposited, which
    ! is required before the caller advances its accounted-mass marker.
    real(c_float), intent(inout) :: state(:,:,:)
    integer, intent(in) :: slot
    real(dp), intent(in) :: emitted_photons(:), cell_volume_code
    real(dp), intent(in) :: length_unit_cm, n_h_unit_cm3
    real(dp), intent(in) :: angular_weights(:)
    real(dp), intent(out) :: deposited_density_code
    integer, intent(out) :: ierr
    ! Optional paired state: actual E = reference(g)*N + shift, in photon
    ! CODE density * eV. Supply BOTH arrays only when subsequent operators
    ! transport the correction. Error 6 rejects a malformed/unrepresentable
    ! energy state atomically; neither photons nor corrections are published.
    real(dp), optional, intent(inout) :: persistent_energy_shift(:,:,:)
    real(dp), optional, intent(in) :: group_mean_energy_ev(:)
    ! Actual source E/Q, independent of the fixed reference used to encode E.
    real(dp), optional, intent(in) :: source_mean_energy_ev(:)
    ! Band two-moment mode injects the representable FP32 packet's energy.
    ! Do not retain an energy-only photon when a source tail rounds to N=0.
    ! Bound the total source-energy rounding, never alter the physical table.
    logical,optional,intent(in)::quantize_source
    integer :: group, allocation_status
    real(dp) :: total_weight, state_limit,source_energy,rounding_energy
    real(dp), allocatable :: group_density(:), directional_density(:,:)
    real(c_float), allocatable :: next_number(:,:)
    real(dp), allocatable :: next_shift(:,:), actual_energy(:)
    real(dp), allocatable :: injection_mean(:)

    ierr = 0
    deposited_density_code = 0.0d0
    if(present(source_mean_energy_ev))then
       if(.not.present(persistent_energy_shift))then
          ierr=6;return
       endif
       if(size(source_mean_energy_ev)/=size(state,2))then
          ierr=6;return
       endif
       if(any(.not.ieee_is_finite(source_mean_energy_ev)).or.any(source_mean_energy_ev<=0d0))then
          ierr=6;return
       endif
    endif
    if(present(quantize_source))then
       if(quantize_source.and..not.present(persistent_energy_shift))then
          ierr=6;return
       endif
    endif
    if (slot < 1 .or. slot > size(state,3)) then
       ierr = 1
       return
    end if
    if (size(emitted_photons) <= 0 .or. size(state,2) /= size(emitted_photons)) then
       ierr = 2
       return
    end if
    if (size(state,1) /= size(angular_weights)) then
       ierr = 3
       return
    end if
    if(present(persistent_energy_shift).neqv.present(group_mean_energy_ev))then
       ierr=6
       return
    endif
    if(present(persistent_energy_shift))then
       injection_mean=group_mean_energy_ev
       if(present(source_mean_energy_ev))injection_mean=source_mean_energy_ev
       if(any(shape(persistent_energy_shift)/=shape(state)).or.size(group_mean_energy_ev)/=size(state,2))then
          ierr=6
          return
       endif
       if(any(.not.ieee_is_finite(group_mean_energy_ev)).or.any(group_mean_energy_ev<=0.0_dp).or. &
            any(.not.ieee_is_finite(persistent_energy_shift(:,:,slot))).or.any(state(:,:,slot)<0.0_c_float))then
          ierr=6
          return
       endif
       allocate(next_number(size(state,1),size(state,2)),next_shift(size(state,1),size(state,2)), &
            actual_energy(size(state,1)),stat=allocation_status)
       if(allocation_status/=0)then
          ierr=6
          return
       endif
       do group=1,size(state,2)
          actual_energy=group_mean_energy_ev(group)*real(state(:,group,slot),dp)+ &
               persistent_energy_shift(:,group,slot)
          if(any(.not.ieee_is_finite(actual_energy)).or.any(actual_energy<0.0_dp).or. &
               any(state(:,group,slot)==0.0_c_float.and.persistent_energy_shift(:,group,slot)/=0.0_dp))then
             ierr=6
             return
          endif
       enddo
    endif
    if (.not. ieee_is_finite(cell_volume_code) .or. cell_volume_code <= 0.0d0 .or. &
         .not. ieee_is_finite(length_unit_cm) .or. length_unit_cm <= 0.0d0 .or. &
         .not. ieee_is_finite(n_h_unit_cm3) .or. n_h_unit_cm3 <= 0.0d0 .or. &
         any(.not. ieee_is_finite(emitted_photons)) .or. &
         any(emitted_photons < 0.0d0) .or. &
         any(.not. ieee_is_finite(angular_weights)) .or. &
         any(angular_weights < 0.0d0) .or. sum(angular_weights) <= 0.0d0) then
       ierr = 5
       return
    end if
    if (.not. all(ieee_is_finite(real(state(:,:,slot),dp)))) then
       ierr = 5
       return
    end if

    total_weight = sum(angular_weights)
    if (.not. ieee_is_finite(total_weight) .or. total_weight <= 0.0d0) then
       ierr = 5
       return
    end if
    allocate(group_density(size(emitted_photons)), &
         directional_density(size(angular_weights),size(emitted_photons)))
    group_density = 0.0d0
    directional_density = 0.0d0
    state_limit = real(huge(0.0_c_float),dp) * 0.5d0

    do group = 1, size(emitted_photons)
       call snrt_agn_photons_to_density_code(emitted_photons(group), &
            cell_volume_code, length_unit_cm, n_h_unit_cm3, group_density(group))
       if (.not. ieee_is_finite(group_density(group)) .or. &
            group_density(group) < 0.0d0) then
          ierr = 5
          deallocate(group_density, directional_density)
          return
       end if
       call snrt_agn_isotropic_packet(group_density(group), angular_weights, &
            directional_density(:,group))
       if (.not. all(ieee_is_finite(directional_density(:,group))) .or. &
            .not. all(ieee_is_finite(real(state(:,group,slot),dp) + &
            directional_density(:,group))) .or. &
            maxval(abs(real(state(:,group,slot),dp) + directional_density(:,group))) > &
            state_limit) then
          ierr = 4
          deallocate(group_density, directional_density)
          return
       end if
    end do

    if(present(persistent_energy_shift))then
       if(.not.ieee_is_finite(sum(group_density)))then
          ierr=6
          return
       endif
       if(present(quantize_source))then
          if(quantize_source)then
             source_energy=sum(injection_mean*group_density)
             rounding_energy=0
             do group=1,size(emitted_photons)
                next_number(:,group)=real(directional_density(:,group),c_float)
                rounding_energy=rounding_energy+injection_mean(group)* &
                     sum(abs(directional_density(:,group)-real(next_number(:,group),dp)))
             enddo
             if(.not.ieee_is_finite(source_energy).or..not.ieee_is_finite(rounding_energy).or. &
                  rounding_energy>8*epsilon(0.0_c_float)*source_energy)then
                ierr=6;return
             endif
             directional_density=real(next_number,dp)
             group_density=sum(directional_density,dim=1)
          endif
       endif
       do group=1,size(emitted_photons)
          ! Keep the established FP32 number arithmetic, including rounding
          ! the increment. Charge the EXACT FP64 directional source to energy.
          next_number(:,group)=state(:,group,slot)+real(directional_density(:,group),c_float)
          ! Algebraically oldN+increment-newN, with the FP32 number difference
          ! formed first so a small source is not lost in a large FP64 oldN.
          next_shift(:,group)=persistent_energy_shift(:,group,slot)+group_mean_energy_ev(group)* &
               ((real(state(:,group,slot),dp)-real(next_number(:,group),dp))+directional_density(:,group))
          if(present(source_mean_energy_ev))then
             next_shift(:,group)=persistent_energy_shift(:,group,slot)+group_mean_energy_ev(group)* &
                  (real(state(:,group,slot),dp)-real(next_number(:,group),dp))+ &
                  injection_mean(group)*directional_density(:,group)
          endif
          actual_energy=group_mean_energy_ev(group)*real(next_number(:,group),dp)+next_shift(:,group)
          if(any(.not.ieee_is_finite(next_shift(:,group))).or.any(.not.ieee_is_finite(actual_energy)).or. &
               any(actual_energy<0.0_dp).or.any(next_number(:,group)==0.0_c_float.and.next_shift(:,group)/=0.0_dp))then
             ierr=6
             return
          endif
       enddo
       ! The complete source is accepted or rejected across ALL groups.
       state(:,:,slot)=next_number
       persistent_energy_shift(:,:,slot)=next_shift
       deposited_density_code=sum(group_density)
       return
    endif

    ! Legacy commit arithmetic is unchanged when paired arguments are absent.
    do group = 1, size(emitted_photons)
       state(:,group,slot) = state(:,group,slot) + &
            real(directional_density(:,group),c_float)
    end do
    deposited_density_code = sum(group_density)
    deallocate(group_density, directional_density)
  end subroutine snrt_agn_deposit_transaction

end module snrt_agn_source
