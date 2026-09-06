! DUST-8 FP64 receiver for the FP32 CUDA species/dust ledger.
!
! This module validates only the trial photon ledger.  It deliberately owns
! no RAMSES state and therefore cannot accidentally commit dust abundance,
! temperature, momentum, or re-emission.  The tolerance is expressed in
! units of FP32 epsilon and is scaled by the largest value in each ledger
! identity; this is the documented boundary between CUDA FP32 and the FP64
! native transaction.
module snrt_dust_transaction
  use, intrinsic :: iso_c_binding, only: c_float
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_parameters, only: dp
  implicit none
  private

  integer, parameter, public :: snrt_dust_transaction_ok = 0
  integer, parameter, public :: snrt_dust_transaction_err_shape = 1
  integer, parameter, public :: snrt_dust_transaction_err_state = 2
  integer, parameter, public :: snrt_dust_transaction_err_closure = 3

  integer, parameter, public :: snrt_dust_transaction_hhe_species = 3
  integer, parameter, public :: snrt_dust_transaction_tolerance_ulps = 64
  real(dp), parameter, public :: snrt_dust_transaction_fp32_epsilon = &
       real(epsilon(1.0),dp)
  real(dp), parameter, public :: snrt_dust_transaction_relative_tolerance = &
       real(snrt_dust_transaction_tolerance_ulps,dp) * &
       snrt_dust_transaction_fp32_epsilon

  public :: snrt_dust_validate_ledgers

contains

  subroutine snrt_dust_validate_ledgers(raw_group, hhe_species_group, &
       dust_group, returned_group, assigned_group, max_relative_error, ierr)
    ! All input arrays are CUDA outputs in c_float.  Arithmetic below is
    ! intentionally promoted to dp before the receiver identities are tested.
    ! hhe_species_group has layout (leaf,group,species), matching the CUDA ABI.
    real(c_float), intent(in) :: raw_group(:,:), hhe_species_group(:,:,:)
    real(c_float), intent(in) :: dust_group(:,:), returned_group(:,:)
    real(c_float), intent(in) :: assigned_group(:,:)
    real(dp), intent(out) :: max_relative_error
    integer, intent(out) :: ierr
    integer :: ileaf, igroup
    real(dp) :: raw, assigned, dust, returned, hhe_sum, scale
    real(dp) :: raw_residual, component_residual

    max_relative_error = 0.0d0
    ierr = snrt_dust_transaction_ok
    if (size(raw_group,1) /= size(dust_group,1) .or. &
         size(raw_group,2) /= size(dust_group,2) .or. &
         size(raw_group,1) /= size(returned_group,1) .or. &
         size(raw_group,2) /= size(returned_group,2) .or. &
         size(raw_group,1) /= size(assigned_group,1) .or. &
         size(raw_group,2) /= size(assigned_group,2) .or. &
         size(hhe_species_group,1) /= size(raw_group,1) .or. &
         size(hhe_species_group,2) /= size(raw_group,2) .or. &
         size(hhe_species_group,3) /= snrt_dust_transaction_hhe_species) then
       ierr = snrt_dust_transaction_err_shape
       return
    end if

    if (size(raw_group) == 0) return
    if (any(.not. ieee_is_finite(raw_group)) .or. &
         any(.not. ieee_is_finite(hhe_species_group)) .or. &
         any(.not. ieee_is_finite(dust_group)) .or. &
         any(.not. ieee_is_finite(returned_group)) .or. &
         any(.not. ieee_is_finite(assigned_group)) .or. &
         any(raw_group < 0.0_c_float) .or. &
         any(hhe_species_group < 0.0_c_float) .or. &
         any(dust_group < 0.0_c_float) .or. &
         any(returned_group < 0.0_c_float) .or. &
         any(assigned_group < 0.0_c_float)) then
       ierr = snrt_dust_transaction_err_state
       return
    end if

    do igroup = 1, size(raw_group,2)
       do ileaf = 1, size(raw_group,1)
          raw = real(raw_group(ileaf,igroup),dp)
          assigned = real(assigned_group(ileaf,igroup),dp)
          dust = real(dust_group(ileaf,igroup),dp)
          returned = real(returned_group(ileaf,igroup),dp)
          hhe_sum = sum(real(hhe_species_group(ileaf,igroup,:),dp))
          scale = max(1.0d-300, abs(raw), abs(assigned), abs(dust), &
               abs(returned), maxval(abs(real(hhe_species_group(ileaf,igroup,:),dp))))
          raw_residual = abs(raw-assigned-returned) / scale
          component_residual = abs(assigned-hhe_sum-dust) / scale
          max_relative_error = max(max_relative_error, raw_residual, component_residual)
       end do
    end do
    if (max_relative_error > snrt_dust_transaction_relative_tolerance) &
         ierr = snrt_dust_transaction_err_closure
  end subroutine snrt_dust_validate_ledgers

end module snrt_dust_transaction
