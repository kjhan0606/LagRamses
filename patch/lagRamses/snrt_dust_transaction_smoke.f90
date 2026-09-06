program snrt_dust_transaction_smoke
  use, intrinsic :: iso_c_binding, only: c_float
  use amr_parameters, only: dp
  use snrt_dust_transaction, only: snrt_dust_validate_ledgers, &
       snrt_dust_transaction_ok, snrt_dust_transaction_err_state, &
       snrt_dust_transaction_err_shape, snrt_dust_transaction_err_closure, &
       snrt_dust_transaction_relative_tolerance
  implicit none

  integer, parameter :: nleaf = 2, ngroup = 3
  real(c_float) :: raw(nleaf,ngroup), hhe(nleaf,ngroup,3)
  real(c_float) :: dust(nleaf,ngroup), returned(nleaf,ngroup), assigned(nleaf,ngroup)
  real(c_float) :: hhe_bad(nleaf,ngroup,2)
  real(dp) :: relative_error
  integer :: ierr

  raw = 0.0_c_float
  hhe = 0.0_c_float
  dust = 0.0_c_float
  returned = 0.0_c_float
  assigned = 0.0_c_float
  hhe(1,1,:) = (/1.0_c_float, 0.25_c_float, 0.25_c_float/)
  dust(1,1) = 1.5_c_float
  returned(1,1) = 1.0_c_float
  raw(1,1) = 4.0_c_float
  hhe(1,2,:) = (/0.1_c_float, 0.2_c_float, 0.1_c_float/)
  dust(1,2) = 0.6_c_float
  returned(1,2) = 1.0_c_float
  raw(1,2) = 2.0_c_float
  hhe(2,1,:) = (/1.0e-8_c_float, 2.0e-8_c_float, 3.0e-8_c_float/)
  dust(2,1) = 4.0e-8_c_float
  returned(2,1) = 1.0e-8_c_float
  raw(2,1) = 1.1e-7_c_float
  assigned = sum(hhe, dim=3) + dust

  call snrt_dust_validate_ledgers(raw, hhe, dust, returned, assigned, &
       relative_error, ierr)
  if (ierr /= snrt_dust_transaction_ok .or. relative_error > &
       real(snrt_dust_transaction_relative_tolerance)) error stop 1
  write(*,'(a,es12.4)') 'SNRT_DUST_LEDGER_VALIDATION_OK relative_error=', relative_error

  assigned(1,1) = assigned(1,1) + 1.0e-3_c_float
  call snrt_dust_validate_ledgers(raw, hhe, dust, returned, assigned, &
       relative_error, ierr)
  if (ierr /= snrt_dust_transaction_err_closure) error stop 2
  assigned(1,1) = assigned(1,1) - 1.0e-3_c_float

  dust(1,1) = -1.0_c_float
  call snrt_dust_validate_ledgers(raw, hhe, dust, returned, assigned, &
       relative_error, ierr)
  if (ierr /= snrt_dust_transaction_err_state) error stop 3
  dust(1,1) = 1.5_c_float

  hhe_bad = 0.0_c_float
  call snrt_dust_validate_ledgers(raw, hhe_bad, dust, returned, assigned, &
       relative_error, ierr)
  if (ierr /= snrt_dust_transaction_err_shape) error stop 4

  write(*,'(a)') 'SNRT_DUST_LEDGER_VALIDATION_NEGATIVE_OK state=1 shape=1 closure=1'
end program snrt_dust_transaction_smoke
