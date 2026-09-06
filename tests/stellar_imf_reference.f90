program stellar_imf_reference
  use stellar_enrichment_config, only: stellar_dp
  use stellar_ssp_sources, only: calculate_imf_normalization, &
       calculate_imf_mass_fraction, evaluate_imf, imf_miller_scalo
  implicit none

  integer :: imf_id, ierr, interval
  real(stellar_dp) :: normalization, fraction, total, local_norm
  real(stellar_dp), parameter :: edges(4) = (/0.08d0, 1.0d0, 10.0d0, 120.0d0/)

  do imf_id = 0, 4
     call emit_case(imf_id, 0.08_stellar_dp, 120.0_stellar_dp)
     call emit_case(imf_id, 1.0_stellar_dp, 120.0_stellar_dp)
  end do

  call calculate_imf_normalization(imf_miller_scalo, edges(1), edges(4), normalization, ierr)
  if (ierr /= 0) error stop 2
  total = 0.0_stellar_dp
  do interval = 1, 3
     call calculate_imf_mass_fraction(imf_miller_scalo, edges(1), edges(4), &
          edges(interval), edges(interval+1), fraction, ierr)
     if (ierr /= 0) error stop 3
     call calculate_imf_normalization(imf_miller_scalo, edges(interval), &
          edges(interval+1), local_norm, ierr)
     if (ierr /= 0) error stop 4
     if (abs(fraction-normalization/local_norm) > 1.0e-13_stellar_dp) error stop 5
     total = total + fraction
  end do
  if (abs(total-1.0_stellar_dp) > 1.0e-13_stellar_dp) error stop 6
  do interval = 2, 3
     if (abs(evaluate_imf(edges(interval)*(1.0d0-1.0d-10), imf_miller_scalo) / &
          evaluate_imf(edges(interval)*(1.0d0+1.0d-10), imf_miller_scalo)-1.0d0) &
          > 1.0d-8) error stop 7
  end do

contains

  subroutine emit_case(imf_id, mass_min, mass_max)
    integer, intent(in) :: imf_id
    real(stellar_dp), intent(in) :: mass_min, mass_max

    call calculate_imf_normalization(imf_id, mass_min, mass_max, &
         normalization, ierr)
    if (ierr /= 0) error stop 1
    write(*,'(I0,1X,ES24.16,1X,ES24.16,1X,ES24.16)') &
         imf_id, mass_min, mass_max, normalization
  end subroutine emit_case
end program stellar_imf_reference
