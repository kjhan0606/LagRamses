program stellar_imf_reference
  use stellar_enrichment_config, only: stellar_dp
  use stellar_ssp_sources, only: calculate_imf_normalization
  implicit none

  integer :: imf_id, ierr
  real(stellar_dp) :: normalization

  do imf_id = 0, 3
     call emit_case(imf_id, 0.08_stellar_dp, 120.0_stellar_dp)
     call emit_case(imf_id, 1.0_stellar_dp, 120.0_stellar_dp)
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
