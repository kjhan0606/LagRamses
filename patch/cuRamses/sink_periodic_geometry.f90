module sink_periodic_geometry
  use amr_parameters, only: dp
  implicit none
  private

  public :: sink_periodic_minimum_image

contains

  pure function sink_periodic_minimum_image(displacement, box_length) result(image)
    real(dp), intent(in) :: displacement, box_length
    real(dp) :: image, half_box, tie_lower, tie_upper

    image = displacement
    half_box = 0.5_dp * box_length

    ! At the antipode, choose the negative image deterministically. Include
    ! exactly one representable neighbor on each side of |L/2| so roundoff in
    ! equivalent periodic coordinates cannot select opposite images.
    if (box_length > 0.0_dp) then
       tie_lower = nearest(half_box, -1.0_dp)
       tie_upper = nearest(half_box,  1.0_dp)
       if (abs(displacement) >= tie_lower .and. abs(displacement) <= tie_upper) then
          if (displacement > 0.0_dp) image = displacement - box_length
          return
       end if
    end if

    ! Preserve the existing single-wrap behavior away from the antipode.
    if (image > half_box) image = image - box_length
    if (image < -half_box) image = image + box_length
  end function sink_periodic_minimum_image

end module sink_periodic_geometry
