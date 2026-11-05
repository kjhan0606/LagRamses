program sink_periodic_geometry_test
  use amr_parameters, only: dp
  use sink_periodic_geometry, only: sink_periodic_minimum_image
  implicit none

  real(dp), parameter :: box_length = 1.0d-3
  real(dp), parameter :: half_box = 0.5_dp * box_length
  integer :: failures

  failures = 0

  call expect('positive exact antipode', half_box, -half_box)
  call expect('negative exact antipode', -half_box, -half_box)
  call expect('positive inner neighbor', nearest(half_box, -1.0_dp), &
       nearest(half_box, -1.0_dp) - box_length)
  call expect('positive outer neighbor', nearest(half_box, 1.0_dp), &
       nearest(half_box, 1.0_dp) - box_length)
  call expect('negative inner neighbor', nearest(-half_box, 1.0_dp), &
       nearest(-half_box, 1.0_dp))
  call expect('negative outer neighbor', nearest(-half_box, -1.0_dp), &
       nearest(-half_box, -1.0_dp))

  ! The first representable values outside the one-ULP tie window use the
  ! original interior/exterior wrap rules.
  call expect('first excluded positive interior', &
       nearest(nearest(half_box, -1.0_dp), -1.0_dp), &
       nearest(nearest(half_box, -1.0_dp), -1.0_dp))
  call expect('first excluded negative interior', &
       nearest(nearest(-half_box, 1.0_dp), 1.0_dp), &
       nearest(nearest(-half_box, 1.0_dp), 1.0_dp))
  call expect('first excluded positive exterior', &
       nearest(nearest(half_box, 1.0_dp), 1.0_dp), &
       nearest(nearest(half_box, 1.0_dp), 1.0_dp) - box_length)
  call expect('first excluded negative exterior', &
       nearest(nearest(-half_box, -1.0_dp), -1.0_dp), &
       nearest(nearest(-half_box, -1.0_dp), -1.0_dp) + box_length)

  call expect('diagnostic job 1271425 positive one-ULP sample', &
       5.0000000000000012d-4, 5.0000000000000012d-4 - box_length)
  call expect('diagnostic job 1271443 positive one-ULP sample', &
       4.9999999999999990d-4, 4.9999999999999990d-4 - box_length)
  call expect('ordinary positive interior', 0.25_dp * box_length, &
       0.25_dp * box_length)
  call expect('ordinary negative interior', -0.25_dp * box_length, &
       -0.25_dp * box_length)
  call expect('ordinary positive exterior', 0.75_dp * box_length, &
       -0.25_dp * box_length)
  call expect('ordinary negative exterior', -0.75_dp * box_length, &
       0.25_dp * box_length)

  if (failures /= 0) then
     write(*,'(A,I0)') 'SINK_PERIODIC_GEOMETRY_FAIL count=', failures
     stop 1
  end if
  write(*,'(A,I0)') 'SINK_PERIODIC_GEOMETRY_PASS cases=', 16

contains

  subroutine expect(label, displacement, expected)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: displacement, expected
    real(dp) :: actual

    actual = sink_periodic_minimum_image(displacement, box_length)
    if (actual /= expected) then
       failures = failures + 1
       write(*,'(A,": displacement=",ES25.16E3," expected=",ES25.16E3," actual=",ES25.16E3)') &
            trim(label), displacement, expected, actual
    end if
  end subroutine expect

end program sink_periodic_geometry_test
