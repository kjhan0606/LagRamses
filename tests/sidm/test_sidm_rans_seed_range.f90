program test_sidm_rans_seed_range
  use random
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer :: test_seeds(5),i,j
  integer,dimension(2,IRandNumSize) :: states_a,states_b
  integer,dimension(IRandNumSize) :: state
  real(kind=8) :: u

  test_seeds=(/1,4095,4096,1042654254,2147471302/)
  do i=1,size(test_seeds)
     call rans(2,test_seeds(i),states_a)
     call rans(2,test_seeds(i),states_b)
     if(any(states_a/=states_b)) error stop 'RANS is not deterministic'
     if(any(states_a<0) .or. any(states_a>4095)) &
          error stop 'RANS returned a limb outside base-4096 range'
     if(all(states_a(1,:)==states_a(2,:))) &
          error stop 'RANS did not split rank streams'
     state=states_a(1,:)
     do j=1,1000
        call ranf(state,u)
        if(.not.ieee_is_finite(u) .or. u<0.0d0 .or. u>=1.0d0) &
             error stop 'RANF returned a value outside [0,1)'
        if(any(state<0) .or. any(state>4095)) &
             error stop 'RANF returned a limb outside base-4096 range'
     end do
  end do
  write(*,'(A)') 'PASS_RANS_HIGH_SEED_BASE4096'
end program test_sidm_rans_seed_range
