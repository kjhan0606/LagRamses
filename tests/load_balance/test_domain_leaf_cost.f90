program test_domain_leaf_cost
  use amr_parameters
  implicit none

  mem_weight_grid=2268
  mem_weight_part=12

  memory_balance=.true.
  cost_weighting=.true.
  call assert_cost('memory empty leaf',domain_leaf_cost(0,0,99_8,1d0),284_8)
  call assert_cost('memory particle leaf',domain_leaf_cost(3,7,99_8,1d0),320_8)

  memory_balance=.false.
  cost_weighting=.true.
  call assert_cost('work empty leaf',domain_leaf_cost(0,0,1_8,1d0),10_8)
  call assert_cost('work particle leaf',domain_leaf_cost(8,0,1_8,1d0),18_8)
  call assert_cost('work SIDM pairs',domain_leaf_cost(8,3,1_8,1d0),42_8)
  call assert_cost('work level mesh scale',domain_leaf_cost(0,0,1_8,2d0),20_8)
  call assert_cost('work subcycles',domain_leaf_cost(8,0,4_8,1d0),72_8)

  cost_weighting=.false.
  call assert_cost('work without subcycles',domain_leaf_cost(8,0,4_8,1d0),18_8)

  remap_thresh=0.05d0
  lb_remap_min_interval=8
  lb_remap_horizon=16
  lb_remap_safety=1.2d0
  call assert_logical('remap blocked by interval', &
       work_remap_is_economic(0.2d0,400d0,80d0,7),.false.)
  call assert_logical('remap repays cost', &
       work_remap_is_economic(0.2d0,400d0,80d0,8),.true.)
  call assert_logical('remap does not repay cost', &
       work_remap_is_economic(0.06d0,100d0,1000d0,8),.false.)
  call assert_logical('first remap allowed without timing history', &
       work_remap_is_economic(0.2d0,0d0,0d0,8),.true.)

  sidm_npart_min=2
  sidm=.false.
  call assert_int('SIDM pairs disabled for DMO',domain_sidm_pair_count(10),0)
  sidm=.true.
  call assert_int('SIDM sampled pairs',domain_sidm_pair_count(10),5)
  call assert_int('SIDM minimum occupancy',domain_sidm_pair_count(1),0)

  write(*,*)'PASS: domain_leaf_cost'

contains

  subroutine assert_cost(label,actual,expected)
    character(len=*),intent(in)::label
  integer(kind=8),intent(in)::actual,expected

    if(actual/=expected)then
       write(*,*)'FAIL: ',trim(label),' expected=',expected,' actual=',actual
       error stop 1
    end if
  end subroutine assert_cost

  subroutine assert_int(label,actual,expected)
    character(len=*),intent(in)::label
    integer,intent(in)::actual,expected

    if(actual/=expected)then
       write(*,*)'FAIL: ',trim(label),' expected=',expected,' actual=',actual
       error stop 1
    end if
  end subroutine assert_int

  subroutine assert_logical(label,actual,expected)
    character(len=*),intent(in)::label
    logical,intent(in)::actual,expected

    if(actual.neqv.expected)then
       write(*,*)'FAIL: ',trim(label),' expected=',expected,' actual=',actual
       error stop 1
    end if
  end subroutine assert_logical

end program test_domain_leaf_cost
