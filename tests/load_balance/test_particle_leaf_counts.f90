program test_particle_leaf_counts
  use amr_parameters
  use amr_commons
  use pm_commons
  implicit none
  integer,dimension(1:twotondim)::npart_leaf,ndm_leaf

  ngridmax=1
  allocate(xg(1,ndim))
  allocate(xp(4,ndim),headp(1),numbp(1),nextp(4))
  allocate(idp(4),ptypep(4))

  xg=0.5d0
  xp=0.25d0
  xp(2,1)=0.75d0
  xp(3,2)=0.75d0
  xp(4,1:3)=0.75d0
  headp(1)=1
  numbp(1)=4
  nextp=(/2,3,4,0/)
  idp=(/1_8,2_8,3_8,4_8/)
  ptypep=PTYPE_DM
  ptypep(3)=PTYPE_STAR

  call count_particles_by_leaf(1,npart_leaf,ndm_leaf)

  call assert_equal('leaf 1 particles',npart_leaf(1),1)
  call assert_equal('leaf 2 particles',npart_leaf(2),1)
  call assert_equal('leaf 3 particles',npart_leaf(3),1)
  call assert_equal('leaf 8 particles',npart_leaf(8),1)
  call assert_equal('leaf 1 DM',ndm_leaf(1),1)
  call assert_equal('leaf 2 DM',ndm_leaf(2),1)
  call assert_equal('leaf 3 excludes star',ndm_leaf(3),0)
  call assert_equal('leaf 8 DM',ndm_leaf(8),1)
  call assert_equal('particle conservation',sum(npart_leaf),4)

  write(*,*)'PASS: particle leaf counts'

contains

  subroutine assert_equal(label,actual,expected)
    character(len=*),intent(in)::label
    integer,intent(in)::actual,expected

    if(actual/=expected)then
       write(*,*)'FAIL: ',trim(label),' expected=',expected,' actual=',actual
       error stop 1
    end if
  end subroutine assert_equal

end program test_particle_leaf_counts
