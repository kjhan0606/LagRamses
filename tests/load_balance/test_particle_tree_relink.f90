program test_particle_tree_relink
  use amr_parameters
  use amr_commons
  use pm_commons
  use amr_index, only: icell_of
  implicit none
  integer::ind,icell,ncellmax
  real(dp),dimension(1:3)::skip_loc

  ngridmax=9
  amr_block_size=ngridmax
  ncoarse=0
  ncellmax=twotondim*amr_block_size
  allocate(xg(ngridmax,ndim),son(ncellmax))
  allocate(xp(5,ndim),headp(ngridmax),tailp(ngridmax),numbp(ngridmax))
  allocate(nextp(5),prevp(5))

  xg=0.5d0
  son=0
  headp=0
  tailp=0
  numbp=0
  do ind=1,twotondim
     icell=icell_of(1,ind)
     son(icell)=ind+1
  end do

  xp=0.25d0
  xp(2,1)=0.75d0
  xp(3,2)=0.75d0
  xp(4,1:ndim)=0.75d0
  headp(1)=1
  tailp(1)=4
  numbp(1)=4
  nextp=(/2,3,4,0,0/)
  prevp=(/0,1,2,3,0/)

  ! Model a child list contributed by a reception parent handled earlier in
  ! the icpu loop.  The direct relink must append instead of overwriting it.
  headp(2)=5
  tailp(2)=5
  numbp(2)=1

  ! Leave child 3 unrefined: particle 3 must remain on the parent while the
  ! other particles are relinked to child grids 1, 2, and 8.
  son(icell_of(1,3))=0
  skip_loc=0d0
  call kill_tree_grid_relink(1,0.5d0,1d0,skip_loc)

  call assert_list('parent',1,(/3/),1)
  call assert_list('child 1 append',2,(/5,1/),2)
  call assert_list('child 2',3,(/2/),1)
  call assert_list('child 8',9,(/4/),1)
  call assert_equal('particle conservation',sum(numbp),5)

  write(*,*)'PASS: particle tree one-pass relink'

contains

  subroutine assert_list(label,igrid,expected,nexpected)
    character(len=*),intent(in)::label
    integer,intent(in)::igrid,nexpected
    integer,dimension(:),intent(in)::expected
    integer::i,ipart,last

    call assert_equal(trim(label)//' count',numbp(igrid),nexpected)
    ipart=headp(igrid)
    last=0
    do i=1,nexpected
       call assert_equal(trim(label)//' particle',ipart,expected(i))
       call assert_equal(trim(label)//' prev link',prevp(ipart),last)
       last=ipart
       ipart=nextp(ipart)
    end do
    call assert_equal(trim(label)//' terminator',ipart,0)
    call assert_equal(trim(label)//' tail',tailp(igrid),last)
  end subroutine assert_list

  subroutine assert_equal(label,actual,expected)
    character(len=*),intent(in)::label
    integer,intent(in)::actual,expected

    if(actual/=expected)then
       write(*,*)'FAIL: ',trim(label),' expected=',expected,' actual=',actual
       error stop 1
    end if
  end subroutine assert_equal

end program test_particle_tree_relink
