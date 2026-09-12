program dark_probe
  use iso_c_binding
  use omp_lib
  use snrt_chimes
  use snrt_thermochemistry,only:snrt_secondary_tables_load_from_environment
  implicit none
  character(500)::path,dir,groups(9)
  character(2)::num
  real(c_double)::elements(11),old(157),next(157),ctl(9),t,elapsed,start,measured(11),q
  real(c_double),parameter::temps(3)=[100d0,1d4,1d6],dens(3)=[1d-4,.2d0,100d0]
  real(c_double)::thread_state(157,32),thread_t(32),thread_elapsed(32)
  integer::thread_status(32)
  integer::i,j,r,s,atomic
  call snrt_secondary_tables_load_from_environment(s)
  if(s/=0)stop 1
  call get_environment_variable('SNRT_CHIMES_MAIN_DATA',path)
  call get_environment_variable('SNRT_CHIMES_GROUP_DIR',dir)
  do i=1,9
    write(num,'(I2.2)')i
    groups(i)=trim(dir)//'/group_'//num//'.hdf5'//c_null_char
  enddo
  s=chimes_initialize(trim(path)//c_null_char,9,groups,500)
  if(s/=0)stop 2
  s=chimes_set_expansion(.01d0)
  if(s/=0)stop 3
  elements=0;elements(1)=1;elements(2)=.079d0
  elements(3)=1d-8;elements(5)=2d-8
  s=chimes_neutral(elements,old)
  if(s/=0)stop 4
  ! Neutral/molecular and ionized carriers in the same valid H inventory.
  old(1)=.01d0;old(2)=.97d0;old(3)=.01d0;old(138)=.01d0
  do j=1,3
    do i=1,2
      ctl=[dens(j),temps(i),272.7d0,1d10,1d20,0d0,1d0,0d0,.01d0]
      atomic=merge(1,0,i==3)
      start=omp_get_wtime()
!$omp parallel do num_threads(4) private(s,t,next,elapsed)
      do r=1,32
        s=chimes_cell_transition_dark(ctl,elements,old,atomic,t,next,elapsed)
        thread_status(r)=s;thread_state(:,r)=next;thread_t(r)=t;thread_elapsed(r)=elapsed
        if(s/=0)then
          print *,'DARK_FAILURE',i,j,s;stop 5
        endif
      enddo
!$omp end parallel do
      s=chimes_cell_transition_dark(ctl,elements,old,atomic,t,next,elapsed)
      if(s/=0)stop 7
      do r=1,32
        if(thread_status(r)/=s.or.thread_t(r)/=t.or.thread_elapsed(r)/=elapsed)stop 8
        if(any(thread_state(:,r)/=next))stop 9
      enddo
      write(*,'(A,2I3,F12.6)')'TIME ',i,j,omp_get_wtime()-start
      s=chimes_budget(next,measured,q)
      if(s/=0.or.abs(q)>1d-8.or.any(abs(measured-elements)>1d-8*max(elements,1d-20)))stop 6
      write(*,'(A,2I3,158ES25.16)')'STATE ',i,j,t,next
    enddo
  enddo
  print *,'DARK_MIXED_THREAD_PARITY_CONSERVATION_PASS'
end program
