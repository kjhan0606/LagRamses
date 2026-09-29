program dark_probe
  use iso_c_binding
  use omp_lib
  use snrt_chimes
  use snrt_thermochemistry,only:snrt_secondary_tables_load_from_environment
  implicit none
  character(500)::path,dir,groups(9),replay_path
  character(8192)::replay_line
  character(64)::target_subdt_text,target_repeats_text,molecular_case
  character(2)::num
  real(c_double)::elements(11),old(157),next(157),ctl(9),t,elapsed,start,measured(11),q
  real(c_double),parameter::temps(3)=[100d0,1d4,1d6],dens(3)=[1d-4,.2d0,100d0]
  real(c_double)::thread_state(157,32),thread_t(32),thread_elapsed(32)
  real(c_double)::target_subdt
  integer::thread_status(32)
  integer::i,j,r,s,atomic,target_status,target_repeats,repeat_index,replay_unit,replay_status,molecular_status
  integer::capture_rank,capture_cell,capture_nd
  real(c_double)::capture_number,capture_energy
  logical::have_meta,have_controls,have_elements,have_species
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
  ! Replay exactly one captured *live* dark input, without assuming that a
  ! manufactured neutral cell represents the integrated simulation state.
  replay_path=''
  call get_environment_variable('SNRT_CHIMES_REPLAY_FILE',replay_path,status=replay_status)
  if(replay_status==0.and.len_trim(replay_path)>0)then
    open(newunit=replay_unit,file=trim(replay_path),status='old',action='read',iostat=s)
    if(s/=0)stop 15
    have_meta=.false.;have_controls=.false.;have_elements=.false.;have_species=.false.
    do
      read(replay_unit,'(A)',iostat=s)replay_line
      if(s<0)exit
      if(s/=0)stop 16
      if(index(replay_line,'SNRT_CHIMES_CAPTURE_META ')==1)then
        read(replay_line(len('SNRT_CHIMES_CAPTURE_META ')+1:),*,iostat=s) &
             capture_rank,capture_cell,capture_nd,capture_number,capture_energy
        if(s/=0)stop 16
        have_meta=.true.
      else if(index(replay_line,'SNRT_CHIMES_CAPTURE_CONTROLS ')==1)then
        read(replay_line(len('SNRT_CHIMES_CAPTURE_CONTROLS ')+1:),*,iostat=s)ctl
        if(s/=0)stop 16
        have_controls=.true.
      else if(index(replay_line,'SNRT_CHIMES_CAPTURE_ELEMENTS ')==1)then
        read(replay_line(len('SNRT_CHIMES_CAPTURE_ELEMENTS ')+1:),*,iostat=s)elements
        if(s/=0)stop 16
        have_elements=.true.
      else if(index(replay_line,'SNRT_CHIMES_CAPTURE_SPECIES ')==1)then
        read(replay_line(len('SNRT_CHIMES_CAPTURE_SPECIES ')+1:),*,iostat=s)old
        if(s/=0)stop 16
        have_species=.true.
      endif
    enddo
    close(replay_unit)
    if(.not.(have_meta.and.have_controls.and.have_elements.and.have_species))stop 17
    if(capture_rank/=1.or.capture_cell<1.or.capture_nd<1)stop 17
    if(capture_number/=0.or.capture_energy/=0)stop 18
    if(ctl(2)>=chimes_molecular_temperature_max())stop 19
    s=chimes_budget(old,measured,q)
    if(s/=0)stop 20
    if(abs(q)>1d-10.or.any(abs(measured-elements)>1d-8*max(elements,1d-20)))stop 20
    start=omp_get_wtime()
    s=chimes_cell_transition_dark(ctl,elements,old,0,t,next,elapsed)
    if(s/=0)then
      print *,'CHIMES_REPLAY_FAILURE status=',s
      stop 21
    endif
    print *,'CHIMES_REPLAY_PASS cell=',capture_cell,' wall_s=',omp_get_wtime()-start, &
         ' T_K=',t,' elapsed_s=',elapsed
    write(*,'(A,158ES25.16)')'REPLAY_STATE ',t,next
    stop
  endif
  elements=0;elements(1)=1;elements(2)=.079d0
  elements(3)=1d-8;elements(5)=2d-8
  s=chimes_neutral(elements,old)
  if(s/=0)stop 4
  ! Exercise the nonzero line-cooling paths as well as an OH-only case.
  molecular_case=''
  call get_environment_variable('SNRT_CHIMES_MOLECULAR_CASE',molecular_case,status=molecular_status)
  if(molecular_status==0)then
    if(trim(molecular_case)/='all'.and.trim(molecular_case)/='oh_only')stop 22
    elements(3)=1d-4;elements(5)=2d-4
    s=chimes_neutral(elements,old)
    if(s/=0)stop 23
    old(138)=.1d0;old(141)=1d-5
    old(2)=1d0-2d0*old(138)-old(141)
    old(24)=elements(5)-old(141)
    if(trim(molecular_case)=='all')then
      old(149)=1d-5;old(142)=1d-5
      old(2)=old(2)-2d0*old(142)
      old(8)=elements(3)-old(149)
      old(24)=old(24)-old(149)-old(142)
    endif
    s=chimes_budget(old,measured,q)
    if(s/=0.or.abs(q)>1d-10.or.any(abs(measured-elements)>1d-8*max(elements,1d-20)))stop 24
    ctl=[100d0,100d0,272.7d0,1d7,1d20,0d0,1d0,0d0,.01d0]
    start=omp_get_wtime()
    s=chimes_cell_transition_dark(ctl,elements,old,0,t,next,elapsed)
    if(s/=0)stop 25
    s=chimes_budget(next,measured,q)
    if(s/=0.or.abs(q)>1d-8.or.any(abs(measured-elements)>1d-8*max(elements,1d-20)))stop 26
    print *,'MOLECULAR_PASS case=',trim(molecular_case),' wall_s=',omp_get_wtime()-start
    write(*,'(A,158ES25.16)')'MOLECULAR_STATE ',t,next
    stop
  endif
  ! Optional one-cell cost probe at the actual M_N material substep interval.
  ! The default multi-cell parity test below is unchanged.
  call get_environment_variable('SNRT_CHIMES_TARGET_SUBDT_S',target_subdt_text,status=target_status)
  if(target_status==0)then
    read(target_subdt_text,*,iostat=s)target_subdt
    if(s/=0)stop 10
    if(target_subdt<=0)stop 10
    target_repeats=1
    call get_environment_variable('SNRT_CHIMES_TARGET_REPEATS',target_repeats_text,status=target_status)
    if(target_status==0)then
      read(target_repeats_text,*,iostat=s)target_repeats
      if(s/=0)stop 14
      if(target_repeats<1.or.target_repeats>256)stop 14
    endif
    elements=0;elements(1)=1d0
    elements(2)=.24d0/(4d0*.76d0)
    elements(3)=1.8523959236974623d-11/(12d0*.76d0)
    elements(4)=4.949817957217135d-12/(14d0*.76d0)
    elements(5)=4.4901016865472994d-11/(16d0*.76d0)
    elements(6)=1.0301037584282488d-11/(20d0*.76d0)
    elements(7)=4.5084894300446575d-12/(24d0*.76d0)
    elements(8)=5.087631701721149d-12/(28d0*.76d0)
    elements(9)=2.261416834318439d-12/(32d0*.76d0)
    elements(10)=4.4786304397262327d-13/(40d0*.76d0)
    elements(11)=9.018767345995897d-12/(56d0*.76d0)
    s=chimes_neutral(elements,old)
    if(s/=0)stop 11
    ctl=[.2d0,165.4d0,272.7d0,target_subdt,3.013357013175163d21,0d0,1d0,0d0,.01d0]
    start=omp_get_wtime()
    do repeat_index=1,target_repeats
      s=chimes_cell_transition_dark(ctl,elements,old,0,t,next,elapsed)
      if(s/=0)then
        print *,'CHIMES_TARGET_FAILURE repeat/status=',repeat_index,s
        stop 12
      endif
    enddo
    s=chimes_budget(next,measured,q)
    if(s/=0.or.abs(q)>1d-8.or.any(abs(measured-elements)>1d-8*max(elements,1d-20)))stop 13
    print *,'CHIMES_TARGET_PASS repeats=',target_repeats,' subdt_s=',target_subdt, &
         ' wall_s=',omp_get_wtime()-start,' T_K=',t,' elapsed_s=',elapsed
    stop
  endif
  ! Neutral/molecular and ionized carriers in the same valid H inventory.
  old(1)=.01d0;old(2)=.97d0;old(3)=.01d0;old(138)=.01d0
  do j=1,3
    do i=1,2
      ctl=[dens(j),temps(i),272.7d0,1d10,1d20,0d0,1d0,0d0,.01d0]
      atomic=merge(1,0,i==3)
      start=omp_get_wtime()
!$omp parallel do num_threads(omp_get_max_threads()) schedule(dynamic,1) private(s,t,next,elapsed)
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
